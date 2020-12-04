### HEADER ###
# HOTSPOT PREDICTION
# description: set up neural network to predict hotspot density along the human proteome
# input: table with amino acids in sliding windows and scores, window matrices
# output: model, metrics, predictions for test data set
# author: HR

import numpy as np
import pandas as pd

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

tf.keras.backend.clear_session()
import horovod.tensorflow.keras as hvd

hvd.init()

from model_helper import print_initialization, open_and_format_matrices, \
    RestoreBestModel, CosineAnnealing, \
    save_training_res, combine_predictions, ensemble_predictions

print_initialization()

### initialize GPU training environment
print('GPU INITIALIZATION')
# pin GPUs (each GPU gets single process)
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)
if gpus:
    tf.config.experimental.set_visible_devices(gpus[hvd.local_rank()], 'GPU')

### HYPERPARAMETERS ###
print('HYPERPARAMETERS')
embeddingDim = 20
windowSize = 25

epochs = 10
batchSize = 32
pseudocounts = 1

no_cycles = int(epochs / 5)

print('number of epochs, adjusted by number of GPUs: ', epochs)
print('batch size, adjusted by number of GPUs: ', batchSize)
print('number of learning rate cycles: ', no_cycles)
print('number of pseudocounts: ', pseudocounts)
print("-------------------------------------------------------------------------")

########## part 1: fit model ##########
### INPUT ###
print('#####')
print('ONE HOT AND AA INDEX ON DIFFERENT RANKS')
print('no extension')
print('#####')

enc = 'oneHOT' if hvd.rank() % 2 == 0 else 'AAindex'

tokens, counts, emb, dist = open_and_format_matrices(group='train', encoding=enc, spec='05',
                                                     extension='',
                                                     windowSize=windowSize, embeddingDim=embeddingDim,
                                                     relative_dist=False)
tokens_test, counts_test, emb_test, dist_test = open_and_format_matrices(group='test', encoding=enc,
                                                                         spec='05',
                                                                         extension='',
                                                                         windowSize=windowSize,
                                                                         embeddingDim=embeddingDim,
                                                                         relative_dist=False)

bias_initializer = np.mean(counts)


### MAIN PART ###
def build_and_compile_model(max_lr, starting_filter, kernel_size, block_size, dense_nodes, num_blocks,
                            include_distance=False):

    num_blocks_list = block_size
    dilation_rate_list = [1] * num_blocks

    # build dense relu layers with batch norm and dropout
    def dense_layer(prev_layer, nodes):
        norm = layers.BatchNormalization(trainable=True)(prev_layer)
        dense = layers.Dense(nodes, activation='relu',
                             kernel_initializer=tf.keras.initializers.HeNormal())(norm)
        return dense

    # activation and batch normalization
    def bn_relu(inp_layer):
        bn = layers.BatchNormalization(trainable=True)(inp_layer)
        relu = layers.LeakyReLU()(bn)
        return relu

    # residual blocks (convolutions)
    def residual_block(inp_layer, downsample, filters, kernel_size, dilation_rate):
        y = layers.Conv2D(filters=filters,
                          kernel_size=kernel_size,
                          strides=(1 if not downsample else 2),
                          padding='same',
                          kernel_initializer=tf.keras.initializers.HeNormal(),
                          data_format='channels_first')(inp_layer)
        y = bn_relu(y)
        y = layers.Conv2D(filters=filters,
                          kernel_size=kernel_size,
                          strides=1,
                          dilation_rate=dilation_rate,
                          padding='same',
                          kernel_initializer=tf.keras.initializers.HeNormal(),
                          data_format='channels_first')(y)
        y = layers.BatchNormalization(trainable=True)(y)

        if downsample:
            inp_layer = layers.Conv2D(filters=filters,
                                      kernel_size=1,
                                      strides=2,
                                      padding='same',
                                      kernel_initializer=tf.keras.initializers.HeNormal(),
                                      data_format='channels_first')(inp_layer)

        out = layers.Add()([inp_layer, y])
        out = layers.LeakyReLU()(out)

        return out

    ## input
    tf.print('model input')
    inp = keras.Input(shape=(1, windowSize, embeddingDim),
                      name='input_window')
    dist = keras.Input(shape=(2),
                       name='distance')

    ## convolutional layers (ResNet)
    tf.print('residual blocks')

    # structure of residual blocks:
    # a) original: weight-BN-ReLU-weight-BN-addition-ReLU --> currently used
    # b) BN after addition: weight-BN-ReLU-weight-addition-BN-ReLU
    # c) ReLU before addition: weight-BN-ReLU-weight-BN-ReLU-addition
    # d) full pre-activation (SpliceAI): BN-ReLU-weight-BN-ReLU-weight-addition

    # initial convolution
    t = layers.BatchNormalization(trainable=True)(inp)
    t = layers.Conv2D(filters=starting_filter,
                      kernel_size=kernel_size,
                      strides=2,
                      padding='same',
                      kernel_initializer=tf.keras.initializers.HeNormal(),
                      data_format='channels_first')(t)
    t = bn_relu(t)

    # residual blocks
    for i in range(len(num_blocks_list)):
        no_blocks = num_blocks_list[i]
        dil_rate = dilation_rate_list[i]

        t_shortcut = layers.Conv2D(filters=starting_filter,
                                   kernel_size=kernel_size,
                                   strides=(1 if i == 0 else 2),
                                   padding='same',
                                   kernel_initializer=tf.keras.initializers.HeNormal(),
                                   data_format='channels_first')(t)

        for j in range(no_blocks):
            t = residual_block(t,
                               downsample=(j == 0 and i != 0),
                               filters=starting_filter,
                               kernel_size=kernel_size,
                               dilation_rate=dil_rate)

        t = layers.Add()([t, t_shortcut])
        starting_filter *= 2

    t = layers.AveragePooling2D(pool_size=4,
                                data_format='channels_first',
                                padding='same')(t)
    flat = layers.Flatten()(t)

    if include_distance:
        flat = layers.Concatenate()([flat, dist])

    ## dense layers
    tf.print('dense layers')
    # fully-connected layer
    dense1 = dense_layer(flat, dense_nodes)

    out_norm = layers.BatchNormalization(trainable=True)(dense1)
    out = layers.Dense(1,
                       activation='linear',
                       kernel_initializer=tf.keras.initializers.GlorotUniform(),
                       bias_initializer=keras.initializers.Constant(value=bias_initializer),
                       name='output')(out_norm)

    ## concatenate to model
    model = keras.Model(inputs=[inp, dist], outputs=out)

    ## compile model
    tf.print('compile model')
    opt = tf.keras.optimizers.Adam(learning_rate=max_lr)
    opt = hvd.DistributedOptimizer(opt)

    model.compile(loss=keras.losses.MeanSquaredError(),
                  optimizer=opt,
                  metrics=['mean_absolute_error', 'mean_absolute_percentage_error', 'accuracy'],
                  experimental_run_tf_function=False)

    # for reproducibility during optimization
    tf.print('......................................................')
    tf.print('MIXTURE OF RESNET AND SPLICEAI')
    tf.print('optimizer: Adam (not distributed)')
    tf.print('loss: mean squared error')
    tf.print('skip-connections between blocks')
    tf.print('channels: first')
    tf.print('activation function: leaky relu/he_normal')
    tf.print('regularization: none')
    tf.print('using batch normalization: yes')
    tf.print('using Dropout layer: no')
    tf.print('......................................................')

    return model


# make prediction
def make_prediction(model, outname):
    pred = model.predict(x=[emb_test, dist_test],
                         batch_size=batchSize,
                         verbose=1 if hvd.rank() == 0 else 0,
                         max_queue_size=64)
    print('counts:')
    print(counts_test)
    print('prediction:')
    pred_counts = np.array(pred.flatten())
    print(pred_counts)

    # merge actual and predicted counts
    prediction = pd.DataFrame({"Accession": tokens_test[:, 0],
                               "window": tokens_test[:, 1],
                               "count": counts_test,
                               "pred_count": pred_counts})

    pd.DataFrame.to_csv(prediction,
                        '/scratch2/hroetsc/Hotspots/results/{}_prediction_rank{}.csv'.format(outname, hvd.rank()),
                        index=False)


def fitness(max_lr, starting_filter, kernel_size, block_size, dense_nodes, num_blocks,
            batch_size,
            include_distance=False):

    tf.keras.backend.clear_session()
    tf.compat.v1.reset_default_graph()

    #### train model
    print('MODEL TRAINING')
    # define callbacks
    callbacks = [RestoreBestModel(),
                 CosineAnnealing(no_cycles=no_cycles, no_epochs=epochs, max_lr=max_lr),
                 hvd.callbacks.BroadcastGlobalVariablesCallback(0)]

    # define number of steps - make sure that no. of steps is the same for all ranks!
    # otherwise, stalled ranks problem might occur
    steps = int(np.ceil(counts.shape[0] / batch_size))
    val_steps = int(np.ceil(counts_test.shape[0] / batch_size))

    # adjust by number of GPUs
    steps = int(np.ceil(steps / hvd.size()))
    val_steps = int(np.ceil(val_steps / hvd.size()))

    ## fit model
    model = build_and_compile_model(max_lr=max_lr, starting_filter=starting_filter, kernel_size=kernel_size,
                                    block_size=block_size, dense_nodes=dense_nodes, num_blocks=num_blocks,
                                    include_distance=include_distance)

    if hvd.rank() == 0:
        model.summary()
        print('train for {}, validate for {} steps per epoch'.format(steps, val_steps))

    fit = model.fit(x=[emb, dist],
                    y=counts,
                    batch_size=batch_size,
                    validation_data=([emb_test, dist_test], counts_test),
                    validation_batch_size=batch_size,
                    steps_per_epoch=steps,
                    validation_steps=val_steps,
                    epochs=epochs,
                    callbacks=callbacks,
                    initial_epoch=0,
                    max_queue_size=64,
                    verbose=2 if hvd.rank() == 0 else 0,
                    shuffle=True)

    loss = min(fit.history['val_loss'])
    loss_str = str(loss).replace(".", "")

    # model metrics
    val = []
    name = list(fit.history.keys())
    for i, elem in enumerate(fit.history.keys()):
        val.append(fit.history[elem])

    m = list(zip(name, val))
    m = pd.DataFrame(m)


    print('MAKE PREDICTION')
    make_prediction(model=model, outname='best_model_opt')

    if hvd.rank() == 0:

        print('#######################################################################################################')
        print('#######################################################################################################')
        print('LOCALS:')
        print(locals())
        print('#######################################################################################################')
        print('VALIDATION LOSS: ', loss)
        print('#######################################################################################################')
        print('#######################################################################################################')

        pd.DataFrame.to_csv(m, '/scratch2/hroetsc/Hotspots/results/opt_ResNet_metrics_{}.txt'.format(loss_str),
                            header=False,
                            index=False)

        combine_predictions(outname='best_model_opt', loss_str=loss_str)

    tf.keras.backend.clear_session()
    tf.compat.v1.reset_default_graph()


# hyperparameters
batch_sizes = [16, 32, 64, 128, 256]

for b in batch_sizes:
    fitness(max_lr=0.001, starting_filter=16, kernel_size=3, block_size=[5,5,5,5],
            dense_nodes=1024, num_blocks=4, batch_size=b, include_distance=True)




tf.keras.backend.clear_session()
