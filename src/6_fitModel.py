### HEADER ###
# HOTSPOT PREDICTION
# description: fit neural network to predict hotspot density along the human proteome
# input: table with amino acids in sliding windows and scores, window matrices
# output: model, metrics, predictions
# author: HR

import os
import numpy as np
import pandas as pd

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import tensorflow.keras.backend as K

K.clear_session()
import horovod.tensorflow.keras as hvd

hvd.init()

from model_helper import print_initialization, open_and_format_matrices, \
    RestoreBestModel, CosineAnnealing, \
    save_training_res

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

spec = '_25aa_mean'

windowSize = 25
embeddingDim = 20
sgt_dim = 400

epochs = 100
pseudocounts = 1  # !!!
no_cycles = int(epochs / 10)

batch_size = 256
max_lr = 0.001
starting_filter = 16
kernel_size = 3
block_size = 5
num_blocks = 4

dense_nodes = 1024  # !!!
noise_sd = 0.3  # !!!

print('number of epochs, adjusted by number of GPUs: ', epochs)
print('batch size, adjusted by number of GPUs: ', batch_size)
print('number of learning rate cycles: ', no_cycles)
print('maximal learning rate in schedule', max_lr)
print('number of pseudocounts: ', pseudocounts)
print("-------------------------------------------------------------------------")

########## part 1: fit model ##########
### INPUT ###

best_model_path = '/scratch2/hroetsc/Hotspots/results/model/best_model_rank{}.h5'.format(hvd.rank())
last_model_path = '/scratch2/hroetsc/Hotspots/results/model/last_model_rank{}.h5'.format(hvd.rank())

print('#####')
# print('ONE HOT AND AA INDEX ON DIFFERENT RANKS')
print('AA INDEX')
print('no extension')
print('#####')

# enc = 'oneHOT' if hvd.rank() % 2 == 0 else 'AAindex'
enc = 'AAindex'

tokens, counts, labels, emb, dist, wholeSeq = open_and_format_matrices(group='train',
                                                                       encoding=enc,
                                                                       spec=spec,
                                                                       extension='',
                                                                       windowSize=windowSize,
                                                                       embeddingDim=embeddingDim,
                                                                       sgt_dim=sgt_dim,
                                                                       relative_dist=False,
                                                                       protein_norm=False,
                                                                       log_counts=True,
                                                                       pseudocounts=pseudocounts)
tokens_test, counts_test, labels_test, emb_test, dist_test, wholeSeq_test = open_and_format_matrices(group='test',
                                                                                                     encoding=enc,
                                                                                                     spec=spec,
                                                                                                     extension='',
                                                                                                     windowSize=windowSize,
                                                                                                     embeddingDim=embeddingDim,
                                                                                                     sgt_dim=sgt_dim,
                                                                                                     relative_dist=False,
                                                                                                     protein_norm=False,
                                                                                                     log_counts=True,
                                                                                                     pseudocounts=pseudocounts)


bias_initializer = np.mean(counts)  # initialise bias with mean of all counts to prevent model from learning the bias
# bias_initializer = 0

### MAIN PART ###
def build_and_compile_model(max_lr, starting_filter, kernel_size, block_size, dense_nodes, num_blocks,
                            include_distance=True, whole_seq=False,
                            additive_noise_input=True, additive_noise_counts=True,
                            weight_constraint=False):
    num_blocks_list = [block_size] * num_blocks
    dilation_rate_list = [1] * num_blocks

    regr_activation = 'relu'

    print(locals())

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
        y = bn_relu(inp_layer)
        y = layers.Conv2D(filters=filters,
                          kernel_size=kernel_size,
                          strides=(1 if not downsample else 2),
                          padding='same',
                          kernel_initializer=tf.keras.initializers.HeNormal(),
                          data_format='channels_first')(y)
        y = bn_relu(y)
        y = layers.Conv2D(filters=filters,
                          kernel_size=kernel_size,
                          strides=1,
                          dilation_rate=dilation_rate,
                          padding='same',
                          kernel_initializer=tf.keras.initializers.HeNormal(),
                          data_format='channels_first')(y)

        if downsample:
            inp_layer = layers.Conv2D(filters=filters,
                                      kernel_size=1,
                                      strides=2,
                                      padding='same',
                                      kernel_initializer=tf.keras.initializers.HeNormal(),
                                      data_format='channels_first')(inp_layer)

        out = layers.Add()([inp_layer, y])

        return out

    ## input
    tf.print('model input')
    inp = keras.Input(shape=(1, windowSize, embeddingDim),
                      name='input_window')
    dist = keras.Input(shape=(2),
                       name='distance')
    sgt = keras.Input(shape=(1, sgt_dim),
                      name='whole_sequence_embedding')

    ## noise layer
    if additive_noise_input:
        inp0 = layers.GaussianNoise(stddev=noise_sd)(inp)
    else:
        inp0 = inp

    ## convolutional layers (ResNet)
    tf.print('residual blocks')
    tf.print('STRUCTURE OF RESIDUAL BLOCK: D')

    # structure of residual blocks:
    # a) original: weight-BN-ReLU-weight-BN-addition-ReLU
    # b) BN after addition: weight-BN-ReLU-weight-addition-BN-ReLU
    # c) ReLU before addition: weight-BN-ReLU-weight-BN-ReLU-addition
    # d) full pre-activation (SpliceAI): BN-ReLU-weight-BN-ReLU-weight-addition

    # initial convolution
    t = layers.BatchNormalization(trainable=True)(inp0)
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

    # concatenate with distance
    if include_distance:
        flat = layers.Concatenate()([flat, dist])

    # add whole sequence embedding
    if whole_seq:
        sgt_dense = layers.Flatten()(sgt)
        sgt_dense = dense_layer(sgt_dense, nodes=256)
        sgt_dense = dense_layer(sgt_dense, nodes=64)
        flat = layers.Concatenate()([flat, sgt_dense])

    ## dense layers
    tf.print('dense layers')
    # fully-connected layer
    dense1 = dense_layer(flat, dense_nodes)

    # regression problem
    regr_norm = layers.BatchNormalization(trainable=True)(dense1)
    if additive_noise_counts:
        regr_norm = layers.GaussianNoise(stddev=noise_sd)(regr_norm)

    regr = layers.Dense(1,
                        activation=regr_activation,
                        kernel_constraint=keras.constraints.MinMaxNorm(min_value=0,
                                                                       max_value=1) if weight_constraint else None,
                        kernel_initializer=tf.keras.initializers.HeNormal(),
                        bias_initializer=keras.initializers.Constant(value=bias_initializer),
                        name='regression')(regr_norm)

    # classification problem
    classif_norm = layers.BatchNormalization(trainable=True)(dense1)
    classif = layers.Dense(1,
                           activation='sigmoid',
                           kernel_initializer=tf.keras.initializers.GlorotUniform(),
                           bias_initializer=tf.keras.initializers.Constant(value=bias_initializer),
                           name='classification')(classif_norm)

    ## concatenate to model
    model = keras.Model(inputs=[inp, dist, sgt], outputs=[regr, classif])

    ## compile model
    tf.print('compile model')
    opt = tf.keras.optimizers.Adam(learning_rate=max_lr)
    opt = hvd.DistributedOptimizer(opt)

    losses = {'regression': keras.losses.MeanSquaredError(),
              'classification': keras.losses.BinaryCrossentropy(label_smoothing=.1)}
    loss_weights = {'regression': 1.0,
                    'classification': 0.1}
    metrics = {'regression': ['mean_absolute_error', 'mean_absolute_percentage_error'],
               'classification': [keras.metrics.AUC(curve='ROC', name='roc_auc'),
                                  keras.metrics.AUC(curve='PR', name='pr_auc'),
                                  keras.metrics.BinaryAccuracy(name='accuracy')]}

    model.compile(loss=losses,
                  loss_weights=loss_weights,
                  optimizer=opt,
                  metrics=metrics,
                  experimental_run_tf_function=False)

    # for reproducibility during optimization
    tf.print('......................................................')
    tf.print('MIXTURE OF RESNET AND SPLICEAI')
    tf.print('optimizer: Adam (distributed)')
    tf.print('loss: mean squared error and binary crossentropy')
    tf.print('regression and classification problem')
    tf.print('skip-connections between blocks')
    tf.print('channels: first')
    tf.print('activation function: leaky relu/he_normal')
    tf.print('regularization: none')
    tf.print('using batch normalization: yes')
    tf.print('using Dropout layer: no')
    tf.print('......................................................')

    return model


#### train model
print('MODEL TRAINING')
# define callbacks
callbacks = [RestoreBestModel(last_model_path),
             CosineAnnealing(no_cycles=no_cycles, no_epochs=epochs, max_lr=max_lr),
             hvd.callbacks.BroadcastGlobalVariablesCallback(0),
             keras.callbacks.ModelCheckpoint(filepath=last_model_path,
                                             verbose=1,
                                             monitor='val_loss',
                                             save_freq='epoch')]

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
                                include_distance=True, whole_seq=False,
                                additive_noise_input=True, additive_noise_counts=True,
                                weight_constraint=False)

# print('CONTINUE TRAINING OF EXISTING MODEL ON NEW DATA')
# model = keras.models.load_model(last_model_path)

if hvd.rank() == 0:
    model.summary()
    print('train for {}, validate for {} steps per epoch'.format(steps, val_steps))

fit = model.fit(x=[emb, dist, wholeSeq],
                y=[counts, labels],
                batch_size=batch_size,
                validation_data=([emb_test, dist_test, wholeSeq_test], [counts_test, labels_test]),
                validation_batch_size=batch_size,
                steps_per_epoch=steps,
                validation_steps=val_steps,
                epochs=epochs,
                callbacks=callbacks,
                initial_epoch=0,
                max_queue_size=64,
                verbose=2 if hvd.rank() == 0 else 0,
                shuffle=True)

### OUTPUT ###
print('SAVE MODEL AND METRICS')
save_training_res(model, fit, best_model_path)

########## part 2: make prediction ##########

print('LOAD MODELS')
last_model = keras.models.load_model(last_model_path)
best_model = keras.models.load_model(best_model_path)


def make_prediction(model, outname):
    print(outname)

    outpath = '/scratch2/hroetsc/Hotspots/results/{}_prediction_rank{}.csv'.format(outname, hvd.rank())
    if os.path.exists(outpath):
        os.remove(outpath)

    pred = model.predict(x=[emb_test, dist_test, wholeSeq_test],
                         batch_size=batch_size,
                         verbose=1 if hvd.rank() == 0 else 0,
                         max_queue_size=64)

    print('counts:')
    print(counts_test)
    print('prediction:')
    pred_counts = np.array(pred[0].flatten())
    print(pred_counts)

    print('labels:')
    print(labels_test)
    print('prediction:')
    pred_labels = np.array(pred[1].flatten())
    print(pred_labels)

    # merge actual and predicted counts
    prediction = pd.DataFrame({"Accession": tokens_test[:, 0],
                               "window": tokens_test[:, 1],
                               "count": counts_test,
                               "pred_count": pred_counts,
                               "class": labels_test,
                               "pred_class": pred_labels})

    pd.DataFrame.to_csv(prediction, outpath, index=False)


print('MAKE PREDICTION')
make_prediction(model=best_model, outname='best_model')
make_prediction(model=last_model, outname='last_model')
