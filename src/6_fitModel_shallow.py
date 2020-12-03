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
scale_counts = False
oversample_large = False

windowSize = 25
embeddingDim = 20
sgt_dim = 400

epochs = 100
pseudocounts = 1   # !!!
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

# oversample large training counts
if oversample_large:
    k1, k2, k3, k4, k5, k6, k7 = np.where(counts > 1)[0], np.where(counts > 1.2)[0], np.where(counts > 1.4)[0], \
                                 np.where(counts > 1.6)[0], np.where(counts > 1.8)[0], np.where(counts > 2)[0], \
                                 np.where(counts > 3)[0]

    ind = np.concatenate([k1, k2, k3, k4, k5, k6, k7])

    # concatenate
    tokens = np.concatenate([tokens, tokens[ind, :]], axis=0)
    counts = np.concatenate([counts, counts[ind]], axis=0)
    labels = np.concatenate([labels, labels[ind]], axis=0)
    emb = np.concatenate([emb, emb[ind, :, :, :]], axis=0)
    dist = np.concatenate([dist, dist[ind, :]], axis=0)
    wholeSeq = np.concatenate([wholeSeq, wholeSeq[ind, :]], axis=0)

    # shuffle
    rnd = np.random.randint(low=0, high=counts.shape[0], size=counts.shape[0])
    tokens, counts, labels, emb, dist, wholeSeq = tokens[rnd, :], counts[rnd], labels[rnd], \
                                                  emb[rnd, :, :, :], dist[rnd, :], wholeSeq[rnd, :]


# scale counts between 0 and 1 for regression
def scaling(x0):
    x1 = (x0 - np.min(x0)) / (np.max(x0) - np.min(x0))
    return x1, np.max(x0), np.min(x0)


def backtransformation(x1, max_x0, min_x0):
    x0 = x1 * (max_x0 - min_x0) + min_x0
    return x0


if scale_counts:
    print('scaling counts between 0 and 1')
    counts, max_counts, min_counts = scaling(counts)
    counts_test, max_counts_test, min_counts_test = scaling(counts_test)

    print('max train count:', max_counts, ', min train count:', min_counts)
    print('max test count:', max_counts_test, ', min test count:', min_counts_test)

bias_initializer = np.mean(counts)  # initialise bias with mean of all counts to prevent model from learning the bias
# bias_initializer = 0


### MAIN PART ###
# custom loss function
def SMAPE(y_true, y_pred):
    numerator = K.sum(K.square(y_pred - y_true))
    denominator = K.sum(y_true + y_pred)
    return numerator / denominator

def build_and_compile_model(max_lr, starting_filter, kernel_size, block_size, dense_nodes, num_blocks,
                            include_distance=True, whole_seq=False,
                            additive_noise_input=True, additive_noise_counts=True,
                            weight_constraint=False):
    num_blocks_list = [block_size] * num_blocks
    dilation_rate_list = [1] * num_blocks

    regr_activation = 'sigmoid' if scale_counts else 'relu'

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


    ## input
    tf.print('model input')
    inp = keras.Input(shape=(1, windowSize, embeddingDim),
                      name='input_window')
    dist = keras.Input(shape=(2),
                       name='distance')
    sgt = keras.Input(shape=(sgt_dim),
                      name='whole_sequence_embedding')

    ## noise layer
    if additive_noise_input:
        inp0 = layers.GaussianNoise(stddev=noise_sd)(inp)
    else:
        inp0 = inp

    ## convolutional layers
    # initial convolution
    t = layers.BatchNormalization(trainable=True)(inp0)
    t = layers.Conv2D(filters=starting_filter,
                      kernel_size=kernel_size,
                      strides=2,
                      padding='same',
                      kernel_initializer=tf.keras.initializers.HeNormal(),
                      data_format='channels_first')(t)
    t = bn_relu(t)

    # second convolution
    t = layers.Conv2D(filters=int(starting_filter*2),
                      kernel_size=kernel_size,
                      strides=2,
                      padding='same',
                      kernel_initializer=tf.keras.initializers.HeNormal(),
                      data_format='channels_first')(t)
    t = bn_relu(t)

    t = layers.AveragePooling2D(pool_size=4,
                                data_format='channels_first',
                                padding='same')(t)
    flat = layers.Flatten()(t)

    # concatenate with distance
    if include_distance:
        flat = layers.Concatenate()([flat, dist])

    # add whole sequence embedding
    if whole_seq:
        sgt_dense = dense_layer(sgt, nodes=256)
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
                    'classification': 0.05 if scale_counts else 0.2}
    metrics = {'regression': ['mean_squared_error', 'mean_absolute_error', 'mean_absolute_percentage_error'],
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
    tf.print('SHALLOW CONV NET WITH 2 CONVOLUTIONS + POOLING')
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
    if scale_counts:
        prediction = pd.DataFrame({"Accession": tokens_test[:, 0],
                                   "window": tokens_test[:, 1],
                                   "count": counts_test,
                                   "pred_count": pred_counts,
                                   "count_transformed": backtransformation(counts_test, max_counts_test,
                                                                           min_counts_test),
                                   "pred_count_transformed": backtransformation(pred_counts, max_counts_test,
                                                                                min_counts_test),
                                   "class": labels_test,
                                   "pred_class": pred_labels})
    else:
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

