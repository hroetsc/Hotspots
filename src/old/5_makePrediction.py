### HEADER ###
# HOTSPOT PREDICTION
# description: predict hotspot density along the human proteome
# input: model (last model and best model)
# output: predictions for test data set
# author: HR

import os
import numpy as np
import pandas as pd

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

tf.keras.backend.clear_session()

import horovod.tensorflow.keras as hvd
hvd.init()

from model_helper import print_initialization, open_and_format_matrices, combine_predictions

print_initialization()

print('GPU INITIALIZATION')
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)
if gpus:
    tf.config.experimental.set_visible_devices(gpus[hvd.local_rank()], 'GPU')


### INPUT ###
embeddingDim = 20
windowSize = 25
batch_size = 128

print('OPEN TEST DATA SET')
enc = 'oneHOT' if hvd.rank() % 2 == 0 else 'AAindex'
tokens_test, counts_test, emb_test, dist_test = open_and_format_matrices(group='test', encoding=enc,
                                                                         spec='50-2',
                                                                         extension='',
                                                                         windowSize=windowSize,
                                                                         embeddingDim=embeddingDim,
                                                                         relative_dist=False)

print('LOAD MODELS')
last_model = keras.models.load_model('/scratch2/hroetsc/Hotspots/results/model/last_model_rank{}.h5'.format(hvd.rank()))
best_model = keras.models.load_model('/scratch2/hroetsc/Hotspots/results/model/best_model_rank{}.h5'.format(hvd.rank()))

### MAIN PART ###
# make prediction
def make_prediction(model, outname):
    print(outname)

    outpath = '/scratch2/hroetsc/Hotspots/results/{}_prediction_rank{}.csv'.format(outname, hvd.rank())
    if os.path.exists(outpath):
        os.remove(outpath)

    pred = model.predict(x=[emb_test, dist_test],
                         batch_size=batch_size,
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

    pd.DataFrame.to_csv(prediction, outpath, index=False)


print('MAKE PREDICTION')
make_prediction(model=best_model, outname='best_model')
make_prediction(model=last_model, outname='last_model')

print('COMBINE PREDICTIONS')
if hvd.rank() == 0:
    combine_predictions(outname='best_model')
    combine_predictions(outname='last_model')
