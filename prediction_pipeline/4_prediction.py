### HEADER ###
# HOTSPOT PREDICTION: VALIDATION
# description: load pre-trained models and predict counts for unknown sequences
# input: table with amino acids in sliding windows and scores, window matrices, whole sequence encodings
# output: predictions
# author: HR

import os
import numpy as np
import pandas as pd

import tensorflow as tf
from tensorflow import keras

tf.keras.backend.clear_session()
import horovod.tensorflow.keras as hvd

hvd.init()

from helper import open_and_format_matrices


### initialize GPU training environment
print('GPU INITIALIZATION')
# pin GPUs (each GPU gets single process)
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)
if gpus:
    tf.config.experimental.set_visible_devices(gpus[hvd.local_rank()], 'GPU')

### HYPERPARAMETERS ###
spec = '_50aa_corona'

windowSize = 50
embeddingDim = 20
sgt_dim = 400
batch_size = 32

### INPUT ###
best_model_path = '/scratch2/hroetsc/Hotspots/results/model/best_model_rank{}.h5'.format(hvd.rank())
last_model_path = '/scratch2/hroetsc/Hotspots/results/model/last_model_rank{}.h5'.format(hvd.rank())

print('#####')
print('ONE HOT AND AA INDEX ON DIFFERENT RANKS')
print('no extension')
print('#####')

enc = 'oneHOT' if hvd.rank() % 2 == 0 else 'AAindex'

tokens, emb, dist, wholeSeq = open_and_format_matrices(encoding=enc,
                                                       spec=spec,
                                                       windowSize=windowSize,
                                                       embeddingDim=embeddingDim,
                                                       sgt_dim=sgt_dim,
                                                       relative_dist=False)


### MAIN PART ###
print('LOAD MODELS')
last_model = keras.models.load_model(last_model_path)
best_model = keras.models.load_model(best_model_path)


def make_prediction(model, outname):
    print(outname)

    outpath = '/scratch2/hroetsc/Hotspots/results/' + outname + '_prediction' + spec + '_rank' + str(hvd.rank()) + '.csv'
    if os.path.exists(outpath):
        os.remove(outpath)

    pred = model.predict(x=[emb, dist, wholeSeq],
                         batch_size=batch_size,
                         verbose=1 if hvd.rank() == 0 else 0,
                         max_queue_size=64)

    print('prediction:')
    pred_counts = np.array(pred.flatten())
    print(pred_counts)

    # merge actual and predicted counts
    prediction = pd.DataFrame({'organism': tokens[:, 0],
                               "accession": tokens[:, 1],
                               'annotation': tokens[:, 2],
                               "window": tokens[:, 3],
                               "pred_count": pred_counts})

    pd.DataFrame.to_csv(prediction, outpath, index=False)


### OUTPUT ###
print('MAKE PREDICTION')
make_prediction(model=best_model, outname='best_model')
make_prediction(model=last_model, outname='last_model')

tf.keras.backend.clear_session()
