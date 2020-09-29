### HEADER ###
# HOTSPOT PREDICTION
# description: helper script for model fitting
# input: -
# output: -
# author: HR

import os
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
import horovod.tensorflow.keras as hvd


########################################################################################################################
# INITIALIZATION
########################################################################################################################

def print_initialization():
    print("-------------------------------------------------------------------------")
    print("ENVIRONMENT VARIABLES")

    jobname = str(os.environ["SLURM_JOB_NAME"])
    print("JOBNAME: ", jobname)
    nodelist = os.environ["SLURM_JOB_NODELIST"]
    print("NODELIST: ", nodelist)
    nodename = os.environ["SLURMD_NODENAME"]
    print("NODENAME: ", nodename)
    num_nodes = int(os.getenv("SLURM_JOB_NUM_NODES"))
    print("NUMBER OF NODES: ", num_nodes)
    print("NUM GPUs AVAILABLE: ", len(tf.config.experimental.list_physical_devices('GPU')))

    print("-------------------------------------------------------------------------")

    print('HOROVOD CONFIGURATION')

    print('number of Horovod processes on current node: ', hvd.local_size())
    print('rank of the current process: ', hvd.rank())  # node
    print('local rank of the current process: ', hvd.local_rank())  # process on node
    print('Is MPI multi-threading supported? ', hvd.mpi_threads_supported())
    print('Is MPI enabled? ', hvd.mpi_enabled())
    print('Is Horovod compiled with MPI? ', hvd.mpi_built())
    print('Is Horovod compiled with NCCL? ', hvd.nccl_built())

    print("-------------------------------------------------------------------------")


########################################################################################################################
# DATA FORMATTING
########################################################################################################################
def format_input(tokensAndCounts, pseudocounts, relative_dist=False):
    tokens = np.array(tokensAndCounts.loc[:, ['Accession', 'window']], dtype='object')

    dist_N = np.array(tokensAndCounts['dist_N'], dtype='int32')
    dist_C = np.array(tokensAndCounts['dist_C'], dtype='int32')

    if relative_dist:
        dist_N = dist_N / np.array(tokensAndCounts['len_protein'], dtype='int32')
        dist_C = dist_C / np.array(tokensAndCounts['len_protein'], dtype='int32')

    dist = np.array([dist_N, dist_C]).transpose()

    counts = np.array(tokensAndCounts['counts'], dtype='float32')
    counts = np.log2((counts + pseudocounts))

    print('number of features: ', counts.shape[0])
    return tokens, counts, dist


def open_and_format_matrices(group, encoding, spec, extension, windowSize, embeddingDim, relative_dist, pseudocounts=1):
    print(group)
    print(encoding)
    print(extension)

    tokensAndCounts = pd.read_csv(str('/scratch2/hroetsc/Hotspots/data/'+extension+'windowTokens_'+ group + 'ing' + spec + '.csv'))
    tokens, counts, dist = format_input(tokensAndCounts, pseudocounts, relative_dist=relative_dist)

    emb_path = str('/scratch2/hroetsc/Hotspots/data/EMBEDDING_' + encoding + '_'+group+spec+'.dat')
    acc_path = str('/scratch2/hroetsc/Hotspots/data/ACCESSIONS_' + encoding + '_'+group+spec+'.dat')

    no_elements = int(windowSize * embeddingDim)
    chunk_size = int(no_elements * 4)

    embMatrix = [None] * tokens.shape[0]
    accMatrix = [None] * tokens.shape[0]
    chunk_pos = 0

    # open weights and accessions binary file
    with open(emb_path, 'rb') as emin, open(acc_path, 'rb') as ain:
        # loop over files to get elements
        for b in range(tokens.shape[0]):
            emin.seek(chunk_pos, 0)
            dt = np.fromfile(emin, dtype='float32', count=no_elements)

            # get current accession (index)
            ain.seek(int(b * 4), 0)
            cnt_acc = int(np.fromfile(ain, dtype='int32', count=1))
            accMatrix[b] = cnt_acc

            # make sure to pass 4D-Tensor to model: (batchSize, depth, height, width)
            dt = dt.reshape((windowSize, embeddingDim))
            # for 2D convolution --> 4D input:
            embMatrix[b] = np.expand_dims(dt, axis=0)

            # increment chunk position
            chunk_pos += chunk_size

        emin.close()
        ain.close()

    # order tokens and count according to order in embedding matrix
    accMatrix = np.array(accMatrix, dtype='int32')
    embMatrix = np.array(embMatrix, dtype='float32')

    tokens = tokens[accMatrix, :]
    counts = counts[accMatrix]
    dist = dist[accMatrix]

    sort = np.argsort(accMatrix)  # sort to be able to concatenate multiple matrices
    tokens = tokens[sort, :]
    counts = counts[sort]
    dist = dist[sort]
    embMatrix = embMatrix[sort, :, :, :]

    print(counts[:10])

    return tokens, counts, embMatrix, dist


########################################################################################################################
# CALLBACKS
########################################################################################################################

class RestoreBestModel(keras.callbacks.Callback):
    def __init__(self):
        super(RestoreBestModel, self).__init__()
        self.best_weights = None  # best weights

    def on_train_begin(self, logs=None):
        self.best = np.Inf  # initialize best as Inf

    def on_epoch_end(self, epoch, logs=None):
        current = logs.get('val_loss')  # get validation loss
        if np.less(current, self.best):
            self.best = current
            self.best_weights = self.model.get_weights()  # record the best weights

    def on_train_end(self, logs=None):
        self.model.save('/scratch2/hroetsc/Hotspots/results/model/last_model_rank{}.h5'.format(hvd.rank()))

        self.model.set_weights(self.best_weights)
        print('RESTORING WEIGHTS FROM VALIDATION LOSS {}'.format(self.best))



# plan for learning rates
class CosineAnnealing(keras.callbacks.Callback):
    def __init__(self, no_cycles, no_epochs, max_lr):
        super(CosineAnnealing, self).__init__()
        self.no_cycles = no_cycles
        self.no_epochs = no_epochs
        self.max_lr = max_lr
        self.lrates = list()
        self.ensemble = 0

    def cos_annealing(self, epoch, no_cycles, no_epochs, max_lr):
        epochs_per_cycle = np.floor(no_epochs / no_cycles)
        cos_arg = (np.pi * (epoch % epochs_per_cycle)) / (epochs_per_cycle)
        return (max_lr / 2) * (np.cos(cos_arg) + 1)

    def on_epoch_begin(self, epoch, logs=None):
        lr = self.cos_annealing(epoch, self.no_cycles, self.no_epochs, self.max_lr)
        tf.keras.backend.set_value(self.model.optimizer.lr, lr)
        self.lrates.append(lr)

    def on_epoch_end(self, epoch, logs=None):
        lr = float(tf.keras.backend.get_value(self.model.optimizer.learning_rate))

        # if (epoch+1) % self.no_cycles == 0:
        #     self.model.save(
        #         '/scratch2/hroetsc/Hotspots/results/model/ensemble/ensemble_rank{}_epoch{}.h5'.format(hvd.rank(), self.ensemble))
        #     self.ensemble += 1
        #     print('saving model for ensemble prediction')

        print("epoch {}: learning rate is {}".format(epoch, lr))



########################################################################################################################
# SAVE METRICS AND PREDICTION
########################################################################################################################

def save_training_res(model, fit):
    print('----- saving training results -----')

    # save entire model
    model.save('/scratch2/hroetsc/Hotspots/results/model/best_model_rank{}.h5'.format(hvd.rank()))

    if hvd.rank() == 0:
        # save weights
        model.save_weights('/scratch2/hroetsc/Hotspots/results/model/weights.h5')

        # save metrics
        val = []
        name = list(fit.history.keys())
        for i, elem in enumerate(fit.history.keys()):
            val.append(fit.history[elem])

        m = list(zip(name, val))
        m = pd.DataFrame(m)
        pd.DataFrame.to_csv(m, '/scratch2/hroetsc/Hotspots/results/model_metrics.txt', header=False, index=False)


def combine_predictions(outname):
    print('----- combining predictions from all ranks -----')

    for i in range(hvd.size()):
        cnt_pred = pd.read_csv('/scratch2/hroetsc/Hotspots/results/{}_prediction_rank{}.csv'.format(outname, i))
        if i == 0:
            res = pd.DataFrame({'Accession': cnt_pred['Accession'],
                                'window': cnt_pred['window'],
                                'count': cnt_pred['count'],
                                'pred_count': cnt_pred['pred_count'] * (1 / hvd.size())})
        else:
            res['pred_count'] += cnt_pred['pred_count'] * (1 / hvd.size())

    pd.DataFrame.to_csv(res,
                        '/scratch2/hroetsc/Hotspots/results/{}_prediction.csv'.format(outname),
                        index=False)


def ensemble_predictions(emb_test, dist_test, epochs, no_cycles, batchSize):
    print('----- using learning rate schedule minima for ensemble prediction -----')

    ranks = hvd.size()
    models_per_rank = int(epochs / no_cycles)

    all_predictions = []

    for i in range(ranks):
        for j in range(models_per_rank):

            print('prediction on rank {}, model no. {}'.format(i, j))

            cnt_model = keras.models.load_model(
                '/scratch2/hroetsc/Hotspots/results/model/ensemble/ensemble_rank{}_epoch{}.h5'.format(i, j))
            cnt_prediction = cnt_model.predict(x=[emb_test, dist_test],
                                               batch_size=batchSize,
                                               verbose=0,
                                               max_queue_size=64)

            all_predictions.append(cnt_prediction.flatten())

    # average of all predictions
    all_predictions = np.array(all_predictions)
    average_prediction = np.mean(all_predictions, axis=0)

    return average_prediction