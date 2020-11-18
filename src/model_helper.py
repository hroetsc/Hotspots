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
def format_input(tokensAndCounts, pseudocounts, relative_dist=False, protein_norm=False, log_counts=False):
    tokens = np.array(tokensAndCounts[['Accession', 'window']])

    dist_N = np.array(tokensAndCounts['dist_N'], dtype='int32')
    dist_C = np.array(tokensAndCounts['dist_C'], dtype='int32')

    if relative_dist:
        dist_N = dist_N / np.array(tokensAndCounts['len_protein'], dtype='int32')
        dist_C = dist_C / np.array(tokensAndCounts['len_protein'], dtype='int32')

    dist = np.array([dist_N, dist_C]).transpose()

    counts = np.array(tokensAndCounts['counts'], dtype='float32')
    labels = np.where(counts > 0, 1, 0)

    if protein_norm:
        print('normalise counts by protein length')
        acc = set(tokens[:, 0])
        acc = [i for i in acc]

        for i in acc:
            cnt = np.where(tokens[:, 0] == i)[0]
            counts[cnt] = counts[cnt] / len(cnt)

    if log_counts:
        counts = np.log2((counts + pseudocounts))
        print('using log-transformed counts')
    else:
        print('using raw counts')

    print('number of features: ', counts.shape[0])

    return tokens, counts, labels, dist


def open_and_format_matrices(group, encoding, spec, extension,
                             windowSize, embeddingDim, sgt_dim,
                             relative_dist, protein_norm, log_counts,
                             pseudocounts=1):
    print(group)
    print(encoding)
    print(spec)
    print(extension)

    # accession, windows and counts
    tokensAndCounts = pd.read_csv(
        str('/scratch2/hroetsc/Hotspots/data/' + extension + 'windowTokens_' + group + 'ing' + spec + '.csv'))
    tokens, counts, labels, dist = format_input(tokensAndCounts,
                                                pseudocounts,
                                                relative_dist=relative_dist,
                                                protein_norm=protein_norm,
                                                log_counts=log_counts)

    # get embeddings
    emb_path = str('/scratch2/hroetsc/Hotspots/data/EMBEDDING_' + encoding + '_' + group + spec + '.dat')
    acc_path = str('/scratch2/hroetsc/Hotspots/data/ACCESSIONS_' + encoding + '_' + group + spec + '.dat')
    whole_seq_path = str('/scratch2/hroetsc/Hotspots/data/WHOLE-SEQ_' + group + spec + '.dat')
    whole_seq_ids = str('/scratch2/hroetsc/Hotspots/data/WHOLE-SEQ_ID_' + group + spec + '.dat')

    no_elements = int(windowSize * embeddingDim)

    embMatrix = [None] * tokens.shape[0]
    accMatrix = [None] * tokens.shape[0]
    wholeSeq = [None] * tokens.shape[0]
    wholeSeq_IDs = [None] * tokens.shape[0]

    # open weights and accessions binary file
    with open(emb_path, 'rb') as emin, open(acc_path, 'rb') as ain, \
            open(whole_seq_path, 'rb') as win, open(whole_seq_ids, 'rb') as win_id:
        # loop over files to get elements
        for b in range(tokens.shape[0]):
            # window embeddings
            emin.seek(int(b * 4 * no_elements), 0)
            dt = np.fromfile(emin, dtype='float32', count=no_elements)
            # make sure to pass 4D-Tensor to model: (batchSize, depth, height, width)
            dt = dt.reshape((windowSize, embeddingDim))
            # for 2D convolution --> 4D input:
            embMatrix[b] = np.expand_dims(dt, axis=0)

            # get current accession (index)
            ain.seek(int(b * 4), 0)
            accMatrix[b] = int(np.fromfile(ain, dtype='int32', count=1))

            # get whole sequence embedding
            win.seek(int(b * 4 * sgt_dim), 0)
            cnt_ws = np.fromfile(win, dtype='float32', count=sgt_dim)
            wholeSeq[b] = cnt_ws

            # get ID of sequence
            win_id.seek(int(b * 4), 0)
            wholeSeq_IDs[b] = int(np.fromfile(win_id, dtype='int32', count=1))

        emin.close()
        ain.close()

    # np arrays
    accMatrix = np.array(accMatrix, dtype='int32')
    embMatrix = np.array(embMatrix, dtype='float32')
    wholeSeq_IDs = np.array(wholeSeq_IDs, dtype='int32')
    wholeSeq = np.array(wholeSeq, dtype='float32')

    # order all data frames according to occurrence in windowTokens
    # (information kept in ids)
    embMatrix = embMatrix[np.argsort(accMatrix), :, :, :]
    wholeSeq = wholeSeq[np.argsort(wholeSeq_IDs), :]

    # should be in the same order as tokens, counts and distance

    print(counts[:10])

    return tokens, counts, labels, embMatrix, dist, wholeSeq


########################################################################################################################
# CALLBACKS
########################################################################################################################

class RestoreBestModel(keras.callbacks.Callback):
    def __init__(self, outpath):
        super(RestoreBestModel, self).__init__()
        self.outpath = outpath
        self.best_weights = None  # best weights

    def on_train_begin(self, logs=None):
        self.best = np.Inf  # initialize best as Inf

    def on_epoch_end(self, epoch, logs=None):
        current = logs.get('val_loss')  # get validation loss
        if np.less(current, self.best):
            self.best = current
            self.best_weights = self.model.get_weights()  # record the best weights

    def on_train_end(self, logs=None):
        self.model.save(self.outpath)

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

def save_training_res(model, fit, outpath):
    print('----- saving training results -----')

    # save entire model
    model.save(outpath)

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
    print(outname)

    outpath = '/scratch2/hroetsc/Hotspots/results/{}_prediction.csv'.format(outname)
    if os.path.exists(outpath):
        os.remove(outpath)

    for i in range(hvd.size()):
        cnt_pred = pd.read_csv('/scratch2/hroetsc/Hotspots/results/{}_prediction_rank{}.csv'.format(outname, i))
        if i == 0:
            res = pd.DataFrame({'Accession': cnt_pred['Accession'],
                                'window': cnt_pred['window'],
                                'count': cnt_pred['count'],
                                'pred_count': cnt_pred['pred_count'] * (1 / hvd.size())})
        else:
            res['pred_count'] += cnt_pred['pred_count'] * (1 / hvd.size())

    pd.DataFrame.to_csv(res, outpath, index=False)


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
