### HEADER ###
# HOTSPOT PREDICTION: VALIDATION
# description: provide helper function for feature generation and data handling
# input: -
# output: -
# author: HR

import pandas as pd
import numpy as np


#######################################################################################################################
# FEATURE ENCODING
#######################################################################################################################

def window_encoding(window, encoding, emb_file, acc_file, embeddingDim, window_size):
    aas = [str for str in [window[1]]][0]
    out = [None] * window_size

    for n, aa in enumerate(aas):
        cnt = encoding[encoding['aa'] == aa].iloc[:, 1:(embeddingDim + 1)]
        out[n] = np.array(cnt, dtype='float32').flatten()

    out = np.asarray(out, dtype='float32')
    idx = np.asarray(window[2], dtype='int32')

    # save accessions and weight matrix as numpy array in binary format
    with open(emb_file, 'ab') as ef, open(acc_file, 'ab') as af:
        out.tofile(ef)
        idx.tofile(af)

        ef.close()
        af.close()


def sgt_extension(window_acc, sgt_emb, ID, outpath, idxpath):
    cnt = sgt_emb[sgt_emb['id'] == window_acc].iloc[:, 1:]
    cnt = np.array(cnt, dtype='float32')

    ID = np.asarray(ID, dtype='int32')

    with open(outpath, 'ab') as of, open(idxpath, 'ab') as idf:
        cnt.tofile(of)
        ID.tofile(idf)

        of.close()
        idf.close()


#######################################################################################################################
# DATA LOADING
#######################################################################################################################

def format_input(tokensAndCounts, relative_dist=False):
    tokens = np.array(tokensAndCounts[['organism', 'accession', 'annotation', 'window']])

    dist_N = np.array(tokensAndCounts['dist_N'], dtype='int32')
    dist_C = np.array(tokensAndCounts['dist_C'], dtype='int32')

    if relative_dist:
        dist_N = dist_N / np.array(tokensAndCounts['len_protein'], dtype='int32')
        dist_C = dist_C / np.array(tokensAndCounts['len_protein'], dtype='int32')

    dist = np.array([dist_N, dist_C]).transpose()

    print('number of features: ', tokens.shape[0])
    return tokens, dist


def open_and_format_matrices(encoding, spec,
                             windowSize, embeddingDim, sgt_dim,
                             relative_dist):
    print(encoding)
    print(spec)

    # accession, windows and counts
    tokensAndCounts = pd.read_csv(str('/scratch2/hroetsc/Hotspots/data/windows' + spec + '.csv'))
    tokens, dist = format_input(tokensAndCounts,
                                relative_dist=relative_dist)

    # get embeddings
    emb_path = str('/scratch2/hroetsc/Hotspots/data/EMBEDDING_' + encoding + spec + '.dat')
    acc_path = str('/scratch2/hroetsc/Hotspots/data/ACCESSIONS_' + encoding + spec + '.dat')
    whole_seq_path = str('/scratch2/hroetsc/Hotspots/data/WHOLE-SEQ' + spec + '.dat')
    whole_seq_ids = str('/scratch2/hroetsc/Hotspots/data/WHOLE-SEQ_ID' + spec + '.dat')

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

    return tokens, embMatrix, dist, wholeSeq
