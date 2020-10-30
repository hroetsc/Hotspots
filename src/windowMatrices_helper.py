### HEADER ###
# HOTSPOT PREDICTION
# description: helper script for parallel calculation of embedding matrices
# input: -
# output: -
# author: HR

import numpy as np

def find_encodings(window, encoding, emb_file, acc_file, embeddingDim, window_size):

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