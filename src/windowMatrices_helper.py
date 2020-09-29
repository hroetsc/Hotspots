### HEADER ###
# HOTSPOT PREDICTION
# description: helper script for parallel calculation of embedding matrices
# input: -
# output: -
# author: HR

import numpy as np

def find_encodings(window, encoding, emb_file, acc_file, embeddingDim, window_size, sgt):

    aas = [str for str in [window[1]]][0]

    if sgt:
        ID = window[0][0]
        encoding = encoding[encoding['Accession'] == ID]

    out = [None] * window_size

    for n, aa in enumerate(aas):
        if sgt:
            cnt = encoding[encoding['aa'] == aa].iloc[:, :embeddingDim]
        else:
            cnt = encoding[encoding['aa'] == aa].iloc[:, 1:(embeddingDim+1)]

        out[n] = np.array(cnt, dtype='float32').flatten()

    out = np.asarray(out, dtype='float32')
    idx = np.asarray(window[2], dtype='int32')

    # save accessions and weight matrix as numpy array in binary format
    with open(emb_file, 'ab') as ef, open(acc_file, 'ab') as af:
        out.tofile(ef)
        idx.tofile(af)

        ef.close()
        af.close()