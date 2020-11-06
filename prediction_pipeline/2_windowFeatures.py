### HEADER ###
# HOTSPOT PREDICTION: VALIDATION
# description: calculate embedding matrices for every window in validation data set
# input: table with amino acids in sliding windows, different encodings
# output: big matrix file in binary format
# author: HR

import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import helper

workers = 16
embeddingDim = 20
window_size = 50

spec = '_50aa_corona'
print(spec)

### MAIN PART ###
pool = mp.Pool(workers)

def generate_matrices(enc, spec):
    print(enc)

    windows = pd.read_csv('validation/data/windows'+spec+'.csv')
    windows['idx'] = windows.reset_index().index
    windows = np.array(windows[['accession', 'window', 'idx']], dtype='object')

    encoding = pd.read_csv('data/ENCODING_'+enc+'.csv')

    emb = str('validation/data/EMBEDDING_' + enc + spec + '.dat')
    acc = str('validation/data/ACCESSIONS_' + enc + spec + '.dat')

    if os.path.exists(emb):
        os.remove(emb)
        os.remove(acc)

    if __name__ == "__main__":
        pool.starmap(helper.window_encoding,
                     [[window, encoding, emb, acc, embeddingDim, window_size] for window in list(windows)])


# window matrices
generate_matrices(enc='AAindex', spec=spec)
generate_matrices(enc='oneHOT', spec=spec)

pool.close()