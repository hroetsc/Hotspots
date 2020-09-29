### HEADER ###
# HOTSPOT PREDICTION
# description: calculate embedding matrices for every window in train and test data set
# input: table with amino acids in sliding windows, different encodings
# output: big matrix files (train and test data) in binary format
# author: HR

import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import windowMatrices_helper

workers = 16
embeddingDim = 20
window_size = 25

### MAIN PART ###

def get_window(tokensAndCounts):
    tokensAndCounts['idx'] = tokensAndCounts.reset_index().index
    windows = np.array(tokensAndCounts.loc[:, ['Accession', 'window', 'idx']], dtype='object')
    return windows

pool = mp.Pool(workers)

def generate_matrices(enc, spec):
    print(enc)

    tokensAndCounts_train = pd.read_csv(str('data/windowTokens_training'+spec+'.csv'))
    tokensAndCounts_test = pd.read_csv(str('data/windowTokens_testing'+spec+'.csv'))

    encoding = pd.read_csv(str('data/ENCODING_'+enc+'.csv'))

    emb_train = str('data/EMBEDDING_'+enc+'_train'+spec+'.dat')
    emb_test = str('data/EMBEDDING_'+enc+'_test'+spec+'.dat')
    acc_train = str('data/ACCESSIONS_'+enc+'_train'+spec+'.dat')
    acc_test = str('data/ACCESSIONS_'+enc+'_test'+spec+'.dat')

    if os.path.exists(str('data/*_'+enc+'_t*'+spec+'.dat')):
        os.remove(str('data/*_'+enc+'_t*'+spec+'.dat'))

    windows_train = get_window(tokensAndCounts_train)
    windows_test = get_window(tokensAndCounts_test)

    if enc == "SGT":
        sgt = True
    else:
        sgt = False

    if __name__ == "__main__":
        print('training')
        pool.starmap( windowMatrices_helper.find_encodings,
                     [[window, encoding, emb_train, acc_train, embeddingDim, window_size, sgt] for window in list(windows_train)] )

    if __name__ == "__main__":
        print('testing')
        pool.starmap( windowMatrices_helper.find_encodings,
                     [[window, encoding, emb_test, acc_test, embeddingDim, window_size, sgt] for window in list(windows_test)] )


generate_matrices(enc='AAindex', spec='50-2')
generate_matrices(enc='oneHOT', spec='50-2')
# generate_matrices(enc='SGT', spec='05')

pool.close()