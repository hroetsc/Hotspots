### HEADER ###
# HOTSPOT PREDICTION
# description: get single-aa embedding using sequence graph transform (SGT) algorithm
# input: proteome
# output: SGT embeddings for all amino acids (protein-specific)
# author: HR

import pandas as pd
import numpy as np
from sgt import SGT

embeddingDim = 20

### INPUT ###
proteome = pd.read_csv('data/proteins_w_hotspots.csv',
                       names=['id', 'sequence'], header=1)

aas = np.array(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])


### MAIN PART ###
proteome['sequence'] = proteome['sequence'].map(list)

# get SGT embedding
sgt = SGT(kappa=5,
          lengthsensitive=False)
sgt_embedding = sgt.fit_transform(corpus=proteome)

# split into aa
for r in range(sgt_embedding.shape[0]):
    if r == 0:
        enc = pd.DataFrame(
            np.array(sgt_embedding.iloc[r, 1:], dtype='float32').reshape((embeddingDim, embeddingDim)))
        enc['Accession'] = np.repeat(sgt_embedding.iloc[r, 0], embeddingDim)
        enc['aa'] = aas

    else:
        cnt_enc = pd.DataFrame(
            np.array(sgt_embedding.iloc[r, 1:], dtype='float32').reshape((embeddingDim, embeddingDim)))
        cnt_enc['Accession'] = np.repeat(sgt_embedding.iloc[r, 0], embeddingDim)
        cnt_enc['aa'] = aas

        enc = pd.concat([enc, cnt_enc])


### OUTPUT ###
pd.DataFrame.to_csv(enc, 'data/ENCODING_SGT.csv', index=False)

