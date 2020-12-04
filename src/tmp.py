
import numpy as np
import pandas as pd

group = 'test'
spec = '_25aa_100-sample'
encoding = 'oneHOT'

tokensAndCounts = pd.read_csv('validation/data/windows_25aa_corona.csv')
pseudocounts = 1
relative_dist = False
protein_norm = False
log_counts = True

emb_path = str('data/EMBEDDING_' + encoding + '_'+group+spec+'.dat')
acc_path = str('data/ACCESSIONS_' + encoding + '_'+group+spec+'.dat')
whole_seq_path = str('data/WHOLE-SEQ_'+group+spec+'.dat')
whole_seq_ids = str('data/WHOLE-SEQ_ID_' + group + spec + '.dat')

windowSize = 25
embeddingDim = 20
sgt_dim = 400

originar_arr = np.array(['A', 'B', 'C', 'D', "E", "F"])
transf_arr = np.array(['B', "D", "F", "A", "C", "E"])
transf_IDs = np.array([1, 3, 5, 0, 2, 4])

no_elements = int(windowSize * embeddingDim)

with open('validation/data/EMBEDDING_AAindex_25aa_corona.dat', 'rb') as emb:
    embMatrix = np.fromfile(emb, dtype='float32')
    emb.seek(0)
    embMatrix_1 = np.fromfile(emb, dtype='float32', count=no_elements).reshape(windowSize, embeddingDim)
    emb.seek(4*no_elements)
    embMatrix_2 = np.fromfile(emb, dtype='float32', count=no_elements).reshape(windowSize, embeddingDim)
    emb.close()


no_windows = int(len(embMatrix) / no_elements)

embMatrix_reshape = embMatrix.reshape(no_windows, windowSize, embeddingDim)
embMatrix_reshape = np.expand_dims(embMatrix_reshape, axis=0)



