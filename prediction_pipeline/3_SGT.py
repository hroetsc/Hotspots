### HEADER ###
# HOTSPOT PREDICTION: VALIDATION
# description: generate whole-sequence embeddings using sequence graph transform algorithm
# input: aggregated viral protein sequences
# output: 400-dimensional embedding for each protein
# author: HR

import os
import pandas as pd
from sgt import SGT
import multiprocessing as mp
import helper

workers = 16
spec = '_50aa_corona'


### INPUT ###
prots = pd.read_csv('validation/data/aggregated_sequences.csv')
windows = pd.read_csv('validation/data/windows'+spec+'.csv')


### MAIN PART ###
print('obtain the SGT embeddings')
# set up sequence graph transform
sgt = SGT(kappa=5,
          lengthsensitive=False,
          mode='multiprocessing',
          processors=workers)

prots['id'] = prots['accession']
prots['sequence'] = prots['sequence'].map(list)

# needs to be run only once (irregardless of window size)
if not os.path.exists('validation/data/sgt.csv'):
    sgt_embedding = sgt.fit_transform(prots)
    pd.DataFrame.to_csv(sgt_embedding,
                        'validation/data/sgt.csv',
                        index=False)
else:
    sgt_embedding = pd.read_csv('validation/data/sgt.csv')


# whole sequence embedding
print('expand the SGT embeddings for all windows')

sgt_path = str('validation/data/WHOLE-SEQ'+spec+'.dat')
sgt_ids = str('validation/data/WHOLE-SEQ_ID'+spec+'.dat')

pool = mp.Pool(workers)
if __name__ == "__main__":
    pool.starmap(helper.sgt_extension,
                 [[windows['accession'][i], sgt_embedding, i, sgt_path, sgt_ids] for i in range(windows.shape[0])])
pool.close()



