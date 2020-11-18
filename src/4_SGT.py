### HEADER ###
# HOTSPOT PREDICTION
# description: generate whole-sequence embeddings using sequence graph transform algorithm
# input: human proteome, train and test protein IDs
# output: 400-dimensional embedding for each protein in train and test data
# author: HR

import os
import pandas as pd
from sgt import SGT
import multiprocessing as mp
import windowMatrices_helper

workers = 12
spec_short = ''
spec_long = '_25aa'


### INPUT ###
human_proteome = pd.read_csv('data/proteins_w_hotspots.csv')
train_IDs = pd.Series(pd.read_csv('data/IDs_train'+spec_short+'.csv')['x'])
test_IDs = pd.Series(pd.read_csv('data/IDs_test'+spec_short+'.csv')['x'])

windows_train = pd.read_csv(str('data/windowTokens_training'+spec_long+'.csv'))
windows_test = pd.read_csv(str('data/windowTokens_testing'+spec_long+'.csv'))


### MAIN PART ###
# set up sequence graph transform
sgt = SGT(kappa=5,
          lengthsensitive=False,
          mode='multiprocessing',
          processors=workers)

def sgt_transform(df):
    df['id'] = df['Accession']
    df['sequence'] = df['seqs'].map(list)

    sgt_embedding = sgt.fit_transform(df)

    return sgt_embedding

# get proteins in train and test data
train = human_proteome[human_proteome['Accession'].isin(train_IDs)]
test = human_proteome[human_proteome['Accession'].isin(test_IDs)]

assert len(train) == len(train_IDs), "train sets are incomplete!"
assert len(test) == len(test_IDs), "test sets are incomplete!"

# apply sgt or open them if already generated
print('obtain the SGT embeddings of the train and test data set')

if not os.path.exists(str('data/sgt_train'+spec_short+'.csv')):
    train_sgt = sgt_transform(df=train)
    test_sgt = sgt_transform(df=test)

    pd.DataFrame.to_csv(train_sgt,
                        str('data/sgt_train' + spec_short + '.csv'),
                        index=False)

    pd.DataFrame.to_csv(test_sgt,
                        str('data/sgt_test' + spec_short + '.csv'),
                        index=False)
else:
    train_sgt = pd.read_csv(str('data/sgt_train' + spec_short + '.csv'))
    test_sgt = pd.read_csv(str('data/sgt_test' + spec_short + '.csv'))


# whole sequence embedding
print('expand the SGT embeddings for all windows')

train_sgt_path = str('data/WHOLE-SEQ_train'+spec_long+'.dat')
test_sgt_path = str('data/WHOLE-SEQ_test'+spec_long+'.dat')

train_sgt_ids = str('data/WHOLE-SEQ_ID_train'+spec_long+'.dat')
test_sgt_ids = str('data/WHOLE-SEQ_ID_test'+spec_long+'.dat')


pool = mp.Pool(workers)
if __name__ == "__main__":
    pool.starmap( windowMatrices_helper.sgt_extension,
                  [[windows_train['Accession'][i], train_sgt, i, train_sgt_path, train_sgt_ids] for i in range(windows_train.shape[0])] )

if __name__ == "__main__":
    pool.starmap( windowMatrices_helper.sgt_extension,
                  [[windows_test['Accession'][i], test_sgt, i, test_sgt_path, test_sgt_ids] for i in range(windows_test.shape[0])] )

pool.close()


