group = 'test'
spec = '_25aa_100-sample'
encoding = 'oneHOT'

tokensAndCounts = pd.read_csv('data/windowTokens_testing_25aa_100-sample.csv')
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