### HEADER ###
# HOTSPOT PREDICTION
# description: generate feature set using AA-Index
# input: list of all AA indices
# output: propensity scales of amino acids (scaled between 0 and 1)
# author: HR

library(plyr)
library(dplyr)
library(stringr)
library(caret)
library(bio3d)
library(seqinr)

embeddingDim = 20

### INPUT ###
data("aaindex")
windows =  read.csv("/media/hanna/Hanna2/DATA/Hotspots/DATA/windowTokens_training05.csv",
                    stringsAsFactors = F, header = T)


### MAIN PART ###
# how many unique amino acids?
AAs = windows$window %>% paste(collapse = "")
AAs = AAs %>% strsplit(coll("")) %>% unlist() %>% as.character() %>% unique() %>% sort()
print(AAs)


# load and transform amino acid indices into data frame
aaindices = matrix(ncol = 20, nrow = length(aaindex))
rownames(aaindices) = names(aaindex)
colnames(aaindices) = aaindex[[1]][["I"]] %>% names() %>% a()

for (k in 1:length(aaindex)) {
  cnt = aaindex[[k]][["I"]]
  aaindices[k, ] = cnt
}

aaindices = aaindices[, order(colnames(aaindices))] %>%
  as.data.frame() %>%
  na.omit() %>%
  t()


# remove highly correlated features
correlation = cor(aaindices)
high_correlation = findCorrelation(correlation, cutoff = .9)
aaindices = aaindices[, -high_correlation]


# PCA --> select principal components as aa encoding
pca = prcomp(aaindices %>% t())
summary(pca)

encoding = pca$rotation
encoding = cbind(rownames(encoding), encoding) %>% as.data.frame()
names(encoding)[1] = "aa"


### OUTPUT ###
write.csv(encoding, "/media/hanna/Hanna2/DATA/Hotspots/DATA/ENCODING_AAindex.csv", row.names = F)


