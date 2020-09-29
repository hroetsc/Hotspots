### HEADER ###
# HOTSPOT PREDICTION
# description: generate feature set one-hot encoding
# input: windows
# output: one-hot encoding for amino acids
# author: HR

library(dplyr)
library(stringr)


### INPUT ###
windows = read.csv("data/windowTokens_training20.csv", stringsAsFactors = F, header = T)


### MAIN PART ###
# how many unique amino acids?
AAs = windows$window %>% paste(collapse = "")
AAs = AAs %>% strsplit(coll("")) %>% unlist() %>% as.character() %>% unique() %>% sort()
print(AAs)

AAS = matrix(ncol = length(AAs), nrow = length(AAs))
diag(AAS) = 1
AAS[is.na(AAS)] = 0

aas = cbind(AAs, AAS %>% as.data.frame())
names(aas)[1] = "aa"


### OUTPUT ###
write.csv(aas, "data/ENCODING_oneHOT.csv", row.names = F)

