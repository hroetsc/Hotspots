### HEADER ###
# HOTSPOT PREDICTION
# description: generate feature set using AA-Index
# input: list of all AA indices
# output: propensity scales of amino acids (scaled between 0 and 1)
# author: HR

library(plyr)
library(dplyr)
library(stringr)

library(bio3d)
library(seqinr)

embeddingDim = 20

### INPUT ###
data("aaindex")
windows = read.csv("data/windowTokens_training20.csv", stringsAsFactors = F, header = T)


### MAIN PART ###
# how many unique amino acids?
AAs = windows$window %>% paste(collapse = "")
AAs = AAs %>% strsplit(coll("")) %>% unlist() %>% as.character() %>% unique() %>% sort()
print(AAs)

# select 20 indices
idx = c("KYTJ820101",  # hydropathy (Kyte-Doolitle)
        "FASG760104",  # pkN
        "FASG760105",  # pKC
        "FAUJ880109",  # no. of hydrogen bond donors
        "KLEP840101",  # net charge
        "DAYM780101",  # amino acid composition
        "CHOC750101",  # average volume of buried residue
        "CHOC760101",  # residue accessible surface area in tripeptide
        "CHAM820101",  # polarizability
        "ARGP820101",  # hydrophobicity index
        "PONP800102",  # average gain in surrounding hydrophobicity
        "RADA880108",  # mean polarity
        "WOLR810101",  # hydration potential
        "ZIMJ680102",  # bulkiness
        "ZIMJ680104",  # isoelectric point
        "TAKK010101",  # side-chain contribution to protein stability
        "EISD860103",  # direction of hydrophobic moment
        "FAUJ880108",  # localized electrical effect
        "HUTJ700101",  # heat capacity
        "HUTJ700102")  # absolute entropy

aas = matrix(ncol = embeddingDim + 1, nrow = embeddingDim) %>% as.data.frame()
names(aas) = c("aa", idx)
aas$aa = AAs

for (i in 2:ncol(aas)) {
  
  cnt.idx = aaindex[[names(aas)[i]]][["I"]]
  names(cnt.idx) = names(cnt.idx) %>% a()
  
  cnt.idx = cnt.idx[order(names(cnt.idx))]
  
  aas[, i] = (cnt.idx - min(cnt.idx)) / (max(cnt.idx) - min(cnt.idx))  # scale
}

# plot
for (j in 1:nrow(aas)) {
  plot(density(aas[j, c(2:ncol(aas))] %>% as.character() %>% as.numeric()))
}


### OUTPUT ###
write.csv(aas, "data/ENCODING_AAindex.csv", row.names = F)


