### HEADER ###
# HOTSPOT PREDICTION: NEURAL ARCHITECTURE SEARCH
# description: subset of proteins with similar characteristics
# input: windowTokens, human proteome, window counts
# output: train and test data set
# author: HR

library(data.table)
library(dplyr)

testSize = .2

### INPUT ###
windowTokens = fread("data/windowTokens25aa_mean.csv") %>%
  as.data.frame()
prots = read.csv("data/proteins_w_hotspots.csv",
                 stringsAsFactors = F)

### MAIN PART ##
# remove isoforms
repeat {
  isoforms = windowTokens$Accession[which(duplicated(windowTokens$window))] %>%
    unique()
  
  if (length(isoforms) == 0) {
    break
    
  } else {
    windowTokens = windowTokens[-which(windowTokens$Accession %in% isoforms), ] 
  }
}

prots = prots[prots$Accession %in% windowTokens$Accession, ]

# mean count
# protein length
# number of peaks

characteristics = data.frame(accession = prots$Accession,
                             min_count = NA,
                             max_count = NA,
                             mean_count = NA,
                             protein_length = nchar(prots$seqs),
                             no_peaks = NA)

pb = txtProgressBar(min = 0, max = nrow(prots), style = 3)
for (i in 1:nrow(prots)) {
  setTxtProgressBar(pb, i)
  cnt = windowTokens[windowTokens$Accession == prots$Accession[i], ]
  
  characteristics$min_count[i] = cnt$counts %>% min()
  characteristics$max_count[i] = cnt$counts %>% max()
  characteristics$mean_count[i] = cnt$counts %>% mean()
  
  slope = rep(NA, nrow(cnt))
  slope[1] = 0
  for (j in 2:nrow(cnt)) {
    slope[j] = cnt$counts[j] - cnt$counts[j-1]
  }
  roundSlope = -(round(log10(abs(mean(slope))))-2)
  slope = round(slope,digits=roundSlope)
  
  turnpoints = rep(FALSE,(nrow(cnt)-1))
  for(j in 1:(length(slope)-1)){
    if(slope[j]>0 & slope[j+1]<0){ turnpoints[j] = TRUE  }
  }
  
  characteristics$no_peaks[i] = which(turnpoints) %>% length()
  
}

# select ca. 2k proteins
summary(characteristics$mean_count)
summary(characteristics$protein_length)
summary(characteristics$no_peaks)

k = which(characteristics$mean_count > 1 & characteristics$protein_length <= 450 & characteristics$no_peaks > 1)

# split in train and test
k_test = k[sample(length(k), ceiling(testSize*length(k)))]
k_train = k[-which(k %in% k_test)]

windowTokens_train = windowTokens[windowTokens$Accession %in% characteristics$accession[k_train], ]
windowTokens_test = windowTokens[windowTokens$Accession %in% characteristics$accession[k_test], ]


### OUTPUT ###
write.csv(characteristics, "NAS/data/sequence_characteristics.csv", row.names = F)

write.csv(windowTokens_train, "NAS/data/windowTokens_training_25aa_mean.csv", row.names = F)
write.csv(windowTokens_test, "NAS/data/windowTokens_testing_25aa_mean.csv", row.names = F)


