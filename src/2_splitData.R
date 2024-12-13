### HEADER ###
# HOTSPOT PREDICTION
# description: split features into training and testing data set
# input: features from 1_featureGeneration.R, protein accessions
# output: training and testing data set
# author: HR

library(plyr)
library(dplyr)
library(stringr)
library(data.table)

testSize = .2

### INPUT ###
windowTokens = fread("data/windowTokens25aa.csv") %>%
  as.data.frame()
accU = windowTokens$Accession %>% unique()

load("data/IEDB+Sarkizova/windowCounts.RData")


### MAIN PART ###
#### clean data: remove all proteins with duplicated windows
repeat {
  isoforms = windowTokens$Accession[which(duplicated(windowTokens$window))] %>%
    unique()
  
  if (length(isoforms) == 0) {
    break
    
  } else {
   windowTokens = windowTokens[-which(windowTokens$Accession %in% isoforms), ] 
  }
}



### characterise count distributions
counts = log2(windowTokens$counts+1)
summary(counts)

zero = which(counts == 0)
length(zero) / length(counts)

plot(density(counts[-zero]))
# lines(density(2^counts[-zero] - 1), type = 'l', col = 'red')
summary(counts[-zero])


### split data sets
accU_noisoforms = windowTokens$Accession %>% unique()
n_test = ceiling(length(accU_noisoforms) * testSize)

testing.acc = accU_noisoforms[sample(length(accU_noisoforms), n_test)]
training.acc = accU_noisoforms[-which(accU_noisoforms %in% testing.acc)]

testing = windowTokens[windowTokens$Accession %in% testing.acc, ]
training = windowTokens[windowTokens$Accession %in% training.acc, ]

# check!
intersect(training$window, testing$window) %>% length()

{
# # sample training data based on average count per protein
# trainCounts = windowCounts[names(windowCounts) %in% training$Accession]
# 
# average_counts = rep(NA, length(trainCounts))
# pb = txtProgressBar(min = 0, max = length(trainCounts), style = 3)
# for (t in 1:length(trainCounts)) {
#   setTxtProgressBar(pb, t)
#   average_counts[t] = trainCounts[[t]] %>% mean()
# }
# average = cbind(names(trainCounts), average_counts) %>% as.data.frame()
# 
# # keep all proteins that have an average count > 0.6
# high_counts = which(average_counts > 0.5)
# # from the ones lower than 0.5, sample as many as are in the > 0.6 data set
# low_counts = which(average_counts <= 0.5)
# low_counts = low_counts[sample(length(low_counts), length(high_counts))]
# 
# keep = which(training$Accession %in% unique(average$V1[c(high_counts, low_counts)]))
# 
# average$average_counts %>% as.character() %>% as.numeric() %>% density() %>% plot()
# average$average_counts[c(low_counts, high_counts)] %>% as.character() %>% as.numeric() %>% density() %>% lines(col = "red")
# 
# training = training[keep, ]
}

# save IDs of proteins in train/test data
train_IDs = training$Accession %>% unique()
test_IDs = testing$Accession %>% unique()


### OUTPUT ###
write.csv(training,
          "data/windowTokens_training_25aa.csv",
          row.names = F)
write.csv(testing,
          "data/windowTokens_testing_25aa.csv",
          row.names = F)

write.csv(train_IDs, "data/IDs_train.csv", row.names = F)
write.csv(test_IDs, "data/IDs_test.csv", row.names = F)




### get same proteins for different window sizes ###

train_IDs = read.csv("data/IDs_train.csv", stringsAsFactors = F) %>% c() %>% unlist()
test_IDs = read.csv("data/IDs_test.csv", stringsAsFactors = F) %>% c() %>% unlist()


# first count vs mean count
windowTokens_firstAA = fread("data/windowTokens25aa_firstAA.csv") %>%
  as.data.frame()
train_firstAA = windowTokens_firstAA[windowTokens_firstAA$Accession %in% train_IDs, ]
test_firstAA = windowTokens_firstAA[windowTokens_firstAA$Accession %in% test_IDs, ]

train = fread("data/windowTokens_training_25aa.csv") %>%
  as.data.frame()
all(train$window == train_firstAA$window)

test = fread("data/windowTokens_testing_25aa.csv") %>%
  as.data.frame()
all(test$window == test_firstAA$window)

write.csv(train_firstAA, "data/windowTokens_training_25aa_firstAA.csv", row.names = F)
write.csv(test_firstAA, "data/windowTokens_testing_25aa_firstAA.csv", row.names = F)

summary(windowTokens_firstAA$counts)
plot(density(windowTokens_firstAA$counts))



# first count vs mean count
windowTokens_mean = fread("data/windowTokens25aa_mean.csv") %>%
  as.data.frame()
train_mean = windowTokens_mean[windowTokens_mean$Accession %in% train_IDs, ]
test_mean = windowTokens_mean[windowTokens_mean$Accession %in% test_IDs, ]

train = fread("data/windowTokens_training_25aa.csv") %>%
  as.data.frame()
all(train$window == train_mean$window)

test = fread("data/windowTokens_testing_25aa.csv") %>%
  as.data.frame()
all(test$window == test_mean$window)

write.csv(train_mean, "data/windowTokens_training_25aa_mean.csv", row.names = F)
write.csv(test_mean, "data/windowTokens_testing_25aa_mean.csv", row.names = F)

summary(windowTokens_mean$counts)
plot(density(windowTokens_mean$counts))




# different window sizes
windows_10aa = fread("data/windowTokens10aa.csv") %>%
  as.data.frame()

train_10aa = windows_10aa[windows_10aa$Accession %in% train_IDs, ]
test_10aa = windows_10aa[windows_10aa$Accession %in% test_IDs, ]

write.csv(train_10aa, "data/windowTokens_training_10aa.csv", row.names = F)
write.csv(test_10aa, "data/windowTokens_testing_10aa.csv", row.names = F)



windows_50aa = fread("data/windowTokens50aa.csv") %>%
  as.data.frame()

train_50aa = windows_50aa[windows_50aa$Accession %in% train_IDs, ]
test_50aa = windows_50aa[windows_50aa$Accession %in% test_IDs, ]

write.csv(train_50aa, "data/windowTokens_training_50aa.csv", row.names = F)
write.csv(test_50aa, "data/windowTokens_testing_50aa.csv", row.names = F)



windows_75aa = fread("data/windowTokens75aa.csv") %>%
  as.data.frame()

train_75aa = windows_75aa[windows_75aa$Accession %in% train_IDs, ]
test_75aa = windows_75aa[windows_75aa$Accession %in% test_IDs, ]

write.csv(train_75aa, "data/windowTokens_training_75aa.csv", row.names = F)
write.csv(test_75aa, "data/windowTokens_testing_75aa.csv", row.names = F)

