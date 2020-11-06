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

### INPUT ###
windowTokens = fread("data/windowTokens25aa.csv") %>%
  as.data.frame()
accU = windowTokens$Accession %>% unique()

load("data/windowCounts.RData")


### MAIN PART ###
### characterise count distributions
counts = log2(windowTokens$counts+1)
summary(counts)

zero = which(counts == 0)
length(zero) / length(counts)

plot(density(counts[-zero]))
# lines(density(2^counts[-zero] - 1), type = 'l', col = 'red')
summary(counts[-zero])


### sample
testSize = 0.05
n_test = ceiling(length(accU) * testSize)

testing.acc = accU[sample(length(accU), n_test)]
training.acc = accU[-which(accU %in% testing.acc)]


#### clean data
# no isoforms of training data in testing data (and vice versa!)
# make sure that all duplicated windows end up in one group
train.preliminary = windowTokens[windowTokens$Accession %in% training.acc, ]
test.preliminary = windowTokens[windowTokens$Accession %in% testing.acc, ]

pb = txtProgressBar(min = 0, max = length(testing.acc), style = 3)
for (t in 1:length(testing.acc)){
  setTxtProgressBar(pb, t)
  
  cnt = test.preliminary[test.preliminary$Accession == testing.acc[t], ]
  k = which(train.preliminary$window %in% cnt$window)
  
  if (length(k) > 0) {
    cnt.training.acc = unique(train.preliminary$Accession[k])
    testing.acc = c(testing.acc,
                    cnt.training.acc)
    training.acc[training.acc == cnt.training.acc] = NA
  }
}

testing.acc = unique(testing.acc)
training.acc = na.omit(training.acc) %>% unique()

# split data sets
training = windowTokens[windowTokens$Accession %in% training.acc, ]
testing = windowTokens[windowTokens$Accession %in% testing.acc, ]

# check!
intersect(training$window, testing$window) %>% length()

# sample training data based on average count per protein
trainCounts = windowCounts[names(windowCounts) %in% training$Accession]

average_counts = rep(NA, length(trainCounts))
pb = txtProgressBar(min = 0, max = length(trainCounts), style = 3)
for (t in 1:length(trainCounts)) {
  setTxtProgressBar(pb, t)
  average_counts[t] = trainCounts[[t]] %>% mean()
}
average = cbind(names(trainCounts), average_counts) %>% as.data.frame()

# keep all proteins that have an average count > 0.6
high_counts = which(average_counts > 0.6)
# from the ones lower than 0.5, sample as many as are in the > 0.6 data set
low_counts = which(average_counts <= 0.6)
low_counts = low_counts[sample(length(low_counts), length(high_counts))]

keep = which(training$Accession %in% unique(average$V1[c(high_counts, low_counts)]))

average$average_counts %>% as.character() %>% as.numeric() %>% density() %>% plot()
average$average_counts[c(low_counts, high_counts)] %>% as.character() %>% as.numeric() %>% density() %>% lines(col = "red")

training = training[keep, ]


# save IDs of proteins in train/test data
train_IDs = training$Accession %>% unique()
test_IDs = testing$Accession %>% unique()


### OUTPUT ###
write.csv(training,
          "data/windowTokens_training_25aa_100-sample.csv",
          row.names = F)
write.csv(testing,
          "data/windowTokens_testing_25aa_100-sample.csv",
          row.names = F)

write.csv(train_IDs, "data/IDs_train100-sample.csv", row.names = F)
write.csv(test_IDs, "data/IDs_test100-sample.csv", row.names = F)


### get same proteins for different window sizes ###
windows_10aa = fread("data/windowTokens10aa.csv") %>%
  as.data.frame()

train_10aa = windows_10aa[windows_10aa$Accession %in% train_IDs, ]
test_10aa = windows_10aa[windows_10aa$Accession %in% test_IDs, ]

write.csv(train_10aa, "data/windowTokens_training_10aa_100-sample.csv", row.names = F)
write.csv(test_10aa, "data/windowTokens_testing_10aa_100-sample.csv", row.names = F)



windows_50aa = fread("data/windowTokens50aa.csv") %>%
  as.data.frame()

train_50aa = windows_50aa[windows_50aa$Accession %in% train_IDs, ]
test_50aa = windows_50aa[windows_50aa$Accession %in% test_IDs, ]

write.csv(train_50aa, "data/windowTokens_training_50aa_100-sample.csv", row.names = F)
write.csv(test_50aa, "data/windowTokens_testing_50aa_100-sample.csv", row.names = F)

