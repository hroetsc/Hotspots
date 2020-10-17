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
windowTokens = fread("data/windowTokens.csv") %>%
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


### remove windows with X, replace U by C
X_idx = str_detect(windowTokens$window, "X")
X_idx[X_idx] %>% length()
windowTokens = windowTokens[X_idx == F, ]

U_idx = str_detect(windowTokens$window, "U")
U_idx[U_idx] %>% length()
windowTokens_U = windowTokens[U_idx, ]

for (i in 1:nrow(windowTokens_U)) {
  windowTokens_U$window[i] = str_replace_all(windowTokens_U$window[i], pattern = "U", replacement = "C") 
}

windowTokens[U_idx, ] = windowTokens_U


### sample
testSize = 0.1

# no isoforms of training data in testing data (and vice versa!)
accU_noisoforms = str_split_fixed(accU, coll("-"), Inf)[, 1] %>% unique()

n_train = ceiling(length(accU_noisoforms) * (1-testSize))
n_test = floor(length(accU_noisoforms) * testSize)

training.acc = accU_noisoforms[sample(length(accU_noisoforms), n_train)]
testing.acc = accU_noisoforms[sample(length(accU_noisoforms), n_test)]


#### clean data
# split data sets
train_idx = which(str_split_fixed(windowTokens$Accession, coll("-"), Inf)[, 1] %>%
                    as.character() %in% training.acc)
training = windowTokens[train_idx, ]

test_idx = which(str_split_fixed(windowTokens$Accession, coll("-"), Inf)[, 1] %>%
                    as.character() %in% testing.acc)
testing = windowTokens[test_idx, ]


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

# do not execute bc positions as additional feature
# remove sliding windows that occur in training data set from testing data set
# intersect(training.tmp$window, testing$window) %>% length()
# rm = which(testing$window %in% training.tmp$window)
# if (length(rm) > 0){ testing.tmp = testing[-rm, ] }


# downsample
# training.5per = training[training$Accession %in% sample(training.acc, ceiling(length(training.acc)*0.05)), ]
# testing.5per = testing[testing$Accession %in% sample(testing.acc, ceiling(length(testing.acc)*0.05)), ]
# 
# training.20per = training[training$Accession %in% sample(training.acc, ceiling(length(training.acc)*0.2)), ]
# testing.20per = testing[testing$Accession %in% sample(testing.acc, ceiling(length(testing.acc)*0.2)), ]
# 
# training.50per = training[training$Accession %in% sample(training.acc, ceiling(length(training.acc)*0.5)), ]
# testing.50per = testing[testing$Accession %in% sample(testing.acc, ceiling(length(testing.acc)*0.5)), ]


### OUTPUT ###
write.csv(training,
          "/media/hanna/Hanna2/DATA/Hotspots/DATA/windowTokens_training100-sample.csv",
          row.names = F)
write.csv(testing,
          "/media/hanna/Hanna2/DATA/Hotspots/DATA/windowTokens_testing100-sample.csv",
          row.names = F)

# write.csv(training.5per, "data/windowTokens_training05.csv", row.names = F)
# write.csv(testing.5per, "data/windowTokens_testing05.csv", row.names = F)
# 
# write.csv(training.20per, "data/windowTokens_training20.csv", row.names = F)
# write.csv(testing.20per, "data/windowTokens_testing20.csv", row.names = F)
# 
# write.csv(training.50per, "data/windowTokens_training50.csv", row.names = F)
# write.csv(testing.50per, "data/windowTokens_testing50.csv", row.names = F)


#### get the second half of the data
training = read.csv("data/windowTokens_training.csv", stringsAsFactors = F)
testing = read.csv("data/windowTokens_testing.csv", stringsAsFactors = F)

training.50per = read.csv("data/windowTokens_training50.csv", stringsAsFactors = F)
testing.50per = read.csv("data/windowTokens_testing50.csv", stringsAsFactors = F)

training.50per.2 = training[!training$Accession %in% training.50per$Accession, ]
testing.50per.2 = testing[!testing$Accession %in% testing.50per$Accession, ]

nrow(training.50per) + nrow(training.50per.2) == nrow(training)
nrow(testing.50per) + nrow(testing.50per.2) == nrow(testing)

write.csv(training.50per.2, "data/windowTokens_training50-2.csv", row.names = F)
write.csv(testing.50per.2, "data/windowTokens_testing50-2.csv", row.names = F)



#### for comparison of different extensions:
# get the same training/testing data sets as in df w/o extension
# to avoid re-calculation of all the windows

training.50per = read.csv("data/windowTokens_training50.csv", stringsAsFactors = F)
testing.50per = read.csv("data/windowTokens_testing50.csv", stringsAsFactors = F)

window = read.csv("data/windowTokens_norm.csv", stringsAsFactors = F)

get_same_data = function(df = "", window = ""){
  
  window_sub = window[which(window$Accession %in% df$Accession), ]
  out = window_sub[order(window_sub$window %in% df$window), ]
  out = out[order(out$Accession %in% window_sub$Accession), ]
  
  if (length(which(out$window != df$window)) > 0) {
    print("ELEMENTS DO NOT MATCH !!!")
  }
  
  return(out)
}

training = get_same_data(df = training.50per, window = window)
testing = get_same_data(df = testing.50per, window = window)

write.csv(training, "data/windowTokens_norm_training50.csv", row.names = F)
write.csv(testing, "data/windowTokens_norm_testing50.csv", row.names = F)

