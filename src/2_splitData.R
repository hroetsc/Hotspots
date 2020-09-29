### HEADER ###
# HOTSPOT PREDICTION
# description: split features into training and testing data set
# input: features from 1_featureGeneration.R, protein accessions
# output: training and testing data set
# author: HR

library(plyr)
library(dplyr)
library(stringr)


### INPUT ###
windowTokens = read.csv("data/windowTokens.csv", stringsAsFactors = F)
accU = windowTokens$Accession %>% unique()

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
testSize = 0.4

n_train = ceiling(length(accU) * (1-testSize))
n_test = floor(length(accU) * testSize)

training.acc = accU[sample(length(accU), n_train)]
testing.acc = accU[sample(length(accU), n_test)]


#### clean data 
# no isoforms of training data in testing data
testing.acc = testing.acc %>%
  as.character()
testing.acc = testing.acc[-which(str_split_fixed(testing.acc, coll("-"), Inf)[,1] %in% training.acc)]

training.acc = training.acc %>%
  as.character()


# split data sets
training = windowTokens[windowTokens$Accession %in% training.acc, ]
testing = windowTokens[windowTokens$Accession %in% testing.acc, ]

# remove sliding windows that occur in training data set from testing data set
intersect(training$tokens, testing$tokens) %>% length()
rm = which(testing$tokens %in% training$tokens)
if (length(rm) > 0){ testing = testing[-rm, ] }


# downsample
training.5per = training[training$Accession %in% sample(training.acc, ceiling(length(training.acc)*0.05)), ]
testing.5per = testing[testing$Accession %in% sample(testing.acc, ceiling(length(testing.acc)*0.05)), ]

training.20per = training[training$Accession %in% sample(training.acc, ceiling(length(training.acc)*0.2)), ]
testing.20per = testing[testing$Accession %in% sample(testing.acc, ceiling(length(testing.acc)*0.2)), ]

training.50per = training[training$Accession %in% sample(training.acc, ceiling(length(training.acc)*0.5)), ]
testing.50per = testing[testing$Accession %in% sample(testing.acc, ceiling(length(testing.acc)*0.5)), ]

### OUTPUT ###
write.csv(training, "data/windowTokens_training.csv", row.names = F)
write.csv(testing, "data/windowTokens_testing.csv", row.names = F)

write.csv(training.5per, "data/windowTokens_training05.csv", row.names = F)
write.csv(testing.5per, "data/windowTokens_testing05.csv", row.names = F)

write.csv(training.20per, "data/windowTokens_training20.csv", row.names = F)
write.csv(testing.20per, "data/windowTokens_testing20.csv", row.names = F)

write.csv(training.50per, "data/windowTokens_training50.csv", row.names = F)
write.csv(testing.50per, "data/windowTokens_testing50.csv", row.names = F)

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

training.5per = read.csv("data/windowTokens_training05.csv", stringsAsFactors = F)
testing.5per = read.csv("data/windowTokens_testing05.csv", stringsAsFactors = F)

window = read.csv("data/NandCext_windowTokens.csv", stringsAsFactors = F)

get_same_data = function(df = "", window = ""){
  
  window_sub = window[which(window$Accession %in% df$Accession), ]
  out = window_sub[order(window_sub$window %in% df$window), ]
  
  if (length(which(out$window != df$window)) > 0) {
    print("ELEMENTS DO NOT MATCH !!!")
  }
  
  return(out)
}

training = get_same_data(df = training.5per, window = window)
testing = get_same_data(df = testing.5per, window = window)

write.csv(training, "data/NandCext_windowTokens_training05.csv", row.names = F)
write.csv(testing, "data/NandCext_windowTokens_testing05.csv", row.names = F)

