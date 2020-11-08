### HEADER ###
# HOTSPOT PREDICTION
# description: find reason for different prediction axes
#               - do windows occur in training data set?
#               - is one of the groups oversampled?
#               - is the protein expression different?
# input: predictions (from 7_performance.R), protein-wise performance (7_analyseProteins.R),
#         windows, mhc ligand sources, protein expressions in human tissues, human proteome
# output: -
# author: HR


library(dplyr)
library(stringr)
library(data.table)


### INPUT ###
train = fread("data/windowTokens_training_50aa_100-sample.csv") %>%
  as.data.frame()
test = fread("data/windowTokens_testing_50aa_100-sample.csv") %>%
  as.data.frame()

overview = read.csv("results/overview.csv", stringsAsFactors = F)

source_antigens = read.csv("HOTSPOTS/DATA/mhc_ligand_table_export_1588248580.csv",
                           stringsAsFactors = F, header = 2)
hp = read.csv("/media/hanna/Hanna2/DATA/GENERAL/proteome_human.csv", stringsAsFactors = F)

protExpr = fread("proteinatlas.tsv") %>% as.data.frame()

### MAIN PART ###
########## identify well and poorly predicted proteins ##########
if(! dir.exists("results/exploratory")) { dir.create("results/exploratory") }

# average count > 0.6 (below that not important)
well.predicted = overview$Accession[overview$PCC > .5] %>% as.character()
poorly.predicted = overview$Accession[overview$PCC < .2] %>% as.character()


########## check for overlap of windows ##########
windows_test.well = test$window[test$Accession %in% well.predicted]
windows_test.poorly = test$window[test$Accession %in% poorly.predicted]

overlap_well = which(windows_test.well %in% train$window) %>% length
overlap_poorly = which(windows_test.poorly %in% train$window) %>% length()

overlap_well / length(windows_test.well)
overlap_poorly / length(windows_test.poorly)

# just to be sure ... check for isoforms
train.acc = str_split_fixed(train$Accession, coll("-"), Inf)[, 1] %>% unique()
test.acc = str_split_fixed(test$Accession, coll("-"), Inf)[, 1] %>% unique()
intersect(train.acc, test.acc) %>% length()  # phew


# fetch some examples
tmp = train$Accession[which(windows_test.well %in% train$window)] %>% unique()
tmp = train$Accession[which(windows_test.well %in% train$window)] %>% unique()

for (i in 1:length(windows_test.well)) {
  k = which(windows_test.well[i] == train$window[train$window %in% windows_test.well])
  if (length(k) > 0) {
    print("..........")
    print(test$Accession[test$window == windows_test.well[i]] %>% unique())
    print(train$Accession[train$window == windows_test.well[i]] %>% unique())
    print("..........")
  }
}

# check how similar counts are
counts.test = c()
counts.train = c()
for (i in 1:20) {
  k = which(windows_test.well[i] == train$window[train$window %in% windows_test.well])
  if (length(k) > 0) {
    counts.test = c(counts.test,
                    test$counts[test$window == windows_test.well[i]])
    counts.train = c(counts.train,
                     train$counts[train$window == windows_test.well[i]])
  }
}


########## check for overstudying ##########
sampling_frequency = data.frame(Accession = str_split_fixed(overview$Accession, coll("-"), Inf)[, 1],
                                no_epitopes = rep(NA, nrow(overview)),
                                no_studies = rep(NA, nrow(overview)),
                                protein_length = rep(NA, nrow(overview)))

pb = txtProgressBar(min = 0, max = nrow(overview), style = 3)
for (i in 1:nrow(overview)) {
  setTxtProgressBar(pb, i)
  cnt = source_antigens[source_antigens$Epitope.8 == overview$Accession[i], ]
  
  if(nrow(cnt) > 0) {
    cnt = cnt[-which(duplicated(cnt$Epitope.2)), ]
    sampling_frequency$no_epitopes[i] = nrow(cnt)
    sampling_frequency$no_studies[i] = cnt$Reference.4 %>% unique() %>% length()
    sampling_frequency$protein_length[i] = hp$seqs[hp$Accession == overview$Accession[i]] %>% nchar()
  }
}

sampling_frequency = na.omit(sampling_frequency)
nrow(sampling_frequency)
# --> might not be reliable

# normalise epitopes by protein length
sampling_frequency$epitope_freq = sampling_frequency$no_epitopes / sampling_frequency$protein_length
cor(sampling_frequency$no_studies, sampling_frequency$epitope_freq, method = "pearson")


# well and poorly predicted proteins
samplingFreg_well = which(sampling_frequency$Accession %in% well.predicted)
samplingFreg_poorly = which(sampling_frequency$Accession %in% poorly.predicted)

lm.well = lm(epitope_freq[samplingFreg_well] ~ no_studies[samplingFreg_well],
             data = sampling_frequency[samplingFreg_well, ])
lm.poorly = lm(epitope_freq[samplingFreg_poorly] ~ no_studies[samplingFreg_poorly],
               data = sampling_frequency[samplingFreg_poorly, ])


# plot
colour = rep(NA, nrow(sampling_frequency))
colour[samplingFreg_well] = "darkgreen"
colour[samplingFreg_poorly] = "firebrick"
colour[is.na(colour)] = "black"

plot(epitope_freq ~ no_studies, data = sampling_frequency,
     pch = 16,
     col = colour)
abline(a = lm.well$coefficients[1],
       b = lm.well$coefficients[2],
       col = "darkgreen")
abline(a = lm.poorly$coefficients[1],
       b = lm.poorly$coefficients[2],
       col = "firebrick")
legend("topright",
       legend = c("well predicted", "poorly predicted"),
       col = c("darkgreen", "firebrick"),
       pch = c(16, 16),
       cex = .8)


########## check for protein expression ##########

Expr_UniProtID = protExpr$Uniprot
Expr_tbl = protExpr[, str_detect(names(protExpr), "Tissue RNA -")]
Rank_tbl = Expr_tbl

for (c in 1:ncol(Expr_tbl)) {
  Rank_tbl[, c] = order(Expr_tbl[, c], decreasing = F)
}

protExpr_mean = data.frame(accession = Expr_UniProtID,
                           mean_rank = rowMeans(Rank_tbl))

protExpr_mean.well = protExpr_mean$mean_rank[protExpr_mean$accession %in% 
                                               str_split_fixed(well.predicted, coll("-"), Inf)[, 1]]
protExpr_mean.poorly = protExpr_mean$mean_rank[protExpr_mean$accession %in% 
                                               str_split_fixed(poorly.predicted, coll("-"), Inf)[, 1]]


plot(density(protExpr_mean.well),
     col = "darkgreen",
     main = "RNA expression of proteins")
lines(density(protExpr_mean.poorly),
     col = "firebrick")
legend("topright",
       legend = c("well predicted", "poorly predicted"),
       col = c("darkgreen", "firebrick"),
       lty = c(1, 1),
       cex = .6)

summary(protExpr_mean.well)
summary(protExpr_mean.poorly)

