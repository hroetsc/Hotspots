### HEADER ###
# HOTSPOT PREDICTION
# description: download results from cluster, plot performance curves, evaluate predictions
# input: -
# output: evaluation results
# author: HR

library(dplyr)
library(stringr)
library(data.table)
library(rhdf5)
library(ggplot2)
library(tidyr)
library(tidymodels)
library(DescTools)
library(zoo)

JOBID = "5636045-2-last-nolog"


### INPUT ###
# download results
if(!dir.exists(paste0("results/", JOBID))) {
  dir.create(paste0("results/", JOBID))
  dir.create(paste0("results/", JOBID, "/model"))
}

system(paste0("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/model_metrics.txt results/", JOBID, "/"))
system(paste0("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/*_model_prediction_rank* results/", JOBID, "/"))
system(paste0("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/*_model_*.h5 results/", JOBID, "/model/"))

metrics = read.table(paste0("results/", JOBID, "/model_metrics.txt"),
                     sep = ",", stringsAsFactors = F)
prediction = read.csv("results/5636045-2/last_model_prediction_rank0.csv",
                      stringsAsFactors = F)



### MAIN PART ###
########## combine predictions of all GPUs ########## 
preds = list.files("results/5636045-2",
                   pattern = "last_model_prediction_rank",
                   full.names = T)

pred_counts_onehot = rep(0, nrow(prediction))
pred_counts_aaindex = rep(0, nrow(prediction))
pred_counts = rep(0, nrow(prediction))

pb = txtProgressBar(min = 0 , max = length(preds), style = 3)
for (p in 1:length(preds)) {
  setTxtProgressBar(pb, p)
  
  rank = str_split(preds[p], "[:punct:]", simplify = T)[, 7]
  rank = str_sub(rank, start = 5, end = nchar(rank)) %>% as.character %>% as.numeric()
  
  cnt_pred = read.csv(preds[p], stringsAsFactors = F)
  
  if (rank %% 2 == 0){
    pred_counts_onehot = pred_counts_onehot + (cnt_pred$pred_count * 1/(length(preds)*0.5))
  } else {
    pred_counts_aaindex = pred_counts_aaindex + (cnt_pred$pred_count * 1/(length(preds)*0.5))
  }

  pred_counts = pred_counts + (cnt_pred$pred_count * 1/length(preds))
}

prediction$pred_count = pred_counts
# prediction$pred_count = pred_counts_onehot
# prediction$pred_count = pred_counts_aaindex

count = prediction$count
pred_count = prediction$pred_count

prediction$count = 2^(prediction$count) - 1
prediction$pred_count = 2^(prediction$pred_count) - 1


########## training metrics ##########

{
  metrics$V2 = str_split_fixed(metrics$V2, coll("["), Inf)[,2]
  metrics$V2 = str_split_fixed(metrics$V2, coll("]"), Inf)[,1]
  
  var = metrics$V1
  val = str_split_fixed(metrics$V2, coll(","), Inf) %>% as.data.frame()
  metrics = cbind(var, val)
  
  metrics = t(metrics) %>% as.data.frame()
  metrics = metrics[-1,]
  
  epochs = as.numeric(seq(1, nrow(metrics)))
  
  rownames(metrics) = epochs
  colnames(metrics) = var
  
  # convert factors into numeric
  for (c in 1:ncol(metrics)){
    metrics[,c] = as.numeric(as.character(metrics[,c]))
  }
}

# plotting function
plotting = function(col1 = "", col2 = "", name = ""){
  
  lower = min(col1, col2) - .2
  upper = max(col1, col2) + .2
  
  
  plot(col1 ~ epochs,
       main = paste0(str_replace_all(name, "_", " "), " for training and validation data set"),
       pch = 20,
       cex = 0.8,
       col = "darkblue",
       ylim = c(lower, upper),
       xlab = "epoch",
       ylab = str_replace_all(name, "_", " "))
  points(col2 ~ epochs,
         pch = 20,
         cex = 0.8,
         col = "firebrick")
  legend("topright", cex = 1,
         legend = c("training", "validation"),
         col = c("darkblue", "firebrick"),
         pch = c(20, 20),
         box.lty = 1,
         pt.cex = 1)
}

{
  pdf(paste0("results/plots/",JOBID,"_metrics.pdf"), width = 12, height = 12)
  par(mfrow = c(3, 2))
  # plot and save
  for (i in 1:(ncol(metrics)/2)){
    plotting(col1 = metrics[,i],
             col2 = metrics[, (i + (ncol(metrics)/2))],
             name = colnames(metrics)[i])
  }
  
  # generalization loss
  gl = 100 * ((metrics$val_loss/min(metrics$val_loss)) - 1)
  plot(log2(gl), col='firebrick',
       ylab = 'log2 generalisation loss',
       pch = 20,
       cex = 0.8,
       xlab = 'epoch',
       main = 'generalisation loss during training')
  dev.off()
}

########## regression ##########

# general
{
  summary(prediction$count) %>% print()
  summary(prediction$pred_count) %>% print()
  
  prediction = na.omit(prediction)
  
  # linear model --> R^2 and adjusted R^2
  pred.lm = lm(pred_count ~ count, data = prediction)
  summary(pred.lm)
  
  
  # correlation coefficients
  pc = cor(prediction$count, prediction$pred_count, method = "pearson")
  sm = cor(prediction$count, prediction$pred_count, method = "spearman")
  
  # mean squared error
  mse = (prediction$count - prediction$pred_count)^2 %>% mean() %>% round(4)
  # root mean squared error
  rmse = sqrt(mse) %>% round(4)
  # mean absolute deviation
  mae = (prediction$count - prediction$pred_count) %>% abs() %>% mean() %>% round(4)
  
  # sumarise
  all.metrics = c(JOBID, summary(pred.lm)$r.squared, pc, mse, rmse, mae)
  names(all.metrics) = c("JOBID", "Rsquared", "PCC", "MSE", "RMSE", "MAE")
  all.metrics
}


start = min(prediction$count, prediction$pred_count) - 0.1
stop = max(prediction$count, prediction$pred_count) + 0.1

# visual
prediction[, c("count", "pred_count")] %>% gather() %>%
  ggplot(aes(x = value, color = key)) +
  geom_density() +
  ggtitle("true and predicted hotspot counts") +
  theme_bw()
ggsave(paste0("results/plots/", JOBID, "_trueVSpredicted-dens.png"), plot = last_plot(),
       device = "png", dpi = "retina")

ggplot(prediction, aes(x = count, y = pred_count)) +
  geom_point(alpha = 0.1, size = 0.1) +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  ggtitle("true and predicted hotspot counts",
          subtitle = paste0("PCC: ", pc %>% round(4), ", R^2: ", summary(pred.lm)$r.squared %>% round(4))) +
  theme_bw()
ggsave(paste0("results/plots/", JOBID, "_trueVSpredicted-scatter.png"), plot = last_plot(),
       device = "png", dpi = "retina")

# append to past metrics
out = "results/performance.csv"
if(file.exists(out)) {
  write(all.metrics, file = out, ncolumns = length(all.metrics),
        append = T, sep = ",")
  
} else {
  
  write(all.metrics, file = out, ncolumns = length(all.metrics),
        append = F, sep = ",")
  
}



########## cluster proteins ##########
# prediction = read.csv("results/5365339-11/last_model_prediction_rank0.csv", stringsAsFactors = F)
# prediction$pred_count = pred_counts

prots = prediction$Accession %>% unique()
dist = list()

pb = txtProgressBar(min = 0, max = length(prots), style = 3)
for (i in 1:length(prots)){
  setTxtProgressBar(pb, i)
  k = which(prediction$Accession == prots[i])
  dist[[i]] = prediction$count[k] %>% as.numeric()
}

# pairwise distances
pairwise.dist = matrix(nrow = length(prots), ncol = length(prots))
colnames(pairwise.dist) = prots
rownames(pairwise.dist) = prots

for (i in 1:length(prots)) {
  setTxtProgressBar(pb, i)
  for (j in 1:length(prots)){
    pairwise.dist[i, j] = ks.test(dist[[i]], dist[[j]], alternative = "two.sided")$statistic
  }
}

write.csv(pairwise.dist, "results/pairwise_ks-statistics_test100-sample.csv", row.names = F)
pairwise.dist = read.csv("results/pairwise_ks-statistics_test100-sample.csv", stringsAsFactors = F)


# clustering
km.results = kmeans(pairwise.dist, centers = 4, iter.max = 20, nstart = 4)
cluster = km.results$cluster %>% as.numeric() %>% as.data.frame()
cluster$Accession = prots
  
names(cluster)[1] = "n_cluster"
prediction.cl = left_join(prediction, cluster) %>% na.omit()

ggplot(prediction.cl, aes(x = count, y = pred_count, col = factor(n_cluster))) +
  geom_point(alpha = 0.1, size = 0.1) +
  scale_color_viridis_d() +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  ggtitle("true and predicted hotspot counts",
          subtitle = paste0("PCC: ", pc %>% round(4), ", R^2: ", summary(pred.lm)$r.squared %>% round(4))) +
  theme_bw()
ggsave(paste0("results/plots/", JOBID, "_trueVSpredicted-scatter-cluster4.png"), plot = last_plot(),
       device = "png", dpi = "retina")

# analyse clusters
cl = prediction.cl$n_cluster %>% unique()
for (c in cl){
  print(".............................................................")
  print(c)
  summary(prediction.cl$count[prediction.cl$n_cluster == c]) %>% print()
  cnt.pc = cor(prediction$count[prediction.cl$n_cluster == c],
               prediction$pred_count[prediction.cl$n_cluster == c], method = "pearson")
  print(cnt.pc)
  print(".............................................................")
}


########## map peptides to protein sequence ##########
window_size = 50

prots = read.csv("data/proteins_w_hotspots.csv", stringsAsFactors = F, header = T)
prots = prots[which(prots$Accession %in% prediction$Accession), ]
# prots = left_join(prots, cluster)

# get counts for every aa position from token counts
countsTrue = list()
countsPred = list()
countsPred.roll = list()

pb = txtProgressBar(min = 0, max = nrow(prots), style = 3)
for (i in 1:nrow(prots)){
  
  setTxtProgressBar(pb, i)
  
  cnt.Data = prediction[prediction$Accession == prots$Accession[i], ]
  cnt.Data = left_join(cnt.Data, prots[i, ])
  
  cnt.True = rep(NA, nchar(prots$seqs[i]))
  cnt.Pred = rep(NA, nchar(prots$seqs[i]))
  
  
  for (j in 1:nrow(cnt.Data)) {
    
    substr = cnt.Data$window[j]
    idx = str_locate(prots$seqs[i], substr) %>% as.numeric()
    
    cnt.True[idx[1]:idx[2]] = rep(cnt.Data$count[j], window_size)
    cnt.Pred[idx[1]:idx[2]] = rep(cnt.Data$pred_count[j], window_size)
    
  }
  
  countsTrue[[i]] = cnt.True
  countsPred[[i]] = cnt.Pred
  countsPred.roll[[i]] = rollmean(cnt.Pred, k = 9)
  
  names(countsTrue)[i] = prots$Accession[i]
  names(countsPred)[i] = prots$Accession[i]
  names(countsPred.roll[[i]]) = prots$Accession[i]
  
}


# plot
pdf(paste0("results/plots/",JOBID,"_Counts_trueVSpred.pdf"), width = 12, height = 8)
par(mfrow = c(2,2))
for (k in 1:length(countsTrue)) {
  
  y_true = countsTrue[[k]]
  y_pred = countsPred[[k]]
  y_pred.roll = countsPred.roll[[k]]
  x = seq(length(y_true))
  # max and min count
  start = min(y_pred %>% na.omit, y_true %>% na.omit()) - .2
  stop = max(y_pred %>% na.omit, y_true %>% na.omit()) + .2
  
  
  # rectangle where true count is 0 and predicted > 0
  idx = which(y_true == 0 & y_pred > 0)
  diff(idx)
  c = which(diff(idx) != 1)
  
  
  plot(x, y_true,
       type = "l",
       col="firebrick",
       main = prots$Accession[k],
       sub = paste0("CLUSTER: ", prots$n_cluster[k]),
       axes = F,
       ylab = "counts",
       xlab = "position",
       ylim = c(start, stop))
  points(x, y_pred,
         type = "l",
         col="darkblue")
  points(seq(1, length(y_pred.roll), ceiling(length(y_pred.roll)/length(y_pred))),
         y_pred.roll,
         type = "l",
         col="black")
  
  # if(length(idx) > 0){
  #   if(length(c) > 0){
  #     pos1 = 1
  #     pos2 = 1
  #     for (j in 1:length(c)+1){
  #       cnt_idx = idx[pos1:c[pos2]]
  #       rect(min(cnt_idx), start, max(cnt_idx), stop, col = rgb(1,1,0, alpha = .5), lwd = 0)
  #       
  #       pos1 = c[pos2] + 1
  #       pos2 = pos2 + 1
  #       
  #       if (j > length(c)){
  #         cnt_idx = idx[pos1:length(idx)]
  #         rect(min(cnt_idx), start, max(cnt_idx), stop, col = rgb(1,1,0, alpha = .5), lwd = 0)
  #       }
  #     }
  #   } else {
  #     rect(min(idx), start, max(idx), stop, col = rgb(1,1,0, alpha = .5), lwd = 0)
  #   }
  # }
  
  axis(1)
  axis(2)
  
}

dev.off()


########## use rolling mean for correlation ##########
pred.rollmean = c()
true.rollmean = c()

pb = txtProgressBar(min = 0, max = length(countsTrue), style = 3)
for (i in 1:length(countsTrue)) {
  setTxtProgressBar(pb, i)
  
  pred.rollmean = c(pred.rollmean,
                    countsPred.roll[[i]])
  t = countsTrue[[i]]
  true.rollmean = c(true.rollmean,
                    t[5:(length(t)-4)])
  
}

prediction.roll = cbind(true.rollmean, pred.rollmean) %>% as.data.frame() %>% na.omit()
names(prediction.roll) = c("count", "pred_count")
pc = cor(prediction.roll$count,
         prediction.roll$pred_count,
         method = "pearson")
pred.lm = lm(pred_count ~ count, data = prediction.roll)

ggplot(prediction.roll, aes(x = count, y = pred_count)) +
  geom_point(alpha = 0.1, size = 0.1) +
  xlim(c(min(prediction.roll$count, prediction.roll$pred_count),
         max(prediction.roll$count, prediction.roll$pred_count))) +
  ylim(c(min(prediction.roll$count, prediction.roll$pred_count),
         max(prediction.roll$count, prediction.roll$pred_count))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  ggtitle("true and predicted hotspot counts",
          subtitle = paste0("PCC: ", pc %>% round(4), ", R^2: ", summary(pred.lm)$r.squared %>% round(4))) +
  theme_bw()
ggsave(paste0("results/plots/", JOBID, "_trueVSpredicted-scatter-rollmean.png"), plot = last_plot(),
       device = "png", dpi = "retina")


########## binarize counts ##########
prediction.binary = prediction.roll
prediction.binary$count[prediction.binary$count > 0] = 1


roc.pr.CURVE = function(prediction = "") {
  
  PRECISION = function(TP = "", FP = "") {
    return(as.numeric(TP) / (as.numeric(TP) + as.numeric(FP)))
  }
  
  RECALL = function(TP = "", P = "") {
    return(as.numeric(TP) / as.numeric(P))
  }
  
  SENSITIVITY = function(TP = "", P = "") {
    return(as.numeric(TP) / as.numeric(P))
  }
  
  SPECIFICITY = function(TN = "", N = "") {
    return(as.numeric(TN) / as.numeric(N))
  }
  
  th_range = c(-Inf,
               seq(min(prediction$pred_count), max(prediction$pred_count), length.out = 300),
               Inf)
  
  sens = rep(NA, length(th_range))
  spec = rep(NA, length(th_range))
  
  prec = rep(NA, length(th_range))
  rec = rep(NA, length(th_range))
  
  mcc = rec = rep(NA, length(th_range))
  
  pb = txtProgressBar(min = 0, max = length(th_range), style = 3)
  for (t in seq_along(th_range)) {
    setTxtProgressBar(pb, t)
    
    cnt_dat = prediction %>% mutate(pred = ifelse(pred_count > th_range[t], 1, 0))
    
    # P, N, TP, TN, FP
    P = cnt_dat[cnt_dat$count == 1, ] %>% nrow()
    N = cnt_dat[cnt_dat$count == 0, ] %>% nrow()
    
    TP = cnt_dat[cnt_dat$pred == 1 & cnt_dat$count == 1, ] %>% nrow()
    TN = cnt_dat[cnt_dat$pred == 0 & cnt_dat$count == 0, ] %>% nrow()
    
    FP = cnt_dat[cnt_dat$pred == 1 & cnt_dat$count == 0, ] %>% nrow()
    FN = cnt_dat[cnt_dat$pred == 0 & cnt_dat$count == 1, ] %>% nrow()
    
    sens[t] = SENSITIVITY(TP, P)
    spec[t] = SPECIFICITY(TN, N)
    
    prec[t] = PRECISION(TP, FP)
    rec[t] = RECALL(TP, P)
  }
  
  curve = data.frame(score = th_range,
                     precision = prec,
                     recall = rec,
                     sensitivity = sens,
                     specificity = spec)
  
  return(curve)
}

curve = roc.pr.CURVE(prediction = prediction.binary)

# AUC
pr.na = which(! is.na(curve$precision | curve$recall))
pr.auc = AUC(curve$recall[pr.na],
             curve$precision[pr.na])
print(paste0("PR-AUC: ", pr.auc))

roc.na = which(! is.na(curve$sensitivity | curve$specificity))
roc.auc = AUC(curve$specificity[roc.na],
              curve$sensitivity[roc.na])
print(paste0("ROC-AUC: ", roc.auc))


theme_set(theme_bw())
roc.curve = curve %>%
  ggplot() +
  geom_path(aes(1 - specificity, sensitivity)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  ggtitle("ROC curve for binarized hotspot counts",
          subtitle = paste0("AUC: ", roc.auc %>% round(4)))
roc.curve

pr.curve = curve %>%
  ggplot() +
  geom_path(aes(recall, precision)) +
  geom_abline(intercept = length(which(prediction.binary$count == 1))/nrow(prediction.binary),
              slope = 0, linetype = "dotted") +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  ggtitle("PR curve for binarized hotspot counts",
          subtitle = paste0("AUC: ", pr.auc %>% round(4)))
pr.curve


### OUTPUT ###
ggsave(paste0("results/plots/", JOBID, "_ROC.png"),
       roc.curve, device = "png", dpi = "retina")
ggsave(paste0("results/plots/", JOBID, "_PR.png"),
       pr.curve, device = "png", dpi = "retina")


######### just to be sure ... ###########

train = fread("data/windowTokens_training100-sample.csv") %>% as.data.frame()
train.acc = str_split_fixed(train$Accession, coll("-"), Inf)[, 1] %>% unique()
pred.acc = str_split_fixed(prediction$Accession, coll("-"), Inf)[, 1] %>% unique()
intersect(train.acc, pred.acc)
