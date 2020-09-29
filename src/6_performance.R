### HEADER ###
# HOTSPOT PREDICTION
# description: download results from cluster, plot performance curves, evaluate predictions
# input: -
# output: evaluation results
# author: HR

library(dplyr)
library(stringr)
library(rhdf5)
library(ggplot2)
library(tidyr)
library(zoo)


setwd("Documents/QuantSysBios/ProtTransEmbedding/Hotspots/AAEmbeddings/")
JOBID = "5365339-8"


### INPUT ###
# download results
system("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/model_metrics.txt results/")
system("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/best_model_prediction_rank0.csv results/")
system("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/last_model_prediction_rank0.csv results/")
system("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/last_model_prediction_rank1.csv results/")
system("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/best_model_prediction.csv results/")
system("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/last_model_prediction.csv results/")

metrics = read.table("results/5365339-8/model_metrics.txt", sep = ",", stringsAsFactors = F)
prediction = read.csv("results/5365339-8/best_model_prediction.csv", stringsAsFactors = F)



### MAIN PART ###

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
  
  
  
  if (min(col1, col2) < 0){
    col1 = log2(col1 + abs(min(col1, col2))+1)
    col2 = log2(col2 + abs(min(col1, col2))+1)
    
  } else {
    col1 = log2(col1 + 1)
    col2 = log2(col2 + 1)
  }
  
  lower = min(col1, col2) - 1
  upper = max(col1, col2) + 1
  
  
  plot(col1 ~ epochs,
       main = paste0(str_replace_all(name, "_", " "), " for training and validation data set"),
       pch = 20,
       cex = 0.8,
       col = "darkblue",
       ylim = c(lower, upper),
       xlab = "epoch",
       ylab = paste0('log2 ', str_replace_all(name, "_", " ")))
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
  
  k = which(prediction$count == 0)
  # plot(density(prediction$pred_count[k]))
  
  k %>% length() / nrow(prediction)
  which(prediction$count == 0 & prediction$pred_count > 0) %>% length() / nrow(prediction)
  
  which(prediction$count == 0 & prediction$pred_count == 0) %>% length()
  
  # prediction$count = log10(2^(prediction$count))
  # prediction$pred_count = log10(2^(prediction$pred_count))
  
  prediction = na.omit(prediction)
  
  prediction$Accession = NULL
  prediction$window = NULL
  
  
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

start = min(prediction) - 0.1
stop = max(prediction) + 0.1

# visual
prediction[, c("count", "pred_count")] %>% gather() %>%
  ggplot(aes(x = value, color = key)) +
  geom_density() +
  ggtitle("true and predicted hotspot counts") +
  theme_bw()
ggsave(paste0("results/plots/", JOBID, "_trueVSpredicted-dens-log10.png"), plot = last_plot(),
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
ggsave(paste0("results/plots/", JOBID, "_trueVSpredicted-scatter-log10.png"), plot = last_plot(),
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



########## map substrates ##########
window_size = 25

prots = read.csv("data/proteins_w_hotspots.csv", stringsAsFactors = F, header = T)
prediction = read.csv("results/5365339-8/best_model_prediction.csv", stringsAsFactors = F)

# prediction$count = 2^(prediction$count) - 1
# prediction$pred_count = 2^(prediction$pred_count) - 1

prots = prots[which(prots$Accession %in% prediction$Accession), ]


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
    # 
    # if (! any(is.na(idx))){
    #   cnt.True[idx[1]:idx[2]] = rep(cnt.Data$count[j], window_size)
    #   cnt.Pred[idx[1]:idx[2]] = rep(cnt.Data$pred_count[j], window_size)
    #   
    # }
    # 
    
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
  # max count
  stop = max(y_pred %>% na.omit, y_true %>% na.omit()) + .2
  
  
  # rectangle where true count is 0 and predicted > 0
  idx = which(y_true == 0 & y_pred > 0)
  diff(idx)
  c = which(diff(idx) != 1)
  
  
  plot(x, y_true,
       type = "l",
       col="firebrick",
       main = prots$Accession[k],
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
  
  if(length(idx) > 0){
    if(length(c) > 0){
      pos1 = 1
      pos2 = 1
      for (j in 1:length(c)+1){
        cnt_idx = idx[pos1:c[pos2]]
        rect(min(cnt_idx), start, max(cnt_idx), stop, col = rgb(1,1,0, alpha = .5), lwd = 0)
        
        pos1 = c[pos2] + 1
        pos2 = pos2 + 1
        
        if (j > length(c)){
          cnt_idx = idx[pos1:length(idx)]
          rect(min(cnt_idx), start, max(cnt_idx), stop, col = rgb(1,1,0, alpha = .5), lwd = 0)
        }
      }
    } else {
      rect(min(idx), start, max(idx), stop, col = rgb(1,1,0, alpha = .5), lwd = 0)
    }
  }
  
  axis(1)
  axis(2)
  
}

dev.off()



