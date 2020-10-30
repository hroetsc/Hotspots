### HEADER ###
# HOTSPOT PREDICTION
# description: fetch and analyse results from hyperparameter optimisation
# input: metrics, predictions
# output: optimal hyperparameters
# author: HR


library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(paletteer)

JOBID = c("5358701", "5365339")
target = "val_loss"  # varibale that should have been minimised

### INPUT ###
system("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/hyperopt/*opt* results/hyperopt/")

metrics = list.files("results/hyperopt", pattern = "metrics", full.names = T)
preds = list.files("results/hyperopt", pattern = "prediction", full.names = T)


### MAIN PART ###
# get losses from file names
targetVar = str_split_fixed(preds, "[:punct:]", Inf)[, 6]
targetVar = paste0("0.", substr(targetVar, 2, nchar(targetVar[1]))) %>%
  as.numeric()


########### linear fit on predictions ########### 
# --> get slope, intercept, PCC, and R^2
linear.fit = data.frame(target = targetVar,
                        slope = rep(NA, length(preds)),
                        intercept = rep(NA, length(preds)),
                        PCC = rep(NA, length(preds)),
                        Rsq = rep(NA, length(preds)))

for (l in 1:nrow(linear.fit)) {
  cnt.pred = read.csv(preds[l], stringsAsFactors = F) %>% na.omit()
  pred.lm = lm(pred_count ~ count, data = cnt.pred)
  
  linear.fit$intercept[l] = pred.lm$coefficients[1]
  linear.fit$slope[l] = pred.lm$coefficients[2]
  linear.fit$Rsq[l] = summary(pred.lm)$r.squared
  linear.fit$PCC[l] = cor(cnt.pred$count, cnt.pred$pred_count, method = "pearson")
  
}

# plot fits
ggplot() + 
  scale_y_continuous(limits = c(-1, 6)) +
  scale_x_continuous(limits = c(-1, 6)) +
  geom_abline(data = linear.fit, aes(slope=slope, intercept=intercept, color=PCC)) +
  scale_color_viridis_c("PCC", direction = -1) +
  ggtitle("linear fit on predicted counts",
          subtitle = "for different hyperparameter combinations") +
  xlab("true count") +
  ylab("predicted count") +
  coord_equal() +
  geom_abline(slope = 1, intercept = 0, lty = "dashed") +
  theme_bw()
ggsave(filename = "results/hyperopt/linear_fit.png", plot = last_plot(),
       device = "png", dpi = "retina")


########### training metrics ###########
for (l in 1:length(metrics)) {
  cnt.metrics = read.csv(metrics[l], stringsAsFactors = F, header = F)
  # format table
  {
    cnt.metrics$V2 = str_split_fixed(cnt.metrics$V2, coll("["), Inf)[,2]
    cnt.metrics$V2 = str_split_fixed(cnt.metrics$V2, coll("]"), Inf)[,1]
    
    var = cnt.metrics$V1
    val = str_split_fixed(cnt.metrics$V2, coll(","), Inf) %>% as.data.frame()
    cnt.metrics = cbind(var, val)
    
    cnt.metrics = t(cnt.metrics) %>% as.data.frame()
    cnt.metrics = cnt.metrics[-1,]
    
    epochs = as.numeric(seq(1, nrow(cnt.metrics)))
    
    rownames(cnt.metrics) = epochs
    colnames(cnt.metrics) = var
    
    # convert factors into numeric
    for (c in 1:ncol(cnt.metrics)){
      cnt.metrics[,c] = as.numeric(as.character(cnt.metrics[,c]))
    }
  }
  cnt.metrics$epochs = seq(nrow(cnt.metrics))
  cnt.metrics$target = rep(targetVar[l], nrow(cnt.metrics))
  
  if (l == 1) {
    all.metrics = cnt.metrics
  } else {
    all.metrics = rbind(all.metrics, cnt.metrics)
  }
}


for (c in 1:(ncol(all.metrics)-2)) {
  cnt.name = colnames(all.metrics)[c] %>% str_replace_all(pattern = "_", " ")
  
  ggplot(data = all.metrics, aes(x = epochs, y = all.metrics[, c], color = target)) +
    geom_point() +
    geom_path(aes(group = target)) +
    scale_color_viridis_c("validation\nloss", direction = 1) +
    ylab(cnt.name) +
    xlab("epoch") +
    ggtitle("training metrics for different hyperparameter combinations",
            subtitle = cnt.name) +
    theme_bw()
  
  ggsave(filename = paste0("results/hyperopt/training_metrics_", colnames(all.metrics)[c], ".png"),
         plot = last_plot(), device = "png", dpi = "retina")
}


########### hyperparameter combinations ########### 

params = c("max_lr", "starting_filter", "block_size", "batch_size", "loss")

if(! dir.exists("results/hyperopt/outfiles")) { dir.create("results/hyperopt/outfiles")}
system(paste0("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/hotspots-", JOBID[1], ".*.out results/hyperopt/outfiles/"))
system(paste0("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/hotspots-", JOBID[2], ".*.out results/hyperopt/outfiles/"))
out.fs = list.files(path = "results/hyperopt/outfiles", full.names = T)

all.params = list()

for (o in 1:length(out.fs)) {
  tmp = read.delim(out.fs[o], header = F)
  k = which(tmp == "LOCALS:") + 1 
  
  if (length(k) > 0){
    p = matrix(ncol = length(params), nrow = length(k))  
    colnames(p) = params
    p = data.frame(p)
    
    for (i in 1:length(k)){
      cnt.param = tmp[k[i], ] %>% as.character() %>%
        str_remove_all(pattern = coll("{")) %>%
        str_remove_all(pattern = coll("}")) %>%
        str_split_fixed(pattern = coll(","), Inf) %>%
        t() %>%
        str_split_fixed(pattern = coll(":"), Inf) %>%
        as.data.frame()
      cnt.param$V1 = str_remove_all(cnt.param$V1, pattern = "'") %>% 
        str_remove_all(pattern = " ") %>% 
        as.character()
      cnt.param$V2 = str_remove_all(cnt.param$V2, pattern = " ") %>% 
        as.character()
      
      for (j in 1:ncol(p)){
        p[i, j] = cnt.param$V2[cnt.param$V1 == colnames(p)[j]] %>% as.character()
      }
    }
    
  } else {
    p = matrix(ncol = length(params), nrow = 1)  
    colnames(p) = params
    p = data.frame(p)
    
  }
  all.params[[o]] = p
}


all.params = plyr::ldply(all.params) %>% as.data.frame() %>% na.omit()
all.params$target = all.params$loss %>% as.numeric()

master = left_join(linear.fit, all.params)

### OUTPUT ###
write.csv(linear.fit, "results/hyperopt/linear_fits.csv", row.names = F)
write.csv(all.params, "results/hyperopt/hyperparam_choices.csv", row.names = F)
write.csv(master, "results/hyperopt/params.csv", row.names = F)



