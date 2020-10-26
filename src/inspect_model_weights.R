### HEADER ###
# HOTSPOT PREDICTION
# description: examine model weights for different stages and GPUs
# input: best and last models
# output: -
# author: HR


library(dplyr)
library(rhdf5)


### INPUT ###
model_path = "results/5542543-1/model/"


### MAIN PART 
# select last hidden layer (1024-d fully connected ReLU)
layers = h5ls(paste0(model_path, "best_model_rank0.h5"))

# load layers
model_onehot.best = h5read(paste0(model_path, "best_model_rank0.h5"), "/model_weights/dense")
model_onehot.last = h5read(paste0(model_path, "last_model_rank0.h5"), "/model_weights/dense")

model_aaindex.best = h5read(paste0(model_path, "best_model_rank1.h5"), "/model_weights/dense")
model_aaindex.last = h5read(paste0(model_path, "last_model_rank1.h5"), "/model_weights/dense")

density(model_onehot.best[["dense"]][["kernel:0"]]) %>% plot(main = "kernels of last hidden layer",
                                                             col = "darkgreen")
  density(model_onehot.last[["dense"]][["kernel:0"]]) %>% lines(col = "green")
density(model_aaindex.best[["dense"]][["kernel:0"]]) %>% lines(col = "firebrick")
density(model_aaindex.last[["dense"]][["kernel:0"]]) %>% lines(col = "red")
  
