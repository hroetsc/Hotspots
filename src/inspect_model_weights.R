### HEADER ###
# HOTSPOT PREDICTION
# description: examine model weights for different stages and GPUs
# input: best and last models
# output: -
# author: HR


library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(rhdf5)
library(ggplot2)


### INPUT ###
models = list.files(path = "results", pattern = ".h5", all.files = T, recursive = T, full.names = T)
best_models = models[str_detect(models, pattern = "best_model")]
last_models = models[str_detect(models, pattern = "last_model")]


### MAIN PART ###
get_weights_and_bias = function(IDs_of_interest, descriptions, layer_of_interest, outname){
  
  layers = h5ls(best_models[str_detect(best_models, IDs_of_interest[1])][1])
  model_weights = h5read(best_models[str_detect(best_models, IDs_of_interest[1])][1],
                         "/model_weights")
  
  kernels = list()
  biases =  list()
  
  for (i in 1:length(IDs_of_interest)) {
    cnt = best_models[str_detect(best_models, IDs_of_interest[i])]
    
    cnt.bias = rep(NA, length(cnt))
    cnt.kernel.dens = list()
    
    for (j in 1:length(cnt)) {
      cnt.rank = h5read(cnt[j], layer_of_interest)
      cnt.bias[j] = cnt.rank$`bias:0`
      cnt.kernel.dens[[j]] = density(cnt.rank$`kernel:0`)
    }
    
    kernels[[i]] = cnt.kernel.dens
    biases[[i]] = cnt.bias
    
  }
  
  names(kernels) = IDs_of_interest
  names(biases) = IDs_of_interest
  
  # plotting
  # biases
  biases_df = plyr::ldply(biases) %>%
    t()
  colnames(biases_df) = descriptions
  biases_df = biases_df[-1, ] %>%
    as.data.frame() %>%
    gather()
  
  ggplot(biases_df, aes(as.factor(key), as.numeric(value), fill = key)) +
    geom_violin(scale = "width", trim = F,
                draw_quantiles = c(0.5)) +
    xlab("description") +
    ylab("regression bias") +
    scale_fill_viridis_d() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30),
          legend.position = "none")
  ggsave(filename = paste0("results/exploratory/best-model-weights_bias_", outname, ".png"), plot = last_plot(),
         device = "png", dpi = "retina")
  
  # kernels
  png(paste0("results/exploratory/best-model-weights_kernel-dens_", outname, ".png"),
      width = 12, height = 8, res = 300, units = "in")
  plot(kernels[[1]][[1]],
       main = "kernel density",
       col = rainbow(length(kernels))[1])
  
  for (i in 1:length(kernels)) {
    for (j in 1:length(kernels[[i]])) {
      
      points(kernels[[i]][[j]],
             type = "l",
             col = rainbow(length(kernels))[i])
      
    }
  }
  legend("topleft",
         legend = descriptions,
         lty = rep(1, length(kernels)),
         col = rainbow(length(kernels)),
         cex = .7,
         bg = "transparent",
         bty = "n")
  dev.off()
}

###### noise effect ###### 
IDs_of_interest = c("5691486-1",  # binary classification
                    "5691486-2",  # bc + noise to input windows
                    "5693780-0",  # bc + noise to counts
                    "5693780-2")  # bc + noise to counts + noise to windows

descriptions = c("binary classification",
                 "bc + noise to input windows",
                 "bc + noise to counts",
                 "bc + noise to counts + noise to windows")
get_weights_and_bias(IDs_of_interest, descriptions, "/model_weights/regression/regression", outname = "noise-effect")


###### using only aa indices ######
IDs_of_interest = c("5693780-3",
                    "5790191-1")

descriptions = c("AA indices and one-hot (50 aa)",
                 "only AA indices (50 aa, relu)")

get_weights_and_bias(IDs_of_interest, descriptions, "/model_weights/regression/regression", outname = "training-data")


###### counts scale ######
IDs_of_interest = c("5693780-3",
                    "5714916-1",
                    "5714916-3",
                    "5761109-0")

descriptions = c("log counts",
                 "raw counts",
                 "scale raw between 0-1 (linear act)",
                 "scale raw between 0-1 (sigmoid act)")

get_weights_and_bias(IDs_of_interest, descriptions, "/model_weights/regression/regression", outname = "counts-scaling")


###### oversampling ######
IDs_of_interest = c("5693780-2",
                    "5790191-2",
                    "5790191-3",
                    "5807441-0")

descriptions = c("no oversampling (25 aa)",
                 "oversampling, fewer n, relu",
                 "oversampling",
                 "oversampling, fewer n, higher noise, relu, SMAPE")

get_weights_and_bias(IDs_of_interest, descriptions, "/model_weights/regression/regression", outname = "oversampling")


###### loss function ######
IDs_of_interest = c("5693780-2",
                    "5807441-0",
                    "5807441-2")
descriptions = c("MSE loss",
                 "SMAPE, oversampling, fewer n, higher noise, relu",
                 "SMSPE, higher noise, relu")

get_weights_and_bias(IDs_of_interest, descriptions, "/model_weights/regression/regression", outname = "loss-function")


