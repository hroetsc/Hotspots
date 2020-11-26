### HEADER ###
# HOTSPOT PREDICTION
# description: characterise and compare different loss functions
# input: predictions (from 7_performance.R)
# output: -
# author: HR


library(dplyr)


ae = abs(prediction$count - prediction$pred_count)

png("results/ae-error_vs_counts.png")
plot(ae~prediction$count,
     pch = 16,
     cex = .1,
     ylab = "absolute error",
     xlab = "true count")
dev.off()

se = ae^2
png("results/se-error_vs_counts.png")
plot(se~prediction$count,
     pch = 16,
     cex = .1,
     ylab = "squared error",
     xlab = "true count")
dev.off()

# distance to mean vs number of data points (binned)
counts = prediction$count
bins = seq(-0.1, max(counts), 0.1)

dist.to.mean = data.frame(bin = bins,
                          dist_to_mean = abs(bins - mean(counts)),
                          no_windows = NA)

for (b in 2:length(bins)){
  dist.to.mean$no_windows[b] = which(counts >= bins[b-1] &
                                       counts <= bins[b]) %>% length()
}

dist.to.mean = na.omit(dist.to.mean)
plot(log10(dist.to.mean$no_windows) ~ dist.to.mean$dist_to_mean)
points(log10(dist.to.mean$no_windows) ~ dist.to.mean$bin, col = "red")


# mean squared error, median squared error, mean absolute percentage error, symmetric mape
abs(prediction$count - prediction$pred_count)^2 %>% mean()
abs(prediction$count - prediction$pred_count)^2 %>% median()
(abs(prediction$count - prediction$pred_count)/(prediction$count+1e-07)) %>% mean()

y = prediction$count
y_hat = prediction$pred_count

mape = mean(abs(y - y_hat)/(y + 1e-12))
smape1 = mean(abs(y_hat - y) / (y + 1e-12 + y_hat))
smape2 = sum(abs(y_hat - y)) / sum(y + y_hat)

which(prediction$pred_count == 0 & prediction$count == 0) %>% length()


# compare smape with mse for imbalanced data sets and random prediction
smape2 = function(y_hat, y){ return(sum(abs(y_hat - y)) / sum(y + y_hat)) }
mse = function(y_hat, y) { return(mean(abs(y_hat - y)^2)) }

# imbalanced data sets
smape_all = rep(NA, 50)
mse_all = rep(NA, 50)


for (i in 1:50){
  y_model = c(seq(0,1, length.out = 100*i),
              seq(1,6, length.out = 100))
  y_pred = runif(min = 0, max = 6, n = (100 + 100*i))
  
  smape_all[i] = smape2(y_model, y_pred)
  mse_all[i] = mse(y_model, y_pred)
  
}

plot(mse_all, type = "l")
x11()
plot(smape_all, type = "l", col = "blue")

# real data ratio
length(y[y <= 1]) / length(y[y >= 1])


# upper limit of predictions
limits = seq(3,6,0.1)
smape_all = rep(NA, length(limits))
mse_all = rep(NA, length(limits))

for (i in 1:length(limits)){
  y_model = c(seq(0,1, length.out = 100*4),
              seq(1,6, length.out = 100))
  y_pred = runif(min = 0, max = limits[i], n = (100 + 100*4))
  
  smape_all[i] = smape2(y_model, y_pred)
  mse_all[i] = mse(y_model, y_pred)
  
}

plot(limits, mse_all, type = "l")
x11()
plot(limits, smape_all, type = "l", col = "blue")

# real data limit
max(y_hat)

sape = abs(prediction$pred_count - prediction$count) / (prediction$count + prediction$pred_count)
png("results/sape-error_vs_counts.png")
plot(sape~prediction$count,
     pch = 16,
     cex = .1,
     ylab = "symmetric absolute percentage error",
     xlab = "true count")
dev.off()
