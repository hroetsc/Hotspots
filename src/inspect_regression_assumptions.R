### HEADER ###
# HOTSPOT PREDICTION
# description: check if least square regression assumptions are fulfilled
# input: predictions (from 7_performance.R)
# output: -
# author: HR


library(ggplot2)

res = (prediction$count = prediction$pred_count)^2


##### normality of residuals #####

png("results/exploratory/qq-plot-residuals.png")
qqnorm(res, pch = 1, frame = F)
qqline(res, col = "steelblue", lwd = 2)
dev.off()

mean(res)


##### homoscedasticity #####
homoscedasticity = data.frame(predicted = prediction$pred_count,
                              residual = res)

ggplot(homoscedasticity, aes(x = predicted, y = res)) +
  geom_point() +
  geom_path()
ggsave(filename = "results/exploratory/homoscedasticity-residuals.png",
       plot = last_plot(), dpi = "retina")



