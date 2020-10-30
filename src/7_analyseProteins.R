### HEADER ###
# HOTSPOT PREDICTION
# description: identify proteins with good/bad predictions and analyse them
# input: predictions, JOBID (loaded and processed in 6_performance.R)
# output: good/bad predicted proteins, plots
# author: HR


library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)


### INPUT ###
prediction = read.csv("results/last_model_prediction.csv", stringsAsFactors = F)
hp = read.csv("/media/hanna/Hanna2/DATA/GENERAL/proteome_human.csv", stringsAsFactors = F)
props = read.csv("/media/hanna/Hanna2/DATA/Seq2Vec/RUNS/HumanProteome/PropMatrix_seqs.csv",
                 stringsAsFactors = F)
props$X = NULL
props = left_join(hp, props)


### MAIN PART ###
prediction = na.omit(prediction)
prots = prediction$Accession %>% unique()

######### general properties #########
overview = data.frame(Accession = prots,
                      R_sq = rep(NA, length(prots)),
                      PCC = rep(NA, length(prots)),
                      len_protein = rep(NA, length(prots)),
                      true_max_count = rep(NA, length(prots)),
                      true_min_count = rep(NA, length(prots)),
                      true_average_count = rep(NA, length(prots)),
                      pred_max_count = rep(NA, length(prots)),
                      pred_min_count = rep(NA, length(prots)),
                      pred_average_count = rep(NA, length(prots)))

pb = txtProgressBar(min = 0, max = nrow(overview), style = 3)
for (i in 1:nrow(overview)) {
  
  setTxtProgressBar(pb, i)
  
  cnt = prediction[prediction$Accession == overview$Accession[i], ]
  cnt.lm = lm(pred_count ~ count, data = cnt)
  
  overview$PCC[i] = cor(cnt$count, cnt$pred_count, method = "pearson")
  overview$R_sq[i] = summary(cnt.lm)$r.squared
  overview$len_protein[i] = nchar(hp$seqs[hp$Accession == overview$Accession[i]])
  
  overview$true_max_count[i] = max(cnt$count)
  overview$true_min_count[i] = min(cnt$count)
  overview$true_average_count[i] = cnt$count %>% mean()
  
  overview$pred_max_count[i] = max(cnt$pred_count)
  overview$pred_min_count[i] = min(cnt$pred_count)
  overview$pred_average_count[i] = cnt$pred_count %>% mean()
  
}
write.csv(overview, "results/overview.csv", row.names = F)

# histogram of protein-level PCC's
png("results/performance_protein_wise.png")
q = quantile(overview$PCC, c(0.05, 0.5, 0.95))
m = overview$PCC %>% mean()
hist(overview$PCC, freq = F, breaks = 100,
     main = "histogram of PCCs of all test proteins",
     xlab = 'PCC')
abline(v = q[1], col = "red")
abline(v = q[2], col = "red")
abline(v = m, col = "blue")
abline(v = q[3], col = "red")
dev.off()


# plot PCC / R^2 vs. protein features
pdf(file = "results/performance_protein_wise.pdf", height = 12, width = 24)
par(mfrow = c(2,3))
plot(overview$PCC ~ overview$max_count, pch = 20, col = "red", ylim = c(0,6))
points(overview$PCC ~ overview$min_count, pch = 20, col = "blue")
plot(overview$PCC ~ log10(overview$average_count), pch = 20, col = "black")
plot(overview$PCC ~ log10(overview$len_protein), pch = 20, col = "black")


plot(overview$R_sq ~ overview$max_count, pch = 20, col = "red", ylim = c(0,6))
points(overview$R_sq ~ overview$min_count, pch = 20, col = "blue")
plot(overview$R_sq ~ log10(overview$average_count), pch = 20, col = "black")
plot(overview$R_sq ~ log10(overview$len_protein), pch = 20, col = "black")
dev.off()


######### scatterplot by features #########
pred_with_features = left_join(prediction, overview)
pred_with_features = left_join(pred_with_features, props)

start = min(pred_with_features$count, pred_with_features$pred_count) - 0.1
stop = max(pred_with_features$count, pred_with_features$pred_count) + 0.1

{
ggplot(pred_with_features, aes(x = count, y = pred_count, col = log10(len_protein))) +
  geom_point(alpha = 0.05, size = 0.1) +
  scale_color_gradient(low = 'darkblue', high = 'orangered') +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  theme_bw()
ggsave("results/scatterplot_protein-length.png", plot = last_plot(),
       device = "png", dpi = "retina")

ggplot(pred_with_features, aes(x = count, y = pred_count, col = average_count)) +
  geom_point(alpha = 0.05, size = 0.1) +
  scale_color_gradient(low = 'darkblue', high = 'orangered') +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  theme_bw()
ggsave("results/scatterplot_average-count.png", plot = last_plot(),
       device = "png", dpi = "retina")

ggplot(pred_with_features, aes(x = count, y = pred_count, col = pI)) +
  geom_point(alpha = 0.05, size = 0.1) +
  scale_color_gradient(low = 'darkblue', high = 'orangered') +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  theme_bw()
ggsave("results/scatterplot_pI.png", plot = last_plot(),
       device = "png", dpi = "retina")

ggplot(pred_with_features, aes(x = count, y = pred_count, col = Hydrophobicity)) +
  geom_point(alpha = 0.05, size = 0.1) +
  scale_color_gradient(low = 'darkblue', high = 'orangered') +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  theme_bw()
ggsave("results/scatterplot_hydrophobicity.png", plot = last_plot(),
       device = "png", dpi = "retina")

ggplot(pred_with_features, aes(x = count, y = pred_count, col = BLOSUM1)) +
  geom_point(alpha = 0.05, size = 0.1) +
  scale_color_gradient(low = 'darkblue', high = 'orangered') +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  theme_bw()
ggsave("results/scatterplot_BLOSUM1.png", plot = last_plot(),
       device = "png", dpi = "retina")

ggplot(pred_with_features, aes(x = count, y = pred_count, col = Polarity)) +
  geom_point(alpha = 0.05, size = 0.1) +
  scale_color_gradient(low = 'darkblue', high = 'orangered') +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  theme_bw()
ggsave("results/scatterplot_polarity.png", plot = last_plot(),
       device = "png", dpi = "retina")

ggplot(pred_with_features, aes(x = count, y = pred_count, col = H.bonding)) +
  geom_point(alpha = 0.05, size = 0.1) +
  scale_color_gradient(low = 'darkblue', high = 'orangered') +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  theme_bw()
ggsave("results/scatterplot_H-bonding.png", plot = last_plot(),
       device = "png", dpi = "retina")
}


########## characterise outliers ##########
high.counts.bad = overview[overview$average_count > 1 & overview$PCC < .3, ]
write.csv(high.counts.bad, "results/badly_predicted_high_counts.csv", row.names = F)

high.counts.well = overview[overview$average_count > 1 & overview$PCC > .7, ]
write.csv(high.counts.well, "results/well_predicted_high_counts.csv", row.names = F)


prediction.categories = prediction
prediction.categories$color = "undefined"
prediction.categories$color[prediction.categories$Accession %in% high.counts$Accession] = "high, bad prediction"
prediction.categories$color[prediction.categories$Accession %in% high.counts.well$Accession] = "high, good prediction"

ggplot(prediction.categories, aes(x = count, y = pred_count, col = color)) +
  geom_point(alpha = 0.05, size = 0.1) +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  theme_bw()
ggsave(paste0("results/plots/", JOBID, "_trueVSpredicted-scatter-categories.png"), plot = last_plot(),
       device = "png", dpi = "retina")


