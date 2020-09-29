### HEADER ###
# HOTSPOT PREDICTION
# description: identify proteins with good/bad predictions and analyse them
# input: predictions
# output: good/bad predicted proteins
# author: HR


library(plyr)
library(dplyr)
library(stringr)
library(GOSemSim)
library(ggplot2)
library(scales)


### INPUT ###
prediction = read.csv("results/best_model_prediction.csv", stringsAsFactors = F)
hp = read.csv("../../files/proteome_human.csv", stringsAsFactors = F)
props = read.csv("../../RUNS/HumanProteome/PropMatrix_seqs.csv", stringsAsFactors = F)
props$X = NULL
props = left_join(hp, props)

### MAIN PART ###
prediction = na.omit(prediction)
prots = prediction$Accession %>% unique()

## general properties
overview = data.frame(Accession = prots,
                      R_sq = rep(NA, length(prots)),
                      PCC = rep(NA, length(prots)),
                      len_protein = rep(NA, length(prots)),
                      max_count = rep(NA, length(prots)),
                      min_count = rep(NA, length(prots)),
                      average_count = rep(NA, length(prots)))

pb = txtProgressBar(min = 0, max = nrow(overview), style = 3)
for (i in 1:nrow(overview)) {
  
  setTxtProgressBar(pb, i)
  
  cnt = prediction[prediction$Accession == overview$Accession[i], ]
  cnt.lm = lm(pred_count ~ count, data = cnt)
  
  overview$PCC[i] = cor(cnt$count, cnt$pred_count, method = "pearson")
  overview$R_sq[i] = summary(cnt.lm)$r.squared
  overview$len_protein[i] = nchar(hp$seqs[hp$Accession == overview$Accession[i]])
  overview$max_count[i] = max(cnt$count)
  overview$min_count[i] = min(cnt$count)
  overview$average_count[i] = cnt$count %>% mean()
  
}


# plot
pdf(file = "results/performance_protein_wise.pdf", height = 12, width = 24)
par(mfrow = c(2,3))
plot(overview$PCC ~ overview$max_count, pch = 20, col = "red", ylim = c(0,6))
points(overview$PCC ~ overview$min_count, pch = 20, col = "blue")
plot(overview$PCC ~ overview$average_count, pch = 20, col = "black")
plot(overview$PCC ~ overview$len_protein, pch = 20, col = "black")


plot(overview$R_sq ~ overview$max_count, pch = 20, col = "red", ylim = c(0,6))
points(overview$R_sq ~ overview$min_count, pch = 20, col = "blue")
plot(overview$R_sq ~ overview$average_count, pch = 20, col = "black")
plot(overview$R_sq ~ overview$len_protein, pch = 20, col = "black")
dev.off()


## scatterplot by features
pred_with_features = left_join(prediction, overview)
pred_with_features = left_join(pred_with_features, props)

start = min(pred_with_features$count, pred_with_features$pred_count) - 0.1
stop = max(pred_with_features$count, pred_with_features$pred_count) + 0.1

{
ggplot(pred_with_features, aes(x = count, y = pred_count, col = len_protein)) +
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

## biological functions
hsGO = godata('org.Hs.eg.db', keytype = "UNIPROT", ont = "BP", computeIC = F)
annot = hsGO@geneAnno

for (i in 1:nrow(overview)) {
  overview$GO_terms[i] = annot[annot$UNIPROT == str_split(overview$accession[i], coll("-"), simplify = T)[, 1], "GO"] %>%
    paste(collapse = " ")
}
  

## characterise outliers
outliers = prediction[prediction$count > 3 & prediction$pred_count < 2, ]
outliers$Accession %>% unique()


### OUTPUT ###




