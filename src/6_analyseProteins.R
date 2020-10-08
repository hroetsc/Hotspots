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


### INPUT ###
prediction = read.csv("results/last_model_prediction.csv", stringsAsFactors = F)
hp = read.csv("../Seq2Vec/files/proteome_human.csv", stringsAsFactors = F)
props = read.csv("../Seq2Vec/RUNS/HumanProteome/PropMatrix_seqs.csv", stringsAsFactors = F)
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
write.csv(overview, "results/performance_protein_wise.csv", row.names = F)

# histogram of protein-level PCC's
q = quantile(overview$PCC, c(0.05, 0.5, 0.95))
m = overview$PCC %>% mean()
hist(overview$PCC, freq = F, breaks = 100,
     main = "histogram of PCCs of all test proteins",
     xlab = 'PCC')
abline(v = q[1], col = "red")
abline(v = q[2], col = "red")
abline(v = m, col = "blue")
abline(v = q[3], col = "red")


# plot
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


## biological functions
hsGO = godata('org.Hs.eg.db', keytype = "UNIPROT", ont = "BP", computeIC = F)
annot = hsGO@geneAnno

for (i in 1:nrow(overview)) {
  overview$GO_terms[i] = annot[annot$UNIPROT == str_split(overview$accession[i], coll("-"), simplify = T)[, 1], "GO"] %>%
    paste(collapse = " ")
}
  

########## characterise outliers ##########
high.counts = overview[overview$average_count > 1 & overview$PCC < .5, ]
write.csv(high.counts, "results/outliers.csv", row.names = F)

high.counts.well = overview[overview$average_count > 1 & overview$PCC > .8, ]
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


##########  source antigens ########## 
source_antigens = read.csv("HOTSPOTS/DATA/mhc_ligand_table_export_1588248580.csv",
                           stringsAsFactors = F, header = 2)

high.counts.acc =  str_split_fixed(high.counts$Accession, coll("-"), Inf)[,1]
k = source_antigens$Epitope.8 %in% high.counts.acc

# cell types
cells = source_antigens$Antigen.Processing.Cells[k] %>% table()
View(cells)

# epitopes
epitope_names = source_antigens$Epitope.6[k] %>% unique()
View(epitope_names)

# haplotype of cells
cell_haplotype = source_antigens$MHC[k] %>% unique()
View(cell_haplotype)
cell_haplotype.tbl = source_antigens$MHC[k] %>% table()
all_cell_haplotype.tbl = source_antigens$MHC %>% table()

# most likely MHC molecule
epitope_haplotype = str_split_fixed(source_antigens$Host.3[k], coll(";"), Inf)[, 1]
View(epitope_haplotype %>% unique())
epitope_haplotype.tbl = epitope_haplotype %>% table()
all_epitope_haplotype.tbl = str_split_fixed(source_antigens$Host.3, coll(";"), Inf)[, 1] %>%
  table()


# disease
disease = source_antigens$In.vivo.Process.1[k] %>% unique()
View(disease)
disease.tbl = source_antigens$In.vivo.Process.1[k] %>% table()
all_disease.tbl = source_antigens$In.vivo.Process.1 %>% table()


## overstudying
sampling_frequency = data.frame(Accession = str_split_fixed(prots, coll("-"), Inf)[, 1],
                                no_epitopes = rep(NA, length(prots)),
                                no_studies = rep(NA, length(prots)),
                                protein_length = rep(NA, length(prots)))

pb = txtProgressBar(min = 0, max = length(prots), style = 3)
for (i in 1:length(prots)) {
  setTxtProgressBar(pb, i)
  cnt = source_antigens[source_antigens$Epitope.8 == prots[i], ]
  
  if(nrow(cnt) > 0) {
    cnt = cnt[-which(duplicated(cnt$Epitope.2)), ]
    sampling_frequency$no_epitopes[i] = nrow(cnt)
    sampling_frequency$no_studies[i] = cnt$Reference.4 %>% unique() %>% length()
    sampling_frequency$protein_length[i] = hp$seqs[hp$Accession == prots[i]] %>% nchar()
  }
}

sampling_frequency = na.omit(sampling_frequency)
sampling_frequency = sampling_frequency[-which(sampling_frequency$no_epitopes == 0), ]

l = which(sampling_frequency$Accession %in% high.counts.acc)
summary(sampling_frequency$no_studies)
summary(sampling_frequency$no_studies[l])
t.test(sampling_frequency$no_studies, sampling_frequency$no_studies[l])

plot(density(sampling_frequency$no_studies))
lines(density(sampling_frequency$no_studies[l]), col = "red")

prediction.studies = prediction
prediction.studies$Accession = str_split_fixed(prediction.studies$Accession, coll("-"), Inf)[, 1]
prediction.studies = left_join(prediction, sampling_frequency) %>% na.omit()

overview.tmp = overview
overview.tmp$Accession = str_split_fixed(overview.tmp$Accession, coll("-"), Inf)[, 1]
overview.studies = left_join(overview.tmp, sampling_frequency)

cor(overview.studies$average_count, overview.studies$no_studies)
ggplot(overview.studies, aes(x = no_studies, y = average_count)) +
  geom_point(alpha = 0.2, size = 0.5) +
  ggtitle("true counts vs no of studies") +
  theme_bw()
ggsave(paste0("results/plots/", JOBID, "_trueVSstudies.png"), plot = last_plot(),
       device = "png", dpi = "retina")


# linking number of epitopes to number of references
pc = cor(sampling_frequency$no_epitopes, sampling_frequency$no_studies, method = "pearson") %>%
  round(4)
plot(sampling_frequency$no_epitopes ~ sampling_frequency$no_studies,
     main = "number of references vs number of unique epitopes per protein",
     sub = paste0("PCC = ", pc),
     xlab = "number of references",
     ylab = "number of unique epitopes")

studies.per.epitope = (sampling_frequency$no_studies / sampling_frequency$no_epitopes) * sampling_frequency$protein_length
hist(studies.per.epitope, breaks = 20,
     main = "average number of references per unique epitope per residue",
     xlab = "(number of references / number of epitopes) * protein length")

# excluding over- and understudied proteins
q = quantile(studies.per.epitope, c(0.05, 0.5, 0.95), na.rm = T)
overstudied = sampling_frequency$Accession[studies.per.epitope >= q[3]] %>%
  as.character() %>%
  na.omit()

understudied = sampling_frequency$Accession[studies.per.epitope <= q[1]] %>%
  as.character() %>%
  na.omit()

acc = str_split_fixed(prediction$Accession, coll("-"), Inf)[, 1]
over = which(acc %in% overstudied)
under = which(acc %in% understudied)

prediction.excl = prediction
prediction.excl$color = "normal"
prediction.excl$color[over] = "overstudied"
prediction.excl$color[under] = "understudied"


ggplot(prediction.excl, aes(x = count, y = pred_count, col = color)) +
  geom_point(alpha = 0.05, size = 0.1) +
  xlim(c(start, stop)) +
  ylim(c(start, stop)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  coord_equal() +
  theme_bw()
ggsave(paste0("results/plots/", JOBID, "_trueVSpredicted-scatter-over-under.png"), plot = last_plot(),
       device = "png", dpi = "retina")

# exclude for correlation
cor(prediction.excl$count, prediction.excl$pred_count, method = "pearson")
cor(prediction.excl$count[-c(over, under)],
    prediction.excl$pred_count[-c(over, under)],
    method = "pearson")
cor(prediction.excl$count[-c(over)],
    prediction.excl$pred_count[-c(over)],
    method = "pearson")
cor(prediction.excl$count[-c(under)],
    prediction.excl$pred_count[-c(under)],
    method = "pearson")

