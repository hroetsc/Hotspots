### HEADER ###
# HOTSPOT PREDICTION
# description: investigate the imapct of different data characteristics on the hotspot count
# input: mhc ligand table, badly predicted high-count proteins
# output: -
# author: HR


library(ggplot2)
library(stringr)
library(dplyr)

JOBID = "5439350-1"

### INPUT ###
source_antigens = read.csv("HOTSPOTS/DATA/mhc_ligand_table_export_1588248580.csv",
                           stringsAsFactors = F, header = 2)

high.counts.bad = read.csv("results/badly_predicted_high_counts.csv",
                           stringsAsFactors = F)


### MAIN PART ###

high.counts.acc =  str_split_fixed(high.counts.bad$Accession, coll("-"), Inf)[,1]
k = source_antigens$Epitope.8 %in% high.counts.acc


############ general features ############
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


############ overstudying ############
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

