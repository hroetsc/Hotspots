### HEADER ###
# HOTSPOT PREDICTION: VALIDATION
# description: identify homologous proteins between different coronaviruses
# input: aggregated protein sequences
# output: homologous proteins
# author: HR

library(dplyr)
library(protr)
library(lattice)
library(Rfast)


### INPUT ###
prots = read.csv("validation/data/aggregated_sequences.csv", stringsAsFactors = F)

viruses = c("SARS-CoV2",
            "HCoV-OC43",
            "HCoV-HKU1",
            "HCoV-NL63",
            "HCoV-229E")

### MAIN PART ###
# pairwise sequence similarity based on local alignment
seq_sim = parSeqSim(split(prots$sequence, seq(nrow(prots))),
                    cores = 16,
                    type = "local",
                    submat = "BLOSUM62")

rownames(seq_sim) = prots$annotation
colnames(seq_sim) = rownames(seq_sim)

sars_idx = which(prots$organism == "SARS-CoV2")
homologs = list()
for (r in sars_idx){
  cnt = seq_sim[r, ]
  
  # remove proteins from same organism
  cnt.org = prots$organism[prots$annotation == rownames(seq_sim)[r]]
  cnt = cnt[-which(names(cnt) %in% prots$annotation[prots$organism == cnt.org])]
  
  homologs[[r]] = c(names(cnt)[cnt == max(cnt)],
                    names(cnt)[cnt == nth(cnt, 2, descending = T)],
                    names(cnt)[cnt == nth(cnt, 3, descending = T)],
                    names(cnt)[cnt == nth(cnt, 4, descending = T)]) %>% unique()
  
  names(homologs)[r] = prots$annotation[r]
}

# go through results by hand
poi_SARS.CoV2 = c("ORF1ab polyprotein [Severe acute respiratory syndrome coronavirus 2]",
                  "surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]",
                  "ORF3a protein [Severe acute respiratory syndrome coronavirus 2]",
                  "envelope protein [Severe acute respiratory syndrome coronavirus 2]",
                  "membrane glycoprotein [Severe acute respiratory syndrome coronavirus 2]",
                  "nucleocapsid phosphoprotein [Severe acute respiratory syndrome coronavirus 2]")
homologs.short = homologs[names(homologs) %in% poi_SARS.CoV2]

poi_HCoV.OC43 = c("ORF1ab polyprotein [Human coronavirus OC43]",
                  "spike surface glycoprotein [Human coronavirus OC43]",
                  "nsp11 [Human coronavirus OC43]",
                  "envelope protein [Human coronavirus OC43]",
                  "membrane protein [Human coronavirus OC43]",
                  "nucleocapsid protein [Human coronavirus OC43]")

poi_HCoV.HKU1 = c("orf1ab polyprotein [Human coronavirus HKU1]",
                  "spike glycoprotein [Human coronavirus HKU1]",
                  "envelope protein [Human coronavirus HKU1]",
                  "envelope protein [Human coronavirus HKU1]",
                  "membrane glycoprotein [Human coronavirus HKU1]",
                  "nucleocapsid phosphoprotein [Human coronavirus HKU1]")

poi_HCoV.NL63 = c("replicase polyprotein 1ab [Human coronavirus NL63]",
                  "Spike protein [Human coronavirus NL63]",
                  "Envelope protein [Human coronavirus NL63]",
                  "Envelope protein [Human coronavirus NL63]",
                  "Membrane protein [Human coronavirus NL63]",
                  "Nucleocapsid protein [Human coronavirus NL63]")

poi_HCoV.229E = c("replicase polyprotein 1ab [Human coronavirus 229E]",
                  "surface glycoprotein [Human coronavirus 229E]",
                  "envelope protein [Human coronavirus 229E]",
                  "envelope protein [Human coronavirus 229E]",
                  "membrane protein [Human coronavirus 229E]",
                  "nucleocapsid protein [Human coronavirus 229E]")


# combine in single df
selected_homologs = data.frame(SARS_CoV2 = poi_SARS.CoV2,
                               HCoV_OC43 = poi_HCoV.OC43,
                               HCoV_HKU1 = poi_HCoV.HKU1,
                               HCoV_NL63 = poi_HCoV.NL63,
                               HCoV_229E = poi_HCoV.229E)

### OUTPUT ###
write.csv(selected_homologs, "validation/results/homologs.csv", row.names = F)

