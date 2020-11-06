### HEADER ###
# HOTSPOT PREDICTION: VALIDATION
# description: format proteins, generate sliding windows, retrieve distance to N and C term
# input: fasta file with protein sequences and accessions
# output: preprocessed features
# author: HR

library(seqinr)
library(plyr)
library(dplyr)
library(stringr)

window_size = 50
spec = "_50aa_corona"


### INPUT ###
viruses = c("SARS-CoV2",
            "HCoV-OC43",
            "HCoV-HKU1",
            "HCoV-NL63",
            "HCoV-229E")

virus_proteomes = list.files(path = "validation/data", pattern = "sequence_", full.names = T)


### MAIN PART ###
# open fasta files, convert into data frame, concatenate
prots = data.frame(organism = NA,
                   accession = NA,
                   annotation = NA,
                   sequence = NA)

for (v in viruses) {
  cnt.proteome = read.fasta(virus_proteomes[str_detect(virus_proteomes, v)],
                            seqtype = "AA")
  
  cnt.df = data.frame(organism = rep(v, length(cnt.proteome)),
                      accession = names(cnt.proteome),
                      annotation = sapply(cnt.proteome, getAnnot),
                      sequence = sapply(cnt.proteome, paste, collapse=""))
  cnt.df$annotation = str_remove_all(cnt.df$annotation, paste0(">", cnt.df$accession, " "))
  
  prots = rbind(prots, cnt.df)
  
  paste0(v, ", number of proteins: ", nrow(cnt.df)) %>% print()
}

prots = na.omit(prots)


# check for unconventional aas
any(str_detect(prots$sequence, "U"))
any(str_detect(prots$sequence, "X"))


# sliding window
windows = list()

pb = txtProgressBar(min = 0, max = nrow(prots), style = 3)
for (i in 1:nrow(prots)) {
  setTxtProgressBar(pb, i)
  
  # split current protein sequence into single amino acids
  cnt.AA = strsplit(prots$sequence[i], coll("")) %>% unlist() %>% as.data.frame()
  names(cnt.AA) = "aa"
  cnt.AA$aa = as.character(cnt.AA$aa)
  cnt.AA$idx = seq(1, nrow(cnt.AA)) %>% as.numeric()
  
  # apply sliding window and check which amino acids are in current window
  # store all window amino acids   for current protein
  # make sure the protein contains more that window size amino acids
  if (nrow(cnt.AA) > window_size ){
    
    windows.cnt.prot = data.frame(organism = rep(prots$organism[i], nrow(cnt.AA) - window_size + 1),
                                  accession = rep(prots$accession[i], nrow(cnt.AA) - window_size + 1),
                                  annotation = rep(prots$annotation[i], nrow(cnt.AA) - window_size + 1),
                                  window = rep(NA, nrow(cnt.AA) - window_size + 1),
                                  dist_N = seq(0, nrow(cnt.AA)-window_size),
                                  dist_C = seq(nrow(cnt.AA)-window_size, 0),
                                  len_protein = rep(nchar(prots$sequence[i]), nrow(cnt.AA) - window_size + 1))
    
    pos1 = 1
    pos2 = pos1 + window_size - 1
   
    while (pos2 <= nrow(cnt.AA)) {
      
      # check which tokens are in current sliding window
      cnt.Wnd = cnt.AA$aa[pos1:pos2] %>% paste(collapse = "")
      windows.cnt.prot$window[pos1] = cnt.Wnd
      
      # increment positions of sliding window
      pos1 = pos1 + 1
      pos2 = pos2 + 1
      
    }
    
    windows[[i]] = windows.cnt.prot
  }
  
}

names(windows) = prots$accession

# convert list into df
windows_df = plyr::ldply(windows) %>% na.omit()
windows_df$.id = NULL


### OUTPUT ###
write.csv(prots, "validation/data/aggregated_sequences.csv", row.names = F)

save(windows, file = paste0("validation/data/windows", spec, ".RData"))
write.csv(windows_df, paste0("validation/data/windows", spec, ".csv"), row.names = F)

