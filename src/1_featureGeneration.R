### HEADER ###
# HOTSPOT PREDICTION
# description: generate feature set by applying a sliding window over a protein sequence
# input: windowCounts (from Juliane)
# output: feature set (counts and tokens for every sequence)
# author: HR

library(plyr)
library(dplyr)
library(stringr)

library(foreach)
library(doParallel)
library(future)
registerDoParallel(availableCores())


window_size = 25
# ext = 25  # extension, if used


### INPUT ###
load("data/IEDB+Sarkizova/accU.RData")
load("data/IEDB+Sarkizova/windowCounts.RData")

# human proteome
prots = read.csv("/media/hanna/Hanna2/DATA/GENERAL/proteome_human.csv",
                 stringsAsFactors = F, header = T)


### MAIN PART ###
names(windowCounts) = accU

# remove proteins with "lonely peptides"
# lonely peptide --> several short disconnected non-overlapping regions with only one count per aa
# conditions: max count must be 1
rm = c()
for (w in 1:length(windowCounts)){
  cnt = windowCounts[[w]]
  if (max(cnt) == 1){
    rm = c(rm, w)
  }

}

length(rm)
windowCounts = windowCounts[-rm]


# keep only proteins with hotspots
prots = prots[which(prots$Accession %in% names(windowCounts)), ]
# remove proteins with X
X_idx = str_detect(prots$seqs, "X")
X_idx[X_idx] %>% length()
prots = prots[X_idx == F, ]

# replace U by C
U_idx = str_detect(prots$seqs, "U")
U_idx[U_idx] %>% length()
prots_U = prots[U_idx, ]

for (i in 1:nrow(prots_U)) {
  prots_U$seqs[i] = str_replace_all(prots_U$seqs[i], pattern = "U", replacement = "C") 
}

prots[U_idx, ] = prots_U


# get windows and average counts per window
get_windows_counts = function(extension = "", outfile = ""){
  
  windowTokens = list()
  
  windowTokens = foreach (i = 1:length(windowCounts)) %dopar% {
    
    cnt.Prot = prots[prots$Accession == names(windowCounts)[i], ][1, ]
    
    cnt.AA = strsplit(cnt.Prot$seqs, coll("")) %>% unlist() %>% as.data.frame()
    names(cnt.AA) = "aa"
    cnt.AA$aa = as.character(cnt.AA$aa)
    cnt.AA$idx = seq(1, nrow(cnt.AA)) %>% as.numeric()
    
    # apply sliding window and check which amino acids are in current window
    # store all window amino acids   for current protein
    if (nrow(cnt.AA) > window_size ){
      wnd.Tokens = data.frame(Accession = rep(cnt.Prot$Accession, nrow(cnt.AA) - window_size + 1),
                              window = rep(NA, nrow(cnt.AA) - window_size + 1),
                              counts = rep(NA, nrow(cnt.AA) - window_size + 1),
                              dist_N = seq(0, nrow(cnt.AA)-window_size),
                              dist_C = seq(nrow(cnt.AA)-window_size, 0),
                              len_protein = rep(nchar(cnt.Prot$seqs), nrow(cnt.AA) - window_size + 1))
      
      pos1 = 1
      pos2 = pos1 + window_size - 1
      
      while (pos2 <= nrow(cnt.AA)) {
        
        # check which tokens are in current sliding window
        cnt.Wnd = cnt.AA$aa[pos1:pos2] %>% paste(collapse = "")
        wnd.Tokens$window[pos1] = cnt.Wnd
        
        
        # different extensions --> refine window that is used to get counts
        if (extension == "none") {
          Wnd.for.counts = cnt.AA$aa[pos1:pos2] %>% paste(collapse = "")
        } else if (extension == "N") {
          Wnd.for.counts = cnt.AA$aa[(pos1+ext):pos2] %>% paste(collapse = "")
        } else if (extension == "C") {
          Wnd.for.counts = cnt.AA$aa[pos1:(pos2-ext)] %>% paste(collapse = "")
        } else if (extension == "NandC") {
          Wnd.for.counts = cnt.AA$aa[(pos1+ext):(pos2-ext)] %>% paste(collapse = "")
        }
        
        
        # get mean counts of current window
        loc = str_locate(cnt.Prot$seqs,
                         Wnd.for.counts) %>% as.numeric()
        
        # centered window
        window_pos1 = loc[1]-floor(window_size/2)
        window_pos2 = loc[2]+floor(window_size/2)
        
        if (window_pos1 < 1) { window_pos1 = 1 }
        if (window_pos2 > nrow(cnt.AA)) { window_pos2 = nrow(cnt.AA) }
        
        wnd.Tokens$counts[pos1] = windowCounts[[i]][window_pos1:window_pos2] %>% mean()
        
        # increment positions of sliding window
        pos1 = pos1 + 1
        pos2 = pos2 + 1
        
      }
      
      windowTokens[[i]] = wnd.Tokens
    }
    
  }
  
  
  ### OUTPUT ###
  save(windowTokens, file = paste0(outfile, ".RData"))
  save(windowCounts, file = "data/windowCounts.RData")
  
  # reduce feature space
  # average counts of neighbouring identical token sets - should not happen anymore
  # plus:
  # flatten and save list in chunks
  # load(outfile)
  
  out = paste0(outfile, ".csv")
  if (file.exists(out)){ system(paste0("rm ", out)) }
  
  
  pb = txtProgressBar(min = 0, max = length(windowTokens), style = 3)
  
  for ( i in 1:length(windowTokens) ) {
    
    setTxtProgressBar(pb, i)
    
    wT = windowTokens[[i]] %>% as.data.frame()
    
    if(file.exists(out)) {
      write.table(wT, out, sep = ",", row.names = F, append = T, col.names = F)
      
    } else {
      write.table(wT, out, sep = ",", row.names = F, append = F, col.names = T)
      
    }
  }
  
  
}


# apply
get_windows_counts(extension = "none",
                   outfile = "data/windowTokens25aa_mean")

# save proteins with hotspots
write.csv(prots, "data/proteins_w_hotspots.csv", row.names = F)


