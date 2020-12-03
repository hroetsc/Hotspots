
library(dplyr)

### INPUT ###

window_size = 25

# human proteome
prots = read.csv("data/proteins_w_hotspots.csv", stringsAsFactors = F)

load("data/IEDB+Sarkizova/windowCounts.RData")
load("data/IEDB+Sarkizova/accU.RData")
names(windowCounts) = accU


### MAIN PART ###
##### compare window counts ##### 

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


# plot proteins
pb = txtProgressBar(min = 0, max = nrow(prots), style = 3)

pdf("data/true-counts_vs_window-counts.pdf", width = 12, height = 8)
par(mfrow = c(2,2))

for (i in 1:nrow(prots)) {
  setTxtProgressBar(pb, i)
  
  cnt.prot = windowCounts[names(windowCounts) == prots$Accession[i]] %>%
    unlist() %>%
    as.numeric()
  
  plot(cnt.prot,
       type = "l",
       main = prots$Accession[i],
       ylab = "count",
       xlab = "position")
  
  mean_counts = c()
  median_counts = c()
  middle_counts = c()
  
  if (length(cnt.prot) >= window_size ){
    pos1 = 1
    pos2 = pos1 + window_size - 1
    
    while (pos2 <= length(cnt.prot)) {
      
      window_pos1 = pos1-floor(window_size/2)
      window_pos2 = pos2+floor(window_size/2)
      
      if (window_pos1 < 1) { window_pos1 = 1 }
      if (window_pos2 > length(cnt.prot)) { window_pos2 = length(cnt.prot) }
      
      mean_counts = c(mean_counts,
                      cnt.prot[window_pos1:window_pos2] %>% mean())
      median_counts = c(median_counts,
                        cnt.prot[window_pos1:window_pos2] %>% median())
      
      # increment positions of sliding window
      pos1 = pos1 + 1
      pos2 = pos2 + 1
    }
    
    points(mean_counts,
           type = "l",
           col = "orange")
    points(median_counts,
           type = "l",
           col = "dodgerblue")
    
  }
  legend("topleft",
         legend = c("true count", "mean window count", "median window count"),
         lty = c(1,1,1,1),
         col = c("black", "orange", "dodgerblue"),
         cex = .6,
         bg = "transparent",
         bty = "n")
  
}

dev.off()



###### different normalisation techniques ######
## use window counts with lonely peptider removed form 1_featureGeneration.R
k = sample(length(windowCounts), 100)

pdf("normalisation_examples.pdf", height = 10)
par(mfrow = c(1,1))
for (i in k){
  cnt = windowCounts[[i]]
  cnt_scale = (cnt - min(cnt)) / (max(cnt) - min(cnt))
  cnt_len = cnt / length(cnt)
  cnt_z = (cnt - mean(cnt)) / sd(cnt)
  
  cnt2 = windowCounts[[i+10]]
  cnt2_scale = (cnt2 - min(cnt2)) / (max(cnt2) - min(cnt2))
  cnt2_len = cnt2 / length(cnt2)
  cnt2_z = (cnt2 - mean(cnt2)) / sd(cnt2)
  
  cnt3 = windowCounts[[i+20]]
  cnt3_scale = (cnt3 - min(cnt3)) / (max(cnt3) - min(cnt3))
  cnt3_len = cnt3 / length(cnt3)
  cnt3_z = (cnt3 - mean(cnt3)) / sd(cnt3)
  
  
  layout(matrix(1:4, 4, 1, byrow = T))
  cnt %>% plot(type = "l",
               main = "no normalisation",
               ylim = c(0, max(cnt, cnt2, cnt3)))
  cnt2 %>% lines(col = "blue")
  cnt3 %>% lines(col = "red")
  
  cnt_len %>% plot(type = "l",
               main = "divided by protein length",
               ylim = c(0, max(cnt_len, cnt2_len, cnt3_len)))
  cnt2_len %>% lines(col = "blue")
  cnt3_len %>% lines(col = "red")
  
  cnt_scale %>% plot(type = "l",
                    main = "scaled between 0 and 1",
                    ylim = c(0, max(cnt_scale, cnt2_scale, cnt3_scale)))
  cnt2_scale %>% lines(col = "blue")
  cnt3_scale %>% lines(col = "red")
  
  cnt_z %>% plot(type = "l",
                 main = "z-transformation",
                 ylim = c(min(cnt_z, cnt2_z, cnt3_z), max(cnt_z, cnt2_z, cnt3_z)))
  cnt2_z %>% lines(col = "blue")
  cnt3_z %>% lines(col = "red")
  
}
dev.off()
