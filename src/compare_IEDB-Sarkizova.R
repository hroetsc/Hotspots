### HEADER ###
# HOTSPOT PREDICTION
# description: compare IEDB and Sarkizova datasets
#               1) overlap of detected proteins
#               2) correlation of counts for overlapping peptides
# input: windowCounts, accU (from Juliane)
# output: correlation among the data, joint data set if comparable
# author: HR

library(dplyr)
library(eulerr)

### INPUT ###
load("data/IEDB/windowCounts.RData")
load("data/IEDB/accU.RData")
iedb.counts = windowCounts
names(iedb.counts) = accU

load("data/Sarkizova/windowCounts.RData")
load("data/Sarkizova/accU.RData")
sarkizova.counts = windowCounts
names(sarkizova.counts) = accU


### MAIN PART ###
### overlap of detected proteins
prot.names = list(IEDB = names(iedb.counts),
                  Sarkizova = names(sarkizova.counts))

eu = euler(prot.names, shape = "ellipse")
png("data/overlap_Sarkizova-IEDB.png")
plot(eu, quantities = T)
dev.off()

overlap = names(sarkizova.counts)[names(sarkizova.counts) %in% names(iedb.counts)]
print(length(overlap))

## correlation of counts
count_cor = data.frame(accession = overlap,
                       pcc = NA)

pb = txtProgressBar(min = 0, max = length(overlap), style = 3)

# profile plots
pdf("data/correlation_Sarkizova-IEDB.pdf", width = 12, height = 8)
par(mfrow = c(2,2))
for (i in 1:length(overlap)) {
  setTxtProgressBar(pb, i)
  
  cnt.iedb = iedb.counts[names(iedb.counts) == overlap[i]] %>% unlist()
  cnt.sarkizova = sarkizova.counts[names(sarkizova.counts) == overlap[i]] %>% unlist()
  
  pc = cor(cnt.iedb, cnt.sarkizova, method = "pearson")
  count_cor$pcc[i] = pc
  
  lower = min(cnt.iedb, cnt.sarkizova)
  upper = max(cnt.iedb, cnt.sarkizova)
  
  plot(cnt.iedb,
       type = "l",
       col = "darkgreen",
       xlab = "position",
       ylab = "count",
       ylim = c(lower, upper),
       main = overlap[i],
       sub = paste0("PCC: ", round(pc, 4)),
       axes = F)
  points(cnt.sarkizova,
         type = "l",
         col = "firebrick")
  legend("topright",
         legend = c("IEDB", "Sarkizova"),
         lty = c(1,1),
         col = c("darkgreen", "firebrick"),
         cex = .6)
  axis(1)
  axis(2)
}
dev.off()

# correlation coefficients
summary(count_cor$pcc)
q = quantile(count_cor$pcc, c(0.1, 0.5, 0.9))

png("data/correlation-coeff_Sarkizova-IEDB.png")
hist(count_cor$pcc,
     main = "correlation between Sarkoviza and IEDB counts",
     breaks = 75,
     xlab = "PCC")
abline(v = q[1], col = "red")
abline(v = q[2], col = "red")
abline(v = q[3], col = "red")
dev.off()

# scatter plots
pdf("data/correlation_Sarkizova-IEDB_scatter.pdf", width = 12, height = 8)
par(mfrow = c(2,2))
for (i in 1:length(overlap)) {
  setTxtProgressBar(pb, i)
  
  cnt.iedb = iedb.counts[names(iedb.counts) == overlap[i]] %>% unlist()
  cnt.sarkizova = sarkizova.counts[names(sarkizova.counts) == overlap[i]] %>% unlist()
  
  pc = cor(cnt.iedb, cnt.sarkizova, method = "pearson")
  lower = min(cnt.iedb, cnt.sarkizova)
  upper = max(cnt.iedb, cnt.sarkizova)
  
  plot(cnt.sarkizova ~ cnt.iedb,
       xlab = "IEDB",
       ylab = "Sarkizova",
       pch = 16,
       ylim = c(lower, upper),
       xlim = c(lower, upper),
       main = overlap[i],
       sub = paste0("PCC: ", round(pc, 4)),
       axes = F)
  abline(a = 0, b = 1)
  axis(1)
  axis(2)
}
dev.off()

### OUTPUT ###


