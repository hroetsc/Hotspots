### HEADER ###
# HOTSPOT PREDICTION
# description: analyse window similarity vs count similarity
# input: windowTokens
# output: -
# author: HR


library(data.table)
library(dplyr)
library(Biostrings)
library(protr)
library(ggplot2)

k = 1e04

### INPUT ###
windowTokens = fread("data/windowTokens25aa_mean.csv") %>%
  as.data.frame()


### MAIN PART ###
# sample k window pairs
sampled_windows = sample(nrow(windowTokens), 2*k)
group1 = sampled_windows[1:k]
group2 = sampled_windows[(k+1):(2*k)]

# fetch count difference and sequence similarity
sim = data.frame(ind1 = group1,
                 ind2 = group2,
                 sequence_similarity = NA,
                 count_difference = NA)

pb = txtProgressBar(min = 0, max = k, style = 3)
for (i in 1:k) {
  setTxtProgressBar(pb, i)
  
  sim$count_difference[i] = abs(windowTokens$counts[sim$ind1[i]] - windowTokens$counts[sim$ind2[i]])
  # sim$sequence_similarity[i] = Biostrings::pairwiseAlignment(windowTokens$window[sim$ind1[i]],
  #                                                            windowTokens$window[sim$ind2[i]],
  #                                                            type = "global",
  #                                                            scoreOnly = T)
  sim$sequence_similarity[i] = parSeqSim(list(s1 = windowTokens$window[sim$ind1[i]],
                                              s2 = windowTokens$window[sim$ind2[i]]),
                                         type = "global",
                                         cores = 8)[1,2]
}

summary(sim$count_difference) 
summary(sim$sequence_similarity)

# plot
x11()
ggplot(data = sim, aes(x = sequence_similarity, y = log2(count_difference))) +
  geom_point(size = 1) + 
  theme_classic()

ggsave(filename = "results/exploratory/count_diff-vs-window_sim.png", plot = last_plot(),
       device = "png", dpi = "retina")

# count difference is independent of sequence similarity




