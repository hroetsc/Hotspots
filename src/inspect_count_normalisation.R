

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
