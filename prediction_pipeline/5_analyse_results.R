### HEADER ###
# HOTSPOT PREDICTION: VALIDATION
# description: analyse prediction results to identify antigens
# input: predictions, aggregated sequences, homologous proteins across viruses
# output: 10-20 antigens of interest (20-40 aa length)
# author: HR, JL

library(dplyr)
library(stringr)
library(zoo)
library(Biostrings)
library(msa)
library(mixtools)

spec = "_50aa_corona2"
window_size = 50


### INPUT ###
prots = read.csv("validation/data/aggregated_sequences.csv", stringsAsFactors = F)
homologs = read.csv("validation/results/homologs.csv", stringsAsFactors = F)

# manually identified epitopes
literature_epitopes = read.csv("validation/data/epitopes_literature.csv", stringsAsFactors = F)
# IEDB epitopes
IEDB_epitopes = read.csv("validation/data/epitope_table_export_1604341084.csv", stringsAsFactors = F)

system("scp -rp hroetsc@transfer.gwdg.de:/usr/users/hroetsc/Hotspots/results/*_model_prediction_*aa_*_rank*.csv validation/results/predictions/")
preds = list.files(path = "validation/results/predictions",
                   pattern = paste0("model_prediction", spec), full.names = T)

viruses = c("SARS-CoV2",
            "HCoV-OC43",
            "HCoV-HKU1",
            "HCoV-NL63",
            "HCoV-229E")

### MAIN PART ###
########## retrieve counts ##########
# aggregate predictions from different ranks, back-transform counts into non-log scale
aggregation = function(model){
  print(model)
  
  cnt.preds = preds[str_detect(preds, model)]
  prediction = read.csv(cnt.preds[1], stringsAsFactors = F)
  
  pred_counts = rep(0, nrow(prediction))
  
  pb = txtProgressBar(min = 0 , max = length(cnt.preds), style = 3)
  for (p in 1:length(cnt.preds)) {
    setTxtProgressBar(pb, p)
    
    cnt_pred = read.csv(cnt.preds[p], stringsAsFactors = F)
    pred_counts = pred_counts + (cnt_pred$pred_count * 1/length(cnt.preds))
    
  }
  
  prediction$pred_count = pred_counts
  prediction$pred_count = 2^(prediction$pred_count) - 1
  
  return(prediction)
}

best.m_prediction = aggregation(model = "best_model")
last.m_prediction = aggregation(model = "last_model")

max_count = max(best.m_prediction$pred_count, last.m_prediction$pred_count) %>% ceiling()
prots = prots[prots$accession %in% unique(best.m_prediction$accession), ] %>% unique()

########## plot profiles ##########
# return list with exact and rolling mean counts
plot_profiles = function(prediction, model) {
  
  cnt.prots = prots[prots$accession %in% unique(prediction$accession), ]
  
  # get counts for each position of the proteins
  counts = list()
  counts.rollmean = list()
  
  for (i in 1:nrow(cnt.prots)) {
    cnt.Data = prediction[prediction$accession == cnt.prots$accession[i], ]
    cnt.counts = rep(NA, nchar(cnt.prots$sequence[i]))
    
    for (j in 1:nrow(cnt.Data)) {
      idx = str_locate(cnt.prots$sequence[i], cnt.Data$window[j]) %>% as.numeric()
      cnt.counts[idx[1] : idx[2]] = rep(cnt.Data$pred_count[j], window_size)
    }
    
    counts[[i]] = cnt.counts
    counts.rollmean[[i]] = rollmean(cnt.counts, k = 9)
    
    names(counts)[i] = cnt.prots$accession[i]
    names(counts.rollmean)[i] = cnt.prots$accession[i]
    
  }
  
  # plot
  pdf(paste0("validation/results/profiles_", model, ".pdf"),
      width = 12, height = 8)
  par(mfrow = c(2,2))
  
  for (i in 1:nrow(cnt.prots)){
    y = counts[[i]]
    y.mean = counts.rollmean[[i]]
    x = seq(length(y))
    
    plot(x, y,
         type = "l",
         col = "darkblue",
         ylim = c(0, max_count),
         axes = F,
         ylab = "counts",
         xlab = "position",
         main = cnt.prots$annotation[i],
         sub = paste0("accession: ", prots$accession[i]))
    
    points(seq(1, length(y.mean), ceiling(length(y.mean)/length(y))),
           y.mean,
           type = "l",
           col="black")
    
    abline(h = 1, col = "darkgreen")
    
    axis(1)
    axis(2)
  }
  dev.off()
  
  # concatenate counts and return
  aggregated_counts = list(counts = counts,
                           counts.rollmean = counts.rollmean)
  return(aggregated_counts)
}

best.m_counts = plot_profiles(prediction = best.m_prediction, model = "best_model")
last.m_counts = plot_profiles(prediction = last.m_prediction, model = "last_model")


get_counts = function(prediction){
  cnt.prots = prots[prots$accession %in% unique(prediction$accession), ]
  
  # get counts for each position of the proteins
  counts = list()
  counts.rollmean = list()
  
  for (i in 1:nrow(cnt.prots)) {
    cnt.Data = prediction[prediction$accession == cnt.prots$accession[i], ]
    cnt.counts = rep(NA, nchar(cnt.prots$sequence[i]))
    
    for (j in 1:nrow(cnt.Data)) {
      idx = str_locate(cnt.prots$sequence[i], cnt.Data$window[j]) %>% as.numeric()
      cnt.counts[idx[1] : idx[2]] = rep(cnt.Data$pred_count[j], window_size)
    }
    
    counts[[i]] = cnt.counts
    counts.rollmean[[i]] = rollmean(cnt.counts, k = 9)
    
    names(counts)[i] = cnt.prots$accession[i]
    names(counts.rollmean)[i] = cnt.prots$accession[i]
    
  }
  
  # concatenate counts and return
  aggregated_counts = list(counts = counts,
                           counts.rollmean = counts.rollmean)
  return(aggregated_counts)
}

# in case profiles do not need to be replotted
best.m_counts = get_counts(prediction = best.m_prediction)
last.m_counts = get_counts(prediction = last.m_prediction)


########## compare last and best model ##########
# include antigens from literature
best_and_last = function() {
  
  cnt.prots = prots[prots$accession %in% unique(best.m_prediction$accession), ]
  
  pdf(paste0("validation/results/profiles_rollmean.pdf"),
      width = 12, height = 8)
  par(mfrow = c(2,2))
  
  for (i in 1:nrow(cnt.prots)){
    y = best.m_counts[["counts"]][[i]]
    
    y.mean_best = best.m_counts[["counts.rollmean"]][[i]]
    y.mean_last = last.m_counts[["counts.rollmean"]][[i]]
    
    
    plot(seq(1, length(y.mean_best), ceiling(length(y.mean_best)/length(y))),
         y.mean_best,
         type = "l",
         col = "black",
         ylim = c(0, max_count),
         axes = F,
         ylab = "counts",
         xlab = "position",
         main = cnt.prots$annotation[i],
         sub = paste0("accession: ", prots$accession[i]))
    
    points(seq(1, length(y.mean_last), ceiling(length(y.mean_last)/length(y))),
           y.mean_last,
           type = "l",
           lty = "dashed",
           col="black")
    
    abline(h = 1, col = "darkgreen")
    
    legend("topright",
           legend = c("best model", "last model"),
           lty = c("solid", "dashed"),
           cex = .8)
    
    axis(1)
    axis(2)
  }
  
  dev.off()
}

best_and_last()


########## compare homologs ##########
# plot titles
homologs$trivial = c("ORF1ab",
                     "surface glycoprotein",
                     "ORF3a",
                     "envelope protein",
                     "membrane glycoprotein",
                     "nucleocapsid phosphoprotein")

get_x_and_y = function(counts.ls, id, return_mean = T){
  y = counts.ls[["counts"]][names(counts.ls[["counts"]]) == id] %>% unlist()
  y.mean = counts.ls[["counts.rollmean"]][names(counts.ls[["counts.rollmean"]]) == id] %>% unlist()
  x = seq(1, length(y.mean), ceiling(length(y.mean)/length(y)))
  
  if (return_mean){
    return(list(x = x,
                y = y.mean))
  } else {
    return(list(x = seq(length(y)),
                y = y))
  }
  
  
}

plot_homologs = function(counts.ls, model){
  
  pdf(paste0("validation/results/homologs_", model, ".pdf"),
      width = 12, height = 8)
  for (i in 1:nrow(homologs)) {
    
    sars.cov2 = prots$accession[prots$annotation == homologs$SARS_CoV2[i]]
    hcov.oc43 = prots$accession[prots$annotation == homologs$HCoV_OC43[i]]
    hcov.hku1 = prots$accession[prots$annotation == homologs$HCoV_HKU1[i]]
    hcov.nl63 = prots$accession[prots$annotation == homologs$HCoV_NL63[i]]
    hcov.229e = prots$accession[prots$annotation == homologs$HCoV_229E[i]]
    
    
    plot(x = get_x_and_y(counts.ls, sars.cov2)[['x']],
         y = get_x_and_y(counts.ls, sars.cov2)[['y']],
         type = "l",
         col = "firebrick",
         ylim = c(0, max_count),
         axes = F,
         ylab = "counts",
         xlab = "position",
         main = homologs$trivial[i])
    
    points(x = get_x_and_y(counts.ls, hcov.oc43)[['x']],
           y = get_x_and_y(counts.ls, hcov.oc43)[['y']],
           type = "l",
           col = "darkslategray")
    
    points(x = get_x_and_y(counts.ls, hcov.hku1)[['x']],
           y = get_x_and_y(counts.ls, hcov.hku1)[['y']],
           type = "l",
           col = "dodgerblue")
    
    points(x = get_x_and_y(counts.ls, hcov.nl63)[['x']],
           y = get_x_and_y(counts.ls, hcov.nl63)[['y']],
           type = "l",
           col = "darkorchid")
    
    points(x = get_x_and_y(counts.ls, hcov.229e)[['x']],
           y = get_x_and_y(counts.ls, hcov.229e)[['y']],
           type = "l",
           col = "aquamarine")
    
    abline(h = 1, col = "darkgreen")
    
    legend("topright",
           legend = viruses,
           lty = rep(1, 5),
           col = c("firebrick", "darkslategray", "dodgerblue", "darkorchid", "aquamarine"),
           cex = 0.8)
    
    axis(1)
    axis(2)
  }
  dev.off()
}

plot_homologs(best.m_counts, "best_model")
plot_homologs(last.m_counts, "last_model")


plot_alignments = function(counts.ls, model) {
  
  pdf(paste0("validation/results/homologs_msa_", model, ".pdf"),
      width = 36, height = 8)
  for (i in 1:nrow(homologs)) {
    sars.cov2.counts = get_x_and_y(counts.ls, prots$accession[prots$annotation == homologs$SARS_CoV2[i]],
                                   return_mean = F)
    hcov.oc43.counts = get_x_and_y(counts.ls, prots$accession[prots$annotation == homologs$HCoV_OC43[i]],
                                   return_mean = F)
    hcov.hku1.counts = get_x_and_y(counts.ls, prots$accession[prots$annotation == homologs$HCoV_HKU1[i]],
                                   return_mean = F)
    hcov.nl63.counts = get_x_and_y(counts.ls, prots$accession[prots$annotation == homologs$HCoV_NL63[i]],
                                   return_mean = F)
    hcov.229e.counts = get_x_and_y(counts.ls, prots$accession[prots$annotation == homologs$HCoV_229E[i]],
                                   return_mean = F)
    
    ls = list(prots$sequence[prots$annotation == homologs$SARS_CoV2[i]] %>% AAString(),
              prots$sequence[prots$annotation == homologs$HCoV_OC43[i]] %>% AAString(),
              prots$sequence[prots$annotation == homologs$HCoV_HKU1[i]] %>% AAString(),
              prots$sequence[prots$annotation == homologs$HCoV_NL63[i]] %>% AAString(),
              prots$sequence[prots$annotation == homologs$HCoV_229E[i]] %>% AAString())
    
    stringset = AAStringSet(ls)
    alignment = msa(stringset) %>% msaConvert("seqinr::alignment")
    
    sars.cov2.df = alignment$seq[1] %>% strsplit("") %>% unlist() %>% as.character %>% as.data.frame()
    hcov.oc43.df = alignment$seq[2] %>% strsplit("") %>% unlist() %>% as.character %>% as.data.frame()
    hcov.hku1.df = alignment$seq[3] %>% strsplit("") %>% unlist() %>% as.character %>% as.data.frame()
    hcov.nl63.df = alignment$seq[4] %>% strsplit("") %>% unlist() %>% as.character %>% as.data.frame()
    hcov.229e.df = alignment$seq[5] %>% strsplit("") %>% unlist() %>% as.character %>% as.data.frame()
    
    counter.sars.cov2 = 1
    counter.hcov.oc43 = 1
    counter.hcov.hku1 = 1
    counter.hcov.nl63 = 1
    counter.hcov.229e = 1
    for (j in 1:nrow(sars.cov2.df)) {
      
      if(j == 1){
        sars.cov2.df$count = NA
        hcov.oc43.df$count = NA
        hcov.hku1.df$count = NA
        hcov.nl63.df$count = NA
        hcov.229e.df$count = NA
      }
      
      if (! sars.cov2.df$.[j] == "-"){
        sars.cov2.df$count[j] = sars.cov2.counts[['y']][counter.sars.cov2]
        counter.sars.cov2 = counter.sars.cov2 + 1
      }
      
      if (! hcov.oc43.df$.[j] == "-"){
        hcov.oc43.df$count[j] = hcov.oc43.counts[['y']][counter.hcov.oc43]
        counter.hcov.oc43 = counter.hcov.oc43 + 1
      }
      
      if (! hcov.hku1.df$.[j] == "-"){
        hcov.hku1.df$count[j] = hcov.hku1.counts[['y']][counter.hcov.hku1]
        counter.hcov.hku1 = counter.hcov.hku1 + 1
      }
      
      if (! hcov.nl63.df$.[j] == "-"){
        hcov.nl63.df$count[j] = hcov.nl63.counts[['y']][counter.hcov.nl63]
        counter.hcov.nl63 = counter.hcov.nl63 + 1
      }
      
      if (! hcov.229e.df$.[j] == "-"){
        hcov.229e.df$count[j] = hcov.229e.counts[['y']][counter.hcov.229e]
        counter.hcov.229e = counter.hcov.229e + 1
      }
    }
    
    plot(sars.cov2.df$count,
         type = "l",
         axes=F,
         col = "firebrick",
         ylim = c(0, max_count),
         ylab = "counts",
         xlab = "aligned sequence",
         main = homologs$trivial[i])
    
    points(hcov.oc43.df$count,
           type = "l",
           col = "darkslategray")
    
    points(hcov.hku1.df$count,
           type = "l",
           col = "dodgerblue")
    
    points(hcov.nl63.df$count,
           type = "l",
           col = "darkorchid")
    
    points(hcov.229e.df$count,
           type = "l",
           col = "aquamarine")
    
    legend("topright",
           legend = viruses,
           lty = rep(1, 5),
           col = c("firebrick", "darkslategray", "dodgerblue", "darkorchid", "aquamarine"),
           cex = 0.6)
    
    axis(1,
         at = seq(nrow(sars.cov2.df))[as.numeric(rownames(sars.cov2.df)) %% 1 == 0],
         labels = sars.cov2.df$.[as.numeric(rownames(sars.cov2.df)) %% 1 == 0])
    axis(2)
  }
  
  dev.off()
}

plot_alignments(best.m_counts, "best_model")
plot_alignments(last.m_counts, "last_model")


########## identify candidate sequences ##########
# gaussian mixture model on rolling mean

identify_candidates = function(counts.ls, model, virus){
  counts.ls = counts.ls[["counts.rollmean"]]
  cnt.prots = prots[prots$organism == virus & prots$accession %in% names(counts.ls), ]
  counts.ls = counts.ls[names(counts.ls) %in% cnt.prots$accession]
  
  out = data.frame(accession = NA,
                   annotation = NA,
                   start = NA,
                   end = NA,
                   antigen_sequence = NA,
                   mean_pred_count = NA)
  
  pb = txtProgressBar(min = 0, max = length(counts.ls), style = 3)
  for (i in 1:length(counts.ls)) {
    setTxtProgressBar(pb, i)
    
    mean_pred_counts = c()
    antigen_sequences = c()
    
    cnt.name = names(counts.ls)[i]
    cnt.seq = cnt.prots$sequence[cnt.prots$accession == cnt.name]
    
    if(min(counts.ls[[i]]) < 0){
      counts.ls[[i]] = counts.ls[[i]] - min(counts.ls[[i]])
    }
    
    # mixture model to identify hotspot regions
    slope = rep(NA, length(counts.ls[[i]]))
    slope[1] = 0
    for (j in 2:length(counts.ls[[i]])) {
      slope[j] = counts.ls[[i]][j] - counts.ls[[i]][j-1]
    }
    roundSlope = -(round(log10(abs(mean(slope))))-2)
    slope = round(slope,digits=roundSlope)
    
    turnpoints = rep(FALSE,(length(counts.ls[[i]])-1))
    for(j in 1:(length(slope)-1)){
      if(slope[j]>0 & slope[j+1]<0){ turnpoints[j] = TRUE  }
    }
    
    numPeaks = length(which(turnpoints==TRUE))
    centrePeaks = which(turnpoints==TRUE)
    
    if(slope[length(slope)]>0 & counts.ls[[i]][length(counts.ls[[i]])]>1){
      numPeaks = numPeaks+1
      centrePeaks = c(centrePeaks,length(counts.ls[[i]]))
    }
    
    if(slope[2]<0 & counts.ls[[i]][1]>1){
      numPeaks = numPeaks+1
      centrePeaks = c(centrePeaks,1)
    }
    
    numInt = round(length(centrePeaks)/4)+1
    Lint = ceiling(length(counts.ls[[i]])/numInt)
    
    allRegions = numeric()
    
    for(p in 1:numInt){
      #  print(numInt-p)
      s = (p-1)*Lint+1
      e = min(p*Lint,length(counts.ls[[i]]))
      N = 100*Lint
      interval = c(s:e)
      inInterval = which(centrePeaks%in%interval)
      if(length(inInterval)>0){
        samples = sample(c(s:e),N,replace=TRUE,prob=counts.ls[[i]][interval])
        
        mixModel = function(){
          tryCatch(
            {
              return(normalmixEM(samples,
                                 mu=centrePeaks[inInterval],
                                 epsilon=10**(-10),
                                 maxit = 2000))
            },
            error = function(error_message){
              message(error_message)
              return(NA)
            }
          )
        }
        
        x = mixModel()
        
        if (! is.na(x)) {
          means = x[[3]]
          sds = x[[4]]
          
          regions = cbind(round(means-2*sds),round(means+2*sds))
          
          k = which(regions<1,arr.ind=TRUE)
          if(length(k)>0){
            regions[k] = 1
          }
          
          k = which(regions>length(counts.ls[[i]]),arr.ind=TRUE)
          if(length(k)>0){
            regions[k] = length(counts.ls[[i]])
          }
          
          
          maxRegion = rep(NA,dim(regions)[1])
          for(j in 1:dim(regions)[1]){
            centre = regions[j,1]+(regions[j,2]-regions[j,1])/2
            maxRegion[j] = max(counts.ls[[i]][round(centre)])
          }
          
          k = which(maxRegion>=mean(counts.ls[[i]]))
          if(length(k)>0){
            allRegions = rbind(allRegions,regions[k,])
            
            for (k2 in 1:length(k)) {
              mean_pred_counts = c(mean_pred_counts,
                                   counts.ls[[i]][regions[k2, 1]: regions[k2, 2]] %>% mean())
              antigen_sequences = c(antigen_sequences,
                                    substr(cnt.seq, start = regions[k2, 1], stop = regions[k2, 2]))
            }
          }
          
        }
        
      }
    }
    
    # append to results table
    if (length(allRegions) > 0){
      cnt = data.frame(accession = rep(cnt.name, nrow(allRegions)),
                       annotation = rep(cnt.prots$annotation[cnt.prots$accession == cnt.name], nrow(allRegions)),
                       start = allRegions[, 1],
                       end = allRegions[, 2],
                       antigen_sequence = antigen_sequences,
                       mean_pred_count = mean_pred_counts)
      
      out = rbind(out,
                  cnt %>% unique())
    }
    
  }
  
  out = na.omit(out)
  
  # filter regions
  # longer than 10 amino acids
  out = out[out$end - out$start >= 10, ]
  
  # high predicted counts
  q = quantile(out$mean_pred_count, 0.8)
  out = out[out$mean_pred_count >= q, ]
  
  write.csv(out, paste0("validation/results/candidates_", model, "_", virus, ".csv"), row.names = F)
  return(out)
}

candidates.best = identify_candidates(counts.ls = best.m_counts,
                                      model = "best_model",
                                      virus = "SARS-CoV2")

candidates.last = identify_candidates(counts.ls = last.m_counts,
                                      model = "last_model",
                                      virus = "SARS-CoV2")


# plot candidate sequences
plot_candidates = function() {
  cnt.prots = prots[prots$organism == "SARS-CoV2" &
                      prots$accession %in% unique(best.m_prediction$accession), ]
  
  pdf(paste0("validation/results/profiles_SARS-CoV2_candidates.pdf"),
      width = 12, height = 8)
  
  for (i in 1:nrow(cnt.prots)){
    
    # profiles
    y = best.m_counts[["counts"]][[i]]
    y.mean_best = best.m_counts[["counts.rollmean"]][[i]]
    y.mean_last = last.m_counts[["counts.rollmean"]][[i]]
    
    plot(seq(1, length(y.mean_best), ceiling(length(y.mean_best)/length(y))),
         y.mean_best,
         type = "l",
         col = "deeppink",
         ylim = c(0, max_count),
         axes = F,
         ylab = "counts",
         xlab = "position",
         main = cnt.prots$annotation[i],
         sub = paste0("accession: ", prots$accession[i]))
    
    points(seq(1, length(y.mean_last), ceiling(length(y.mean_last)/length(y))),
           y.mean_last,
           type = "l",
           col="darkturquoise")
    
    abline(h = 1, col = "darkgreen")
    
    # candidates as rectangles
    cnt.candidates_best = candidates.best[candidates.best$accession == cnt.prots$accession[i], ]
    cnt.candidates_last = candidates.last[candidates.last$accession == cnt.prots$accession[i], ]
    
    if(nrow(cnt.candidates_best) > 0) {
      for (b in 1:nrow(cnt.candidates_best)) {
        rect(xleft = cnt.candidates_best$start[b], ybottom = 0,
             xright = cnt.candidates_best$end[b], ytop = max_count,
             col = rgb(255/255,20/255,147/255, alpha = 0.4),
             lwd = 0)
      }
    }
    
    if(nrow(cnt.candidates_last) > 0) {
      for (l in 1:nrow(cnt.candidates_last)) {
        rect(xleft = cnt.candidates_last$start[l], ybottom = 0,
             xright = cnt.candidates_last$end[l], ytop = max_count,
             col = rgb(0,206/255,209/255, alpha = 0.4),
             lwd = 0)
      }
    }
    
    legend("topright",
           legend = c("best model", "last model"),
           col = c("deeppink", "darkturquoise"),
           lty = c(1,1),
           cex = .8)
    
    axis(1)
    axis(2)
  }
  
  dev.off()
}

# antigen sequences are not correct!!! (???)
plot_candidates()

# re-open
candidates.best = read.csv("validation/results/candidates_best_model_SARS-CoV2.csv",
                           stringsAsFactors = F)
candidates.last = read.csv("validation/results/candidates_last_model_SARS-CoV2.csv",
                           stringsAsFactors = F)

########## compare with epitopes in literature ##########
# preprocessing
literature_epitopes$epitope = str_remove_all(literature_epitopes$epitope, " ")

for (e in 1:nrow(literature_epitopes)) {
  cnt.ref = prots$annotation[str_detect(prots$sequence,
                                        literature_epitopes$epitope[e])]
  if (length(cnt.ref) > 1) {
    cnt.ref = cnt.ref[cnt.ref == literature_epitopes$annotation[e]]
  }

  if (length(cnt.ref) > 0) {
    literature_epitopes$correct_annotation[e] = cnt.ref
    literature_epitopes[e, c("start", "end")] = str_locate(prots$sequence[prots$annotation == cnt.ref],
                                                           literature_epitopes$epitope[e])

    } else { literature_epitopes$correct_annotation[e] = NA }

}

literature_epitopes = na.omit(literature_epitopes)
# check for which epitopes annotation was wrong / incomplete
# literature_epitopes[literature_epitopes$annotation != literature_epitopes$correct_annotation, ]

epitopes = IEDB_epitopes[, c("Epitope.2", "Epitope.5", "Epitope.6", "Epitope.9")]
epitopes = epitopes[-1, ]
colnames(epitopes) = c("epitope", "start", "end", "correct_annotation")
epitopes$start = as.numeric(epitopes$start)
epitopes$end = as.numeric(epitopes$end)

all_epitopes = rbind(epitopes, literature_epitopes[, c("epitope", "start", "end", "correct_annotation")]) %>%
  unique()
names(all_epitopes) = c("epitope", "start", "end", "annotation")

plot_literature_epitopes = function(literature_epitopes){
  cnt.prots = prots[prots$organism == "SARS-CoV2", ]
  
  pdf(paste0("validation/results/profiles_literature-epitopes_SARS-CoV2.pdf"),
      width = 12, height = 8)
  for (i in 1:nrow(cnt.prots)){
    
    cnt.epitopes = literature_epitopes[toupper(literature_epitopes$annotation) == toupper(cnt.prots$annotation[i]), ]
    
    y = best.m_counts[["counts"]][[i]]
    y.mean_best = best.m_counts[["counts.rollmean"]][[i]]
    y.mean_last = last.m_counts[["counts.rollmean"]][[i]]
    
    
    plot(seq(1, length(y.mean_best), ceiling(length(y.mean_best)/length(y))),
         y.mean_best,
         type = "l",
         col = "black",
         ylim = c(0, max_count),
         axes = F,
         ylab = "counts",
         xlab = "position",
         main = cnt.prots$annotation[i],
         sub = paste0("accession: ", prots$accession[i]))
    
    points(seq(1, length(y.mean_last), ceiling(length(y.mean_last)/length(y))),
           y.mean_last,
           type = "l",
           lty = "dashed",
           col="black")
    
    abline(h = 1, col = "darkgreen")
    
    if (nrow(cnt.epitopes) > 0){
      cnt.epitopes = cnt.epitopes[order(cnt.epitopes$start), ]
      y0 = 0.1
      for (p in 1:nrow(cnt.epitopes)){
        
        if (length(which(cnt.epitopes$start[1:p] <= cnt.epitopes$start[p] & 
                   cnt.epitopes$end[1:p] >= cnt.epitopes$start[p])) > 1) {
          
          cnt.intercept = length(which(cnt.epitopes$start[1:p] <= cnt.epitopes$start[p] & 
                                         cnt.epitopes$end[1:p] >= cnt.epitopes$start[p]))*0.1
          
        } else {cnt.intercept = y0}
        
        segments(x0 = cnt.epitopes$start[p], y0 = cnt.intercept,
                 x1 = cnt.epitopes$end[p], y1 = cnt.intercept,
                 col = "firebrick", lwd = 1.5)
      }
    }
    
    legend("topright",
           legend = c("best model", "last model"),
           lty = c("solid", "dashed"),
           cex = .8)
    
    axis(1)
    axis(2)
  }
  
  dev.off()
}

plot_literature_epitopes(literature_epitopes = all_epitopes)



########## identify candidates ##########
p_50.1.last = read.csv("validation/results/candidates_last_model_SARS-CoV2.csv", stringsAsFactors = F)
p_50.1.best = read.csv("validation/results/candidates_best_model_SARS-CoV2.csv", stringsAsFactors = F)
p_50.2.last = read.csv("validation/results/low_lr_50aa/candidates_last_model_SARS-CoV2.csv", stringsAsFactors = F)
p_50.2.best = read.csv("validation/results/low_lr_50aa/candidates_best_model_SARS-CoV2.csv", stringsAsFactors = F)

p.50 = rbind(p_50.1.last, p_50.1.best, p_50.2.last, p_50.2.best)
p.50[duplicated(p.50$antigen_sequence), ]

p.50 = split.data.frame(p.50, p.50$annotation)

# select regions by hand
candidates = read.csv("validation/candidates.csv", stringsAsFactors = F)
candidates = candidates[candidates$in_final, ]
candidates$length = candidates$end - candidates$start + 1
candidates$sequence = NA

for (i in 1:nrow(candidates)){
  candidates$sequence[i] = substr(prots$sequence[prots$annotation == candidates$annotation[i]],
                                  start = candidates$start[i],
                                  stop = candidates$end[i])
}


# open predictions with protein embedding
spec = "_50aa_corona_"
preds = list.files(path = "validation/results/predictions",
                   pattern = paste0("model_prediction", spec), full.names = T)
best.m_prediction_protEmb = aggregation("best_model")
last.m_prediction_protEmb = aggregation("last_model")

best.m_counts_protEmb = get_counts(best.m_prediction_protEmb)
last.m_counts_protEmb = get_counts(last.m_prediction_protEmb)


# plot
plot_final = function(only_immunogenic = F){
  
  # tmp
  max_count = 2
  
  if (only_immunogenic){
    all_epitopes = literature_epitopes[literature_epitopes$source %in% c("Peng_NatureImmunology_2020",
                                                                         "Ferretti_CellImmunity_2020"),
                                       c("epitope", "start", "end", "correct_annotation")]
    names(all_epitopes) = c("epitope", "start", "end", "annotation")
    fname = "validation/results/profiles_candidates+literature(immunogenic)_SARS-CoV2.pdf"
    
  } else {
    fname = "validation/results/profiles_candidates+literature_SARS-CoV2.pdf"
  }
  
  
  cnt.prots = prots[prots$organism == "SARS-CoV2", ]
  cnt.prots = cnt.prots[-which(duplicated(cnt.prots$sequence)), ]
  
  pdf(paste0(fname),
      width = 12, height = 8)
  for (i in 1:nrow(cnt.prots)){
    
    cnt.epitopes = all_epitopes[toupper(all_epitopes$annotation) == toupper(cnt.prots$annotation[i]), ]
    cnt.candidates = candidates[candidates$annotation == cnt.prots$annotation[i], ]
    
    y = best.m_counts[["counts"]][[i]]
    y.mean_best = best.m_counts[["counts.rollmean"]][[i]]
    y.mean_last = last.m_counts[["counts.rollmean"]][[i]]
    
    y.mean_best_protEmb = best.m_counts_protEmb[["counts.rollmean"]][[i]]
    y.mean_last_protEmb = last.m_counts_protEmb[["counts.rollmean"]][[i]]
    
    # no protein embedding
    plot(seq(1, length(y.mean_best), ceiling(length(y.mean_best)/length(y))),
         y.mean_best,
         type = "l",
         col = "deeppink",
         ylim = c(0, max_count),
         axes = F,
         ylab = "counts",
         xlab = "position",
         main = cnt.prots$annotation[i],
         sub = paste0("accession: ", prots$accession[i]))
    
    points(seq(1, length(y.mean_last), ceiling(length(y.mean_last)/length(y))),
           y.mean_last,
           type = "l",
           lty = "dashed",
           col= "deeppink")
    
    # protein embedding
    points(seq(1, length(y.mean_best_protEmb), ceiling(length(y.mean_best_protEmb)/length(y))),
         y.mean_best_protEmb,
         type = "l",
         col = "darkviolet",
         ylim = c(0, max_count),
         axes = F,
         ylab = "counts",
         xlab = "position",
         main = cnt.prots$annotation[i],
         sub = paste0("accession: ", prots$accession[i]))
    
    points(seq(1, length(y.mean_last_protEmb), ceiling(length(y.mean_last_protEmb)/length(y))),
           y.mean_last_protEmb,
           type = "l",
           lty = "dashed",
           col= "darkviolet")
    
    # "threshold"
    abline(h = 1, col = "black")
    
    # known epitopes
    if (nrow(cnt.epitopes) > 0){
      cnt.epitopes = cnt.epitopes[order(cnt.epitopes$start), ]
      y0 = 0.1
      for (p in 1:nrow(cnt.epitopes)){
        
        if (length(which(cnt.epitopes$start[1:p] <= cnt.epitopes$start[p] & 
                         cnt.epitopes$end[1:p] >= cnt.epitopes$start[p])) > 1) {
          
          cnt.intercept = length(which(cnt.epitopes$start[1:p] <= cnt.epitopes$start[p] & 
                                         cnt.epitopes$end[1:p] >= cnt.epitopes$start[p]))*0.1
          
        } else {cnt.intercept = y0}
        
        segments(x0 = cnt.epitopes$start[p], y0 = cnt.intercept,
                 x1 = cnt.epitopes$end[p], y1 = cnt.intercept,
                 col = "cyan3", lwd = 1.5)
      }
      
    }
    
    # substrate candidates
    if (nrow(cnt.candidates) > 0) {
      for (c in 1:nrow(cnt.candidates)) {
        segments(x0 = cnt.candidates$start[c], y0 = .05,
                 x1 = cnt.candidates$end[c], y1 = .05,
                 col = "chartreuse3", lwd = 1.5)
      }
    }
    
    # legend
    legend("topright",
           legend = c("best model", "last model",
                      "no protein embedding", "protein embedding",
                      "epitopes", "substrates"),
           lty = c("solid", "dashed",
                   "solid", "solid",
                   "solid", "solid"),
           col = c("black", "black",
                   "deeppink", "darkviolet",
                   "cyan3", "chartreuse3"),
           cex = .8)
    
    axis(1)
    axis(2)
  }
  
  dev.off()
}

plot_final()
plot_final(only_immunogenic = T)


candidates$in_final = NULL
write.csv(candidates, "validation/candidates_final.csv", row.names = F)


