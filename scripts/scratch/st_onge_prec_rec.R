prec_rec_vs_stonge <- function(gi_data,
                                control_name = "NoDrug",
                                condition_name = "MMS",
                                fdr_prefix = "FDR.neutral_xy",
                                z_prefix = "Z_GIS_xy",
                                gi_prefix = "GIS_xy",
                                fdr_cutoff = 0.01,
                                metr1 = 'prec',
                                metr2 = 'rec',
                                metr1_name = 'Precision',
                                metr2_name = 'Recall',
                                col1 = 'darkblue',
                                col2 = 'darkred',
                                xlims = c(-0.1, 0.1),
                                ylims = c(0,100),
                                xlab = 'GIS Cutoff',
                                ylab1 = 'Precision - St.Onge Validation Rate (%)',
                                ylab2 = 'Recall - St.Onge Recovery Rate (%)',
                                draw_cutoffs = T,
                                legend_pos = c(0,30),
                                cutoffs_drawn = NULL) {
  #Filter
  gi_data_filtered <-
    dplyr::filter(gi_data,
                  SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
  
  gi_data_filtered <-
    dplyr::filter(gi_data_filtered,
                  Remove_by_Chromosomal_distance_or_SameGene == 'no')
  
  
  #Joint performance
  labels_pos <- c()
  labels_neg <- c()
  scores_cond <- c()
  
  for (condition in c(control_name, condition_name)) {
    if (condition == control_name) {
      st_onge_class <- 'SOJ_Class_NoMMS'
    }
    if (condition == condition_name) {
      st_onge_class <- 'SOJ_Class_MMS'
    }
    
    fdr_column <- paste(c(fdr_prefix, condition), collapse = '.')
    z_column <- paste(c(z_prefix, condition), collapse = '.')
    gi_column <- paste(c(gi_prefix, condition), collapse = '.')
    
    labels_pos <-
      c(labels_pos, gi_data_filtered[, st_onge_class] == 'ALLEVIATING')
    labels_neg <-
      c(labels_neg, gi_data_filtered[, st_onge_class] == 'AGGRAVATING')
    
    scores_cond <-
      c(scores_cond,
        gi_data_filtered[, gi_column] * as.numeric(gi_data_filtered[, fdr_column] < fdr_cutoff))

  }
  
  
  
  pos_metr1 <-
    ROCR::performance(ROCR::prediction(scores_cond, labels_pos), metr1)
  neg_metr1 <-
    ROCR::performance(ROCR::prediction(-scores_cond, labels_neg), metr1)
  
  pos_metr2 <-
    ROCR::performance(ROCR::prediction(scores_cond, labels_pos), metr2)
  neg_metr2 <-
    ROCR::performance(ROCR::prediction(-scores_cond, labels_neg), metr2)
  
  
  
  pos_cutoff <- pos_metr1@x.values[[1]]
  pos_perf1 <- pos_metr1@y.values[[1]]*100
  pos_perf2 <- pos_metr2@y.values[[1]]*100
  
  
  
  
  neg_cutoff <- (-1)*neg_metr1@x.values[[1]]
  neg_perf1 <- neg_metr1@y.values[[1]]*100
  neg_perf2 <- neg_metr2@y.values[[1]]*100
  
  pars <- list('neg'=list(neg_cutoff,neg_perf1,neg_perf2),
               'pos'=list(pos_cutoff,pos_perf1,pos_perf2))
  signs <- c(`<`,`>`)
  titles <- c('Negative Interactions','Positive Interactions')
  
  par(mfrow=c(1,2))
  i <- 0
  for(par in pars){
    par <- par
    i <- i + 1
    cutoff <- par[[1]]
    perf1 <- par[[2]]
    perf2 <- par[[3]]
    
    par(las = 3)
    par(mar = c(5,5,2,5.5))
    par(xpd = F)
    plot(
      cutoff[signs[[i]](cutoff,0)],
      perf1[signs[[i]](cutoff,0)],
      type = 'l',
      xlim = c(0, xlims[i]),
      ylim = ylims,
      ylab = ylab1,
      xlab = xlab,
      col = col1,
      lwd = 2,
      main = titles[i]
    )
    axis(side = 4)
    
    mtext(side = 4, line = 3, ylab2)
    lines(
      cutoff[signs[[i]](cutoff,0)],
      perf2[signs[[i]](cutoff,0)],
      xlim = c(0, xlims[i]),
      ylim = ylims,
      col = col2,
      lwd = 2
    )
    legend(legend_pos[1],
           legend_pos[2],
           legend = c('Precision', 'Recall'),
           col = c(col1, col2),
           pch = 15,
           cex = 1.1)
    if (!is.null(cutoffs_drawn)) {
      abline(v = cutoffs_drawn[i], lty = 3)
      
      p1 <- perf1[which.min(abs(cutoff - cutoffs_drawn[i]))]
      p2 <- perf2[which.min(abs(cutoff - cutoffs_drawn[i]))]
      
      p1_text <- paste(c(metr1_name,' = ',format(p1,digits=2),'%'),collapse='')
      p2_text <- paste(c(metr2_name,' = ',format(p2,digits=2),'%'),collapse='')
      
      text(cutoffs_drawn[i],p1,p1_text,pos=1)
      text(cutoffs_drawn[i],p2,p2_text,pos=1)
    }
  }
  #neg_cutoff <- neg_perf@x.values[[1]]
  #neg_precision <- neg_perf@y.values[[1]]
  
  #lines(pos_cutoff[pos_cutoff > 0],
  #      (pos_precision * 100)[pos_cutoff > 0],
  #      lwd = 1,
  #      col = 'black')
  #lines(-neg_cutoff[neg_cutoff > 0],
  #      (neg_precision * 100)[neg_cutoff > 0],
  #      lwd = 1,
  #      col = 'black')
  
  #best_neg <- neg_cutoff[which.max(neg_precision)]
  #best_pos <- pos_cutoff[which.max(pos_precision)]
  
  #print(best_neg)
  #labels_pos <-
  #  gi_data_filtered[, 'SOJ_Class_NoMMS'] == 'ALLEVIATING'
  #labels_neg <-
  #  gi_data_filtered[, 'SOJ_Class_NoMMS'] == 'AGGRAVATING'
  #
  #  labels_pos <- c(labels_pos,gi_data_filtered[, 'SOJ_Class_MMS'] == 'ALLEVIATING')
  #  labels_neg <- c(labels_neg,gi_data_filtered[, 'SOJ_Class_MMS'] == 'AGGRAVATING')
  
  
  #if (draw_cutoffs == T) {
  #  if (is.null(cutoffs_drawn)) {
  #     abline(v = -best_neg,
  #            lty = 3,
  #            lwd = 0.7)
  #     abline(v = best_pos, lty = 3, lwd = 0.7)
  #   } else{
  #     abline(v = -cutoffs_drawn[1],
  #            lty = 3,
  #            lwd = 0.7)
  #     abline(v = cutoffs_drawn[2],
  #            lty = 3,
  #            lwd = 0.7)
  #   }
  # }
  # legend(
  #   xlims[1],
  #   35,
  #   legend = c(control_name, condition_name, 'Combined'),
  #   fill = c('blue', 'red', 'black')
  # )
  # 
  # return(c(best_neg, best_pos))
}
