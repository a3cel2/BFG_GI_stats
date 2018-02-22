#devtools::use_package('ROCR')
#devtools::use_package('dplyr')


#Used for calculating AUC
#' Trapezoid method for integration
trapezoid_integration <- function(x, y) {
  area_sum <- 0
  for (i in 2:length(x)) {
    base <- x[i] - x[i - 1]
    mid_height <- mean(c(y[i], y[i - 1]))
    area_sum <- area_sum + base * mid_height
  }
  return(area_sum)
}



##Legacy function
#' #' Plots St. Onge et al validation as a function of the internal FDR estimates
#' #'
#' #' @param gi_data input genetic interaction table
#' #' @param control_name condition which corresponds to 'NoDrug' in the St onge data
#' #' @param condition_name condition which corresponds to 'MMS' in the St onge data
#' #' @param fdr_prefix pasted with control_name and condition_name to find the FDR columns
#' #' @param z_prefix pasted with control_name and condition_name to find the Z columns
#' #' @param fdr_cutoff strains above this FDR are filtered out when drawing the plot
#' #' @param xlims plot x axis limits
#' #' @param gi_prefix pasted with control_name and condition_name to find the GIS columns
#' #' @param metr performance metric string, must be a valid ROCR metric
#' #' @param ylab plot y axis limits
#' #' @param draw_cutoffs boolean, whether to draw lines at cutoffs 
#' #' @param gi_cutoff_pos 
#' #' @param gi_cutoff_neg 
#' #' @param cutoffs_drawn 
#' #' @param xlab plot x axis label
#' precision_vs_stonge <- function(gi_data,
#'                                 control_name = "NoDrug",
#'                                 condition_name = "MMS",
#'                                 fdr_prefix = "FDR.neutral_xy",
#'                                 z_prefix = "Z_GIS_xy",
#'                                 gi_prefix = "GIS_xy",
#'                                 fdr_cutoff = 0.05,
#'                                 metr = 'prec',
#'                                 xlims = c(-5, 5),
#'                                 xlab = expression('GIS'),
#'                                 ylab = 'St.Onge Validation Rate (%)',
#'                                 draw_cutoffs = T,
#'                                 gi_cutoff_pos = 0.05,
#'                                 gi_cutoff_neg = -0.07,
#'                                 cutoffs_drawn = NULL) {
#'   #Individual performance
#'   for (condition in c(control_name, condition_name)) {
#'     if (condition == control_name) {
#'       st_onge_class <- 'SOJ_Class_NoMMS'
#'     }
#'     if (condition == condition_name) {
#'       st_onge_class <- 'SOJ_Class_MMS'
#'     }
#'     
#'     gi_data_filtered <-
#'       dplyr::filter(gi_data,
#'                     SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
#'     
#'     gi_data_filtered <-
#'       dplyr::filter(gi_data_filtered,
#'                     Remove_by_Chromosomal_distance_or_SameGene == 'no')
#'     
#'     
#'     
#'     fdr_column <- paste(c(fdr_prefix, condition), collapse = '.')
#'     z_column <- paste(c(z_prefix, condition), collapse = '.')
#'     gi_column <- paste(c(gi_prefix, condition), collapse = '.')
#'     
#'     labels_pos <-
#'       gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
#'     labels_neg <-
#'       gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
#'     
#'     scores_cond <-
#'       gi_data_filtered[, gi_column] * as.numeric(gi_data_filtered[, fdr_column] < fdr_cutoff)
#' 
#'     
#'     pos_perf <-
#'       ROCR::performance(ROCR::prediction(scores_cond, labels_pos), metr)
#'     neg_perf <-
#'       ROCR::performance(ROCR::prediction(-scores_cond, labels_neg), metr)
#'     
#'     pos_cutoff <- pos_perf@x.values[[1]]
#'     pos_precision <- pos_perf@y.values[[1]]
#'     
#'     neg_cutoff <- neg_perf@x.values[[1]]
#'     neg_precision <- neg_perf@y.values[[1]]
#'     
#'     if (condition == control_name) {
#'       if (is.null(xlims)) {
#'         xlims <-
#'           c(-1 * neg_perf@x.values[[1]][2], pos_perf@x.values[[1]][2])
#'       }
#'       par(mar = c(4.5, 5, 3, 1))
#'       par(las = 1)
#'       plot(
#'         pos_cutoff[pos_cutoff > 0],
#'         (pos_precision * 100)[pos_cutoff > 0],
#'         type = 'l',
#'         ylab = ylab,
#'         xlab = xlab,
#'         main = '',
#'         #GI Precision vs St. Onge',
#'         lwd = 1,
#'         xlim = xlims,
#'         ylim = c(0, 100),
#'         col = 'blue'
#'       )
#'       lines(-neg_cutoff[neg_cutoff > 0],
#'             (neg_precision * 100)[neg_cutoff > 0],
#'             lwd = 1,
#'             col = 'blue')
#'     } else{
#'       lines(pos_cutoff[pos_cutoff > 0],
#'             (pos_precision * 100)[pos_cutoff > 0],
#'             lwd = 1,
#'             col = 'red')
#'       lines(-neg_cutoff[neg_cutoff > 0],
#'             (neg_precision * 100)[neg_cutoff > 0],
#'             lwd = 1,
#'             col = 'red')
#'     }
#'   }
#'   
#'   #Joint performance
#'   labels_pos <- c()
#'   labels_neg <- c()
#'   scores_cond <- c()
#'   
#'   for (condition in c(control_name, condition_name)) {
#'     if (condition == control_name) {
#'       st_onge_class <- 'SOJ_Class_NoMMS'
#'     }
#'     if (condition == condition_name) {
#'       st_onge_class <- 'SOJ_Class_MMS'
#'     }
#'     
#'     fdr_column <- paste(c(fdr_prefix, condition), collapse = '.')
#'     z_column <- paste(c(z_prefix, condition), collapse = '.')
#'     gi_column <- paste(c(gi_prefix, condition), collapse = '.')
#'     
#'     labels_pos <-
#'       c(labels_pos, gi_data_filtered[, st_onge_class] == 'ALLEVIATING')
#'     labels_neg <-
#'       c(labels_neg, gi_data_filtered[, st_onge_class] == 'AGGRAVATING')
#'     
#'     scores_cond <-
#'       c(scores_cond,
#'         gi_data_filtered[, gi_column] * as.numeric(gi_data_filtered[, fdr_column] < fdr_cutoff))
#'     #sign(gi_data_filtered[, z_column]) * -log10(as.numeric(gi_data_filtered[, fdr_column]))# *
#'     #as.numeric(
#'     # gi_data_filtered[, gi_column] > gi_cutoff_pos |
#'     #    gi_data_filtered[, gi_column] < gi_cutoff_neg
#'     #)
#'     #)
#'   }
#'   
#'   
#'   
#'   pos_perf <-
#'     ROCR::performance(ROCR::prediction(scores_cond, labels_pos), metr)
#'   neg_perf <-
#'     ROCR::performance(ROCR::prediction(-scores_cond, labels_neg), metr)
#'   
#'   pos_cutoff <- pos_perf@x.values[[1]]
#'   pos_precision <- pos_perf@y.values[[1]]
#'   
#'   neg_cutoff <- neg_perf@x.values[[1]]
#'   neg_precision <- neg_perf@y.values[[1]]
#'   
#'   lines(pos_cutoff[pos_cutoff > 0],
#'         (pos_precision * 100)[pos_cutoff > 0],
#'         lwd = 1,
#'         col = 'black')
#'   lines(-neg_cutoff[neg_cutoff > 0],
#'         (neg_precision * 100)[neg_cutoff > 0],
#'         lwd = 1,
#'         col = 'black')
#'   
#'   best_neg <- neg_cutoff[which.max(neg_precision)]
#'   best_pos <- pos_cutoff[which.max(pos_precision)]
#'   
#'   #print(best_neg)
#'   #labels_pos <-
#'   #  gi_data_filtered[, 'SOJ_Class_NoMMS'] == 'ALLEVIATING'
#'   #labels_neg <-
#'   #  gi_data_filtered[, 'SOJ_Class_NoMMS'] == 'AGGRAVATING'
#'   #
#'   #  labels_pos <- c(labels_pos,gi_data_filtered[, 'SOJ_Class_MMS'] == 'ALLEVIATING')
#'   #  labels_neg <- c(labels_neg,gi_data_filtered[, 'SOJ_Class_MMS'] == 'AGGRAVATING')
#'   
#'   
#'   if (draw_cutoffs == T) {
#'     if (is.null(cutoffs_drawn)) {
#'       abline(v = -best_neg,
#'              lty = 3,
#'              lwd = 0.7)
#'       abline(v = best_pos, lty = 3, lwd = 0.7)
#'     } else{
#'       abline(v = cutoffs_drawn[1],
#'              lty = 3,
#'              lwd = 0.7)
#'       abline(v = cutoffs_drawn[2],
#'              lty = 3,
#'              lwd = 0.7)
#'     }
#'   }
#'   legend(
#'     xlims[1],
#'     35,
#'     legend = c(control_name, condition_name, 'Combined'),
#'     fill = c('blue', 'red', 'black')
#'   )
#'   
#'   return(c(best_neg, best_pos))
#' }



#' Plots FPR-Sensitivty curve versus the St. Onge data
#'
#' @param gi_data input genetic interaction table
#' @param control_name condition which corresponds to 'NoMMS' in the St onge data
#' @param condition_name condition which corresponds to 'MMS' in the St onge data
#' @param fdr_prefix pasted with control_name and condition_name to find the FDR columns
#' @param z_prefix pasted with control_name and condition_name to find the Z columns
#' @param old_data if true, uses the raw GI scores, rather than FDR scores to do AUC
#' @param gi_prefix pasted with control_name and condition_name to find the raw GI columns
#' @param neg_col colour to draw performance curve for negative interactions
#' @param pos_col colour to draw performance curve for positive interactions
#' @param lwd line width in the plot
st_onge_auc_plot <- function(gi_data,
                             control_name = "NoDrug",
                             condition_name = "MMS",
                             fdr_prefix = "FDR.neutral_xy",
                             z_prefix = "Z_GIS_xy",
                             old_data = F,
                             gi_prefix = "GIS_xy",
                             neg_col = rgb(230 / 255, 155 / 255, 34 / 255),
                             pos_col = rgb(90 / 255, 179 / 255, 228 / 255),
                             lwd = 3) {
  conditions <- c(control_name, condition_name)
  
  for (condition in c(conditions, 'both')) {
    if (condition != 'both') {
      condname <- condition
    }
    else{
      condname <- ''
    }
    
    gi_data_filtered <-
      dplyr::filter(gi_data,
                    SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
    
    if (condition != 'both') {
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
        gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
      labels_neg <-
        gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
      
      if (old_data == F) {
        scores_cond <- gi_data_filtered[, z_column]
      } else if (old_data == T) {
        scores_cond <- gi_data_filtered[, gi_column]
      }
    } else{
      labels_pos <- c()
      labels_neg <- c()
      for (st_onge_class in c('SOJ_Class_NoMMS', 'SOJ_Class_MMS')) {
        labels_pos <-
          c(labels_pos, gi_data_filtered[, st_onge_class] == 'ALLEVIATING')
        labels_neg <-
          c(labels_neg, gi_data_filtered[, st_onge_class] == 'AGGRAVATING')
      }
      
      scores_cond <- c()
      if (old_data == F) {
        for (condition in conditions) {
          z_column <- paste(c(z_prefix, condition), collapse = '.')
          scores_cond <- c(scores_cond, gi_data_filtered[, z_column])
        }
      } else if (old_data == T) {
        for (condition in conditions) {
          gi_column <- paste(c(gi_prefix, condition), collapse = '.')
          scores_cond <-
            c(scores_cond, gi_data_filtered[, gi_column])
        }
      }
    }
    
    perf_pos <-
      ROCR::performance(ROCR::prediction(scores_cond, labels_pos), 'sens', 'fpr')
    perf_neg <-
      ROCR::performance(ROCR::prediction(-scores_cond, labels_neg), 'sens', 'fpr')
    
    par(las = 1)
    plot(
      perf_pos@x.values[[1]] * 100,
      perf_pos@y.values[[1]] * 100,
      lwd = lwd,
      col = pos_col,
      main = condname,
      xlab = 'False Positive Rate (%)',
      ylab = 'Sensitivity (%)',
      type = 'l',
      cex.lab = 1.9,
      cex.axis = 1.7
    )
    
    lines(c(0, 100),
          c(0, 100),
          col = 'grey80',
          lty = 5,
          lwd = lwd / 2)
    
    lines(
      perf_neg@x.values[[1]] * 100,
      perf_neg@y.values[[1]] * 100,
      lwd = lwd,
      col = neg_col
    )
    
    auc_neg <-
      trapezoid_integration(perf_neg@x.values[[1]], perf_neg@y.values[[1]])
    auc_pos <-
      trapezoid_integration(perf_pos@x.values[[1]], perf_pos@y.values[[1]])
    
    auc_neg <- format(auc_neg, digits = 2)
    auc_pos <- format(auc_pos, digits = 2)
    
    
    t1_neg <- bquote(phantom(bold("Negative GIs"))*" - "*.(auc_neg))
    t2_neg <- bquote(bold("Negative GIs")*phantom(" - "*.(auc_neg)))
    
    t1_pos <- bquote(phantom(bold("Positive GIs"))*" - "*.(auc_pos))
    t2_pos <- bquote(bold("Positive GIs")*phantom(" - "*.(auc_pos)))
    

    text(10, 40, expression(underline('AUC')), col = 'black', cex = 1.7, adj = c(0,1))
    
    text(10, 30, t1_neg, col = 'black', cex = 1.7, adj = c(0,1))
    text(10, 30, t2_neg, col = neg_col, cex = 1.7, adj = c(0,1))
    
    text(10, 20, t1_pos, col = 'black', cex = 1.7, adj = c(0,1))
    text(10, 20, t2_pos, col = pos_col, cex = 1.7, adj = c(0,1))
    

  }
  
  
  
}



#' Creates a scatterplot of genetic interaction scores vs St. Onge
#'
#' @param gi_data input genetic interaction table
#' @param control_name condition which corresponds to 'NoMMS' in the St onge data
#' @param condition_name condition which corresponds to 'MMS' in the St onge data
#' @param point_col point colour
#' @param gi_prefix pasted with control_name and condition_name to find the raw GI columns
st_onge_scatterplot <- function(gi_data,
                                control_name = "NoDrug",
                                condition_name = "MMS",
                                gi_prefix = "GIS_xy",
                                point_col = rgb(0, 0, 0, 0.5)) {
  gi_data <-
    dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')
  
  
  for (condition in c(control_name, condition_name)) {
    if (condition == control_name) {
      st_onge_e <- 'SOJ_E_NoMMS'
      condtext <- 'No Drug'
    }
    if (condition == condition_name) {
      st_onge_e <- 'SOJ_E_MMS'
      condtext <- 'MMS'
    }
    gi_column <- paste(c(gi_prefix, condition), collapse = '.')
    
    par(las = 1)
    par(mar = c(4, 5, 1, 1))
    plot(
      gi_data[, st_onge_e],
      gi_data[, gi_column],
      pch = 16,
      main = condtext,
      cex = 0.8,
      col = point_col,
      xlab = '',
      ylab = '',
      axes = 1,
      bty = 'l',
      xaxt = 'n',
      yaxt = 'n'
    )
    par(las = 3)
    axis(side = 1)
    par(las = 1)
    axis(side = 2)
    
    par(las = 3)
    mtext(expression(bold(GIS) ~'          '),
          2,
          line = 3,
          cex = 1.5)
    
    mtext('           (BFG-GI)',
          2,
          line = 3.2,
          cex = 1)
    par(las = 1)
    mtext(expression(epsilon ~ "         "), 1, line = 2.8 , cex = 2)
    mtext(expression( ~(St.*Onge)), 1, line = 2.9)
    #c('Genetic Interaction Score','(this study)'))
    cor_text <- sprintf('r = %s', format(cor(gi_data[, st_onge_e],
                                             gi_data[, gi_column], use = 'pair'), digits = 2))
    
    text(min(gi_data[, st_onge_e], na.rm = T),
         max(gi_data[, gi_column], na.rm = T),
         cor_text,
         adj = c(0, 1),
         cex = 1.3)
  }
}



#Legacy scripts not used to generate final figure, compares old and new GI data
#on tpr-fpr AUC and scatterplot correlation wih St. Onge data

# st_onge_comparison_plot <- function(gi_data_old,
#                                     gi_data_new,
#                                     control_name = "GIS_ij.DMSO",
#                                     condition_name = "GIS_ij.MMS"){
#   #par(mfrow = c(1, 2))
#
#   gi_data_list <- list(gi_data_new, gi_data_old)
#   names(gi_data_list) <- c('new', 'old')
#
#
#   #Compare AUC
#   #gi_scores <- gi_data_new[, grep('^GIS', colnames(gi_data_new))]
#   for (condition in c(control_name, condition_name)) {
#     for (i in 1:length(gi_data_list)) {
#       gi_data <- gi_data_list[[i]]
#       gi_data_filtered <-
#         dplyr::filter(gi_data,
#                       SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
#
#       #gi_scores <- grep('^GIS', colnames(gi_data))]
#
#
#       if (condition == control_name) {
#         st_onge_class <- 'SOJ_Class_NoMMS'
#       }
#       if (condition == condition_name) {
#         st_onge_class <- 'SOJ_Class_MMS'
#       }
#
#       labels_pos <- gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
#       labels_neg <- gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
#       z_scores_cond_pos <-
#         gi_data_filtered[, condition]#, collapse = '')]
#       z_scores_cond_neg <- -z_scores_cond_pos
#
#
#       perf_pos <-
#         ROCR::performance(ROCR::prediction(z_scores_cond_pos, labels_pos), 'sens', 'fpr')
#       perf_neg <-
#         ROCR::performance(ROCR::prediction(z_scores_cond_neg, labels_neg), 'sens', 'fpr')
#
#       ROCR::plot(
#         perf_pos,
#         lwd = 2,
#         col = 'blue',
#         main = paste(names(gi_data_list)[i],condition)
#       )
#       lines(perf_neg@x.values[[1]],
#             perf_neg@y.values[[1]],
#             lwd = 2,
#             col = 'red')
#
#       auc_neg <-
#         trapezoid_integration(perf_neg@x.values[[1]], perf_neg@y.values[[1]])
#       auc_pos <-
#         trapezoid_integration(perf_pos@x.values[[1]], perf_pos@y.values[[1]])
#
#       text(0.8, 0.6, sprintf('AUC neg = %s', format(auc_neg, digits = 2)), col =
#              'red')
#       text(0.8, 0.4, sprintf('AUC pos = %s', format(auc_pos, digits = 2)), col =
#              'blue')
#
#     }
#   }
#
#   #Correlate with GI scores
#   for (condition in c(control_name,condition_name)) {
#     for (i in 1:length(gi_data_list)) {
#       gi_data <- gi_data_list[[i]]
#       if (condition == control_name) {
#         st_onge_e <- 'SOJ_E_NoMMS'
#       }
#       if (condition == condition_name) {
#         st_onge_e <- 'SOJ_E_MMS'
#       }
#       plot(gi_data[, st_onge_e], gi_data[, condition], pch = 16,
#            main= paste(names(gi_data_list)[i],condition), col = rgb(0,0,0,0.3),
#            xlab = 'St Onge E',ylab ='GIS')
#       lines(smooth.spline(gi_data[is.finite(gi_data[, st_onge_e]), st_onge_e],
#                           gi_data[is.finite(gi_data[, st_onge_e]), condition],
#                           penalty=10), col = 'blue')
#       abline(c(0,1),lwd=2,col='red')
#       correl <- cor(gi_data[, st_onge_e], gi_data[, condition],use='pair')
#
#       text(0,0,sprintf('R = %s',format(correl,digits=2)),col='red',cex=1.5)
#       #stop()
#     }
#
#
#   }
#
# }


# 
# precision_recall_vs_stonge <- function(gi_data,
#                                        control_name = "NoDrug",
#                                        condition_name = "MMS",
#                                        fdr_prefix = "FDR.neutral_xy",
#                                        z_prefix = "Z_GIS_xy",
#                                        gi_prefix = "GIS_xy",
#                                        fdr_cutoff = 0.05,
#                                        metr = 'prec',
#                                        xlims = c(-5, 5),
#                                        xlab = expression('-Log'[10] * '(FDR'['neutral'] *
#                                                            ')'),
#                                        ylab = 'St.Onge Validation Rate (%)',
#                                        draw_cutoffs = T,
#                                        gi_cutoff_pos = 0.05,
#                                        gi_cutoff_neg = -0.07,
#                                        cutoffs_drawn = NULL) {
#   #Individual performance
#   for (condition in c(control_name, condition_name)) {
#     if (condition == control_name) {
#       st_onge_class <- 'SOJ_Class_NoMMS'
#     }
#     if (condition == condition_name) {
#       st_onge_class <- 'SOJ_Class_MMS'
#     }
#     
#     gi_data_filtered <-
#       dplyr::filter(gi_data,
#                     SOJ_Class_NoMMS %in% c('NEUTRAL', 'AGGRAVATING', 'ALLEVIATING'))
#     
#     gi_data_filtered <-
#       dplyr::filter(gi_data_filtered,
#                     Remove_by_Chromosomal_distance_or_SameGene == 'no')
#     
#     
#     
#     fdr_column <- paste(c(fdr_prefix, condition), collapse = '.')
#     z_column <- paste(c(z_prefix, condition), collapse = '.')
#     gi_column <- paste(c(gi_prefix, condition), collapse = '.')
#     
#     labels_pos <-
#       gi_data_filtered[, st_onge_class] == 'ALLEVIATING'
#     labels_neg <-
#       gi_data_filtered[, st_onge_class] == 'AGGRAVATING'
#     
#     scores_cond <-
#       gi_data_filtered[, gi_column] * as.numeric(gi_data_filtered[, fdr_column] < fdr_cutoff)
#     #gi_data_filtered[, z_column]
#     #as.numeric(gi_data_filtered[, gi_column] > gi_cutoff_pos |
#     #             gi_data_filtered[, gi_column] < gi_cutoff_neg)
#     
#     for (metr in c('prec','rec')) {
#       pos_perf <-
#         ROCR::performance(ROCR::prediction(scores_cond, labels_pos), metr)
#       neg_perf <-
#         ROCR::performance(ROCR::prediction(-scores_cond, labels_neg), metr)
#       
#       pos_cutoff <- pos_perf@x.values[[1]]
#       pos_precision <- pos_perf@y.values[[1]]
#       
#       neg_cutoff <- neg_perf@x.values[[1]]
#       neg_precision <- neg_perf@y.values[[1]]
#       
#       if (condition == control_name) {
#         if (is.null(xlims)) {
#           xlims <-
#             c(-1 * neg_perf@x.values[[1]][2], pos_perf@x.values[[1]][2])
#         }
#         par(mar = c(4.5, 5, 3, 1))
#         par(las = 1)
#         plotting_func <- ifelse(metr=='prec',`plot`,`points`)
#         plotting_func(
#           pos_cutoff[pos_cutoff > 0],
#           (pos_precision * 100)[pos_cutoff > 0],
#           type = 'l',
#           ylab = ylab,
#           xlab = xlab,
#           main = '',
#           #GI Precision vs St. Onge',
#           lwd = 1,
#           xlim = xlims,
#           ylim = c(0, 100),
#           col = 'blue'
#         )
#         lines(-neg_cutoff[neg_cutoff > 0],
#               (neg_precision * 100)[neg_cutoff > 0],
#               lwd = 1,
#               col = 'blue')
#       } else{
#         lines(pos_cutoff[pos_cutoff > 0],
#               (pos_precision * 100)[pos_cutoff > 0],
#               lwd = 1,
#               col = 'red')
#         lines(-neg_cutoff[neg_cutoff > 0],
#               (neg_precision * 100)[neg_cutoff > 0],
#               lwd = 1,
#               col = 'red')
#       }
#     }
#   }
#   
#   #Joint performance
#   labels_pos <- c()
#   labels_neg <- c()
#   scores_cond <- c()
#   
#   for (condition in c(control_name, condition_name)) {
#     if (condition == control_name) {
#       st_onge_class <- 'SOJ_Class_NoMMS'
#     }
#     if (condition == condition_name) {
#       st_onge_class <- 'SOJ_Class_MMS'
#     }
#     
#     fdr_column <- paste(c(fdr_prefix, condition), collapse = '.')
#     z_column <- paste(c(z_prefix, condition), collapse = '.')
#     gi_column <- paste(c(gi_prefix, condition), collapse = '.')
#     
#     labels_pos <-
#       c(labels_pos, gi_data_filtered[, st_onge_class] == 'ALLEVIATING')
#     labels_neg <-
#       c(labels_neg, gi_data_filtered[, st_onge_class] == 'AGGRAVATING')
#     
#     scores_cond <-
#       c(scores_cond,
#         gi_data_filtered[, gi_column] * as.numeric(gi_data_filtered[, fdr_column] < fdr_cutoff))
#     #sign(gi_data_filtered[, z_column]) * -log10(as.numeric(gi_data_filtered[, fdr_column]))# *
#     #as.numeric(
#     # gi_data_filtered[, gi_column] > gi_cutoff_pos |
#     #    gi_data_filtered[, gi_column] < gi_cutoff_neg
#     #)
#     #)
#   }
#   
#   
#   for (metr in c('prec', 'rec')) {
#     pos_perf <-
#       ROCR::performance(ROCR::prediction(scores_cond, labels_pos), metr)
#     neg_perf <-
#       ROCR::performance(ROCR::prediction(-scores_cond, labels_neg), metr)
#     
#     pos_cutoff <- pos_perf@x.values[[1]]
#     pos_precision <- pos_perf@y.values[[1]]
#     
#     neg_cutoff <- neg_perf@x.values[[1]]
#     neg_precision <- neg_perf@y.values[[1]]
#     
#     lines(pos_cutoff[pos_cutoff > 0],
#           (pos_precision * 100)[pos_cutoff > 0],
#           lwd = 1,
#           col = 'black')
#     lines(-neg_cutoff[neg_cutoff > 0],
#           (neg_precision * 100)[neg_cutoff > 0],
#           lwd = 1,
#           col = 'black')
#   }
#   
#   
#   best_neg <- neg_cutoff[which.max(neg_precision)]
#   best_pos <- pos_cutoff[which.max(pos_precision)]
#   
#   #print(best_neg)
#   #labels_pos <-
#   #  gi_data_filtered[, 'SOJ_Class_NoMMS'] == 'ALLEVIATING'
#   #labels_neg <-
#   #  gi_data_filtered[, 'SOJ_Class_NoMMS'] == 'AGGRAVATING'
#   #
#   #  labels_pos <- c(labels_pos,gi_data_filtered[, 'SOJ_Class_MMS'] == 'ALLEVIATING')
#   #  labels_neg <- c(labels_neg,gi_data_filtered[, 'SOJ_Class_MMS'] == 'AGGRAVATING')
#   
#   
#   if (draw_cutoffs == T) {
#     if (is.null(cutoffs_drawn)) {
#       abline(v = -best_neg,
#              lty = 3,
#              lwd = 0.7)
#       abline(v = best_pos, lty = 3, lwd = 0.7)
#     } else{
#       abline(v = -cutoffs_drawn[1],
#              lty = 3,
#              lwd = 0.7)
#       abline(v = cutoffs_drawn[2],
#              lty = 3,
#              lwd = 0.7)
#     }
#   }
#   legend(
#     xlims[1],
#     35,
#     legend = c(control_name, condition_name, 'Combined'),
#     fill = c('blue', 'red', 'black')
#   )
#   
#   return(c(best_neg, best_pos))
# }



##Legacy function
#' #' Plots St. Onge et al validation as a function of the internal FDR estimates
#' #'
#' #' @param gi_data input genetic interaction table
#' #' @param control_name condition which corresponds to 'NoDrug' in the St onge data
#' #' @param condition_name condition which corresponds to 'MMS' in the St onge data
#' #' @param fdr_prefix 
#' #' @param z_prefix 
#' #' @param fdr_cutoff 
#' #' @param xlims plot x axis limits
#' #' @param gi_prefix p
#' #' @param metr performance metric string, must be a valid ROCR metric
#' #' @param ylab plot y axis limits
#' #' @param draw_cutoffs boolean, whether to draw lines at cutoffs 
#' #' @param gi_cutoff_pos 
#' #' @param gi_cutoff_neg 
#' #' @param cutoffs_drawn 
#' #' @param xlab plot x axis label


#' Plot two ROCR metrics simultaneously as a function of GIS cutoff
#' (defaults to precision and recall)
#'
#' @param gi_data input genetic interaction table
#' @param control_name condition name in gi_data which corresponds to 'NoDrug' in the St onge data
#' @param condition_name condition name in gi_data which corresponds to 'MMS' in the St onge data
#' @param fdr_prefix pasted with control_name and condition_name to find the FDR columns
#' @param z_prefix pasted with control_name and condition_name to find the Z columns
#' @param gi_prefix pasted with control_name and condition_name to find the GIS columns
#' @param fdr_cutoff strains above this FDR are filtered out when drawing the plot
#' @param metr1 first metric to be plotted as function of GIS. Must be valid ROCR performance string. defaults to 'prec'
#' @param metr2 second metric to be plotted as function of GIS. Must be valid ROCR performance string. defaults to 'rec'
#' @param metr1_name plot name for first metric
#' @param metr2_name plot name for second metric
#' @param col1 line colour for drawing first metric performance
#' @param col2 line colour for drawing second metric performance
#' @param xlims x limits for plot
#' @param ylims y limits for plot
#' @param xlab x label for plot
#' @param ylab1 left y label for plot
#' @param ylab2 right y label for plot
#' @param draw_cutoffs boolean, whether draw lines at GIS cutoffs in cutoffs_drawn
#' @param cutoffs_drawn where to draw lines and summarize perforamnce
#' @param legend_pos where to position the legend
#'
#' @return NULL, makes a plot
#'
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
                               col1 = '#C28202',
                               col2 = '#7C32C2',
                               xlims = c(-0.1, 0.1),
                               ylims = c(0,100),
                               xlab = 'GIS Cutoff',
                               ylab1 = 'Precision (%)',
                               ylab2 = 'Recall (%)',
                               draw_cutoffs = T,
                               legend_pos = c(0,30),
                               cutoffs_drawn = c(-0.07,0.07)) {
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
    
    par(mar = c(5,4.5,2,4.5))
    par(xpd = F)
    plot(
      cutoff[signs[[i]](cutoff,0)],
      perf1[signs[[i]](cutoff,0)],
      type = 'l',
      xlim = c(0, xlims[i]),
      ylim = ylims,
      ylab = '',
      xlab = '',
      col = col1,
      lwd = 2,
      main = titles[i],
      axes = F,
      cex.main = 1.3
    )
    par(las = 1)
    axis(side = 4)
    axis(side = 2)
    par(las = 3)
    axis(side = 1)
    
    mtext(side = 4, line = 2.5, ylab2,cex = 1.3)
    mtext(side = 2, line = 2.5, ylab1,cex = 1.3)
    par(las = 1)
    mtext(side = 1, line = 3.5, xlab,cex = 1.3)
    
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
           legend = c(metr1_name, metr2_name),
           col = c(col1, col2),
           pch = 15,
           cex = 1.1)
    if (!is.null(cutoffs_drawn)) {
      abline(v = cutoffs_drawn[i], lty = 3)
      
      p1 <- perf1[which.min(abs(cutoff - cutoffs_drawn[i]))]
      p2 <- perf2[which.min(abs(cutoff - cutoffs_drawn[i]))]
      
      p1_text <- paste(c(metr1_name,': ',format(p1,digits=2),'%'),collapse='')
      p2_text <- paste(c(metr2_name,': ',format(p2,digits=2),'%'),collapse='')
      
      text(cutoffs_drawn[i] - 0.05*cutoffs_drawn[i],p1 + 0.05*p1,p1_text,adj=c(1,0),cex=1.2)
      text(cutoffs_drawn[i] - 0.05*cutoffs_drawn[i],p2 - 0.05*p2,p2_text,adj=c(1,1),cex=1.2)
    }
  }
}
