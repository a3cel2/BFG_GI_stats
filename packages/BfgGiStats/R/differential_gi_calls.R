#' Differential Genetic Interaction analysis
#'
#' @param gi_data processed genetic interaction data
#' @param fdr_cutoff false discovery rate cutoff for calling significant GI changes
#' @param require_sign_change require different GI classifications in order to call a differential interaction; True or False
#' @param nn_pair_type null distribution" of GI scores used to
#' compute FDR.  Either 'broad' for all interactions with neutral
#' pairs or 'narrow' for only neutral-neutral pairs
#' @param make_plots make histogram of Z scores while executing?
#'
#' @return a data frame with differential genetic interaction comparisons, conditions pairs sorted in alphabetical order
differential_gi_analysis <- function(gi_data,
                                     fdr_cutoff = 0.05,
                                     delta_gi_cutoff = 0,
                                     require_sign_change = T,
                                     nn_pair_type = 'broad',
                                     make_plots = F){
  gi_data <-
    dplyr::filter(gi_data, Remove_by_Chromosomal_distance_or_SameGene == 'no')
  
  
  if(nn_pair_type == 'broad'){
    #Can define all non DNA_repair - DNA_repair as neutral
    nn_pairs <-
      gi_data$Type_of_gene_x != 'DNA_repair' | gi_data$Type_of_gene_y != 'DNA_repair'
  }else if(nn_pair_type == 'narrow'){
    #Or just neutral-neutral pairs
    nn_pairs <-
      gi_data$Type_of_gene_x == 'Neutral' &
      gi_data$Type_of_gene_y == 'Neutral'
  }
  
  ddr_pairs <-
    gi_data$Type_of_gene_x == 'DNA_repair' &
    gi_data$Type_of_gene_y == 'DNA_repair'
  
  
  z_cols <- grep('Class',grep('Z_GIS',colnames(gi_data),val=T),invert=T,val=T)
  
  
  conditions <- sapply(strsplit(grep('^GIS',colnames(gi_data),val=T),split='\\.'),function(x){x[2]})
  
  
  ret_list <- list()
  ret_df <- c()
  
  for(condition1 in conditions) {
    ret_list[[condition1]] <- list()
    for (condition2 in conditions) {
      if (condition1  > condition2) {
        
        
        
        z_name1 <- paste(c('Z_GIS_xy', condition1), collapse = '.')
        z_name2 <- paste(c('Z_GIS_xy', condition2), collapse = '.')
        
        gi_name1 <- paste(c('GIS_xy', condition1), collapse = '.')
        gi_name2 <- paste(c('GIS_xy', condition2), collapse = '.')
        
        gi_err_name1 <- paste(c('SE_GIS_xy', condition1), collapse = '.')
        gi_err_name2 <- paste(c('SE_GIS_xy', condition2), collapse = '.')
        
        
        z_class_name1 <- paste(c(z_name1,'Class'),collapse='_')
        z_class_name2 <- paste(c(z_name2,'Class'),collapse='_')
        
        #nn_pairs <- (gi_data[,z_class_name1] == gi_data[,z_class_name2])
        
        
        #nn_pairs <- (gi_data[,z_class_name1] == 'NEUTRAL' & gi_data[,z_class_name2] == 'NEUTRAL')# | (gi_data$Type_of_gene_i != 'DNA_repair' | gi_data$Type_of_gene_y != 'DNA_repair')
        nn_pairs <- (gi_data$Type_of_gene_x != 'DNA_repair' | gi_data$Type_of_gene_y != 'DNA_repair')
        
        nn_z_scores_cond <-
          (gi_data[nn_pairs, gi_name1] - gi_data[nn_pairs, gi_name2])/sqrt(gi_data[nn_pairs, gi_err_name1]^2 + gi_data[nn_pairs, gi_err_name2]^2)
        nn_gi_scores_cond <-
          gi_data[nn_pairs, gi_name1] - gi_data[nn_pairs, gi_name2]
        
         
        non_nn_z_scores_cond <-
          (gi_data[!nn_pairs, gi_name1] - gi_data[!nn_pairs, gi_name2])/sqrt(gi_data[!nn_pairs, gi_err_name1]^2 + gi_data[!nn_pairs, gi_err_name2]^2)
        non_nn_gi_scores_cond <-
          gi_data[!nn_pairs, gi_name1] - gi_data[!nn_pairs, gi_name2]
        
        
        all_z_scores_cond <-
          (gi_data[, gi_name1] - gi_data[, gi_name2])/sqrt(gi_data[, gi_err_name1]^2 + gi_data[, gi_err_name2]^2)
        all_gi_scores_cond <-
          gi_data[, gi_name1] - gi_data[, gi_name2]
        
        
        
        nn_scores_cond <- nn_z_scores_cond
        non_nn_scores_cond <- non_nn_z_scores_cond
        all_scores_cond <- all_z_scores_cond
        
        ##If plots of process desired
        if(make_plots){
          first_hist_test <-
            hist(
              all_z_scores_cond,
              breaks = 100,
              plot = F
            )
          second_hist_test <-
            hist(
              nn_z_scores_cond,
              breaks = first_hist_test$breaks,
              plot = F
            )
          
          hist(
            non_nn_z_scores_cond,
            breaks = first_hist_test$breaks,
            freq = F,
            border = NA,
            col = rgb(1, 0, 0, 0.5),
            xlab = bquote(Delta~ Z[GIS] ~ .(condition1) - ~.(condition2)),#expression(Delta~Z[GIS]),
            ylab = 'Density',
            main = '',#paste(c(condition1,condition2),collapse=' - '),
            ylim = c(0, max(
              c(first_hist_test$density, second_hist_test$density)
            ))
          )
          hist(
            nn_z_scores_cond,
            breaks = first_hist_test$breaks,
            freq = F,
            border = NA,
            col = rgb(0, 0, 0, 0.7),
            add = T
          )
          mu <- mean(nn_z_scores_cond)
          sigma <- sd(nn_z_scores_cond)
          
          ranges <- c(non_nn_z_scores_cond,nn_z_scores_cond)
          line_range <- seq(min(ranges),max(ranges),by = 0.05)
          
          lines(line_range, dnorm(line_range, mean = mu, sd = sigma), lwd = 0.5)
        }
        
        
        #Hacky loop
        unequal_tests <- c(`>=`, `<=`)
        precision_list <- lapply(unequal_tests, function(unequal_test) {
          mu <- mean(nn_scores_cond)
          sigma <- sd(nn_scores_cond)
          sapply(1:length(all_scores_cond), function(i) {
            observed <-
              sum(unequal_test(all_scores_cond, all_scores_cond[i]))
            if (all.equal(unequal_test, `<=`) == T) {
              use_lower_tail <- T
            } else{
              use_lower_tail <- F
            }

            p_val <-
              pnorm(
                all_scores_cond[i],
                mean = mu,
                sd = sigma,
                lower.tail = use_lower_tail
              )
            
            #Two-tailed p values
            return(min(c(p_val,1 - p_val))*2)
          })
        })
        p_values_pos <- precision_list[[1]]
        p_values_neg <- precision_list[[2]]
        comb_p_values <- cbind(p_values_pos, p_values_pos)
        
        #Return positive or negative p-value based on the nominal sign
        #of the interaction
        p_values <- sapply(1:nrow(comb_p_values), function(i) {
          if (all_scores_cond[i] < 0) {
            return(p_values_neg[i])
          } else{
            return(p_values_pos[i])
          }
        })
        
       fdr_scores <- qvalue(p_values)$q
        
        if(require_sign_change){
          sig <- (fdr_scores < fdr_cutoff) & (gi_data[,z_class_name1] != gi_data[,z_class_name2]) & (abs(all_gi_scores_cond) >= delta_gi_cutoff)
        }else{
          sig <- (fdr_scores < fdr_cutoff) & (abs(all_gi_scores_cond) >= delta_gi_cutoff)
        }
        
        #Used to report only for DDR pairs, now reporting all
        ddr_data <- gi_data[,]#gi_data[ddr_pairs,]
        
        sig_names <- apply(ddr_data[which(sig),1:2],1,function(x){paste(x,collapse='_')})
        
        ret_list[[condition1]][[condition2]] <- sig_names
        
        
        
        retvec <- cbind(ddr_data[which(sig),c("Gene_x","Gene_y")],
                        rep(condition1,sum(sig)),
                        rep(condition2,sum(sig)),
                        ddr_data[which(sig),c("Type_of_gene_x","Type_of_gene_y")],
                        ddr_data[which(sig),c(z_name1,z_name2)],
                        ddr_data[which(sig),c(gi_name1,gi_name2)],
                        ddr_data[which(sig),c(z_class_name1,z_class_name2)],
                        all_z_scores_cond[which(sig)],
                        all_gi_scores_cond[which(sig)],
                        fdr_scores[which(sig)])
        
        colnames(retvec) <-
          c(
            'Gene_x',
            'Gene_y',
            'Condition1',
            'Condition2',
            'Type_of_gene_x',
            'Type_of_gene_y',
            'Z_Condition1',
            'Z_Condition2',
            'GI_Condition1',
            'GI_Condition2',
            'Class_Condition1',
            'Class_Condition2',
            'DeltaZ',
            'DeltaGI',
            'DeltaZ_FDR'
          )
        
        ret_df <- rbind(ret_df,retvec)
      }
    }
  }
  return(ret_df)
}

# A histogram function for diagnosis, no longer used
# differential_calls_histogram <- function(differential_calls){
#   orig_hist <-
#     hist(
#       differential_calls$DeltaZ,
#       breaks = 300,
#       freq = F,
#       border = NA,
#       col = rgb(0, 0, 0, 0.7),
#       xlim = c(-20, 20),
#       main = '',
#       xlab = '∆Z',
#       ylim=c(0,0.3)
#     )
#   
#   hist(
#     as.matrix(
#       filter(
#         differential_calls,
#         Class_Condition1 == "Expected" &
#           Class_Condition2 == "Expected"
#       )$DeltaZ
#     ),
#     breaks = orig_hist$breaks,
#     add = T,
#     col = rgb(1, 0, 0, 0.5),
#     freq = FALSE,
#     border = NA
#   )
#   
#   hist(
#     as.matrix(
#       filter(
#         differential_calls,
#         Type_of_gene_x == 'Neutral' |
#           Type_of_gene_y == 'Neutral'
#       )$DeltaZ
#     ),
#     add = T,
#     col = rgb(0, 0, 1, 0.5),
#     freq = F,
#     border = NA,
#     breaks = orig_hist$breaks
#   )
#   legend('topleft',legend=c('All pairs','non-DDR Pairs','neutral GI pairs'),col=c('grey30','blue','red'),pch=15)
# }