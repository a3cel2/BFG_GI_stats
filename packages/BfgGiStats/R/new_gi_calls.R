devtools::use_package('dplyr')
#devtools::use_package('ExtDist')

#' Adds p-value columns to the genetic interaction data
#'
#' @param gi_data input genetic interaction table
#' @param nn_pair_type "null distribution" of GI scores used to
#' compute p-values. Either 'broad' for all interactions with neutral
#' pairs or 'narrow' for only neutral-neutral pairs
#' @param z_grep_pattern grep pattern to find Z_GIS columns
#' @param gi_grep_pattern grep pattern to find GIS columns
#' @param make_plots boolean, whether to make a histogram of p-value distribution
#'
#' @return gi_data with p values added
add_p_values <- function(gi_data,
                         nn_pair_type = 'broad',
                         z_grep_pattern = '^Z_GIS',
                         gi_grep_pattern = '^GIS',
                         make_plots = F) {
  gi_data_full <- gi_data
  
  gi_data <-
    dplyr::filter(gi_data,
                  Remove_by_Chromosomal_distance_or_SameGene == 'no')
  
  ddr_pairs <-
    gi_data$Type_of_gene_x == 'DNA_repair' &
    gi_data$Type_of_gene_y == 'DNA_repair'
  
  #Neutral-Neutral pairs
  if (nn_pair_type == 'broad') {
    #Can define all non DNA_repair - DNA_repair as neutral
    nn_pairs <-
      gi_data$Type_of_gene_x != 'DNA_repair' &
      gi_data$Type_of_gene_y != 'DNA_repair'
  } else if (nn_pair_type == 'narrow') {
    #Or just neutral-neutral pairs
    nn_pairs <-
      gi_data$Type_of_gene_x == 'Neutral' &
      gi_data$Type_of_gene_y == 'Neutral'
  }
  
  gi_cols <- grep(gi_grep_pattern, colnames(gi_data), val = T)
  
  z_cols <-
    grep(
      'Class',
      grep(z_grep_pattern, colnames(gi_data), val = T),
      invert = T,
      val = T
    )
  
  z_scores <- gi_data[, z_cols]
  z_scores_full <- gi_data_full[, z_cols]
  
  nn_scores <- z_scores[nn_pairs, ]
  non_nn_scores <- z_scores[ddr_pairs, ]
  all_scores <- z_scores_full
  
  
  p_value_matrix <- c()
  
  for (condition in colnames(z_scores)) {
    nn_scores_cond <-  nn_scores[, condition]
    non_nn_scores_cond <- non_nn_scores[, condition]
    all_scores_cond_full <- z_scores_full[, condition]
    
    #Hacky way to make p-values for both positive and negative
    #in a loop, full intention now deprecated
    unequal_tests <- c(`>=`, `<=`)
    
    precision_list <- lapply(unequal_tests, function(unequal_test) {
      mu <- mean(nn_scores_cond)
      sigma <- sd(nn_scores_cond)
      
      sapply(1:length(all_scores_cond_full), function(i) {
        if (all.equal(unequal_test, `<=`) == T) {
          use_lower_tail <- T
        } else{
          use_lower_tail <- F
        }
        
        #Using normal distribuion to calculate p value
        p_value <-
          pnorm(
            all_scores_cond_full[i],
            mean = mu,
            sd = sigma,
            lower.tail = use_lower_tail
          )
        
        #Two-tailed
        p_value <- min(c(p_value, 1 - p_value)) * 2
        return(p_value)
      })
    })
    
    if (make_plots) {
      mu <- median(nn_scores_cond)
      sigma <- sd(nn_scores_cond)
      
      
      mu_non_nn <- median(non_nn_scores_cond)
      sigma_non_nn <- sd(non_nn_scores_cond)
      common_breaks <-
        hist(non_nn_scores_cond, breaks = 100, plot = F)$breaks
      
      
      y_limits1 <-
        hist(nn_scores_cond, breaks = common_breaks, plot = F)$density
      y_limits2 <-
        hist(non_nn_scores_cond, breaks = common_breaks, plot = F)$density
      
      
      cond_name <- strsplit(condition, split = '\\.')[[1]][2]
      
      line_range <- (-200:200) / 10
      hist(
        non_nn_scores_cond,
        breaks = common_breaks,
        main = bquote(Z[GIS] ~  distribution ~ .(cond_name)),
        xlab = expression(Z[GIS]),
        freq = F,
        add = F,
        axes = F,
        col = rgb(1, 0, 0, 0.7),
        border = NA,
        ylim = c(0, max(c(
          y_limits1, y_limits2
        ))),
        cex.lab = 1.5
      )
      
      
      hist(
        nn_scores_cond,
        breaks = common_breaks,
        main = condition,
        freq = F,
        add = T,
        border = NA,
        col = rgb(0, 0, 0, 0.7)
      )
      lines(line_range, dnorm(line_range, mean = mu, sd = sigma))
      
      
      box(bty = 'l')
      par(las= 3)
      axis(side = 1)
      par(las= 1)
      axis(side = 2)
    }
    p_values_pos <- precision_list[[1]]
    p_values_neg <- precision_list[[2]]
    
    comb_p_values <- cbind(p_values_pos, p_values_neg)
    
    #Return positive or negative p-values based on the nominal sign
    #of the interaction (loop deprecated, but still keeping here)
    p_values <- sapply(1:nrow(comb_p_values), function(i) {
      if (z_scores_full[i, condition] < 0) {
        return(p_values_neg[i])
      } else{
        return(p_values_pos[i])
      }
    })
    
    p_value_matrix <- cbind(p_value_matrix, p_values)
  }
  
  #Add the appropriate columns at the end and fill in the data
  drugs <- sapply(colnames(z_scores), function(name) {
    strsplit(name, split = '\\.')[[1]][2]
  })
  p_value_names <- sapply(drugs, function(drug) {
    paste(c('P.neutral_xy.', drug), collapse = '')
  })
  
  colnames(p_value_matrix) <- p_value_names
  
  gi_data_full[, p_value_names] <- p_value_matrix
  
  return(gi_data_full)
  
}

#' Update genetic interaction calls based on either Z metrics or internal FDR
#'
#' @param gi_data input genetic interaction table
#' @param z_score_cols which columns correspond to Z scores - default is to use grep to find automatically
#' @param fdr_cols which columns correspond to internal FDRs - default is to use grep to find automatically
#' @param gi_cols which columns correspond to GI scores - default is to use grep to find automatically
#' @param z_cutoff_neg Z_GIS  cutoff for negative interactions
#' @param z_cutoff_pos  Z_GIS cutoff for positive interactions
#' @param fdr_cutoff_pos FDR cutoff for positive interactions
#' @param fdr_cutoff_neg FDR cutoff for negative interactions
#' @param gi_cutoff_pos GIS cutoff for positive interactions
#' @param gi_cutoff_neg GIS cutoff for negative interactions
#
#' @return a version of gi_data with GI calls updated
update_calls <- function(gi_data,
                         z_score_cols = NULL,
                         fdr_cols = NULL,
                         gi_cols = NULL,
                         fdr_cutoff_pos = 0.05,
                         fdr_cutoff_neg = 0.05,
                         #use_z = F,
                         z_cutoff_neg = 0,
                         z_cutoff_pos = 0,
                         gi_cutoff_neg = 0,
                         gi_cutoff_pos = 0) {
  if (is.null(z_score_cols)) {
    z_score_cols <-
      grep('Class',
           grep('^Z_GIS', colnames(gi_data), val = T),
           invert = T,
           val = T)
  }
  if (is.null(gi_cols)) {
    gi_score_cols <-
      grep('^GIS', colnames(gi_data), val = T)
  }
  
  #Add columns for making classifications
  z_scores <- gi_data[, z_score_cols]
  colnames(z_scores) <-
    sapply(colnames(z_scores), function(name) {
      sprintf('%s_Class', name)
    })
  
  gi_data <- cbind(gi_data, z_scores)
  
  z_class_cols <-
    grep('Class',
         grep('^Z_GIS', colnames(gi_data), val = T),
         invert = F,
         val = T)
  
  
  if (is.null(fdr_cols)) {
    fdr_cols <- grep('^FDR', colnames(gi_data), val = T)
  }
  for (i in 1:length(z_score_cols)) {
    gi_data[, z_class_cols[i]] <- "Expected"
    
    neg_crit <- gi_data[, fdr_cols[i]] <= fdr_cutoff_neg &
      gi_data[, z_score_cols[i]] < z_cutoff_neg &
      gi_data[, gi_score_cols[i]] < gi_cutoff_neg
    
    pos_crit <- gi_data[, fdr_cols[i]] <= fdr_cutoff_pos &
      gi_data[, z_score_cols[i]] > z_cutoff_pos &
      gi_data[, gi_score_cols[i]] > gi_cutoff_pos
    
    gi_data[neg_crit, z_class_cols[i]] <- "Negative"
    gi_data[pos_crit, z_class_cols[i]] <- "Positive"
  }
  
  return(gi_data)
}
