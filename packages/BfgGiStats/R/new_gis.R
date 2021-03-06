devtools::use_package('dplyr')


strain_averaging_function <- function(x){
  return(mean(x))
}

#' Update BFG-GI genetic interaction data
#'
#' @param gi_data input genetic interaction table
#' @param pseudocount integer, pseudocount to add in the calculations
#' @param count_data_grep_pattern grep pattern used to find count columns
#' @param g_wt_vec how many generations each pool grew, in the same order as the count columns
#'
#' @return an updated gi_data
update_gis <- function(gi_data,
                       #Old parameter, to use only well-measured pairs, no longer needed
                       #well_measured_cutoff = 100,
                       #This defines the number of generations each pool grew
                       pseudocount = 0.5,
                       count_data_grep_pattern = '^C_',
                       g_wt_vec) {
  #Define special pairs
  same_same <-
    sapply(gi_data$Barcode_x, function(x) {
      strsplit(x, split = '_')[[1]][1]
    }) == sapply(gi_data$Barcode_y, function(x) {
      strsplit(x, split = '_')[[1]][1]
    })
  
  nn_pairs <-
    gi_data$Type_of_gene_x == 'Neutral' &
    gi_data$Type_of_gene_y == 'Neutral' &
    gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
  
  non_nn_unlinked_pairs <-
    gi_data$Type_of_gene_x != 'Neutral' &
    gi_data$Type_of_gene_y != 'Neutral' &
    gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
  
  #Add a pseudocount of 1
  count_data <- gi_data[, grep(count_data_grep_pattern, colnames(gi_data))] + pseudocount
  
  condition_sums <- apply(count_data, 2, sum)
  freq_data <- count_data / condition_sums
  
  
  #See methods for rationale of these metrics
  f_xy_data <- count_data / condition_sums
  r_xy_data <- (f_xy_data / f_xy_data[, 1])[, 2:ncol(f_xy_data)]
  g_xy_data <- t(t(log2(r_xy_data)) + g_wt_vec)
  
  #Count error estimated by the poisson distribution
  r_xy_error <-
    r_xy_data * condition_sums[2:length(condition_sums)]/condition_sums[1] *
    sqrt((sqrt(count_data[, 2:ncol(f_xy_data)]) /count_data[, 2:ncol(f_xy_data)])^2 +
           (sqrt(count_data[, 1]) / count_data[, 1])^2)
  
  #Error propagation to log
  g_xy_error <- abs(r_xy_error/(r_xy_data*log(2)))
  
  #wt fitness estimate and error
  g_xy_wt <- apply(g_xy_data[nn_pairs,],2,strain_averaging_function)
  
  g_wt_error <- apply(g_xy_data[nn_pairs,],2,function(x){sd(x)})
  
  #Normalize fitness by wildtype
  w_xy_data <- t(t(g_xy_data) / g_xy_wt)
  
  barcodes <- unique(unlist(gi_data[, 1:2]))

  
  #Using mean now instead of median because
  #otherwise standard deviation doesn't make sense
  
  
  w_xy_single_barcodes <- t(sapply(barcodes, function(barcode) {
    gene <- barcode
    #gene <- strsplit(barcode,split='_')[[1]][1]
    #gene <- paste(c(gene,'donor'),collapse='_')
    
    criteria <-
      #gi_data[, 1] == barcode &
      (grepl(gene, gi_data[,1]) &
      gi_data$Type_of_gene_y == 'Neutral') |
      #gi_data[, 2] == barcode &
      (grepl(gene, gi_data[,2]) &
      gi_data$Type_of_gene_x == 'Neutral')
    criteria <-
      criteria &
      gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
    
    #Used to filter by counts, no need with error model
    #criteria <-
    #  criteria &
    #  gi_data$C_xy.HetDipl >= well_measured_cutoff
    
    return(apply(w_xy_data[criteria, ], 2, strain_averaging_function))
    #return(apply(w_xy_data[criteria, ], 2, median))
    #return(w_xy_data[criteria, ])
  }))
  
  g_xy_single_barcodes <- t(sapply(barcodes, function(barcode) {
    gene <- barcode
    #gene <- strsplit(barcode,split='_')[[1]][1]
    #gene <- paste(c(gene,'donor'),collapse='_')
    
    criteria <-
      (grepl(gene, gi_data[,1]) &
         gi_data$Type_of_gene_y == 'Neutral') |
      (grepl(gene, gi_data[,2]) &
         gi_data$Type_of_gene_x == 'Neutral')
      #gi_data[, 1] == barcode &
      #gi_data$Type_of_gene_y == 'Neutral' |
      #gi_data[, 2] == barcode & gi_data$Type_of_gene_x == 'Neutral'
    criteria <-
      criteria &
      gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
    #Used to filter by counts, no need with error model
    #criteria <-
    #  criteria &
    #  gi_data$C_xy.HetDipl >= well_measured_cutoff
    
    return(apply(g_xy_data[criteria, ], 2, strain_averaging_function))
  }))
  
  
  #Using standard error of mean as error
  g_xy_error_single_barcodes <- t(sapply(barcodes, function(barcode) {
    gene <- barcode
    #gene <- strsplit(barcode,split='_')[[1]][1]
    #gene <- paste(c(gene,'donor'),collapse='_')
    
    criteria <-
      (grepl(gene, gi_data[,1]) &
         gi_data$Type_of_gene_y == 'Neutral') |
      (grepl(gene, gi_data[,2]) &
         gi_data$Type_of_gene_x == 'Neutral')
    #criteria <-
    #  gi_data[, 1] == barcode &
    #  gi_data$Type_of_gene_y == 'Neutral' |
    #  gi_data[, 2] == barcode & gi_data$Type_of_gene_x == 'Neutral'
    criteria <-
      criteria &
      gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
    #Used to filter by counts, no need with error model
    #criteria <-
    #  criteria &
    #  gi_data$C_xy.HetDipl >= well_measured_cutoff
    return(apply(g_xy_data[criteria, ], 2, function(x){sd(x)}))#sqrt(length(x))}))
  }))
  
  
  bc1 <- gi_data$Barcode_x
  bc2 <- gi_data$Barcode_y
  
  gis <- t(sapply(1:nrow(w_xy_data), function(i) {
    w_x <- w_xy_single_barcodes[bc1[i], ]
    w_y <- w_xy_single_barcodes[bc2[i], ]
    
    w_xy <- w_xy_data[i, ]
    
    #Quick fix
    #Set fitness to 0 if estimated less than 0
    w_x[w_x < 0] <- 0
    w_y[w_y < 0] <- 0
    w_xy[w_xy < 0] <- 0
    
    gis <- (w_xy) - (w_x * w_y)
    
    #Alternative log-model of genetic interactions, not used
    #gis <- log2((w_xy + 0.01)/((w_x * w_y) + 0.01))
    
    #For debugging
    #if(sum(is.na(gis) > 0)){
    #  print(w_xy)
    #  print(w_x)
    #  print(w_y)
    #  #print(c(w_xy,w_x,w_y))
    #}
    
    return(gis)
  }))
  
  
  
  gi_uncertainty <- t(sapply(1:nrow(w_xy_data), function(i) {
    
    g_xy <- g_xy_data[i, ]
    g_x <- g_xy_single_barcodes[bc1[i], ]
    g_y <- g_xy_single_barcodes[bc2[i], ]
    g_wt <- g_xy_wt
    
    
    g_xy_error <- g_xy_error[i, ]
    g_x_error <- g_xy_error_single_barcodes[bc1[i], ]
    g_y_error <- g_xy_error_single_barcodes[bc2[i], ]
    g_wt_error <- g_wt_error
    
    gis <- (g_xy*g_wt)/(g_x*g_y)
    
    
    # numerator <- g_xy * g_wt
    # numerator_error <- numerator*sqrt((g_xy_error/g_xy)^2 + (g_wt_error/g_wt)^2)
    # 
    # denominator <- g_x * g_y
    # denominator_error <- denominator*sqrt((g_x_error/g_x)^2 + (g_y_error/g_y)^2)
    # 
    # ratio <- numerator / denominator
    # 
    # ratio_error <- ratio * sqrt((numerator_error/numerator)^2 + (denominator_error/denominator)^2)
    # 
    # log_ratio_error <- ratio_error/(ratio*log(2))
    # 
    # return(log_ratio_error)
    
    numerator <- g_xy - g_x - g_y

    gis <- numerator/g_wt

    numerator_error <- sqrt(g_xy_error^2 + g_x_error^2 + g_y_error^2)

    gi_error <- abs(gis)*sqrt((numerator_error/numerator)^2 + (g_wt_error/g_wt)^2)

    return(gi_error)
    
    
    
  }))
  
  
  
  
  
  #gis_global <<- gis
  #stop()
  
  #Modify data to replace old definition with new definitions
  colnames(gis) <- sapply(colnames(gis),function(name){
    gsub('^C_','GIS_',name)
  })
  
  
  gi_data <- cbind(gi_data,gis)
  
  
  z_scores <- gis / gi_uncertainty
  
  
  w_x_data <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_x <- g_xy_single_barcodes[bc1[i], ]
    g_wt <- g_xy_wt
    w_x <- g_x/g_wt
    w_x[w_x < 0] <- 0
    
    return(w_x)
  }))
  
  w_y_data <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_y <- g_xy_single_barcodes[bc2[i], ]
    g_wt <- g_xy_wt
    w_y <- g_y/g_wt
    w_y[w_y < 0] <- 0
    
    return(w_y)
  }))
  
  
  w_x_error <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_x <- g_xy_single_barcodes[bc1[i], ]
    g_wt <- g_xy_wt
    g_x_error <- g_xy_error_single_barcodes[bc1[i], ]
    g_wt_error <- g_wt_error
    
    w_x <- g_x/g_wt
    err <- abs(w_x)*sqrt((g_x_error/g_x)^2 + (g_wt_error/g_wt)^2)
    
    return(err)
  }))
  
  w_y_error <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_y <- g_xy_single_barcodes[bc2[i], ]
    g_wt <- g_xy_wt
    g_y_error <- g_xy_error_single_barcodes[bc2[i], ]
    g_wt_error <- g_wt_error
    
    
    
    
    w_y <- g_y/g_wt
    err <- abs(w_y)*sqrt((g_y_error/g_y)^2 + (g_wt_error/g_wt)^2)
    return(err)
  }))
  
  w_xy_error <- t(sapply(1:nrow(w_xy_data), function(i) {
    g_xy <- g_xy_data[i, ]
    g_wt <- g_xy_wt
    
    g_xy_error <- g_xy_error[i, ]
    g_wt_error <- g_wt_error
    
    
    
    w_xy <- g_xy/g_wt
    
    err <- abs(w_xy)*sqrt((g_xy_error/g_xy)^2 + (g_wt_error/g_wt)^2)
    
    #err <- abs(w_y)*sqrt((g_y_error/g_y)^2 + (g_wt_error/g_wt)^2)
    return(err)
  }))
  
  
  
  print(colnames(w_xy_error))
  # Set Some column names
  colnames(w_x_data) <- colnames(w_x_data) %>% sapply(function(name){
    gsub('^C_xy.','W_x.', name)
  })
  colnames(w_y_data) <- colnames(w_y_data) %>% sapply(function(name){
    gsub('^C_xy.','W_y.', name)
  })
  colnames(w_xy_data) <- colnames(w_xy_data) %>% sapply(function(name){
    gsub('^C_xy.','W_xy.', name)
  })
  
  colnames(w_x_error) <- colnames(w_x_error) %>% sapply(function(name){
    gsub('^C_xy.','W_x_SE.', name)
  })
  colnames(w_y_error) <- colnames(w_y_error) %>% sapply(function(name){
    gsub('^C_xy.','W_y_SE.', name)
  })
  colnames(w_xy_error) <- colnames(w_x_error) %>% sapply(function(name){
    gsub('^W_x_SE.','W_xy_SE.', name)
  })
  
  
  
  
  
  
  #print(colnames(w_x_data))
  #
  #print(colnames(w_y_))
  print(colnames(w_xy_error))
  
  
  
  #colnames()
  
  
  colnames(z_scores) <- sapply(colnames(z_scores),function(name){
    gsub('^GIS_','Z_GIS_',name)
  })
  
  #Not returned in latest version
  colnames(gi_uncertainty) <- colnames(gis)
  colnames(gi_uncertainty) <- sapply(colnames(gi_uncertainty),function(name){
    gsub('^GIS_','SE_GIS_',name)
  })
  
  
  #Correct fitness
  w_xy_data[w_xy_data < 0] <- 0
  
  
  gi_data <- cbind(gi_data,
                   
                   w_x_data,
                   w_y_data,
                   w_xy_data,
                   
                   w_x_error,
                   w_y_error,
                   w_xy_error,
                   
                   gi_uncertainty,
                   z_scores)
  
  #Old way of updating the data
  #gi_data[, grep('^GI', colnames(gi_data))] <-
  #  gis
  #z_cols <-
  #  grep('Class',
  #       grep('^Z_GIS', colnames(gi_data), val = T),
  #       invert = T,
  #       val = T)
  #gi_data[, z_cols] <- gis / gi_uncertainty
  
  return(gi_data)
}


update_error_model_by_replicates <- function(gi_data_combined){
  
  
  #This is all to make sure we keep same pairs in both
  r1_data <- dplyr::filter(gi_data_combined, Sample == 'R1')
  r2_data <- dplyr::filter(gi_data_combined, Sample == 'R2')
  
  pairs_r1 <- apply(r1_data[,c('Barcode_x','Barcode_y')],1,function(x){paste(x,collapse='_')})
  pairs_r2 <- apply(r2_data[,c('Barcode_x','Barcode_y')],1,function(x){paste(x,collapse='_')})
  
  r1_data <- cbind(r1_data,pairs_r1)
  r2_data <- cbind(r2_data,pairs_r2)
  colnames(r1_data)[ncol(r1_data)] <- 'pair'
  colnames(r2_data)[ncol(r2_data)] <- 'pair'
  
  merged_data <- merge(r1_data,r2_data,by='pair')
  
  r1_unmerged_data <- merged_data[,grep('.x$',colnames(merged_data),val=T)]
  r2_unmerged_data <- merged_data[,grep('.y$',colnames(merged_data),val=T)]
  
  colnames(r1_unmerged_data) <- colnames(gi_data_combined)
  colnames(r2_unmerged_data) <- colnames(gi_data_combined)
  
  
  conditions <- sapply(grep('^GIS',colnames(r1_data),val=T),function(x){strsplit(x,split='\\.')[[1]][2]})
  
  #Update double mutant error
  for(condition in conditions){
    dm_name <- paste(c('W_xy',condition),collapse='.')
    err_name <- paste(c('W_xy_SE',condition),collapse='.')
    
    r1_unmerged_data[,err_name] <- median(abs(r1_unmerged_data[,dm_name] - r2_unmerged_data[,dm_name]))
    r2_unmerged_data[,err_name] <- r1_unmerged_data[,err_name]
  }
  
  #Update Z values
  r1_old <- r1_unmerged_data
  r2_old <- r2_unmerged_data
  
  for(condition in conditions){
    
    w_x_col <- paste(c('W_x',condition),collapse='.')
    w_y_col <- paste(c('W_y',condition),collapse='.')
    
    w_x_err_col <- paste(c('W_x_SE',condition),collapse='.')
    w_y_err_col <- paste(c('W_y_SE',condition),collapse='.')
    
    w_xy_err_col <- paste(c('W_xy_SE',condition),collapse='.')
    
    gis_col <- paste(c('GIS_xy',condition),collapse='.')
    gis_err_col <- paste(c('SE_GIS_xy',condition),collapse='.')
    
    z_col <- paste(c('Z_GIS_xy',condition),collapse='.')
    
    
    #
    w_x_r1 <- r1_unmerged_data[ ,w_x_col]
    w_y_r1 <- r1_unmerged_data[ ,w_y_col]
    w_x_err_r1 <- r1_unmerged_data[ ,w_x_err_col]
    w_y_err_r1 <- r1_unmerged_data[ ,w_y_err_col]
    w_xy_err_r1 <- r1_unmerged_data[ ,w_xy_err_col]
    
    
    w_x_r2 <- r2_unmerged_data[,w_x_col]
    w_y_r2 <- r2_unmerged_data[,w_y_col]
    w_x_err_r2 <- r2_unmerged_data[ ,w_x_err_col]
    w_y_err_r2 <- r2_unmerged_data[ ,w_y_err_col]
    w_xy_err_r2 <- r2_unmerged_data[ ,w_xy_err_col]
    
    w_x_w_y_err_r1 <- abs(w_x_r1*w_y_r1)*sqrt((w_x_err_r1/w_x_r1)^2 + (w_y_err_r1/w_y_r1)^2)
    w_x_w_y_err_r2 <- abs(w_x_r2*w_y_r2)*sqrt((w_x_err_r2/w_x_r2)^2 + (w_y_err_r2/w_y_r2)^2)
    
    updated_gis_error_r1 <- sqrt(w_x_w_y_err_r1^2 + w_xy_err_r1^2)
    updated_gis_error_r2 <- sqrt(w_x_w_y_err_r2^2 + w_xy_err_r2^2)
    
    #updated_gis_error_r1 <- median(abs(r1_unmerged_data[,gis_col] - r2_unmerged_data[,gis_col]))
    #updated_gis_error_r2 <- median(abs(r1_unmerged_data[,gis_col] - r2_unmerged_data[,gis_col]))
    
    r1_unmerged_data[,gis_err_col] <- updated_gis_error_r1
    r2_unmerged_data[,gis_err_col] <- updated_gis_error_r2
    
    r1_unmerged_data[,z_col] <- r1_unmerged_data[,gis_col]/r1_unmerged_data[,gis_err_col]
    r2_unmerged_data[,z_col] <- r2_unmerged_data[,gis_col]/r2_unmerged_data[,gis_err_col]
    
  }
  
  return(rbind(r1_unmerged_data,r2_unmerged_data))
}