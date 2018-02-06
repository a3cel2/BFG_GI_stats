devtools::use_package('dplyr')

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
                       pseudocount = 1e-200,
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
  #g_xy_wt <- apply(g_xy_data[nn_pairs,],2,mean)
  g_xy_wt <- apply(g_xy_data[nn_pairs,],2,median)
  g_wt_error <- apply(g_xy_data[nn_pairs,],2,function(x){sd(x)})#sqrt(length(x))})
  
  #Normalize fitness by wildtype
  w_xy_data <- t(t(g_xy_data) / g_xy_wt)
  
  barcodes <- unique(unlist(gi_data[, 1:2]))
  
  g_xy_single_barcodes <- sapply(barcodes, function(barcode) {
    criteria <-
      gi_data[, 1] == barcode &
      gi_data$Type_of_gene_y == 'Neutral' |
      gi_data[, 2] == barcode & gi_data$Type_of_gene_x == 'Neutral'
    criteria <-
      criteria &
      gi_data$Remove_by_Chromosomal_distance_or_SameGene == 'no'
    #Used to filter by counts, no need with error model
    #criteria <-
    #  criteria &
    #  gi_data$C_xy.HetDipl >= well_measured_cutoff
    
    return(g_xy_data[criteria, ])
  })
  
  
  
  bc1 <- gi_data$Barcode_x
  bc2 <- gi_data$Barcode_y
  
  gis_full <- lapply(1:nrow(g_xy_data), function(i) {
    print(i)
    g_x <- g_xy_single_barcodes[[bc1[i]]]
    g_y <- g_xy_single_barcodes[[bc2[i]]]
    g_wt <- g_xy_data[nn_pairs,]
    
    g_xy <- g_xy_data[i, , drop =F]
    
    #Quick fix
    #Set fitness to 0 if less than 0
    g_x[g_x < 0] <- 0
    g_y[g_y < 0] <- 0
    g_xy[g_xy < 0] <- 0
    
    
    
    gis <- sapply(1:ncol(g_xy_data),function(j){
      val_comb <- expand.grid(g_x[, j], g_y[, j], g_xy[, j], g_wt[, j])
      return(apply(val_comb,1,function(x){

       w_xy <- x[3]/x[4] 
       gis <- x[3]/x[4] - (x[1]*x[2])/(x[4]^2)
       w_x <- x[1]/x[4]
       w_y <- x[2]/x[4]

       return(list(c(gis,w_x,w_y,w_xy)))
      }))
      
      # return(sapply(1:1000,function(x){
      #   #Rep is to stop 'sample' from treating it as single number
      #   g_x_samp <- sample(rep(g_x[,j],2),1)
      #   g_y_samp <- sample(rep(g_y[,j],2),1)
      #   g_xy_samp <- sample(rep(g_xy[,j],2),1)
      #   g_wt_samp <- sample(rep(g_wt[,j],2),1)
      #   
      #   w_x <- g_x_samp/g_wt_samp
      #   w_y <- g_y_samp/g_wt_samp
      #   w_xy <- g_xy_samp/g_wt_samp
      #   gis <- w_xy - w_x*w_y
      #   
      #   return(list(c(gis,w_x,w_y,w_xy)))
      #   
      # }))
    })
    
    return(t(apply(gis,2,function(gi){
      
      gi_n <- sapply(gi,function(entry){
        entry[1]
      })
      
      w_x <- sapply(gi,function(entry){
        entry[2]
      })
      
      w_y <- sapply(gi,function(entry){
        entry[3]
      })
      
      w_xy <- sapply(gi,function(entry){
        entry[4]
      })
      
      gi_nominal <- median(gi_n)
      gi_error <- quantile(gi_n,probs=c(pnorm(1),pnorm(-1)))
      gi_error <- abs(gi_error[2] - gi_error[1])/2
      
      gi_p <- min(c(sum(gi_n > 0), sum(gi_n < 0)))/length(gi_n)
      gi_p <- 2*gi_p
      
      w_x_nominal <- median(w_x)
      w_y_nominal <- median(w_y)
      w_xy_nominal <- median(w_xy)
      
      w_x_error <- quantile(w_x,probs=c(pnorm(1),pnorm(-1)))
      w_y_error <- quantile(w_y,probs=c(pnorm(1),pnorm(-1)))
      w_xy_error <- quantile(w_xy,probs=c(pnorm(1),pnorm(-1)))
      
      w_x_error <- abs(w_x_error[2] - w_x_error[1])/2
      w_y_error <- abs(w_y_error[2] - w_y_error[1])/2
      w_xy_error <- abs(w_xy_error[2] - w_xy_error[1])/2
      
      return(c(gi_nominal,
               gi_p,
                    w_x_nominal,
                    w_y_nominal,
                    w_xy_nominal,
                    gi_error,
                    w_x_error,
                    w_y_error,
                    w_xy_error))
      
    })))
    
    
  })
  
  
  
  gis <- t(sapply(gis_full,function(gi){gi[,1]}))
  gi_uncertainty <- t(sapply(gis_full,function(gi){gi[,6]}))
  
  w_x_data <- t(sapply(gis_full,function(gi){gi[,3]}))
  w_y_data <- t(sapply(gis_full,function(gi){gi[,4]}))
  
  w_x_error <- t(sapply(gis_full,function(gi){gi[,7]}))
  w_y_error <- t(sapply(gis_full,function(gi){gi[,7]}))
  
  w_xy_error <- t(sapply(gis_full,function(gi){gi[,8]}))
  
  z_scores <- gis / gi_uncertainty
  
  colnames(gis) <- sapply(colnames(g_xy_data),function(name){
    gsub('^C_','GIS_',name)
  })
  
  
  gi_data <- cbind(gi_data,gis)
  
  
  
  # Set Some column names
  colnames(w_x_data) <- colnames(g_xy_data) %>% sapply(function(name){
    gsub('^C_xy.','W_x.', name)
  })
  colnames(w_y_data) <- colnames(g_xy_data) %>% sapply(function(name){
    gsub('^C_xy.','W_y.', name)
  })
  colnames(w_xy_data) <- colnames(g_xy_data) %>% sapply(function(name){
    gsub('^C_xy.','W_xy.', name)
  })
  
  colnames(w_x_error) <- colnames(g_xy_data) %>% sapply(function(name){
    gsub('^C_xy.','W_x_SE.', name)
  })
  colnames(w_y_error) <- colnames(g_xy_data) %>% sapply(function(name){
    gsub('^C_xy.','W_y_SE.', name)
  })
  colnames(w_xy_error) <- colnames(g_xy_data) %>% sapply(function(name){
    gsub('^C_xy.','W_xy_SE.', name)
  })
  
  colnames(z_scores) <- sapply(colnames(g_xy_data),function(name){
    gsub('^C_xy.','Z_GIS.',name)
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

