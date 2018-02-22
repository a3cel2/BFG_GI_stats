library(reshape2)



gi_matr_by_barcode <- function(gi_data,
                               gis_grep_pattern = '^GIS',
                               barcode_name_1 = 'Barcode_x',
                               barcode_name_2 = 'Barcode_y'){
  gi_cols <- grep(gis_grep_pattern, colnames(gi_data), val = T)
  
  genes_x <-
    unique(sapply(gi_data[, 1], function(x) {
      strsplit(x, split = '_')[[1]][1]
    }))
  
  genes_y <-
    unique(sapply(gi_data[, 2], function(y) {
      strsplit(y, split = '_')[[1]][1]
    }))
  
  genes <- unique(c(genes_x, genes_y))
  
  all_data <- c()
  all_donor_data <- c()
  all_recipient_data <- c()
  for (gi_col in gi_cols) {
    test_data <- gi_data[, c(barcode_name_1, barcode_name_2, gi_col)]
    wide_data <-
      dcast(test_data, Barcode_x ~ Barcode_y, fun.aggregate = mean)
    
    
    donor_data <- t(sapply(1:nrow(wide_data), function(i) {
      return(sapply(genes, function(gene) {
        positions <- grep(gene, colnames(wide_data))
        mean_gi <- mean(as.numeric(wide_data[i, positions]))
        return(mean_gi)
      }))
    }))
    
    rownames(donor_data) <- wide_data[, 1]
    
    recipient_data <- t(sapply(2:ncol(wide_data), function(i) {
      return(sapply(genes, function(gene) {
        positions <- grep(gene, wide_data[, 1])
        #retval <-
        #name(retval <- paste(c))
        return(mean(wide_data[positions, i]))
      }))
    }))
    
    #stop()
    rownames(recipient_data) <- colnames(wide_data)[2:ncol(wide_data)]
    comb_data <- rbind(donor_data, recipient_data)
    
    
    colnames(comb_data) <-
      sapply(colnames(comb_data), function(column_name) {
        return(paste(c(column_name, gi_col), collapse = '.'))
      })
    
    
    
    
    #rownames(wide_data) <-
    #  sapply(wide_data[,1], function(row_name) {
    #    return(paste(c(row_name, gi_col), collapse = '.'))
    #  })
    
    
    donor_data <- wide_data[,2:ncol(wide_data)]
    rownames(donor_data) <- wide_data[,1]
    recipient_data <- t(donor_data)
    
    
    colnames(donor_data) <-
      sapply(colnames(donor_data), function(column_name) {
        return(paste(c(column_name, gi_col), collapse = '.'))
      })
    
    colnames(recipient_data) <-
      sapply(colnames(recipient_data), function(column_name) {
        return(paste(c(column_name, gi_col), collapse = '.'))
      })
  
    #rownames(wide_data_new) <- rownames(wide_data)
    
    all_donor_data <- cbind(all_donor_data,as.matrix(donor_data))
    all_recipient_data <- cbind(all_recipient_data,as.matrix(recipient_data))
    
    #paste(c(colnames(comb_data),gi_col),collapse = '.')
    
    all_data <- cbind(all_data, comb_data)
    
    
    #stop()
    #
    
  }
  return(all_data)
}
#stop()

# cor_structure <- sapply(genes,function(gene){
#   to_search <- paste(c(gene,'_'),collapse='')
#   gene_data <- all_data[grep(to_search,rownames(all_data)),]
#   
#   if(length(gene_data) > 1){
#     cors <- cor(t(gene_data),use='pair')
#     return(min(cors[upper.tri(cors)],na.rm=T))
#   }
#   return(NA)
#   
# })
# 
# 




barcode_deviance_mean_rmsd <- function(all_data) {
  rmsds <- sapply(1:nrow(all_data), function(i) {
    print(i)
    barcode_gis <- all_data[i, ]
    barcode_name <- rownames(all_data)[i]
    
    gene <- strsplit(barcode_name, split = '_')[[1]][1]
    gene <- sprintf('%s_',gene)
    
    barcode_with_gene <- grep(gene, rownames(all_data), val = T)
    gis_with_gene <- all_data[barcode_with_gene, ]
    
    #What GIs this barcode will be compared to
    comparison_indeces <-
      grep(barcode_name, rownames(gis_with_gene), invert = T)
    
    if (length(comparison_indeces) > 0) {
      comparison_gis <<- gis_with_gene[comparison_indeces, , drop = F]
      
      #What GIs this barcode has
      barcode_gi <<- gis_with_gene[barcode_name,]
      
      #print(gene)
      if(gene == 'RAD5'){
        print(gene)
        stop()
      }
      
      
      comparison_rmsds <-
        apply(comparison_gis, 1, function(comparison_gi) {
          
          #if_comparable_gi <-
          #  !(is.na(comparison_gi) | is.na(barcode_gi))
          
          #comparison_gi <- comparison_gi[if_comparable_gi]
          #barcode_gi <- barcode_gi[if_comparable_gi]
          
          #rmsd <-
          #  sqrt(sum((comparison_gi - barcode_gi) ^ 2) / sum(if_comparable_gi))
          #return(cor(comparison_gi,barcode_gi))
          
          #return(rmsd)
          
          diff <- abs(comparison_gi - barcode_gi)
          return(diff)
        })
      
      return(apply(comparison_rmsds,1,function(comparison_rmsd){
        return(min(comparison_rmsd,na.rm=T))
      }))
      
      #return(sapply(1:nrow(comparison_rmsds),function(i){
      #  min(comparison_rmsds[i,],na.rm=T)#/abs(barcode_gi[i] + 0.1)
      #}))
      #return(apply(comparison_rmsds,1,function(x){x/barcode_gi}))
      #return(min(comparison_rmsds,na.rm=T))
    } else {
      return(NA)
    }
  })
  
  names(rmsds) <- rownames(all_data)
  return(rmsds)
  
}

gi_data_min_diffs <- function(gi_data,
                           gis_grep_pattern = '^GIS',
                           barcode_1_name = 'Barcode_x',
                           barcode_2_name = 'Barcode_y') {
  genes_x <-
    sapply(gi_data[, barcode_1_name], function(x) {
      strsplit(x, split = '_')[[1]][1]
    })
  
  genes_y <-
    sapply(gi_data[, barcode_2_name], function(y) {
      strsplit(y, split = '_')[[1]][1]
    })
  
  
  gi_cols <- grep(gis_grep_pattern, colnames(gi_data), val = T)
  
  
  
  min_diffs_by_gi <- sapply(1:nrow(gi_data), function(i) {
    #print(i)
    gene_x <- genes_x[i]
    gene_y <- genes_y[i]
    
    grep_gene_x <- paste(c(gene_x, '_'), collapse = '')
    grep_gene_y <- paste(c(gene_y, '_'), collapse = '')
    
    crit <-
      (grepl(grep_gene_x, gi_data[ , barcode_1_name]) &
         grepl(grep_gene_y, gi_data[ , barcode_2_name]))
    
    #stop()
    crit <-
      crit |
      (grepl(grep_gene_x, gi_data[ , barcode_2_name]) &
         grepl(grep_gene_y, gi_data[ , barcode_1_name]))
    
    
    
    
    if (sum(crit) < 2) {
      crit_data <- gi_data[which(crit), gi_cols]
      crit_data[!is.na(crit_data)] <- NA
      return(crit_data)
    }
    
    else{
      relevant_data <- gi_data[which(crit), ]
      this_crit <-
        grepl(gi_data[i, 1], relevant_data[,barcode_1_name]) &
        grepl(gi_data[i, 2], relevant_data[,barcode_2_name])
      
      this_data <- relevant_data[which(this_crit), gi_cols]
      other_data <- relevant_data[which(!this_crit), gi_cols, drop = F]
      
      diff_vec <- t(apply(other_data, 1, function(x) {
        abs(x - as.numeric(this_data))
      }))
      
      return(apply(diff_vec, 2, min))
    }
    
  })
  return(t(min_diffs_by_gi))
}

stop()

##MAIN##
barcode_gi_matr <- gi_matr_by_barcode(gi_data)
barcode_gi_matr_bkup <- barcode_gi_matr
#outlier_names <- 
for(i in 1:10){
  
  barcode_deviance_rmsd_vec <- barcode_deviance_mean_rmsd(barcode_gi_matr)
  
  
  plot(barcode_deviance_rmsd_vec)
  outlier_name <- which.max(barcode_deviance_rmsd_vec)
  print(names(outlier_name))
  print(max(barcode_deviance_rmsd_vec,na.rm=T))
  barcode_gi_matr <- barcode_gi_matr[-outlier_name,]
  
  
}
plot(barcode_deviance_rmsd_vec)


  #apply(all_data,1,function(barcode_gis){
  #barcode_gene <- names(barcode_gis)
  #return(barcode_gene)
#})
