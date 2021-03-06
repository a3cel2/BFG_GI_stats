devtools::use_package('metap')


#' Combine data from multiple barcodes for each gene
#'
#' @param gi_data barcode-wise genetic interaction data
#' @param type_column_grep regex to find 'Type of gene' columns
#' @param gi_column_grep regex to fid GIS columns
#' @param gi_err_column_grep regex to find GIS error columns
#' @param p_val_column_grep regex to find p-value columns
#' @param z_column_grep regex to find Z_GIS and Z_GIS class columns
#' @param w_x_column_grep regex to find gene x fitness columns
#' @param w_x_err_column_grep regex to find gene x fitness error columns
#' @param w_y_column_grep regex to find gene y fitness columns
#' @param w_y_err_column_grep regex to find gene y fitness error columns
#' @param w_xy_column_grep regex to find double mutant fitness columns
#' @param w_xy_err_column_grep regex to find double mutant fitness error columns
#' @param count_column_grep regex to find count columns
#'
#' @return a version of gi_data with data from multiple barcodes combined for each gene
average_gi_data_by_gene <-
  function(gi_data,
           type_column_grep = 'Type_of_gene',
           gi_column_grep = '^GIS',
           gi_err_column_grep = '^SE_GIS',
           fdr_column_grep = '^P.neutral',
           z_column_grep = '^Z_GIS',
           w_x_column_grep = '^W_x\\.',
           w_x_err_column_grep = '^W_x_SE\\.',
           w_y_column_grep = '^W_y\\.',
           w_y_err_column_grep = '^W_y_SE\\.',
           w_xy_column_grep = '^W_xy\\.',
           w_xy_err_column_grep = '^W_xy_SE\\.',
           count_column_grep = '^C_xy') {
    Gene1 <-
      sapply(gi_data$Barcode_x, function(name) {
        paste(strsplit(name, split = '_')[[1]][1],collapse='_')
      })
    
    Gene2 <-
      sapply(gi_data$Barcode_y, function(name) {
        paste(strsplit(name, split = '_')[[1]][1],collapse='_')
      })
    
    gene_columns <- cbind(Gene1, Gene2)
    
    #Alphabetically re-order genes for easier queries
    alpha_order <-
      t(apply(gene_columns, 1, function(x) {
        sort(x, index.return = T)$ix
      }))
    
    gene_columns <- t(sapply(1:nrow(gi_data), function(i) {
      gene_columns[i, alpha_order[i,]]
    }))
    
    type_columns <- grep(type_column_grep, colnames(gi_data))
    
    gi_data[, type_columns] <-
      as.data.frame(t(sapply(1:nrow(gi_data), function(i) {
        dat <- gi_data[i, type_columns]
        return(dat[alpha_order[i,]])
      })))
    
    gi_data[, type_columns[1]] <- unlist(gi_data[, type_columns[1]])
    gi_data[, type_columns[2]] <- unlist(gi_data[, type_columns[2]])
    
    Gene1 <- gene_columns[, 1]
    Gene2 <- gene_columns[, 2]
    
    preserved_colnames <- colnames(gi_data)
    preserved_colnames <-
      preserved_colnames[grep('Barcode', preserved_colnames, invert = T)]
    preserved_colnames <- c('Gene1', 'Gene2', preserved_colnames)
    
    gi_data_grp <- cbind(gi_data, cbind(Gene1, Gene2))
    
    split_gi_data_grp <- split(gi_data_grp,
                               list(Gene1,
                                    Gene2))
    
    
    ret_df <- c()
    for (i in 1:length(split_gi_data_grp)) {
      testset <- split_gi_data_grp[[i]]
      
      #print(testset)
      
      if (nrow(testset) > 0) {
        ret_vec <- testset[1,]
        
        gis_columns <- grep(gi_column_grep, colnames(gi_data))
        err_columns <- grep(gi_err_column_grep, colnames(testset))
        z_columns <-
          grep(
            'Class',
            grep(z_column_grep, colnames(gi_data), val = T),
            val = T,
            invert = T
          )
        z_class_columns <-
          grep(
            'Class',
            grep(z_column_grep, colnames(gi_data), val = T),
            val = T,
            invert = F
          )
        fdr_columns <- grep(fdr_column_grep, colnames(gi_data))
        
        
        
        new_gis <- sapply(1:length(gis_columns), function(i) {
          var_weights <- 1 / (testset[, err_columns[i]] ^ 2)
          return(weighted.mean(x = testset[, gis_columns[i]], w = var_weights))
        })
        
        new_errs <- sapply(1:length(gis_columns), function(i) {
          var_weights <- 1 / (testset[, err_columns[i]] ^ 2)
          sum_var <- sum(var_weights)
          
          return(sqrt(sum((
            testset[, err_columns[i]] * (var_weights / sum_var)
          ) ^ 2)))
          
        })
        
        new_zs <- new_gis / new_errs
        
        new_fdrs <- sapply(1:length(fdr_columns), function(i) {
          var_weights <- 1 / (testset[, err_columns[i]] ^ 2)
          
          p_vals <- testset[, fdr_columns[i]]
          
          
          #metap::sumz exits with one p value instead of just returning the value...
          #so just return it directly
          if (length(p_vals) == 1) {
            return(p_vals)
          }
          
          p_positive <- p_vals
          p_negative <- p_vals
          
          gi_scores <- testset[, gis_columns[i]]
          nominal_gi <- mean(testset[, gis_columns[i]])
          
          if (nominal_gi < 0) {
            p_vals[gi_scores > 0] <- 1 - p_vals[gi_scores > 0]
          } else{
            p_vals[gi_scores < 0] <- 1 - p_vals[gi_scores < 0]
          }
          
          #Next two lines stops metap::sumz from giving error
          #p-values should never be 0, but can be because of float precision
          p_vals[p_vals < .Machine$double.xmin] <- .Machine$double.xmin
          
          #Affects no cases, but added for good measure
          if (max(p_vals) == 1) {
            return(1)
          }
          
          meta_p_val <-
            metap::sumz(p_vals, weights = var_weights)$p[1]
          
          #Two-tailed meta-analysis
          return(min(meta_p_val, 1 - meta_p_val) * 2)
          
        })
        
        new_zs <- new_gis / new_errs
        
        ret_vec[gis_columns] <- new_gis
        ret_vec[z_columns] <- new_zs
        ret_vec[fdr_columns] <- new_fdrs
        ret_vec[err_columns] <- new_errs
        
        #Update w metrics
        w_x_columns <- grep(w_x_column_grep, colnames(gi_data))
        w_x_err_columns <-
          grep(w_x_err_column_grep, colnames(gi_data))
        
        w_y_columns <- grep(w_y_column_grep, colnames(gi_data))
        w_y_err_columns <-
          grep(w_y_err_column_grep, colnames(gi_data))
        
        w_xy_columns <- grep(w_xy_column_grep, colnames(gi_data))
        w_xy_err_columns <-
          grep(w_xy_err_column_grep, colnames(gi_data))
        
        
        
        #Can probably loop this somehow
        new_w_x <- sapply(1:length(w_x_columns), function(i) {
          var_weights <- 1 / (testset[, err_columns[i]] ^ 2)
          return(weighted.mean(x = testset[, w_x_columns[i]], w = var_weights))
        })
        new_w_x_errs <-
          sapply(1:length(w_x_err_columns), function(i) {
            var_weights <- 1 / (testset[, err_columns[i]] ^ 2)
            sum_var <- sum(var_weights)
            return(sqrt(sum((
              testset[, w_x_err_columns[i]] * (var_weights / sum_var)
            ) ^ 2)))
          })
        
        new_w_y <- sapply(1:length(w_y_columns), function(i) {
          var_weights <- 1 / (testset[, err_columns[i]] ^ 2)
          return(weighted.mean(x = testset[, w_y_columns[i]], w = var_weights))
        })
        new_w_y_errs <-
          sapply(1:length(w_y_err_columns), function(i) {
            var_weights <- 1 / (testset[, err_columns[i]] ^ 2)
            sum_var <- sum(var_weights)
            return(sqrt(sum((
              testset[, w_y_err_columns[i]] * (var_weights / sum_var)
            ) ^ 2)))
          })
        
        new_w_xy <- sapply(1:length(w_xy_columns), function(i) {
          var_weights <- 1 / (testset[, err_columns[i]] ^ 2)
          return(weighted.mean(x = testset[, w_xy_columns[i]], w = var_weights))
        })
        new_w_xy_errs <-
          sapply(1:length(w_xy_err_columns), function(i) {
            var_weights <- 1 / (testset[, err_columns[i]] ^ 2)
            sum_var <- sum(var_weights)
            return(sqrt(sum((
              testset[, w_xy_err_columns[i]] * (var_weights / sum_var)
            ) ^ 2)))
          })
        
        ret_vec[w_x_columns] <- new_w_x
        ret_vec[w_x_err_columns] <- new_w_x_errs
        ret_vec[w_y_columns] <- new_w_y
        ret_vec[w_y_err_columns] <- new_w_y_errs
        ret_vec[w_xy_columns] <- new_w_xy
        ret_vec[w_xy_err_columns] <- new_w_xy_errs
        
        
        ret_df <- rbind(ret_df, ret_vec)
      }
    }
    
    gi_data_grp <- ret_df
    
    gi_data_grp <- gi_data_grp[, preserved_colnames]
    
    #Averaging count data doesn't really make sense
    #So best to remove
    gi_data_grp <-
      gi_data_grp[, grep(count_column_grep, colnames(gi_data_grp), invert = T)]
    #Samples are averaged together in this step
    gi_data_grp <- gi_data_grp[, grep("Sample", colnames(gi_data_grp), invert =T)]
    
    
    colnames(gi_data_grp)[1:2] <- c('Gene_x', 'Gene_y')
    
    
    gi_data_grp[, 1] <- as.character(gi_data_grp[, 1])
    gi_data_grp[, 2] <- as.character(gi_data_grp[, 2])
    
    return(gi_data_grp)
  }

# Not used, was experimenting with alternate filtering method
# 
# find_discordant_pairs <-
#   function(gi_data,
#            type_column_grep = 'Type_of_gene',
#            gi_column_grep = '^GIS',
#            gi_err_column_grep = '^SE_GIS',
#            fdr_column_grep = '^P.neutral',
#            z_column_grep = '^Z_GIS',
#            w_x_column_grep = '^W_x\\.',
#            w_x_err_column_grep = '^W_x_SE\\.',
#            w_y_column_grep = '^W_y\\.',
#            w_y_err_column_grep = '^W_y_SE\\.',
#            w_xy_column_grep = '^W_xy\\.',
#            w_xy_err_column_grep = '^W_xy_SE\\.',
#            count_column_grep = '^C_xy') {
#     Gene1 <-
#       sapply(gi_data$Gene_x, function(name) {
#         paste(strsplit(name, split = '_')[[1]][1],collapse='_')
#       })
#     
#     Gene2 <-
#       sapply(gi_data$Gene_y, function(name) {
#         paste(strsplit(name, split = '_')[[1]][1],collapse='_')
#       })
#     
#     gene_columns <- cbind(Gene1, Gene2)
#     
#     #Alphabetically re-order genes for easier queries
#     alpha_order <-
#       t(apply(gene_columns, 1, function(x) {
#         sort(x, index.return = T)$ix
#       }))
#     
#     gene_columns <- t(sapply(1:nrow(gi_data), function(i) {
#       gene_columns[i, alpha_order[i,]]
#     }))
#     
#     type_columns <- grep(type_column_grep, colnames(gi_data))
#     
#     gi_data[, type_columns] <-
#       as.data.frame(t(sapply(1:nrow(gi_data), function(i) {
#         dat <- gi_data[i, type_columns]
#         return(dat[alpha_order[i,]])
#       })))
#     
#     gi_data[, type_columns[1]] <- unlist(gi_data[, type_columns[1]])
#     gi_data[, type_columns[2]] <- unlist(gi_data[, type_columns[2]])
#     
#     Gene1 <- gene_columns[, 1]
#     Gene2 <- gene_columns[, 2]
#     
#     preserved_colnames <- colnames(gi_data)
#     preserved_colnames <-
#       preserved_colnames[grep('Barcode', preserved_colnames, invert = T)]
#     preserved_colnames <- c('Gene1', 'Gene2', preserved_colnames)
#     
#     gi_data_grp <- cbind(gi_data, cbind(Gene1, Gene2))
#     
#     split_gi_data_grp <- split(gi_data_grp,
#                                list(Gene1,
#                                     Gene2))
#     
#     pee_vals <- sapply(split_gi_data_grp,function(grp){
#       if(nrow(grp) == 2){
#         print(grp)
#         g1 <- as.numeric(grp[1,grep(gi_column_grep,colnames(grp))])
#         g2 <- as.numeric(grp[2,grep(gi_column_grep,colnames(grp))])
#         
#         return(abs(g1 - g2))
#         #return(t.test(g1,g2,paired=T)$p.val)
#         #t.test(grp[1,grep(gi_column_grep,colnames(grp))],grp[2,grep(gi_column_grep,colnames(grp))],paired=T)$p.val
#       }
#     })
#     
#     return(pee_vals)
#   }