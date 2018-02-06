##################################################################
##################################################################
##################################################################
########                                                  ########
########                                                  ########
######## Pipeline to calculate Genetic Interaction Scores ########
######## from Barcode Fusion Genetics experiments         ########
######## and stastical evaluation of condition-dependent  ########
######## genetic interactions                             ########         
########                                                  ########
########                                                  ########
########                                                  ########
##################################################################
##################################################################


require(Cairo)
require(qvalue)
require(dplyr)


#Needs to be called from /scripts directory for this to run
this.dir <- getwd()
setwd(this.dir)

devtools::load_all('../packages/BfgGiStats')
devtools::document('../packages/BfgGiStats')

setwd('../data')

#Read in genetic interaction data from two tables, using sum column from second table
gi_data_original <-
  read.table('table_s1_input.tsv',
             head = T,
             stringsAsFactors = F)
#Remove count data from first table
gi_data_original <- gi_data_original[, -grep('^C_', colnames(gi_data_original))]



#Can change this to look at only one of the two technical replicates, analyzes the sum by default
#Change  to 'R1' or 'R2' to look at one of the two only
to_analyze <- c('R1','R2')



#Update linkage removal criteria to 75kb
#gi_data_original[gi_data_original$Chromosomal_distance_bp <= 75000 &
#          !is.na(gi_data_original$Chromosomal_distance_bp), 'Remove_by_Chromosomal_distance_or_SameGene'] <-
#  'YES'
#Remove non well-measured strains




#For troubleshooting/checking
#gi_data_original <- gi_data

#Change genetic interactions and significance calls
doubling_vec <- list(R1 =  c('NoDrug' = 12.24,
                             'DMSO' = 8.2,
                             'MMS' = 7.8,
                             '4NQO' = 7.6,
                             'BLMC' = 7.12,
                             'ZEOC' = 6.24,
                             'HYDX' = 7.92,
                             'DXRB' = 7.36,
                             'CMPT' = 7.68,
                             'CSPL' = 8.48),
                      R2 = c('NoDrug' = 13,
                             'DMSO' = 8.47,
                             'MMS' = 7.88,
                             '4NQO' = 7.4,
                             'BLMC' = 6.76,
                             'ZEOC' = 6.32,
                             'HYDX' = 7.6,
                             'DXRB' = 6.72,
                             'CMPT' = 7.72,
                             'CSPL' = 8.4))



gi_data_combined<- c()
for(i in 1:length(to_analyze)){
  analyze_now <- to_analyze[i]
  
  grep_pattern1 <- sprintf('^C_xy.*%s', analyze_now)
  grep_pattern2 <- sprintf('.%s', analyze_now)
  
  gi_data_2 <-
    read.table('Cxy_Table_WithReplicates.txt',
               head = T,
               stringsAsFactors = F)
  gi_data <-
    cbind(gi_data_original, gi_data_2[, grep(grep_pattern1, colnames(gi_data_2))])
  colnames(gi_data) <- gsub(grep_pattern2, '', colnames(gi_data))
  
  
  well_measured <- gi_data[, grep('HetDipl', colnames(gi_data))] >= 30
  gi_data <- gi_data[well_measured,]
  
  gi_data <- gi_data[grep('MMS4_donor',gi_data$Barcode_x,invert=T),]
  gi_data <- gi_data[grep('MUS81_donor',gi_data$Barcode_x,invert=T),]
  #gi_data <- gi_data[grep('MMS4_recipient',gi_data$Barcode_y,invert=T),]
  #gi_data <- gi_data[grep('MUS81_recipient',gi_data$Barcode_y,invert=T),]
  
  
  new_gis <- update_gis(gi_data,
              #These values are empirically determined and have to be
              #in the same order as in gi_data count columns
              g_wt_vec = doubling_vec[[to_analyze[i]]])
  
  
  
  new_gis <- cbind(new_gis,rep(to_analyze[i],nrow(new_gis)))
  colnames(new_gis)[ncol(new_gis)] <- 'Sample'
  
  new_gis <- new_gis[,c(1:2,ncol(new_gis),3:(ncol(new_gis) - 1))]
  
  gi_data_combined <- rbind(gi_data_combined,new_gis)
}



#Instead of a poisson error model for double mutant fitness, add a global model for each condition
gi_data <- update_error_model_by_replicates(gi_data_combined)
gi_data <- gi_data[, -grep('CMPT', colnames(gi_data))]


#Analyze linkage patterns to justify above removal criteria
setwd(this.dir)
dir.create('../results', showWarnings = FALSE)
setwd('../results')
Cairo::CairoPDF(file = 'GIS_NoDrug_vs_distance.pdf',
                width = 4,
                height = 4)
par(mar=c(4.5,4.5,1,1))
plot(log10(gi_data$Chromosomal_distance_bp + 1),
     gi_data$GIS_xy.NoDrug,
     xlab = expression(Log[10](Chromosomal~dist.~+1)),
     ylab = expression(GIS[xy]~(no~Drug)),
     pch = 16,
     col=rgb(0,0,0,0.3))
abline(v= log10(75001),col='red',lty=3,lwd=2)
dev.off()

#Plot z distribution
setwd(this.dir)
setwd('../results')
Cairo::CairoPDF(file = 'Z_distribution.pdf',
                width = 10,
                height = 4)
par(mfrow = c(2, 5))
par(mar = c(5, 4, 2, 0))
gi_data <- add_p_values(
  gi_data,
  z_grep_pattern = "^Z_GIS",
  nn_pair_type = 'broad',
  make_plots = T
)
dev.off()



#Preserve at this step for differential gi_calls
gi_data_old <- gi_data

## Reviewer plot, have to temporarily allow non well-measured for this plot to work
# setwd(this.dir)
# setwd('../results')
# cut <- 30
# Cairo::CairoPDF(file = 'gis_well_measured_vs_non.pdf',
#                 width = 4,
#                 height = 3.5)
# par(mar=c(4,4,3,1))
# plot(density(as.matrix(
#   gi_data_old %>% dplyr::filter(C_xy.HetDipl >= cut, Chromosomal_distance_bp == 0) %>% dplyr::select(grep('^GIS_xy.', colnames(gi_data_old)))
# )), xlim = c(-1.1, 0.5),lwd=2,col='red',xlab='GIS',main='GIS of Same-Same Pairs')
# lines(density(as.matrix(
#   gi_data_old %>% dplyr::filter(C_xy.HetDipl < cut, Chromosomal_distance_bp == 0) %>% dplyr::select(grep('^GIS_xy.', colnames(gi_data_old)))
# )),
# col='blue',
# lwd='2'
# )
# legend(-0.1,3,
#        legend=c(expression(C[xy]>=30),
#                 expression(C[xy]<30)),
#        col=c('red','blue'),
#        pch=15)
# dev.off()




#Score comparison scatterplot - by barcode
Cairo::CairoPDF(file = 'st_onge_scatterplot_barcodewise.pdf',
                width = 7.5,
                height = 3.5)
par(mfrow = c(1, 2))
st_onge_scatterplot(gi_data)
dev.off()


########
#Calculate a per-gene GI score, Z score, and p value
########


gi_data <- average_gi_data_by_gene(gi_data)


#Make sure p-value distribution looks sane
Cairo::CairoPDF(file = 'gene_averaged_p_value_histogram.pdf',
                width = 5.5,
                height = 4.5)
hist(as.matrix(gi_data[, grep('^P.neutral', colnames(gi_data))]),
     main = '',
     xlab = 'p-value',
     col = 'grey30')
dev.off()



#Calculate FDR from per-gene p values
gi_data[, grep('^P.neutral', colnames(gi_data))] <-
  apply(gi_data[, grep('^P.neutral', colnames(gi_data))], 2, function(x) {
    qvalue(x)$q
  })

fdr_colnames <-
  sapply(grep('P.neut', colnames(gi_data), val = T), function(x) {
    gsub('^P.', 'FDR.', x)
  })

colnames(gi_data)[grep('^P.neutral', colnames(gi_data))] <- fdr_colnames

#Precision at the 0.05,0.07 GIs cutoffs
setwd(this.dir)
setwd('../results')
Cairo::CairoPDF(file = 'gi_prec_vs_st_onge.pdf',
                width = 4.5,
                height = 4)
precision_vs_stonge(
  gi_data,
  fdr_cutoff = 0.01,
  xlims = c(-0.12, 0.12),
  cutoffs_drawn = c(0.07,0.05)
)
dev.off()


#Performance heatmaps in precision and recall as a function of FDR and GIS cutoff
setwd(this.dir)
setwd('../results')

performance_data <- gi_data
performance_data[, grep('^FDR', colnames(performance_data))] <-
  -log10(performance_data[, grep('^FDR', colnames(performance_data))]) *
  sign(performance_data[, grep('^GI', colnames(performance_data))])


prec_vs_st_onge <- make_performance_matrix(performance_data,metr='prec',y_cutoff_vec = seq(0,0.1,length.out=50),
                                           x_cutoff_vec = seq(0,8,length.out=50))
rec_vs_st_onge <- make_performance_matrix(performance_data,metr='rec',y_cutoff_vec = seq(0,0.1,length.out=50),
                                          x_cutoff_vec = seq(0,8,length.out=50))

perf_data <- list(prec_vs_st_onge,rec_vs_st_onge)
names(perf_data) <- c('prec','rec')
directions <- c('positive_interactions','negative_interactions')

for(i in 1:length(perf_data)){
  for(j in 1:length(directions)){
    
    if(names(perf_data)[i] == 'prec'){
      legend_title <- 'Precision'
      min_val <- 0.3
      max_val <- 0.9
    }else if(names(perf_data)[i] == 'rec'){
      legend_title <- 'Recall'
      min_val <- 0.2
      max_val <- 0.8
    }
    
    if(directions[j] == 'positive_interactions'){
      main_title <- '(positive GIs)'
    }else if(directions[j] == 'negative_interactions'){
      main_title <- '(negative GIs)'
    }
    
    
    filename <- sprintf('%s vs St. Onge %s.pdf',legend_title,main_title)
    Cairo::CairoPDF(file = filename,
                    width = 5,
                    height = 4)
    performance_heatmap(
      perf_data[[i]][[directions[j]]],
      xlab = '-Log10(FDR) Threshold',
      ylab = 'GIS Threshold',
      min_val = min_val,
      max_val = max_val,
      highlight_max_vals = F,
      add_legend = T,
      legend_title = legend_title,
      main = sprintf('%s vs St. Onge %s',legend_title,main_title)
    )
    dev.off()
    
  }
}


#Make AUC plot
setwd(this.dir)
setwd('../results')

Cairo::CairoPDF(file = 'auc_vs_st_onge_new.pdf',
                width = 11,
                height = 4)
par(mfrow = c(1, 3))
st_onge_auc_plot(gi_data)
dev.off()

#stop()

#Update calls based on 1% FDR cutoffs
#Effect size cutoffs are 
gi_data <- update_calls(
  gi_data,
  fdr_cutoff_pos = 0.01,#0.01,
  gi_cutoff_pos = 0.05,#.05,
  fdr_cutoff_neg = 0.01,#0.01,#0.01,
  gi_cutoff_neg = -0.08#.05
)


#Make gene-wise scatterplot
#Manuscript uses barcode-wise plot
Cairo::CairoPDF(file = 'st_onge_scatterplot_genewise.pdf',
                width = 7.5,
                height = 3.5)
par(mfrow = c(1, 2))
st_onge_scatterplot(gi_data)
dev.off()


###Write some data
setwd(this.dir)
setwd('../data')
dir.create('output', showWarnings = FALSE)
setwd('output')

write.table(
  gi_data_old,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s1.tsv'
)
write.table(
  gi_data,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s1_genewise.tsv'
)

##Make/write differential calls
##Make differential calls, reporting everything

setwd(this.dir)
setwd('../results')
Cairo::CairoPDF(file = 'delta_z_distribution.pdf',
                width = 12,
                height = 10)
par(mfrow = c(6, 6))
par(mar = c(4, 4, 1, 1))
differential_calls <- differential_gi_analysis(
  gi_data,
  #Report everything
  fdr_cutoff = 1.1,
  require_sign_change = F,
  make_plots = T
)
dev.off()

setwd(this.dir)
setwd('../data/output')
write.table(
  differential_calls,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s2_all.tsv'
)


#Filter for significance and sign change
differential_calls_sig <- differential_gi_analysis(
  gi_data,
  fdr_cutoff = 0.01,
  delta_gi_cutoff = 0.1,
  require_sign_change = T
)

#Can compare what happens if no sign change enforced
differential_calls_sig_no_sign <- differential_gi_analysis(
  gi_data,
  fdr_cutoff = 0.01,
  delta_gi_cutoff = 0.1,
  require_sign_change = F
)



#Write some data
setwd(this.dir)

setwd('../data')
dir.create('output', showWarnings = FALSE)
setwd('output')

write.table(
  differential_calls_sig,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s2_significant.tsv'
)


#Looks at frequency of differetial calls by gene
setwd(this.dir)
setwd('../results')
ddr_data <- dplyr::filter(gi_data, Type_of_gene_x != "Neutral", Type_of_gene_y != "Neutral")
genes <- unique(c(ddr_data$Gene_x,ddr_data$Gene_y))
differential_count <-
  sapply(genes, function(gene) {
    return(nrow(
      dplyr::filter(differential_calls_sig, Gene_x == gene |
               Gene_y == gene)
    ))
  })
Cairo::CairoPDF(file = 'differential_inters_by_gene.pdf',
                width = 4,
                height = 4)
hist(
  differential_count,
  breaks = 20,
  col = rgb(0, 0, 0, 0.7),
  xlab = 'Number of Differential Pairs',
  ylab = 'Frequency (number of genes)',
  main = ''
)
dev.off() 


#Look at only sign reversals
reversal_count <-
  sapply(genes, function(gene) {
    return(nrow(
      dplyr::filter(differential_calls_sig, Gene_x == gene |
               Gene_y == gene, Class_Condition1 != "Expected", Class_Condition2 != "Expected")
    ))
  })

Cairo::CairoPDF(file = 'sign_reversals_by_gene.pdf',
                width = 4,
                height = 4)
hist(
  reversal_count,
  breaks = 20,
  col = rgb(0, 0, 0, 0.7),
  xlab = 'Number of Sign-Reversed Pairs',
  ylab = 'Frequency (number of genes)',
  main = ''
)
dev.off()

#For re-running
setwd(this.dir)