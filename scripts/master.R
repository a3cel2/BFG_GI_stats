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



require(qvalue)
require(dplyr)

#Cairo must be installed for most accurate rendering of plots
#pdf base function will work for 
if('Cairo' %in% rownames(installed.packages())){
  pdf_function <- Cairo::CairoPDF
}else{
  pdf_function <- pdf
}


#Needs to be called from /scripts directory for this to run
this.dir <- getwd()
setwd(this.dir)

devtools::load_all('../packages/BfgGiStats')
devtools::document('../packages/BfgGiStats')

setwd('../data')

#Read in genetic interaction data from two tables, using sum column from second table
gi_data_original <-
  read.table('table_s2.txt',
             head = T,
             stringsAsFactors = F)



#Can change this to look at only one of the two technical replicates, analyzes both by default
#Change  to 'R1' or 'R2' to look at one of the two only
to_analyze <- c('R1','R2')



#Update linkage removal criteria to 75kb
gi_data_original[gi_data_original$Chromosomal_distance_bp <= 75000 &
          !is.na(gi_data_original$Chromosomal_distance_bp), 'Remove_by_Chromosomal_distance_or_SameGene'] <-
  'YES'

#Lists the doubling times for each pool and technical replicate
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
  grep_pattern2 <- sprintf('^C_xy', analyze_now)
  
  
  gi_data <-
    gi_data_original[, c(
      grep(
        grep_pattern2,
        colnames(gi_data_original),
        val = T,
        invert = T
      ),
      grep(grep_pattern1, colnames(gi_data_original), val = T)
    )]
  
  
  #Remove non well-measured strains
  #Comment out the next two lines to allow comparison of same-same pairs in well-measured vs non
  well_measured <- gi_data[, grep('HetDipl', colnames(gi_data))] >= 30
  gi_data <- gi_data[well_measured, ]
  
  new_colnames <- colnames(gi_data)[grep(grep_pattern1,colnames(gi_data))]
  new_colnames <- sapply(new_colnames,function(name){gsub(paste(c('_',analyze_now),collapse=''),'',name)})
  colnames(gi_data)[grep(grep_pattern1,colnames(gi_data))] <- new_colnames
  
  
  new_gis <- update_gis(gi_data,
              #These values are empirically determined and have to be
              #in the same order as in gi_data count columns
              g_wt_vec = doubling_vec[[to_analyze[i]]])
  
  
  
  new_gis <- cbind(new_gis,rep(to_analyze[i],nrow(new_gis)))
  colnames(new_gis)[ncol(new_gis)] <- 'Sample'
  
  new_gis <- new_gis[,c(1:2,ncol(new_gis),3:(ncol(new_gis) - 1))]
  
  gi_data_combined <- rbind(gi_data_combined,new_gis)
}


#Estimate double mutant error globally
gi_data <- update_error_model_by_replicates(gi_data_combined)


#Analyze linkage patterns to illustrate new removal criteria
setwd(this.dir)
neutral_gi_data <- filter(gi_data,Type_of_gene_x == 'Neutral' | Type_of_gene_y == 'Neutral')
dir.create('../results', showWarnings = FALSE)
setwd('../results')
pdf_function(file = 'GIS_NoDrug_vs_distance.pdf',
                width = 4.5,
                height = 3.5)
par(mar=c(4.5,5,1,2))
gi_cols <- grep('^GIS',colnames(neutral_gi_data))
plot(log10(rep(neutral_gi_data$Chromosomal_distance_bp,length(gi_cols)) + 1),
     as.matrix(neutral_gi_data[,gi_cols]),
     xlab = expression(Log[10](Chromosomal~dist.~+1)),
     ylab = expression(GIS[xy]~(neutral~pairs)),
     pch = 16,
     col=rgb(0,0,0,0.3),
     las = 1,
     bty = 'l',
     cex.lab = 1.7)
abline(v= log10(75001),col='red',lty=3,lwd=2)
dev.off()

#Make Fig 2 panel
setwd(this.dir)
setwd('../results')
pdf_function(file = 'GIS_by_linkage.pdf',
                width = 4.5,
                height = 4)
gis_vs_linkage_hist(gi_data)
dev.off()



#Plot z distribution
setwd(this.dir)
setwd('../results')
pdf_function(file = 'Z_distribution.pdf',
                width = 10,
                height = 4)
par(mfrow = c(2, 5))
par(mar = c(5, 4.1, 2, 0))
gi_data <- add_p_values(
  gi_data,
  z_grep_pattern = "^Z_GIS",
  nn_pair_type = 'broad',
  make_plots = T
)
dev.off()

#Filter out camptothecin
gi_data <- gi_data[, -grep('CMPT', colnames(gi_data))]

#Preserve at this step for differential gi_calls
gi_data_old <- gi_data

###########################
## Have to temporarily allow non well-measured above for this plot to work
###########################
# setwd(this.dir)
# setwd('../results')
# cut <- 30
# pdf_function(file = 'gis_well_measured_vs_non.pdf',
#                 width = 4.5,
#                 height = 3.5)
# par(mar=c(4,4.5,3,1))
# 
# same_same_gis_well_measured <- as.matrix(
#   gi_data_old %>% dplyr::filter(C_xy.HetDipl >= cut, Chromosomal_distance_bp == 0) %>% dplyr::select(grep('^GIS_xy.', colnames(gi_data_old)))
# )
# 
# same_same_gis_non_well_measured <- as.matrix(
#   gi_data_old %>% dplyr::filter(C_xy.HetDipl < cut, Chromosomal_distance_bp == 0) %>% dplyr::select(grep('^GIS_xy.', colnames(gi_data_old)))
# )
# 
# h1 <- hist(c(same_same_gis_well_measured,same_same_gis_non_well_measured),breaks=100,plot=F)
# h2 <- hist(c(same_same_gis_well_measured),breaks=100,plot=F)
# h3 <- hist(c(same_same_gis_non_well_measured),breaks=100,plot=F)
# 
# hist(
#   same_same_gis_non_well_measured,
#   breaks = h1$breaks,
#   ylim = c(0, max(c(
#     h2$density, h3$density
#   ))),
#   freq = F,
#   border = NA,
#   col = rgb(1,0,0,0.5),
#   xlab='GIS',
#   main = 'GIS of Same-Gene Pairs',
#   cex.lab = 1.7
# )
# hist(
#   same_same_gis_well_measured,
#   breaks = h1$breaks,
#   ylim = c(0, max(c(
#     h2$density, h3$density
#   ))),
#   freq = F,
#   border = NA,
#   col = rgb(0,0,1,0.5),
#   main = '',
#   add = T
# )
# legend(-0.65,3.5,
#        legend=c(expression(C[xy]>=30),
#                 expression(C[xy]<30)),
#        col=c('blue','red'),
#        pch=15)
# dev.off()



#Score comparison scatterplot - by barcode
pdf_function(file = 'st_onge_scatterplot_barcodewise.pdf',
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
pdf_function(file = 'gene_averaged_p_value_histogram.pdf',
                width = 5.5,
                height = 4.5)
hist(as.matrix(gi_data[, grep('^P.neutral', colnames(gi_data))]),
     main = '',
     xlab = 'p-value',
     ylab = 'Number of Interactions',
     col = 'grey30',
     las = 1,
     cex.lab = 1.5)
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

#Precision/recall at the 0.05,0.07 GIs cutoffs
setwd(this.dir)
setwd('../results')
pdf_function(file = 'prec_rec_vs_st_onge.pdf',
                width = 9.4,
                height = 3.8)
prec_rec_vs_stonge(gi_data, fdr_cutoff = 0.01,cutoffs_drawn = c(-0.075,0.075))
dev.off()

#Make AUC plot
setwd(this.dir)
setwd('../results')

pdf_function(file = 'auc_vs_st_onge_new.pdf',
                width = 12,
                height = 4)
par(mfrow = c(1, 3))
par(mar=c(5,5,5,5))
st_onge_auc_plot(gi_data)
dev.off()


#Update calls based on 1% FDR cutoffs
#Effect size cutoffs are 
gi_data <- update_calls(
  gi_data,
  fdr_cutoff_pos = 0.01,
  gi_cutoff_pos = 0.075,
  fdr_cutoff_neg = 0.01,
  gi_cutoff_neg = -0.075
)


#Make gene-wise scatterplot
#Manuscript uses barcode-wise plot
pdf_function(file = 'st_onge_scatterplot_genewise.pdf',
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
  file = 'table_s3.txt'
)
write.table(
  gi_data,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s4.txt'
)

##Make/write differential calls
##Make differential calls, reporting everything

setwd(this.dir)
setwd('../results')
pdf_function(file = 'delta_z_distribution.pdf',
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

#Make an overall histogram for publication
setwd(this.dir)
setwd('../results')
pdf_function(file = 'delta_z_distribution_overall.pdf',
                width = 4,
                height = 3)
delta_z_histogram(differential_calls)
dev.off()

setwd(this.dir)
setwd('../data/output')
write.table(
  differential_calls,
  row.names = F,
  quote = F,
  sep = '\t',
  file = 'table_s5.tsv'
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
  file = 'table_s6.txt'
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
pdf_function(file = 'differential_inters_by_gene.pdf',
                width = 4,
                height = 4)
hist(
  differential_count,
  breaks = 20,
  col = rgb(0, 0, 0, 0.7),
  xlab = 'Number of Significant Differential Pairs',
  ylab = 'Number of genes',
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

pdf_function(file = 'sign_reversals_by_gene.pdf',
                width = 4,
                height = 4)
hist(
  reversal_count,
  breaks = 20,
  col = rgb(0, 0, 0, 0.7),
  xlab = 'Number of Significant Sign-Reversed Pairs',
  ylab = 'Number of genes',
  main = ''
)
dev.off()

#Reset directory
setwd(this.dir)