gis_vs_linkage_hist <- function(gi_data,
                                col_same_same = rgb(91 / 255, 179 / 255, 228 / 255, 0.7),
                                col_linked = rgb(230 / 255, 154 / 255, 36 / 255, 0.7),
                                col_unlinked = rgb(0, 153 / 255, 128 / 255, 0.7)){
  gi_columns <- grep('^GIS',colnames(gi_data),val=T)
  
  same_same_data <- unlist(dplyr::filter(gi_data,Chromosomal_distance_bp == 0)[,gi_columns])
  linked_data <- unlist(dplyr::filter(gi_data,Chromosomal_distance_bp > 0, Chromosomal_distance_bp < 75000)[,gi_columns])
  unlinked_data <- unlist(dplyr::filter(gi_data,is.na(Chromosomal_distance_bp))[,gi_columns])
  
  #Set up histograms ahead of time to calculate proper breaks
  h1 <- hist(c(same_same_data,linked_data,unlinked_data),breaks = 40, plot = F)
  h2 <- hist(same_same_data,breaks=h1$breaks)
  h3 <- hist(linked_data,breaks=h1$breaks)
  h4 <- hist(unlinked_data,breaks=h1$breaks)
  
  
  
  par(las = 1)
  par(mar=c(4,4.5,1,1))
  plot(
    h2,
    col = col_same_same,
    freq = F,
    ylim = c(0, max(c(
      h2$density, h3$density, h4$density
    ))),
    border = NA,
    xlab = expression(bold('GIS')),
    main = '',
    cex.lab = 1.5,
    cex.axis = 1.5
  )
  box(bty = 'l')
  plot(
    h3,
    col = col_linked,
    freq = F,
    add = T,
    border = NA
  )
  plot(
    h4,
    col = col_unlinked,
    freq = F,
    add = T,
    border = NA
  )
  
  
  legend(
    min(h1$breaks),
    max(c(h2$density, h3$density, h4$density)),
    legend = c('Same Gene Pairs',
               'Linked Gene Pairs',
               'Unlinked Gene Pairs'),
    cex = 1.5,
    col = c(col_same_same, col_linked, col_unlinked),
    pch = 15
  )
}







