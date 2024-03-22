plot_gene_density = function(d, genes){
  
  n = length(genes)
  genes_presence = matrix(0, nrow = 1, ncol = dim(d)[2])
  genes_pixels = list()
  df_tot = NULL
  
  for (i in 1:n){
    
    genes_presence = d[["Spatial"]]$counts[genes[i], ]
    genes_pixels[[i]] = colnames(d)[genes_presence>0]  # barcodes
    x = d@images[[1]]@coordinates[genes_pixels[[i]], 2]
    y = d@images[[1]]@coordinates[genes_pixels[[i]], 3]
    c = d[["Spatial"]]$counts[genes[i], genes_pixels[[i]]]
    df = data.frame(x = x, y = y, c = c, gene = rep(genes[i], length(x)))
    df$barcode = rownames(df)
    rownames(df) = NULL
    
    df_tot = rbind(df_tot, df)
    
  }
  
  xmin = min(d@images[[1]]@coordinates[, 2])
  xmax = max(d@images[[1]]@coordinates[, 2])
  ymin = min(d@images[[1]]@coordinates[, 3])
  ymax = max(d@images[[1]]@coordinates[, 3])
  
  
  g = ggplot(df_tot, aes(x=x, y=y) ) +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, bins=100) + 
    geom_point(size=.2) +
    facet_wrap(vars(gene), nrow = 1) + 
    scale_fill_viridis_c()
  
  return (g)
  
}


# example:
# plot_gene_density(data, c("rpoB", "MAB-0545"))
# data: seurat_object
# second element: vector of genes for which to plot the spatial density