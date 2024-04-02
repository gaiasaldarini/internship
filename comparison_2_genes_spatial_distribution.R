divide_spatial = function (coord, x_seq, y_seq){
  
  tot_coord = coord
  tot_coord$quadrante = rep(0, dim(tot_coord)[1])
  rmin = min(tot_coord[,1])
  rmax = max(tot_coord[,1])
  cmin = min(tot_coord[,2])
  cmax = max(tot_coord[,2])
  row_index = as.integer(seq(rmin, rmax, length.out = x_seq+1))
  col_index = as.integer(seq(cmin, cmax, length.out = y_seq+1))
  quad = 1
  for (r in 1:x_seq){
    for (c in 1:y_seq){
      tot_coord$try = (tot_coord$row >= row_index[r]) * (tot_coord$row <= row_index[r+1]) * (tot_coord$col >= col_index[c]) * (tot_coord$col <= col_index[c+1])
      
      tot_coord[tot_coord$try == 1, 3] = quad
      quad = quad+1
    }
  }
  
  vect = data.frame(barcodes = rownames(tot_coord), quadrante = tot_coord$quadrante)
  return (vect)
  
}

# use:
# coord = data_tot@images[[1]]@coordinates[, c(2,3)]
# info_quad = divide_spatial(coord, 5, 5)
# numero di divisioni su asse x e y 
# coord is dataframe with row_index, column_index as columns
# a row for each barcode
# returns a vector of length equal to the number of spots: the value corresponds to the box 
# the spot belongs to 


conte_per_quadrante = function(d, num_quad, gene_x){
  
  cc = NULL
  
  for (i in num_quad){
    data_now = d@assays$Spatial$counts[gene_x, d$quadrante == i]
    cc[i] = sum(data_now)
    names(cc)[i] = paste0("quadrante", i)
  }
  
  cc = na.omit(cc)
  return (cc)
}


cfr_2genes_distances_pvalues = function(d, genes){
  
  gene1_quad = conte_per_quadrante(d, sort(unique(d$quadrante)), genes[1])
  gene2_quad = conte_per_quadrante(d, sort(unique(d$quadrante)), genes[2])
  tot_quad = as.numeric(table(d$quadrante))
  confronto_df = data.frame(gene1_quad, gene2_quad, tot = tot_quad)
  colnames(confronto_df) = c(genes[1], genes[2], "tot")
  
  pvalues = NULL
  for (q in 1:dim(confronto_df)[1]){
    t = matrix(0, ncol=2, nrow=2)
    colnames(t) = c(genes[1], genes[2])
    rownames(t) = c("myquad", "otherquads")
    t[1,1] = confronto_df[q, 1]
    t[1,2] = confronto_df[q, 2]
    t[2,1] = sum(confronto_df[,1])
    t[2,2] = sum(confronto_df[,2])
    ftest = fisher.test(t) 
    pvalues[q] = ftest$p.value}
  
  pvalues_adjusted = p.adjust(pvalues, method =c("BH"))
  p_val = data.frame(p_values = pvalues, p_adjusted = pvalues_adjusted, sig_padj = pvalues_adjusted< 0.05, row.names = rownames(confronto_df))
  
  return (p_val)
  
}


# example:
# pvals = cfr_2genes_distances_pvalues(data, c("rpoB", "MAB-0545"))
# data: seurat object
# second element: vector of 2 genes names for which to compare the gene expression distribution 
# across the different boxes
# returns a dataframe with pvalues, pvalues adjusted and significance of pvalues adjusted at 0.05 threshold
# for the fisher test for each box
