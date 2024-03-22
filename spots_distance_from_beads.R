spots_distance_from_beads = function(beads, points, seurat_data){
  points = points[! points %in% beads] # tolgo le beads stesse
  
  distance_df = NULL
  names_distance = NULL
  
  coord = seurat_data@images[[1]]@coordinates[, c(2,3)]
  
  beads_coord = coord[beads, ]
  rmin = min(beads_coord[,1])-2
  rmax = max(beads_coord[,1])+2
  cmin = min(beads_coord[,2])-2
  cmax = max(beads_coord[,2])+2
  
  for (i in 1:length(points)){
    r = coord[points[i], 1]
    c = coord[points[i], 2]
    
    conditions = (r<rmin | r>rmax)+(c< cmin | c>cmax)  # 1.a neighborhood delle beads, punti da togliere, per cui non calcolo distanze
    if (conditions>0  ){
      list_dist = list()
      for (j in 1:length(beads)){
        r1 = beads_coord[j,1]
        c1 = beads_coord[j,2]
        
        dist_1 = sqrt((r-r1)^2+(c-c1)^2)
        list_dist[[j]] = dist_1 
      }
      
      d_min = min(unlist(list_dist))
      
      distance_df = rbind(distance_df, d_min)
      names_distance = rbind(names_distance, points[i])
    }
  }
  d_df = data.frame(distance = distance_df, row.names = names_distance)
  
  seurat_data$distance = rep(0, length(seurat_data$orig.ident))
  seurat_data$distance[names_distance] = d_df$distance
  
  p = SpatialFeaturePlot(seurat_data, features = "distance", alpha = c(0.6, 0.9))
  
  return (list(distance_matrix = d_df, plot_distance = p))
  
}


# example use:
# a = spots_distance_from_beads(beads, rpoB_pixels, data)                
# beads = vector of barcodes 
# points = vector of barcodes
# data = seurat object
# a is a list of:
# - $distance_matrix
# - $p: spatial dim plot with spots colour proportional to distance from closest beads