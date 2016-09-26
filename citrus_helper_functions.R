extract_citrus_clusters <- function(RDSfile, anno, cluster_id=NULL) {
  
  require(tools)
  short_name <- basename(file_path_sans_ext(RDSfile)) 
  rds <- readRDS(RDSfile)
  if (is.null(cluster_id)) {
  # no cluster_id supplied, get the glmnet cv.min clusters from the RDS  
    if ( class(rds$conditionRegressionResults[[short_name]]$glmnet$differentialFeatures$cv.min) == "character") { #only 1 cluster ?!
      clusters <- as.numeric(rds$conditionRegressionResults[[short_name]]$glmnet$differentialFeatures$cv.min[2])
    } else {
      clusters <- rds$conditionRegressionResults[[short_name]]$glmnet$differentialFeatures$cv.min$clusters
    }
    expr <- as.data.frame(rds$citrus.combinedFCSSet$data)
    cat("supplied Citrus RDS file contains", length(clusters), "clusters significant via glmnet:", cvmin_clusters )
    cat("\n extracting expression values...")
    events <- members[clusters]
    names(events) <- clusters
    
  } else {
  # cluster_id supplied, just extract what the user asked for
    events <- members[cluster_id]
    names(events) <- cluster_id  
  }
  
  mat <- list()
  
  for (i in 1:length(events)){
    idxs <- events[[i]]
    x <- expr[idxs, chnls]
    x <- expr[idxs, c(chnls,"fileEventNumber", "fileId")]
    x <- merge(x, pd, by="fileId")
    x <- x[,c(chnls, "ID")]
    
    # use MFI instead of raw
    x <- aggregate(.~ID, data=x, median, na.rm=TRUE)
    x <- merge(x, pd, by="ID")
    
    cluster_name <- names(events)[i]
    mat[[i]] <- x
    names(mat)[i] <- names(events)[i]
    
    # name rows with PTIDs, colnames with friendly marker names
    rownames(mat[[i]]) = mat[[i]]$ID
    setnames(mat[[i]], as.vector(markers), names(markers))
  }
  # `mat` is a list of matrices, named for the Citrus cluster IDs
  return(mat)
}


