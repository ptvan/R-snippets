library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(dbscan)
library(sp)
library(pheatmap)

# plots the output of a tsne row on ICS flow data, highlights user-specified window, and return 
# cell events within that window

extract_tsne_rect <- function(dat, x_min, x_max, y_min, y_max, antigen, facet) {
  fct <- as.symbol(facet)

  dat <- subset(dat, stim==antigen)
  tmp <- subset(dat, x>x_min&x<x_max&y>y_min&y<y_max)
  
  plotCall <- substitute(ggplot(subset(dat, stim==antigen), aes(x = x, y = y)) + 
                           geom_point(aes(colour = degree), size = 1, alpha = 0.8) +
                           geom_rect(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, 
                                     color="red",alpha=0) +
                           labs(title = "tSNE comparison") +
                           facet_grid(fg ~ .) + 
                           scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) 
                           
                         ,
                         list(fg = fct)
  )
  p <- eval(plotCall)
  print(p)
  return(tmp)
}


cluster_with_dbscan <- function(dat, epsilon, minpts = 5 ) {
  out <- dbscan(dat[,.(x,y)], eps = epsilon)
  dat[,"dbscan_cluster"] <- out$cluster
  return(dat)
}
