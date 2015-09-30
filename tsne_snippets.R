require(Rtsne)
require(ggplot2)

# plots the output of a tsne row on ICS flow data, highlights user-specified window, and return 
# cell events within that window

extract_tsne_region <- function(dat, x_min, x_max, y_min, y_max, stim, marker) {
  mrkr <- as.symbol(marker)
  tmp <- subset(dat, stim==stim&x>x_min&x<x_max&y>y_min&y<y_max)
  plotCall <- substitute(ggplot(dat, aes(x = x, y = y)) +
    geom_point(aes(colour = cl), size = 1.5, alpha = 0.5) +
     facet_grid(neut_status ~ .) + 
    scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) +
    theme_bw() +
    geom_rect(xmin=x_min,xmax=x_max,ymin=y_min,ymax=y_max, color="red", alpha=0) + 
    ggtitle()
    
  ,
  list(cl = mrkr)
  )
  p <- eval(plotCall)
  
  print(p)
  return(tmp)
}

