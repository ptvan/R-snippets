library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(dbscan)
library(sp)
library(pheatmap)

# plots the output of a tsne run on ICS flow data, highlights user-specified window, and return 
# cell events within that window. Also calls tsne_proportion_boxplot() to perform t.test() on
# extracted events

extract_tsne_rect <- function(dat, x_min, x_max, y_min, y_max, antigen, marker, facet, title) {
  facet <- as.symbol(facet)
  marker <- as.symbol(marker)
  
  # `marker` argument can be an actual marker (eg. "IFNg") or "degree" 
  
  dat <- subset(dat, stim==antigen)
  tmp <- subset(dat, x>x_min&x<x_max&y>y_min&y<y_max)
  
  plotCall <- substitute(ggplot(dat, aes(x = x, y = y)) + 
                           geom_point(aes(colour = mkr), size = 1, alpha = 0.8) +
                           geom_rect(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, 
                                     color="red",alpha=0) +
                           labs(title = title) +
                           facet_grid(fct ~ .) + 
                           scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) 
                         ,
                         list(fct = facet, mkr = marker)
  )
  
  p1 <- eval(plotCall)
  p2 <- tsne_proportion_boxplots(dat, x_min, x_max, y_min, y_max, antigen, marker, facet)
  p3 <- grid.rect(gp=gpar(col="white"))
  pl <- grid.arrange(p1, arrangeGrob(p2,p3, nrow=1), heights=c(1,0.5))
  print(pl)
  return(tmp)
}

# performs a t.test() on extracted events
tsne_proportion_boxplots <- function(dat, x_min, x_max, y_min, y_max, antigen, marker, facet) {
  m <- unique(dat[,.(ptid,neut_status,controller_status)])
  totals <- dat[, length(x), by=ptid]
  setnames(totals, "V1", "total")
  
  dat <- subset(dat, stim==antigen)
  dat <- subset(dat, x>x_min&x<x_max&y>y_min&y<y_max)
  
  non_h <- subset(dat, neut_status =="non")
  neut_h <- subset(dat, neut_status =="neut")
  
  non_s <- non_h[, length(x), by=ptid]
  setnames(non_s, "V1", "in_hull")
  neut_s <- neut_h[, length(x), by=ptid]
  setnames(neut_s, "V1", "in_hull")
  
  subs <- rbind(non_s, neut_s)
  
  m <- merge(m, subs, by="ptid")
  m <- merge(m, totals, by="ptid")
  m[,proportion := in_hull/total]
  m[,log_proportion := log10(proportion)]
  
  p <- t.test(log_proportion ~ neut_status, data=m)
  
  pl <- ggplot(m, aes(x=neut_status, y=log(in_hull))) + 
    geom_boxplot(outlier.size = 3) +
    geom_point(aes(col=ptid), size=3) +
    theme(legend.position="none") +
    labs(title=paste0("Selected cell proportion, p=", round(p$p.value, 3)), x="") 
  
  print(m)
  return(pl)
}




cluster_with_dbscan <- function(dat, epsilon, minpts = 5 ) {
  out <- dbscan(dat[,.(x,y)], eps = epsilon)
  dat[,"dbscan_cluster"] <- out$cluster
  return(dat)
}
