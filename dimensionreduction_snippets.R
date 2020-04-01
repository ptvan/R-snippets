library(keras) # for autoencoders
library(diffusionMap)
library(tidyverse)
library(magrittr)
library(ggplot2)

# reference: Dimensionality Reduction: A Comparative Review, van der Maaten et al 2009 
# https://lvdmaaten.github.io/publications/

# HTRU2 from the UCI ML repository
# https://archive.ics.uci.edu/ml/datasets/HTRU2
dat <- read.csv("~/working/datasets/HTRU.csv")
colnames(dat) <- c(paste0("var", 1:(ncol(dat)-1)), "class")
dat <- dat[sample(nrow(dat),5000),]


##########################
# classical (metric) MDS
##########################
# needs a distance structure
# k is number of dimensions
d <- dist(dat)
mds <- cmdscale(d, eig=TRUE, k=2)
plot(mds$points[,1], mds$points[,2], xlab="Coordinate 1", ylab="Coordinate 2",
     main="metric MDS", type="n")
text(mds$points[,1], mds$points[,2], labels = row.names(dat), cex=.7) 

################# 
# non-metric MDS
################# 
library(MASS)
nmds <- isoMDS(d, k=2)
plot(nmds$points[,1], nmds$points[,2], xlab="Coordinate 1", ylab="Coordinate 2",
     main="non-metric MDS", type="n")
text(nmds$points[,1], nmds$points[,2], labels = row.names(dat), cex=.7) 


###########
# T-SNE
###########
library(Rtsne)
tsne_out <- Rtsne(dat[,c(1:8)],perplexity=30,theta=0.0) 

tsne_out <- tsne_out$Y %>%
          as.data.frame() %>% 
          set_colnames(c("X","Y")) %>%
          inset("class", value=dat$class) %>%
          mutate(class = as.factor(class))

#######
# UMAP
#######
library(uwot)
umap_out <- umap(dat, init = "spca") %>%
          as.data.frame() %>% 
          set_colnames(c("X","Y")) %>%
          inset("class", value=dat$class) %>%
          mutate(class = as.factor(class))

ggplot(umap_out) +
      aes(x=X, y=Y, col=class) +
      geom_point()

############
# KERNEL PCA
############
library(kernlab) 
kpca_out <- kpca(dat, features=2) 
rotated(kpca_out)

############
# SAMMON
###########
library(MASS)
sammon_out <- sammon(dist(dat[,-180]))
plot(sammon_out$points, type="n")
text(sammon_out$points, labels = as.character(1:nrow(dat[,-180])))

