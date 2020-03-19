library(uwot)
library(Rtsne)
library(kernlab) # for kernelpca
library(keras) # for autoencoders
library(MASS) # for Sammon mapping
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

# T-SNE
tsne_out <- Rtsne(dat[,c(1:8)],perplexity=30,theta=0.0) 

tsne_out <- tsne_out$Y %>%
          as.data.frame() %>% 
          set_colnames(c("X","Y")) %>%
          inset("class", value=dat$class) %>%
          mutate(class = as.factor(class))

# UMAP
umap_out <- umap(dat, init = "spca") %>%
          as.data.frame() %>% 
          set_colnames(c("X","Y")) %>%
          inset("class", value=dat$class) %>%
          mutate(class = as.factor(class))

ggplot(umap_out) +
      aes(x=X, y=Y, col=class) +
      geom_point()


# KERNEL PCA
kpca_out <- kpca(dat, features=2) 
rotated(kpca_out)



