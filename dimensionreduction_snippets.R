library(uwot)
library(Rtsne)
library(kernlab)
library(keras) # for autoencoders
library(MASS) # for Sammon mapping
library(diffusionMap)
library(tidyverse)
library(magrittr)
library(ggplot2)

# reference: Dimensionality Reduction: A Comparative Review, van der Maaten et al 2009 
# https://lvdmaaten.github.io/publications/

# seizure dataset from the UCI ML repository
# https://archive.ics.uci.edu/ml/datasets/Epileptic+Seizure+Recognition

dat <- read.csv("~/working/datasets/seizure/data.csv")
colnames(dat) <- gsub("X", "var", colnames(dat))

umap_out <- umap(dat, init = "spca") %>%
          as.data.frame() %>% 
          set_colnames(c("X","Y")) %>%
          inset("class", value=dat$y) %>%
          mutate(class = as.factor(class))
          

ggplot(umap_out) +
      aes(x=X, y=Y, col=class) +
      geom_point()