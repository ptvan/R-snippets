library(ape)
library(treeio)
library(tidytree)
library(ggtree)

setwd("~/working/workflows/")

mytree <- read.newick("multiDNA.phy.treefile")

ggtree(mytree, layout="circular") +
  geom_tiplab(size=3, color="blue")
