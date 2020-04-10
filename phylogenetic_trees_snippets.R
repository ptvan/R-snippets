library(ape)
library(treeio)
library(tidytree)
library(ggtree)
library(Biostrings)

setwd("~/working/workflows/")
dnaData <- readDNAStringSet("multiDNA.phy")

# one sequence is 1 base short, which causes Hamming to fail
# so we drop its sequence
bad <- names(dnaData)[5]
dnaData <- dnaData[-5]

write.dna(dnaData, "out.fa", "fasta")
tipseq_dist <- stringDist(dnaData, method = "hamming")

# also drop the corresponding tip from the tree
mytree <- read.newick("multiDNA.phy.treefile")
mytree <- drop.tip(mytree, mytree$tip.label[which(grepl("Thermoproteus", mytree$tip.label))])


ggtree(mytree, layout="circular") +
  geom_tiplab(size=3, color="blue")

p <- ggtree(mytree, branch.length='none') + 
  geom_tiplab() + xlim_tree(5.5) +
  geom_cladelabel(1, "first clade") 

# plot the tree along the multiple sequence alignment
msaplot(p, "out.fa")