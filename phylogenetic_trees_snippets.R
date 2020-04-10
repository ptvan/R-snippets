library(ape)
library(treeio)
library(tidytree)
library(ggtree)
library(Biostrings)

setwd("~/working/workflows/")

# read in the aligned sequences
COVIDseqs <- readDNAStringSet("COVID19multi.phy")
names(COVIDseqs) <- gsub("\\|","_", names(COVIDseqs))
write.dna(COVIDseqs, "COVID19MSA.fa", "fasta")

# calculate Hamming distance
tipseq_dist <- stringDist(COVIDseqs, method = "hamming")

COVIDtree <- read.newick("COVID19multi.phy.treefile")

# drop a duplicate tip from the tree
# COVIDtree <- drop.tip(COVIDtree, COVIDtree$tip.label[which(grepl("MT042773", mytree$tip.label))])

#
# ggtree(COVIDtree, layout="circular") +
#   geom_tiplab(size=3, color="blue")

# plot the tree along the multiple sequence alignment
p <- ggtree(COVIDtree, branch.length='none') + 
  geom_tiplab() + xlim_tree(5.5) 

msaplot(p, "COVID19MSA.fa")
