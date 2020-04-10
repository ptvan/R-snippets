library(ape)
library(treeio)
library(tidytree)
library(ggtree)
library(Biostrings)
library(stringr)

setwd("~/working/workflows/")

## read in the aligned sequences
COVIDseqs <- readDNAStringSet("COVID19multi.phy")
meta <- data.frame(matrix(unlist(str_split(names(COVIDseqs), pattern = "\\|")), nrow=length(COVIDseqs), byrow=TRUE))
names(COVIDseqs) <- gsub("\\|","_", names(COVIDseqs))
write.dna(COVIDseqs, "COVID19MSA.fa", "fasta")
## calculate Hamming distance
# tipseq_dist <- stringDist(COVIDseqs, method = "hamming")


COVIDtree <- read.newick("COVID19multi.phy.treefile")
## drop a duplicate tip from the tree
# COVIDtree <- drop.tip(COVIDtree, COVIDtree$tip.label[which(grepl("MT042773", mytree$tip.label))])

## COVIDtree is a 'phylo' objects, which has some useful methods
# plot(COVIDtree)
# plot(COVIDtree, type="unrooted")
# nodelabels()
# tiplabels()
# rerooted <- root(COVIDtree)
# all.equal(rerooted, COVIDtree)

## plotting using ggtree
# ggtree(COVIDtree, layout="circular") +
#   geom_tiplab(size=3, color="blue")

ggtree(COVIDtree) +
  geom_tiplab() +
  theme_tree2() +
  xlim_tree(3) +
  geom_treescale(x=0, y=45, width=1, color='red')

## plot the tree along the multiple sequence alignment
p <- ggtree(COVIDtree, branch.length='none') + 
  geom_tiplab() + xlim_tree(5.5) 
msaplot(p, "COVID19MSA.fa", offset=10)
