library(ape)
library(treeio)
library(tidytree)
library(ggtree)
library(Biostrings)
library(stringr)

setwd("~/working/workflows/")

## read in the aligned sequences
# and create metadata
COVIDseqs <- readDNAStringSet("COVID19multi.phy")
meta <- data.frame(matrix(unlist(str_split(names(COVIDseqs), pattern = "\\|")), nrow=length(COVIDseqs), byrow=TRUE))
colnames(meta) <- c("genBank", "date", "location")
names(COVIDseqs) <- gsub("\\|","_", names(COVIDseqs))
meta$label <- names(COVIDseqs)

# write out the dates for treetime
dates <- meta[,c("label","date")]
colnames(dates) <- c("accession","date")
dates$date <- gsub("_", "-", dates$date)
write.csv(dates, file="COVID19dates.csv", row.names = FALSE, quote=FALSE)

# write the aligned sequences out as a FASTA
write.dna(COVIDseqs, "COVID19MSA.fa", "fasta")

## calculate Hamming distance
# tipseq_dist <- stringDist(COVIDseqs, method = "hamming")

## read IQ-TREE output
phylo <- read.newick("COVID19multi.phy.treefile")

## drop a duplicate tip from the tree
# phylo <- drop.tip(phylo, phylo$tip.label[which(grepl("MT042773", mytree$tip.label))])

## phylo is a S3 'phylo' objects, which has some useful methods
# plot(phylo)
# plot(phylo, type="unrooted")
# nodelabels()
# tiplabels()
# rerooted <- root(phylo)
# all.equal(rerooted, phylo)

## tidytree implements a S4 'treedata' object
# since it is internally a tibble, we just full_join to merge metadata 
tidyTree <- full_join(as.treedata(phylo), meta, by="label")

## plotting using ggtree
# ggtree(tidyTree, layout="circular") +
# geom_tiplab(size=3, color="blue")

## with timescale
ggtree(tidyTree) +
  geom_tiplab() +
  theme_tree2() +
  xlim_tree(3) +
  geom_treescale(x=0, y=45, width=1, color='red')

## with metadata
ggtree(tidyTree, branch.length='none') +
  geom_tiplab(aes(label=genBank))

## plot the tree along the multiple sequence alignment
p <- ggtree(tidyTree, branch.length='none') + 
  geom_tiplab() + 
  xlim_tree(5.5) 
msaplot(p, "COVID19MSA.fa", offset=10)

## read the treetime output
# time tree
timed <- as.treedata(read.nexus("2020-04-10_treetime/timetree.nexus"))
ggtree(timed) +
  geom_tiplab() 

# divergence tree
div <- as.treedata(read.nexus("2020-04-10_treetime/divergence_tree.nexus"))
ggtree(div) +
  geom_tiplab() 
