# Make igraph from adjacency matrix
library(WGCNA)
library(igraph)
library(ggraph2)
library(ggnetwork)

######################
# GENERATE DUMMY DATA
######################

# get a bunch of human gene names
genes <- unlist(lookUp(as.character(1:50000), 'org.Hs.eg', 'SYMBOL')) 
genes <- sample(unique(genes[!is.na(genes)]))
names(genes) <- NULL

# generate some dummy subjects, p001 to p050
ptids <- paste0("p",str_pad(as.character(1:50), 3, "left", "0"))

# generate dummy RNASeq counts
counts <- matrix(sample(1:1e6, length(genes)*length(ptids)), nrow=length(genes), ncol=length(ptids)*2)

# each subject has "MEDIA" and "STIM" samples
samples <- c(paste0(ptids,"_MEDIA"), paste0(ptids,"_STIM"))
colnames(counts) <- samples
rownames(counts) <- genes

# create ExpressionSet
dat <- ExpressionSet(assayData=as.matrix(counts))

######################
# CREATE ADJACENCY MATRIX
######################

adjMat  <- adjacency( t(dat[which(rownames(dat)%in%geneList$gene),]))
diag(adjMat) <- 0
g <- graph_from_adjacency_matrix(adjMat, weighted=T, mode=c("undirected"))

# see what the correlation looks like for our genes and filter
# hist(edge_attr(g, "weight"))

g2 <- subgraph.edges(g, E(g)[weight>0.1])
V(g2)$logFC <- geneList[match(vertex_attr(g2, "name"), geneList$gene),]$logFC
V(g2)$geneName <- V(g2)$name

######################
# PLOTTING
######################
pu <- colorRampPalette(c("purple","mediumorchid","purple4"))(20)
ye <- colorRampPalette(c("lemonchiffon", "khaki","yellow"))(20)

ggplot(ggnetwork(g2, layout="fruchtermanreingold"), aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_edges(aes(size = weight), color = "grey80") +
  geom_nodes(aes(color=logFC), size = 10) +
  scale_color_gradientn(colours=c(pu,"black", ye), na.value = "grey98", limits = c(-3, 3)) +
  geom_nodetext(aes(label=geneName)) +
  labs(title="WGCNA network of DEGs (pttype x stim), adjacency > 0.1") +
  theme_blank()

# using Partial Correlation Coefficient with Information Theory (PCIT)
library(PCIT)
corMat <- cor(t(vDat$E[which(rownames(vDat$E)%in%geneList$gene),]))
pcitMemoryRequirement(nrow(corMat), units="MB")
pcit(corMat)

# filter unconnected vertices
net <- induced_subgraph(net, v=which(igraph::degree(g=net, v=V(net))>1))

############################
# CREATE NETWORK FROM "SEED" GENES
############################

# create a network by expanding some "seed" genes with their KEGG-annotated pathway co-members
myGenes <- c("TLR1"
               ,"TLR2"
               ,"TLR4"
               ,"TLR6"
               ,"NOD2"
)

# hsa2path <- keggLink("hsa", "pathway")
# path2hsa <- keggLink("pathway", "hsa")

gn2entrez <- ldply(mget(myGenes, envir=org.Hs.egSYMBOL2EG, ifnotfound=NA), data.frame)
colnames(gn2entrez) <- c("geneName", "entrezNum")
gn2entrez <- gn2entrez[complete.cases(gn2entrez),]

entrez2kegg <- ldply(mget(gn2entrez$entrezNum, envir=org.Hs.egPATH, ifnotfound=NA), data.frame)
colnames(entrez2kegg) <- c("entrezNum", "keggPathway")
myPathways <- paste0("hsa", as.character(unique(entrez2kegg$keggPathway)))

keggMatches <- NULL

for (i in 1:length(myPathways)){
  s <- keggGet(myPathways[i])[[1]]$GENE
  s <- s[grep(";", s)]
  s <- gsub("; .+$", "\\1", s)
  keggMatches <- c(keggMatches, s)
}
