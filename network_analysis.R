make_gsea_overlap_igraph <- function(expressionList, cameraMat, geneSets, verbose=FALSE){
  # Phu T. Van, FHCRC 2017, w/ substantial help and input from C.Murie & V.Voillet
  
  # this function takes a bioConductor EList (eg. voom-transformed RNASeq counts),
  # GSEA output matrix from limma's camera() and a user-specified geneSet
  # and create an igraph graph where the nodes are gene lists and the edge weight
  # is how many expressed genes are common between the gene lists
  
  # example usage:
  # vDat <- voom(exprs(eDat), design=designMat, plot=FALSE, lib.size=libNorm)
  # res <- camera(vDat, setsIndices, design=designMat, contrast=cons[i], sort=TRUE)
  # geneSets <- geneIds(getGmt(gmtFile))
  # n <- make_gsea_igraph(vDat, res, geneSets)
  # plot(n, edge.width=E(n)$overlap*.1, vertex.color="white", vertex.label.cex=0.5)
  
  require(igraph)
  require(stringr)
  
  # get the geneset categories, and the expressed genes
  geneSetCats <- rownames(cameraMat)
  expressedGenes <- rownames(expressionList$E)
  
  # create unconnected graph with geneSet categories as vertices
  net <- make_empty_graph(n = nrow(cameraMat), directed=FALSE)
  
  # add vertex attributes
  V(net)$label <- geneSetCats
  V(net)$totalGeneCount <- unlist(lapply(geneSets[rownames(cameraMat)], length))
  V(net)$sigGeneCount <- cameraMat$NGenes
  V(net)$direction <- cameraMat$Direction
  V(net)$FDR <- cameraMat$FDR
  V(net)$pval <- cameraMat$PValue
  
  # create list of edges
  tmp <- list()
  for (i in 1:length(geneSetCats)){
    if (verbose){
      cat("working on", geneSetCats[i], ":", length(unlist(geneSets[geneSetCats[i]])), "total genes,"
          ,length(which(expressedGenes%in%unlist(geneSets[geneSetCats[i]]))), "of which were expressed \n"
      )}
    tmp[[i]] <- expressedGenes[which(expressedGenes%in%unlist(geneSets[geneSetCats[i]]))]
    names(tmp)[i] <- geneSetCats[i]
  }
  nms <- combn( names(tmp) , 2 , FUN = paste0 , collapse = "xxx" , simplify = FALSE )
  ll <- combn( tmp , 2 , simplify = FALSE )
  edges <- lapply( ll , function(x)  intersect( x[[1]] , x[[2]] )  )
  names(edges) <- nms
  edges <- edges[lapply(edges, length)>0]
  
  # add edges from edge list to  graph
  for (i in 1:length(edges)){
    nodes <- unlist(str_split(names(edges)[i],"xxx"))
    net <- add_edges(net, c(which(nodes[1] == V(net)$label), which(nodes[2] == V(net)$label))) 
    net <- set_edge_attr(net, "overlap", index=i, value=length(edges[[i]]))
  }
  
  return(net)
}  

make_STRING_igraph <- function(inputMatrix, STRINGdbObj){
  # Phu T. Van, FHCRC 2017, w/ substantial help and input from C.Murie & V.Voillet
  
  # this function takes an input matrix having one column named "gene" containing genes of interest and a stringDB object
  # it returns an igraph graph of the STRING network containing those genes
  # the node attributes are "name" (STRING id) and "geneName" (original gene names from input matrix)
  # genes that do not have STRING annotation are omitted
  
  # example usage : 
  # DEGout <- topTable(fitBayes, number=nrow(vDat), coef="pttype", sort="P")
  # DEGout$gene <- rownames(DEGout)
  # string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory="/STRINGdb/HomoSapiens/")
  # n <- make_STRING_igraph(DEGout, string_db)
  # plot(n, vertex.label=V(n)$geneName)
  
  require(STRINGdb)
  require(igraph)
  
  # get gene names <-> STRING_id mapping
  smap <- STRINGdbObj$map(inputMatrix, "gene", removeUnmappedRows = TRUE )
  
  # get network of all STRING edges (HUGE!!!)
  stringNet <- STRINGdbObj$get_graph()
  
  # trim the network down to just our genes and potential edges
  smap <- subset(smap, STRING_id%in%V(stringNet)$name)
  g <- induced_subgraph(stringNet, v=smap$STRING_id)
  
  # add the gene names to the vertexes since STRING ids are non-informative
  V(g)$geneName <- smap[match(vertex_attr(g, "name"), smap$STRING_id),]$gene
  
  return(g)
}

# Make igraph from adjacency matrix
# library(WGCNA)
# library(igraph)
# 
# adjMat  <- adjacency( t(vDat$E[which(rownames(vDat$E)%in%geneList$gene),]))
# diag(adjMat) <- 0
# g <- graph_from_adjacency_matrix(adjMat, weighted=T, mode=c("undirected"))
# 
# # see what the correlation looks like for our genes and filter
# # hist(edge_attr(g, "weight"))
# g2 <- subgraph.edges(g, E(g)[weight>0.1])
# V(g2)$logFC <- geneList[match(vertex_attr(g2, "name"), geneList$gene),]$logFC
# V(g2)$geneName <- V(g2)$name
# 
# 
# # plotting 
# pu <- colorRampPalette(c("purple","mediumorchid","purple4"))(20)                      
# ye <- colorRampPalette(c("lemonchiffon", "khaki","yellow"))(20)
# 
# ggplot(ggnetwork(g2, layout="fruchtermanreingold"), aes(x=x,y=y,xend=xend,yend=yend)) +
#   geom_edges(aes(size = weight), color = "grey80") +
#   geom_nodes(aes(color=logFC), size = 10) +
#   scale_color_gradientn(colours=c(pu,"black", ye), na.value = "grey98", limits = c(-3, 3)) +
#   geom_nodetext(aes(label=geneName)) +
#   labs(title="WGCNA network of DEGs (pttype x stim), adjacency > 0.1") +
#   theme_blank()
# 
# library(PCIT)
# corMat <- cor(t(vDat$E[which(rownames(vDat$E)%in%geneList$gene),]))
# pcitMemoryRequirement(nrow(corMat), units="MB")
# pcit(corMat)

# filter unconnected vertices
# net <- induced_subgraph(net, v=which(igraph::degree(g=net, v=V(net))>1))


