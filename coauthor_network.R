library(rentrez)
library(igraph)
library(ggnetwork)
library(ggraph2)
library(visNetwork)

# search PubMed for PMIDs of my papers, exclude spurious hits
hits <- entrez_search(db="pubmed", term="Phu Van")
titles <- extract_from_esummary(entrez_summary("pubmed", hits$ids), "title")
titles <- titles[!grepl("leaf sprays", titles)]
pmids <- names(titles)

# get the authors of all my papers, harmonize names
authors <- extract_from_esummary(entrez_summary("pubmed", pmids), "authors")
authors[1,]$`29868771`[1] <- "Van PT"
authors[1,]$`24935033`[5] <- "Minden JS"
authors[1,]$`19536208`[14] <- "Martin DB"
authors <- authors[1,]

# make one graph for publications
# label to be displayed can be either PMID or title
g <- make_empty_graph(n = length(pmids), directed = FALSE) 
V(g)$type <- "publication"
V(g)$PMID <- pmids
V(g)$title <- titles
V(g)$display <- titles

# make another graph for authors
# here the label to be displayed is the author name
g2 <- make_empty_graph(n = length(unique(unlist(authors))), directed = FALSE) 
V(g2)$type <- "author"
V(g2)$author <- unique(unlist(authors))
V(g2)$display <- unique(unlist(authors))

# merge the two graphs
g <- disjoint_union(g, g2)

# added the edges
for (i in 1:length(authors)){
  pub <- names(authors)[i]
  cat(pub, ":")
  for (j in 1:length(authors[[pub]])){
    # cat(authors[[i]][j], "\n")
    g <- add_edges(g, c(which(vertex_attr(g, "author") == authors[[i]][j]), which(vertex_attr(g, "PMID") == pub )))
  }
}

# diagnostic plot to make sure the network contains what we added
plot(g, vertex.label=V(g)$display)

# pretty plot using ggnetwork and ggraph
ggplot(ggnetwork(g, layout="fruchtermanreingold"), aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_edges() +
  geom_nodes(aes(color=type),size = 10) +
  geom_nodetext(aes(label=display), size=3) +
  theme_blank()

# for interactive graphs
visIgraph(g) %>%
  visNodes()

