geneCompare <- function (dat, geneList, grp1name="group1", grp2name="group2", grp1idx, grp2idx, plot=FALSE){
  # compares the mean expression between two groups of genes in a given 
  # expression matrix, optionally makes boxplots of the comparison, faceted by gene
  # used in a pipeline, so checks are pretty rudimentary
  
  # example: out <- summarizeGenelist(expressionMatrix,
  #                                   c("LAP3","NCOA7","IFIT2","STAT1","TRIM21"), 
  #                                   grp1name = "control", 
  #                                   grp2name = "experimental", 
  #                                   which(grepl("control", colnames(expressionMatrix))), 
  #                                   which(grepl("experimental", colnames(expressionMatrix)))  
  #                                   )

  require(ggplot2)
  require(data.table)
  
  # check that at least some of data has genes in the specified geneList
  if (length(intersect(geneList, rownames(dat))) > 0){
      
    dat <- subset(dat, rownames(dat) %in% geneList)
      # check that all columns have group assignments and the groups don't overlap
      if (length(grp1idx) + length(grp2idx) == ncol(dat) && !identical(grp1idx, grp2idx)){
        rows <- rownames(dat)
        cols <- colnames(dat)
        idxs <- cols
        idxs[grp1idx] <- grp1name
        idxs[grp2idx] <- grp2name
        idx <- data.frame(cbind(cols,idxs))
        d <- data.table(merge(melt(dat), idx, by.x="Var2", by.y="cols"))
        setnames(d, c("Var1","Var2","value", "idxs"), c("geneName", "inputColumn","expr","group"))
        
        tab <- merge(subset(d, group==grp1name)[,list(group1expr=mean(expr)),by=geneName],
                     subset(d, group==grp2name)[,list(group2expr=mean(expr)),by=geneName],
                     by="geneName")
        setnames(tab, c("group1expr", "group2expr"), c(grp1name,grp2name))
        tab$higher <- tab[,..grp2name] > tab[,..grp1name]
        tab$higher <- gsub("FALSE", grp1name, tab$higher)
        tab$higher <- gsub("TRUE", grp2name, tab$higher)
        
        if(plot) {
          ggplot(d, aes(x=group, y=expr, col=group)) +
            geom_boxplot(outlier.shape = NA) +
            facet_wrap(~ geneName, scales="free") +
            theme(legend.position = "none") +
            labs(title=paste0(length(rows), " genes in list"))
        }
        return(tab)
        
    } else { stop ("must specify one unique index for each column of data !!!") }
    
    colnames(d) <- c("gene", "group", "expression")
  } else { stop ("geneList does not match ANY columns in data !!!") }
}


categoryCompare <- function(dat, setsIndices, categories, grp1name="group1", grp2name="group2", grp1idx, grp2idx){
  
  require(data.table)
  
  if (length(intersect(geneList, rownames(dat))) > 0){
    tab <- data.table(cbind(categories, rep(1, length(categories)), rep(0, length(categories))))
    setnames(tab, c("categories", "V2","V3"), c("category", "totalGenes", "genesUp"))
    tab$totalGenes <- as.numeric(tab$totalGenes)
    tab$genesUp <- as.numeric(tab$genesUp)
    
    for (i in 1:length(categories)){
      cat <- categories[i]
      cidxs <- unlist(setsIndices[cat], use.names = F)
      tab[tab$category==cat]$totalGenes <- length(cidxs)
      d <- dat[cidxs,]
      
      if (length(grp1idx) + length(grp2idx) == ncol(d) && !identical(grp1idx, grp2idx)){
        rows <- rownames(d)
        cols <- colnames(d)
        idxs <- cols
        idxs[grp1idx] <- grp1name
        idxs[grp2idx] <- grp2name
        idx <- data.frame(cbind(cols,idxs))
        d <- data.table(merge(melt(d), idx, by.x="Var2", by.y="cols"))
        setnames(d, c("Var1","Var2","value", "idxs"), c("geneName", "inputColumn","expr","group"))
        
        row <- merge(subset(d, group==grp1name)[,list(group1expr=mean(expr)),by=geneName],
                     subset(d, group==grp2name)[,list(group2expr=mean(expr)),by=geneName],
                     by="geneName")
        tab[tab$category==cat]$genesUp <- length(which(row$group2expr - row$group1expr < 0))
      }
      
    }
    tab$proportionUp <- round(tab$genesUp / tab$totalGenes, 3)
    setnames(tab, c("genesUp"), c(paste0("genesUpIn", group2)))
    return(tab)
    
  } else { stop ("geneList does not match ANY columns in data !!!") }
    
  
  
}