library(Biobase)
library(data.table)
library(ggplot2)
library(gridExtra)
library(limma)
library(edgeR)
library(DESeq2)
library(xtable)
library(GSEABase)
library(DT)
library(plyr)     
library(stringr)  
library(plotly)
library(cowplot)
library(biomaRt)

workDir <- "~/RNASeq/data"
setwd(workDir)

# read in counts and annotation CSV
anno <- read.table("TB_Annotation.csv", sep=",", header = T)
dat <- read.table("TBData.txt", sep="\t", header=T)

# create ExpressionSet
dat <- ExpressionSet(assayData=as.matrix(dat))

# clean and format annotation
# here we have two contrasts: "TB vs. no TB" and "TSTneg vs. TSTpos"
anno <- anno[,c(2,3,4,7,8)]
colnames(anno) <- c("sample_name", "TB_status", "TST_status", "subject", "sample_type")
anno$TB_status <- factor(anno$TB_status)
anno$TST_status <- factor(anno$TST_status)
anno$sample_type <- factor(anno$sample_type)

# create model matrix and corresponding labels for contrasts
mmatrix <- model.matrix(~0+TST_status, data=anno)
cons <- colnames(mmatrix)[-1] 
labs <- c("pttype", "gender", "pttype*gender")


# normalize using voom
normy <- calcNormFactors(dat)
libNorm <- colSums(exprs(dat))*normy
v <- voom(exprs(dat), design=mmatrix, plot=FALSE, lib.size=libNorm)
ranCor <- duplicateCorrelation(v, design=mmatrix, block=anno$subject)$consensus.correlation

# look up gene symbols' chromosomal location on ENSEMBL
# explicitly set which server we use, since mirrors can go down for maintenance

ensemblMart <- useMart("ensembl"
              ,host = "www.ensembl.org"
              ,ensemblRedirect = FALSE)
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensemblMart)

chrom <- c(1:22,"X","Y")
attrs <- c("hgnc_symbol", "external_gene_name", "chromosome_name")

# NOTE: this query can be a bit slow, so likely we'll want to save results
# rather than running it everytime
chrMap <- getBM(attributes=attrs
                ,filters=c("chromosome_name")
                ,values=list(chrom)
                ,mart=ensembl)


expressedGenes <- rownames(v$E)

# get genes on Y chromosome, remove them from data matrix ...
chrYgenes <- intersect(expressedGenes, subset(chrMap, chromosome_name=="Y")$external_gene_name)
eDatnoY <- ExpressionSet(assayData=as.matrix(counts[!rownames(counts) %in% chrYgenes,]))

# ... and see how the data separates without sex genes
plotMDS(vDatnoY, col=rainbow(length(unique(colTB)))[as.numeric(colTB)],  main="TB status, chrY genes removed", pch=1)
legend("bottomleft", c("noTB","TB"), pch=1, col=c("#FF0000FF", "#00FFFFFF"))


# make the contrasts matrix
aovCon <- makeContrasts(status=(negTST - posTST),
                        levels=mmatrix)

### Differentially Expressed Genes (DEG) analysis
# fit different models
fit1 <- lmFit(v, mmatrix, block=anno$subject, correlation=ranCor)
fit2 <- contrasts.fit(fit1, aovCon)
fit2 <- eBayes(fit2, trend=FALSE)

# generate the table of top DEGs from the linear model fit
allOut <- list()
for(i in 1:ncol(aovCon)) {
  allOut[[i]] <- topTable(fit2, number=nrow(v), coef=i, sort="P")
}


### Gene Set Enrichment Analysis (GSEA)
# load GMT files (http://software.broadinstitute.org/gsea/msigdb/collections.jsp)
# GMT files can be concatenated together
FDRCut <- 0.2

gmtFile <- "c2c7.concatenated.v6.0.symbols.gmt"
minGeneSetCut <- 5

geneSet <- getGmt(gmtFile)
geneIds <- geneIds(geneSet)
setsIndices <- ids2indices(geneIds, rownames(vDat$E))
setsIndices <- setsIndices[sapply(setsIndices, length) > minGeneSetCut]

combo <- comboTab <- list()
comboGSEA <- list()
for(i in 1:length(cons)) {
  res <- camera(vDat, setsIndices, design=designMat, contrast=cons[i], sort=TRUE)
  indo <- res$FDR <= FDRCut
  if(sum(indo) > 0) {
    combo[[labs[i]]] <- res[indo, c("Direction","PValue","FDR")]
    comboTab[[labs[[i]]]] <- datatable(res[indo, c("Direction","PValue","FDR")], caption=labs[i]
                                       , extensions = 'Buttons'
                                       , options = list(dom = 'Bfrtip',
                                                        buttons = c('csv', 'excel'))
    )
  }
  comboGSEA[[i]] <- res
}