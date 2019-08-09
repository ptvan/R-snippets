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
library(org.Hs.eg.db)
library(topGO)

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
attrs <- c("hgnc_symbol", "external_gene_name", "chromosome_name", "gene_biotype")


# NOTE: this query can be a bit slow, so likely we'll want to save results
# rather than running it everytime
chrMap <- getBM(attributes=attrs
                ,filters=c("chromosome_name")
                ,values=list(chrom)
                ,mart=ensembl)

# get only protein coding genes
proteinCodingGenes <- subset(chrMap, gene_biotype=="protein_coding")$external_gene_name

## alternatively, we can fetch the gene names from a GTF file, eg. from UCSC hg38
library(refGenome)
hg38gtf <- ensemblGenome()
read.gtf(hg38gtf, "UCSC_hg38_genes.gtf"
         , useBasedir = FALSE)

allGenes <- unique(hg38gtf@ev$gtf$gene_name)

exon <- extractFeature(hg38gtf, "exon")
exonicGenes <- unique(exon@ev$gtf$gene_name)

# mapping mouse <-> human genes via ENSEMBL gene_id
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 

# extra orthology columns so we can filter the mouse genes if necessary
attrs <- c("ensembl_gene_id"
           ,"mmusculus_homolog_ensembl_gene"
           ,"mmusculus_homolog_perc_id_r1"
           ,"mmusculus_homolog_orthology_type"
           ,"mmusculus_homolog_subtype"
           ,"mmusculus_homolog_perc_id"
           ,"mmusculus_homolog_associated_gene_name"
) 

mouse <- getBM( attrs
                ,filters="with_mmusculus_homolog"
                ,values =TRUE
                , mart = human
                , bmHeader=FALSE)

# two separate getBM() calls since BioMart doesn't allow
# gene- and transcript-level queries in the same call
attrs2 <- c("ensembl_gene_id"
            ,"hgnc_symbol"
            ,"external_gene_name")

E2GN <- getBM( attrs2
               ,filters="with_mmusculus_homolog"
               ,values =TRUE
               , mart = human
               , bmHeader=FALSE)

map <- merge(mouse, E2GN, by="ensembl_gene_id")
map$mmusculus_homolog_associated_gene_name <- toupper(map$mmusculus_homolog_associated_gene_name)
map[c("hgnc_symbol","mmusculus_homolog_associated_gene_name")]


# these are the genes actually in our data matrix
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


### Gene Set Enrichment Analysis (GSEA) using CAMERA
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

### GeneOntology analysis using topGO
# which in turn requires annotation from org.Hs.eg.db
# read in a table of genes, logFC, AveExpr, t, P.Value, adj.P.Val and B 
# from limma's topTable()
favGenes <- readRDS("DEG_FDR020.Rds")

# topGO expects an named vector of p.values, and a function to select top genes
# to create a GOdata object
GOinput <- favGenes$P.Value
names(GOinput) <- rownames(favGenes)
selection <- function(allScore){ return(allScore < 0.05)}
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")

GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=GOinput,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=10)

# run test on on the GOdata object, tests are listed by topGO:whichTests()
# which includes Fisher, Kolmogoov-Smirnov (KS) and others
# these tests are for the GO enrichment, and *not* for differential expression 
GOoutput <- runTest(GOdata
                    , algorithm = "classic"
                    , statistic = "ks")

GOtable <- GenTable(GOdata, KS=GOoutput, orderBy="KS", topNodes=20)
