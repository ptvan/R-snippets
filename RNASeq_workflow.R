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

workDir <- "~/RNASeq/data"
setwd(workDir)

# read in counts and annotation CSV
anno <- read.table("TB_Annotation.csv", sep=",", header = T)
dat <- read.table("TBData.txt", sep="\t", header=T)

# create ExpressionSet
dat <- ExpressionSet(assayData=as.matrix(dat))

# clean and format annotation
anno <- anno[,c(2,3,4,7,8)]
colnames(anno) <- c("sample_name", "TB_status", "TST_status", "subject", "sample_type")
anno$TB_status <- factor(anno$TB_status)
anno$TST_status <- factor(anno$TST_status)
anno$sample_type <- factor(anno$sample_type)

# create model matrix
mmatrix <- model.matrix(~0+TST_status, data=anno)

# normalize using voom
normy <- calcNormFactors(dat)
libNorm <- colSums(exprs(dat))*normy
v <- voom(exprs(dat), design=mmatrix, plot=FALSE, lib.size=libNorm)
ranCor <- duplicateCorrelation(v, design=mmatrix, block=anno$subject)$consensus.correlation

# fit different models
fit1 <- lmFit(v, mmatrix, block=anno$subject, correlation=ranCor)
