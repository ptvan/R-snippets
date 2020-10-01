############################################
# generating multiple linked omics datasets
############################################
library(InterSIM)

# proportion of samples in the clusters, length = K
prop <- c(0.20,0.30,0.27,0.23) 

# cluster mean shift, can be set for each dataset
effect <- 5

sim.data <- InterSIM(n.sample=500
                     , cluster.sample.prop = prop
                     , delta.methyl=effect
                     , delta.expr=effect
                     , delta.protein=effect
                     , p.DMP=0.2
                     , p.DEG=NULL
                     , p.DEP=NULL
                     , sigma.methyl=NULL
                     , sigma.expr=NULL
                     , sigma.protein=NULL
                     , cor.methyl.expr=NULL
                     , cor.expr.protein=NULL
                     , do.plot=TRUE
                     , sample.cluster=TRUE
                     , feature.cluster=TRUE)
names(sim.data)

#########################
# using the MOFA package
#########################
# NOTE: need to update below to MOFA v2
library(MOFA)

# simulate some data & create MOFA object
N <- 100
D1 <- 250
D2 <- 500

view1 <- matrix(rnorm(N*D1),nrow=D1, ncol=N)
view2 <- matrix(rnorm(N*D2),nrow=D2, ncol=N)
colnames(view1) <- paste0("sample_",1:N)
colnames(view2) <- paste0("sample_",1:N)
rownames(view1) <- paste0("feature",1:D1,"_view",1)
rownames(view2) <- paste0("feature",1:D2,"_view",2)

MOFAobject <- createMOFAobject(list("view1" = view1, "view2" = view2))
print(MOFAobject)
plotDataOverview(MOFAobject)

# get options
data_opts <- getDefaultDataOptions()
model_opts <- getDefaultModelOptions(MOFAobject)
train_opts <- getDefaultTrainOptions()

# build & train MOFA object
MOFAobject <- prepareMOFA(
  object = MOFAobject,
  DataOptions  = data_opts,
  ModelOptions =  model_opts,
  TrainOptions =  train_opts
)

outfile = "MOFAtest.hdf5"
MOFAobject.trained <- runMOFA(MOFAobject, outfile)


############################
# using the mixOmics package
############################
library(mixOmics)

# load example data
data(breast.TCGA)

X <- list(mRNA = breast.TCGA$data.train$mrna,
          miRNA = breast.TCGA$data.train$mirna,
          protein = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype
list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))

# run Projection to Latent Structure with sparse Discriminant Analysis
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)

# plot individuals and variables
plotIndiv(MyResult.diablo
          ,ind.names = FALSE
          ,legend=TRUE, cex=c(1,2,3)
          # title = 'BRCA with DIABLO'
          )

plotVar(MyResult.diablo
        , var.names = c(FALSE, FALSE, TRUE)
        ,legend=TRUE, pch=c(16,16,1)
        )

# relevance network
network(MyResult.diablo, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'),
        cutoff = 0.6, save = 'png', name.save = 'DIABLOnetwork')

# circos of correlation between data types
circosPlot(MyResult.diablo
           , line=TRUE
           , cutoff=0.7)

# heatmap of a single component
cimDiablo(MyResult.diablo
          , color.blocks = c('darkorchid', 'brown1', 'lightgreen')
          , comp = 1
          , margin=c(8,20)
          , legend.position = "right")

# loadings
plotLoadings(MyResult.diablo
             , comp = 1, contrib = "max")



############################
# using parallelized glm
############################
# many thanks to Chad Young for inspiring the code below
library(glmnet)
library(doParallel)
library(foreach)

# parallelize
numCores <- detectCores()
registerDoParallel(numCores)

# simulate a multi-omics dataset
N <- 50
D1 <- 500
D2 <- 500

genes <- matrix(rnorm(N*D1),nrow=D1, ncol=N)
proteins <- matrix(rnorm(N*D2),nrow=D2, ncol=N)
colnames(genes) <- paste0("sample_",1:N)
colnames(proteins) <- paste0("sample_",1:N)
rownames(genes) <- paste0("orf",1:D1,"gene")
rownames(proteins) <- paste0("orf",1:D2,"protein")
allData <- rbind(genes, proteins)
set.seed(100)
m <- as.data.frame(cbind(colnames(genes), sample(c(rep("HIV", N/2), rep("healthy", N/2)))))
colnames(m) <- c("sample","status")

# covariate
m$status <- as.factor(m$status)

# single gene/protein as predictor
results <- foreach (i = 1:nrow(allData),
                    .combine = bind_rows,
                    .packages = c("glmnet", "dplyr", "tidyr")) %dopar% {
              orf <- rownames(allData)[i]
              fit <- glm(allData[i,] ~ status, data=m)
              coefs <- coef(summary(fit))
              coefs <- coefs %>%
                as_tibble() %>%
                mutate(predictor = orf,
                       var = rownames(coefs),
                       aic = summary(fit)$aic,
                       bic = BIC(fit)) %>%
                filter(var!="(Intercept)") %>%
                rename(estimate = Estimate,
                       se = `Std. Error`,
                       tval = `t value`,
                       pval = `Pr(>|t|)`)
}


###################################
# using the iClusterPlus package
###################################
library(iClusterPlus)
library(GenomicRanges)
library(gplots)

# load bundled data, which includes expression, mutation and copy number
data(variation.hg18.v10.nov.2010)
data(gbm)

# calculate mutation rates then filter the boolean mutation matrix
mut.rate <- apply(gbm.mut,2,mean)
gbm.mut <- gbm.mut[,which(mut.rate > 0.02)]

# remove redundant copy number regions
gbm.cn <- CNregions(seg=gbm.seg
                    , epsilon=0
                    , adaptive=FALSE
                    , rmCNV=TRUE
                    , cnv=variation.hg18.v10.nov.2010[,3:5]
                    , frac.overlap=0.5
                    , rmSmallseg=TRUE
                    , nProbes=5)

# sort the datasets and make sure all samples are in same order
gbm.cn <- gbm.cn[order(rownames(gbm.cn)),]
all(rownames(gbm.cn) == rownames(gbm.exp))
all(rownames(gbm.cn) == rownames(gbm.mut))

# run the integrated clustering with K=2
set.seed(123)
fit.single <- iClusterPlus(dt1=gbm.mut
                           , dt2=gbm.cn
                           , dt3=gbm.exp
                           , type=c("binomial","gaussian","gaussian")
                           , lambda=c(0.04,0.61,0.90)
                           , K=2
                           , maxiter=10)

# more realistically, we would choose K using tune.iclusterPlus
# in this case K=1 to K=3
date()
for(k in 4:6){
    cv2.fit <- tune.iClusterPlus(cpus=8
                                 , dt1=gbm.mut
                                 , dt2=gbm.cn
                                 , dt3=gbm.exp
                                 , type=c("binomial","gaussian","gaussian")
                                 , K=k
                                 , n.lambda=101
                                 , scale.lambda=c(0.05,1,1)
                                 , maxiter=20)
      save(cv2.fit, file=paste("cv2.fit.k",k,".Rdata",sep=""))
}
date()

output2 <- list()
files <- grep("cv2.fit",dir())

for(i in 1:length(files)){
  load(dir()[files[i]])
  output2[[i]]=cv2.fit
}

nLambda <- nrow(output2[[1]]$lambda)
nK <- length(output2)
BIC <- getBIC(output2)
devR <- getDevR(output2)
minBICid <- apply(BIC,2,which.min)
devRatMinBIC <- rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] <- devR[minBICid[i],i]
}

# plot %variance explained vs. K to select K using elbow
plot(1:(nK+1),c(0,devRatMinBIC)
     , type="b"
     , xlab="Number of clusters (K+1)"
     , ylab="%Explained Variation")

# choose clusters from the best K
clusters <- getClusters(output2)
rownames(clusters) <- rownames(gbm.exp)
colnames(clusters) <- paste("K=",2:(length(output2)+1),sep="")
#write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
k=2
best.cluster <- clusters[,k]
best.fit <- output2[[k]]$fit[[which.min(BIC[,k])]]

# extract the features
features <- alist()
features[[1]] <- colnames(gbm.mut)
features[[2]] <- colnames(gbm.cn)
features[[3]] <- colnames(gbm.exp)
sigfeatures <- alist()
for(i in 1:3){
  rowsum <- apply(abs(best.fit$beta[[i]]),1, sum)
  upper <- quantile(rowsum,prob=0.75)
  sigfeatures[[i]] <- (features[[i]])[which(rowsum > upper)] 
}
names(sigfeatures) <- c("mutation","copy number","expression")

# plot a heatmap
col.scheme = alist()
col.scheme[[1]] <- bw.col
col.scheme[[2]] <- bluered(256)
col.scheme[[3]] <- bluered(256)

chr <- unlist(strsplit(colnames(gbm.cn),"\\."))
chr <- chr[seq(1,length(chr),by=2)]
chr <- gsub("chr","",chr)
chr <- as.numeric(chr)

plotHeatmap(fit=best.fit
            , datasets=list(gbm.mut , gbm.cn, gbm.exp)
            , type=c("binomial","gaussian","gaussian")
            # , col.scheme = col.scheme
            , row.order = c(F,F,T)
            , chr = chr
            , plot.chr = c(F,T,F)
            , sparse=c(T,F,T)
            , cap=c(F,T,F))

###################################
# using mocluster package
###################################
library(mogsa)



#########################
# using omicade4 package
#########################
library(omicade4)
