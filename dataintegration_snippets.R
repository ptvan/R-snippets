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
