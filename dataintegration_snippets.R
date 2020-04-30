#########################
# using the MOFA package
#########################
# NOTE: need to update below to MOFA v2
library(MOFA)

# simulate some data & create MOFA object
N <- 100
D1 <- 250
D2 <- 500

view1 = matrix(rnorm(N*D1),nrow=D1, ncol=N)
view2 = matrix(rnorm(N*D2),nrow=D2, ncol=N)
colnames(view1) <- colnames(view2) <- paste0("sample_",1:N)
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


