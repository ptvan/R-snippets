#########################
# using the MOFA package
#########################
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


