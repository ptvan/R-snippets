library(tidyverse)
library(magrittr)
library(ggplot2)

# reference: Dimensionality Reduction: A Comparative Review, van der Maaten et al 2009 
# https://lvdmaaten.github.io/publications/

# HTRU2 from the UCI ML repository
# https://archive.ics.uci.edu/ml/datasets/HTRU2
dat <- read.csv("~/working/datasets/HTRU.csv")
colnames(dat) <- c(paste0("var", 1:(ncol(dat) - 1)), "class")
dat <- dat[sample(nrow(dat),5000),]

##########################
# classical (metric) MDS
##########################
# needs a distance structure
# k is number of dimensions
d <- dist(dat)
mds <- cmdscale(d, eig = TRUE, k = 2)
plot(mds$points[,1], mds$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2",
     main = "metric MDS", type = "n")
text(mds$points[,1], mds$points[,2], labels = row.names(dat), cex = 0.7) 

################# 
# non-parametric MDS
################# 
library(MASS)
nmds <- isoMDS(d, k = 2)
plot(nmds$points[,1], nmds$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2",
     main = "non-metric MDS", type = "n")
text(nmds$points[,1], nmds$points[,2], labels = row.names(dat), cex = 0.7) 


###########
# T-SNE
###########
library(Rtsne)
tsne_out <- Rtsne(dat[,c(1:8)], perplexity = 30,theta = 0.0) 

tsne_out <- tsne_out$Y %>%
          as.data.frame() %>% 
          set_colnames(c("X","Y")) %>%
          inset("class", value = dat$class) %>%
          mutate(class = as.factor(class))

#######
# UMAP
#######
library(uwot)
umap_out <- umap(dat, init = "spca") %>%
          as.data.frame() %>% 
          set_colnames(c("X","Y")) %>%
          inset("class", value = dat$class) %>%
          mutate(class = as.factor(class))

ggplot(umap_out) +
      aes(x = X, y = Y, col = class) +
      geom_point()

############
# KERNEL PCA
############
library(kernlab) 
kpca_out <- kpca(dat, features = 2) 
rotated(kpca_out)


############
# SAMMON
###########
library(MASS)
sammon_out <- sammon(dist(dat[,-180]))
plot(sammon_out$points, type = "n")
text(sammon_out$points, labels = as.character(1:nrow(dat[,-180])))


################
# AUTO-ENCODERS
################
library(keras) 
x_train <- as.matrix(dat[,c(1:7)])

# build neuralnet, 2 bottleneck units ~ 2D
model <- keras_model_sequential()
model %>%
  layer_dense(units = 4, activation = "tanh", input_shape = ncol(x_train)) %>%
  layer_dense(units = 2, activation = "tanh", name = "bottleneck") %>%
  layer_dense(units = 4, activation = "tanh") %>%
  layer_dense(units = ncol(x_train))

model %>% compile(
  loss = "mean_squared_error", 
  optimizer = "adam"
)

summary(model)

# train
model %>% fit(
  x = x_train, 
  y = x_train, 
  epochs = 2000,
  verbose = 1
)

# evaluate neuralnet
mse.ae2 <- evaluate(model, x_train, x_train)

# extract bottleneck layer
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, x_train)

# plot principal components
ggplot(data.frame(PC1 = intermediate_output[,1], PC2 = intermediate_output[,2])
       , aes(x = PC1, y = PC2, col = dat$class)) + 
        geom_point()


#################
# DIFFUSION MAPS
#################
library(diffusionMap)

# works on a `dist` object
d <- dist(dat[,c(1:8)])
dmap <- diffuse(D, eps.val = 0.1) 
plot(dmap)


####################### 
# Self-organizing Maps
####################### 
library(kohonen)
set.seed(100)
grd <- somgrid(xdim = 10, ydim = 10, topo = "hexagonal")
sommodel <- som(data, grd)
plot(sommodel, type = "mapping", pchs = 19)
plot(sommodel, type = "codes", pchs = 19)
plot(sommodel, type = "changes")
plot(sommodel, type = "counts")

################################### 
# Nonnegative Matrix Factorization
################################### 
library(NMF)

# example ExpressionSet with Sample phenotypic column removed
data(esGolub)
esGolub
esGolub <- esGolub[1:200,]
esGolub$Sample <- NULL

# run NMF with rank 3
nmf_out <- nmf(esGolub, 3)

# fit the output
nmf_fitted <- fitted(nmf_out)

# get some quality metrics
summary(nmf_fitted, target = esGolub)

#####################################
# Linear Discriminant Analysis (LDA)
#####################################
# see the corresponding section in classificationregression_snippets.R

