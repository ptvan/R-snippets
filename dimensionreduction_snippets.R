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
# non-metric MDS
################# 
library(MASS)
nmds <- isoMDS(d, k = 2)
plot(nmds$points[,1], nmds$points[,2], xlab = "Coordinate 1", ylab = "Coordinate 2",
     main = "non-metric MDS", type = "n")
text(nmds$points[,1], nmds$points[,2], labels = row.names(dat), cex = .7) 


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
dmap <- diffuse(D, eps.val = .1) 
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

####################### 
# Linear Discriminate Analysis (LDA)
####################### 
library(klaR)
library(ggord)

indices <- sample(2, nrow(iris),
              replace = TRUE,
              prob = c(0.6, 0.4))

training <- iris[indices == 1,]
testing <- iris[indices == 2,]

# train
lda_train <- lda(Species~., training)

# test
lda_prediction <- predict(lda_train, training)

# histogram of LDA1, good separation
ldahist(data = lda_prediction$x[,1], g = training$Species)

# histogram of LDA2, separation is not as good
ldahist(data = lda_prediction$x[,2], g = training$Species)

# additional diagnostic plots
ggord(lda_train, training$Species, ylim = c(-10, 10))
partimat(Species~., data = training, method = "lda")

# check predicted vs. observed counts
p1 <- predict(lda_train, training)$class
training_table <- table(Predicted = p1, Actual = training$Species)

p2 <- predict(lda_train, testing)$class
testing_table <- table(Predicted = p2, Actual = testing$Species)

# calculate accuracy
sum(diag(training_table))/sum(training_table)
