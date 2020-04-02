library(mlbench) # provides data generation funcsions

## load 2D Swiss rolls from mlbench pkg, discard class info
truth <- mlbench.spirals(500, 1, 0.025)
data <- truth$x

###########
#  K-means
###########
# needs K
km <- kmeans(data, centers = 2)


#############
# Mean Shift
############
library(meanShiftR)
# meanshiftR supports multicore acceleration and KDE (second arg)
# faster than LPCM and older implementation `MeanShift``
ms <- meanShift(data, data,
          algorithm="KDTREE", 
          nNeighbor=8, 
          parameters=c(5,7.1) )


###########
#  DBSCAN
###########
# unsupervised
library(dbscan)
db <- dbscan(data, eps=5)

######################
# Spectral Clustering
######################
## using the kernlab package
# NOTE: with higher stddevs, this doesn't work so well (eg. 0.05)
library(kernlab)
sc <- specc(data, centers=2)
plot(data, col=sc, pch=4)            # estimated classes (x)
points(data, col=truth$classes, pch=5)


##########################
# Gaussian Mixture Models
##########################
library(ClusterR) 
gauss2D <- mlbench.2dnormals(n=10000, cl=2, r=3)
g <- GMM(gauss2D, 2, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 10,
    em_iter = 10, verbose = F)   

plot(gauss2D)
points(g$centroids,col="red", pch=16)

pr <- predict_GMM(gauss2D, g$centroids, g$covariance_matrices, g$weights)    
cols <- gsub("1","green",(gsub("0","blue",pr$cluster_labels)))
plot(gauss2D)
points(gauss2D,col=cols)


#######################################
# Agglomerative hierarchical clustering
#######################################
library(cluster) # provides agnes()
ag <- agnes(Ionosphere[,c(3:34)])
cl <- cutree(ag, 2)


############################## 
# Divisive analysis clustering
############################## 
di <- agnes(Ionosphere[,c(3:34)])
