##########################
# Spectral Clustering
##########################
library(mlbench)

## load 2D Swiss rolls from mlbench pkg, discard class info
truth <- mlbench.spirals(500, 1, 0.025)
data <- truth$x

## using the kernlab package
# NOTE: with higher stddevs, this doesn't work so well (eg. 0.05  )
library(kernlab)
sc <- specc(data, centers=2)
plot(data, col=sc, pch=4)            # estimated classes (x)
points(data, col=truth$classes, pch=5)

