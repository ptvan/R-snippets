library(mclust)
library(lubridate)
library(dplyr)
library(mvtnorm)
library(mixtools)

# Thanks to Chad Young for useful mclust discussions

#############
# UNIVARIATE 
############
# read in our iPhone daily step counts
raw <- read.csv("datasets/iphone_health/stepsData.csv")

# day-level data
steps <- raw %>%
  mutate(startDate = factor(as.Date(startDate))) %>%
  dplyr::select(date = startDate, stepsWalked) %>%
  group_by(date) %>%
  summarize(stepsWalked = sum(stepsWalked))

## mixtools
# no mu supplied, data is binned, then K is picked
mixauto <- normalmixEM(steps$stepsWalked)
plot(mixauto, density=TRUE)

# K is supplied
mix4 <- normalmixEM(steps$stepsWalked, k=4)
plot(mix4, density=TRUE)

## mclust
# do density estimation, let mclust pick our distribution count for us
dens <- densityMclust(steps$stepsWalked)

# look at dianogstics
plot(dens, what = "BIC")
summary(dens$BIC)
summary(dens, parameters = TRUE)
plot(dens, what = "density", data=steps$stepsWalked)

# extract parameters of the Gaussians
params <- dens$parameters
nDistros <- length(params$pro)

distros <- data.frame(matrix(0,4,4))
colnames(distros) <- c("id","n", "mean","sd")
distros$id <- 1:4

# re-sample the Gaussians
N <- nrow(steps)
for (i in 1:nDistros){
  distros[distros$id==i,]$n <- floor(N * params$pro[i])
  distros[distros$id==i,]$mean <- params$mean[i]
  distros[distros$id==i,]$sd <- floor(sqrt(params$variance$sigmasq[i]))
}
set.seed(1)
densities <- apply(distros[, -1], 1, function(x) rnorm(n = x[1], mean = x[2], sd = x[2]))

# plot out the estimated densities
plot(0, type = "n", xlim = c(min(unlist(lapply(densities, min))),max(unlist(lapply(densities, max)))), ylim = c(0, 1.5e-4))
for (d in 1:length(densities)){
  lines(density(densities[[d]]), lty = d)
}
legend("topright", legend=1:length(densities), lty=1:length(densities))


##############
# MULTIVARIATE
##############

# probably not a great idea since running and biking are temporally independent...
biking <- read.csv("cyclingData.csv", header=T)
biking <- biking %>%
  group_by(startDate) %>%
  summarize(kcalRun = sum(kcalBurned),milesCycled = sum(milesCycled)) %>%
  mutate(startDate = as.Date(startDate))  %>%
  dplyr::filter(milesCycled < 20) %>%
  dplyr::select(c(startDate,milesCycled))

health <- merge(biking, steps, by="startDate", all=TRUE)
health[is.na(health$milesCycled),]$milesCycled <- 0
health <- dplyr::select(health, c(milesCycled, stepsWalked))

mdens <- densityMclust(health, classification=TRUE)

# retrieve model summary: model type, log-likelihood, BIC, etc.
summary(mdens)

# cluster membership
mdens$classification

# plot the multivariate density
plot(mdens, what="density", type="persp")

# retrieve the covariance matrix
covMat <- cov(health)

# resample
rmvnorm(n = nrow(health), mean=mdens$parameters$mean, sigma=mdens$parameters$variance$sigma[,,8])
