library(mclust)
library(lubridate)
library(tidyr)
library(dplyr)
# Thank you for Chad Young for useful discussions

# read in our iPhone daily step counts
steps <- read.csv("stepsData.csv")

# day-level data
steps <- steps %>%
  mutate(startDate = as.Date(startDate)) %>%
  group_by(startDate) %>%
  summarize(stepsWalked = sum(stepsWalked))

dens <- densityMclust(steps$stepsWalked)
plot(dens, what = "BIC")
summary(dens$BIC)
summary(dens, parameters = TRUE)
plot(dens, what = "density", data=steps$stepsWalked)

params <- dens$parameters
nDistros <- length(params$pro)

distros <- data.frame(matrix(0,4,4))
colnames(distros) <- c("id","n", "mean","sd")
distros$id <- 1:4

N <- nrow(steps)
for (i in 1:nDistros){
  distros[distros$id==i,]$n <- floor(N * params$pro[i])
  distros[distros$id==i,]$mean <- params$mean[i]
  distros[distros$id==i,]$sd <- floor(sqrt(params$variance$sigmasq[i]))
}

densities <- apply(distros[, -1], 1, function(x) rnorm(n = x[1], mean = x[2], sd = x[2]))


plot(0, type = "n", xlim = c(min(unlist(lapply(densities, min))),max(unlist(lapply(densities, max)))), ylim = c(0, 1.5e-4))
for (d in 1:length(densities)){
  lines(density(densities[[d]]), lty = d)
}
legend("topright", legend=1:length(densities), lty=1:length(densities))