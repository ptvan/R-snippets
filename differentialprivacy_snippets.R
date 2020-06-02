library(diffpriv)
target <- function(X) mean(X)
mech <- DPMechLaplace(target = target)