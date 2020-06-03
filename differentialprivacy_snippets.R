library(diffpriv)

target <- function(X) mean(X)
mech <- DPMechLaplace(target = target)

distr <- function(n) rnorm(n)
mech <- sensitivitySampler(mech, oracle = distr, n = 5, gamma = 0.1)