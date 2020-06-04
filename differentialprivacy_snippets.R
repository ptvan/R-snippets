library(diffpriv)

target <- function(X) mean(X)
mech <- DPMechLaplace(target = target)

distr <- function(n) rnorm(n)
mech <- sensitivitySampler(mech, oracle = distr, n = 5, gamma = 0.1)

# length is sensitivitySampler() n
X <- c(4.11,-0.1,-0.14, 1.516, 2.43) 
r <- releaseResponse(mech, privacyParams = DPParamsEps(epsilon = 1), X = X)