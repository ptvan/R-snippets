library(diffpriv)
library(tidycensus)
library(tidyverse)

target <- function(X) mean(X)
mech <- DPMechLaplace(target = target)

distr <- function(n) rnorm(n)
mech <- sensitivitySampler(mech, oracle = distr, n = 5, gamma = 0.1)

# length is sensitivitySampler() n
X <- c(4.11,-0.1,-0.14, 1.516, 2.43) 
r <- releaseResponse(mech, privacyParams = DPParamsEps(epsilon = 1), X = X)

# US Census 2020 uses DP
census_api_key("MY_API_KEY")

age20 <- get_decennial(geography = "state", 
                       variables = "P13_001N", 
                       year = 2020,
                       sumfile = "dhc")