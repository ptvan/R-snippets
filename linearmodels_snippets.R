library(glmm)
library(lme4)

data(Dyestuff)
# using restricted maximum likelihood
fm1 <- lmer(Yield ~ 1 + (1 | Batch), data = Dyestuff, REML=FALSE)
summary(fm1)

# using log-likelihood
fm2 <- lmer(Yield ~ 1 + (1 | Batch), data = Dyestuff, REML=TRUE)
summary(fm2)