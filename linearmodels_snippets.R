library(glmm)
library(lme4)

data("Dyestuff")
fm1 <- lmer(Yield ~ 1 + (1 | Batch), Dyestuff)