library(glmm)
library(lme4)
library(glmmADMB)

data(Dyestuff)
# using restricted maximum likelihood
fm1 <- lmer(Yield ~ 1 + (1 | Batch), data = Dyestuff, REML=FALSE)
summary(fm1)

# using log-likelihood
fm2 <- lmer(Yield ~ 1 + (1 | Batch), data = Dyestuff, REML=TRUE)
summary(fm2)

# using MCML
glmm(Yield ~ Batch, data = Dyestuff, varcomps.names = "Yield", family.glmm = bernoulli.glmm)

# handling zero-inflated data
data(Owls)
Owls <- transform(Owls
                  ,Nest=reorder(Nest,NegPerChick)
                  ,logBroodSize=log(BroodSize)
                  ,NCalls=SiblingNegotiation)

zipoiss <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+offset(logBroodSize)+(1|Nest)
                    ,data=Owls
                    ,zeroInflation=TRUE
                    ,family="poisson")

hurdle <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+BroodSize+(1|Nest)
                         ,data=subset(Owls,NCalls>0)
                         ,family="truncnbinom1")
