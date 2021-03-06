library(glmm)
library(lme4)
library(glmmADMB)
library(brms)
library(MASS) # for Box-Cox transformation
library(car) # for VIF, Durbin-Watson test

# Box-Cox transformation
boxcox(lm(dist~speed,data=cars),lambda=seq(0,1,by=.1))

############################
# MODEL DIAGNOSTICS
##########################
# fitted vs. residuals plot
lmod <- lm(mpg ~ disp + hp + wt + drat, data=mtcars)
plot(fitted(lmod), residuals(lmod), xlab="Fitted", ylab="Residuals")
abline(h=0)
qqPlot(lmod)

# evaluate multi-collinearity by looking at variance inflation factors
vif(lmod)
sqrt(vif(lmod)) > 2

# evaluate independence of errors
durbinWatsonTest(lmod)

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

# hurdle model
hurdle <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+BroodSize+(1|Nest)
                         ,data=subset(Owls,NCalls>0)
                         ,family="truncnbinom1")

# using BRMS
bf1 <- brm(bf(NCalls ~ ArrivalTime), 
            data = Owls, family = gaussian())

