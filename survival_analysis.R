library(survival)
library(ggfortify)
library(ranger)
library(dplyr)

# credit: https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

# note: score = Karnofsky performance score
data(veteran)
head(veteran)

#### Kaplan-Meier Survival Curve
km <- with(veteran, Surv(time, status))

km_fit <- survfit(Surv(time, status) ~ 1, data=veteran)
summary(km_fit, times = c(1,30,60,90*(1:10)))

# stratifying subjects by age & rerunning the model
vet <- mutate(veteran, AG = ifelse((age < 60), "LT60", "OV60"),
              AG = factor(AG),
              trt = factor(trt,labels=c("standard","test")),
              prior = factor(prior,labels=c("N0","Yes")))
km_AG_fit <- survfit(Surv(time, status) ~ AG, data=vet)

# running the G-rho Flemming-Harrington test on the stratified data
survdiff(Surv(time, status) ~ AG, data=vet)

# calculating concordance
# NOTE: only handles a single covariate
survConcordance(status ~ karno
                , data = vet)

#### Cox Proportional Hazard model
# NOTE: Cox assumes that covariates do not vary with time
# fitting all covars
cox <- coxph(Surv(time, status) ~ trt + celltype + 
               karno + diagtime + age + prior 
               , data = vet)
summary(cox)
cox_fit <- survfit(cox)
autoplot(cox_fit)


#### Aalen's additive regression model
# time-varying covariates are handled
# NOTE: aareg() performs the fit so we don't need to call survfit()
aa_fit <- aareg(Surv(time, status) ~ trt + celltype +
                 karno + diagtime + age + prior 
                , data = vet)
autoplot(aa_fit)


#### Random-Forest survival using ranger on age-stratified data
r_fit <- ranger(Surv(time, status) ~ trt + celltype + 
                  karno + diagtime + age + prior,
                data = vet,
                mtry = 4,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE)

# getting the variable importance
vi <- data.frame(sort(round(r_fit$variable.importance, 4), decreasing = TRUE))
names(vi) <- "importance"

# getting Harrell's c-index, basically 
# n_concordance_pairs / (n_conconrdent_pairs + n_discordant_pairs)
# more details: https://statisticaloddsandends.wordpress.com/2019/10/26/what-is-harrells-c-index/
cat("Prediction Error = 1 - Harrell's c-index = ", r_fit$prediction.error)

# Average the survival models
death_times <- r_fit$unique.death.times 
surv_prob <- data.frame(r_fit$survival)
avg_prob <- sapply(surv_prob,mean)