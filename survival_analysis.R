library(survival)
library(ggfortify)
library(ranger)
library(dplyr)

# note: score = Karnofsky performance score
data(veteran)
head(veteran)

#### Kaplan-Meier Survival Curve
km <- with(veteran, Surv(time, status))

km_fit <- survfit(Surv(time, status) ~ 1, data=veteran)
summary(km_fit, times = c(1,30,60,90*(1:10)))

# splitting subjects by age & rerunning the model
vet <- mutate(veteran, AG = ifelse((age < 60), "LT60", "OV60"),
              AG = factor(AG),
              trt = factor(trt,labels=c("standard","test")),
              prior = factor(prior,labels=c("N0","Yes")))
km_AG_fit <- survfit(Surv(time, status) ~ AG, data=vet)

#### Cox Proportional Hazard model
# fitting all covars
# NOTE: Cox assumes that covariates do not vary with time
cox <- coxph(Surv(time, status) ~ trt + celltype + 
               karno + diagtime + age + prior 
               , data = vet)
summary(cox)
cox_fit <- survfit(cox)
autoplot(cox_fit)


#### Random-Forest survival using ranger on age-split data
r_fit <- ranger(Surv(time, status) ~ trt + celltype + 
                  karno + diagtime + age + prior,
                data = vet,
                mtry = 4,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE)

# Average the survival models
death_times <- r_fit$unique.death.times 
surv_prob <- data.frame(r_fit$survival)
avg_prob <- sapply(surv_prob,mean)