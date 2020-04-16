library(survival)
library(ggfortify)
library(ranger)
library(dplyr)

data(veteran)

head(veteran)
# note: score = Karnofsky performance score

## Kaplan-Meier Survival Curve
km <- with(veteran, Surv(time, status))

km_fit <- survfit(Surv(time, status) ~ 1, data=veteran)
summary(km_fit, times = c(1,30,60,90*(1:10)))