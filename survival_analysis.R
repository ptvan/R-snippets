library(survival)
library(ggfortify)
library(ranger)
library(dplyr)

data(veteran)

head(veteran)
# note: score = Karnofsky performance score

## Kaplan-Meier Survival Curve
km <- with(veteran, Surv(time, status))
