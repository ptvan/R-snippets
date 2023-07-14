library(boot)
library(nptest)
library(broom)

set.seed(1)

# bootstrap a point statistic (trimmed mean in this case)
x <- rnorm(10000)
trimmedmean <- function(x, d, trim=0) {
  return(mean(x[d], trim/length(x)))
}

mean_reps <-  boot(x, trimmedmean, R = 1000, trim = 5)
mean_reps$t0

# bootstrap a linear model
rsq_function <- function(formula, data, indices) {
  d <- data[indices,] 
  fit <- lm(formula, data = d)
  return(summary(fit)$r.square)
}

lm_reps <- boot(mtcars, 
     R = 200, 
     statistic = rsq_function,
     formula = mpg~disp
     )

lm_reps$t0
