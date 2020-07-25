library(stats)
library(optimization)
library(GenSA)

## Himmelblau's Function
hi <- function(x,y) { (x^2 + y -11)^2 +
    (x + y^2 - 7)^2 }

## Rosenbrock function
ro <- function(x){
  100*(x[2]-x[1]^2)^2+(1-x[1])^2
}

### Nelder-Mead
nm_out <- stats::optim(fn = ro, par = c(10, 10), method = "Nelder-Mead")

### Simulated Annealing
sa_out <- optim_sa(fun = ro
         , start = (c(-0.60,  -0.030))
         , trace = TRUE
         ,lower = c(-5, -5)
         , upper = c(5, 5)
         , control = list(t0 = 500, nlimit = 50, r = 0.85, rf = 3, ac_acc = 0.1, dyn_rf = TRUE))