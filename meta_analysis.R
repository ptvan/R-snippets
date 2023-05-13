library(metafor)

## loading example data
data("dat.bcg")
print(dat.bcg)

## run effect-size calculation
dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg
              , data = dat.bcg
              , append = TRUE)

## fitting a random-effects model
res <- rma(ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat,  measure = "RR")
confint(res)

## make a forest plot
forest(res, slab = paste(dat$author, dat$year, sep = ", ")
       , xlim = c(-16, 6)
       , at = log(c(0.05, 0.25, 1, 4))
       , atransf = exp
       , ilab = cbind(dat$tpos, dat$tneg, dat$cpos, dat$cneg)
       , ilab.xpos = c(-9.5, -8, -6, -4.5), cex = 0.75)

op <- par(cex = 0.75, font = 2)
text(c(-9.5, -8, -6, -4.5), 15, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75, -5.25), 16, c("Vaccinated", "Control"))
text(-16, 15, "Author(s) and Year", pos = 4)
text(6, 15, "Relative Risk [95% CI]", pos = 2)

par(op)

## fitting a mixed-effects model using year and latitude of study's location
res <- rma(yi, vi, mods = cbind(ablat, year), data = dat)

## make a QQ plot
qqnorm(res)

## make a funnel plot
funnel(res, refline = 0)

## predict new values
predict(res, newmods = cbind(seq(from = 10, to = 60, by = 10), 1970)
        , transf = exp, addx = TRUE)

