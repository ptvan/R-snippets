library(metafor)

## loading example data
data("dat.bcg")
print(dat.bcg)

## run effect-size calculation
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg
              ,data=dat.bcg
              ,append=TRUE)

## escalc also takes formula input, which requires reformatting
k <- length(dat.bcg$trial)
dat.fm <- data.frame(study = factor(rep(1:k, each = 4)))
dat.fm$grp <- factor(rep(c("T", "T", "C", "C"), k), levels = c("T", "C"))
dat.fm$out <- factor(rep(c("+", "-", "+", "-"), k), levels = c("+", "-"))
dat.fm$freq <- with(dat.bcg, c(rbind(tpos, tneg, cpos,cneg)))
# NOTE: this throws an error with metafor 2.4 ?!
escalc(out ~ grp | study, weights = freq, data = dat.fm, measure = "RR")

## fitting a random-effects model
res <- rma(ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat,  measure = "RR")
confint(res)

## plot output using the package's forest()
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

predict(res, newmods = cbind(seq(from = 10, to = 60, by = 10), 1970)
        , transf = exp, addx = TRUE)