library(metafor)

# loading example data
data("dat.bcg")
print(dat.bcg)

# run effect-size calculation
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg
              ,data=dat.bcg
              ,append=TRUE)

# escalc also takes formula input, which requires reformatting
k <- length(dat.bcg$trial)
dat.fm <- data.frame(study = factor(rep(1:k, each = 4)))
dat.fm$grp <- factor(rep(c("T", "T", "C", "C"), k), levels = c("T", "C"))
dat.fm$out <- factor(rep(c("+", "-", "+", "-"), k), levels = c("+", "-"))
dat.fm$freq <- with(dat.bcg, c(rbind(tpos, tneg, cpos,cneg)))
# NOTE: this throws an error with metafor 2.4 ?!
escalc(out ~ grp | study, weights = freq, data = dat.fm, measure = "RR")
