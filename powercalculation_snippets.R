library(pwr)

# t-tests
pwr.t2n.test(n1 = 30, n2= 30, d = 0.2, sig.level =0.05)

# correlation
pwr.r.test(n = 30, r = 0.5, sig.level = 0.05) 

# ANOVA
pwr.anova.test(k = 3, n = 30, f = 0.1, sig.level = 0.05)

# Proportions
pwr.2p.test(h = 0.2, n = 30, sig.level = 0.05) 
