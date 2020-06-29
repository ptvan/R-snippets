library(pwr)

# t-tests
pwr.t2n.test(n1 = 30, n2= 30, d = 0.2, sig.level =0.05)

# correlation
pwr.r.test(n = 30, r = 0.5, sig.level = 0.05) 
