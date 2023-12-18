library(MASS)
library(tseries)
library(nortest)
library(car)

dat1 <- runif(1000)
dat2 <- runif(1000)
dat3 <- runif(1000)

df <- reshape2::melt(cbind(dat1, dat2, dat3))[,c("Var2","value")]
colnames(df) <- c("group","value")

### Shapiro-Wilk test for normality
shapiro.test(dat1)

### Anderson-Darling test for normality
ad.test(dat1)

### Student's t-test
t.test(dat1, dat2, alternative = "two.sided", var.equal = TRUE)

## Kruskal-Wallis test
kruskal.test(dat1, dat2)

## Kendall's Tau
cor.test(dat1, dat2, method="kendall")

## ANOVA
anovaOut <- aov(value ~ group, data=df)
summary(anovaOut)

mod <- lm(conformity ~ fcategory*partner.status, data=Moore,
          contrasts=list(fcategory=contr.sum, partner.status=contr.sum))
Anova(mod, type="II", test.statistic ="Wald")

### Pearson's Chi-Squared
tbl <- table(survey$Smoke, survey$Exer) 
chisq.test(tbl) 

# correction for 2x2 contingency table, ignored if input > 2x2
tbl2x2 <- tbl[c(1:2),c(1:2)]
chisq.test(tbl2x2, correct=TRUE)

### Mann-Whitney U test
wilcox.test(dat1, dat2, alternative = "two.sided")

### test for proportions
prop.test(cbind(dat1, dat2), alternative="two.sided")

# using counts with multiple testing correction
counts_tbl <- data.frame(totals = c(76413234, 86464210, 79836106) , cases = c(36, 19, 22) )
rownames(counts_tbl) <- c("A", "B", "C")

A_vs_B <- prop.test(counts_tbl$cases[c(1,2)],
          counts_tbl$totals[c(1,2)])

B_vs_C <- prop.test(counts_tbl$cases[c(2,3)],
                    counts_tbl$totals[c(2,3)])

A_vs_C <- prop.test(counts_tbl$cases[c(1,3)],
                    counts_tbl$totals[c(1,3)])

pvals <- c(A_vs_B$p.value, B_vs_C$p.value, A_vs_C$p.value)
names(pvals) <- c("A_vs_B","B_vs_C","A_vs_C")

pvals_adjusted <- p.adjust(pvals, method="fdr")

# using Fisher's exact test
# NOTE: prop.test is generally preferred since p-vals will be very similar,
# and fisher.test is limited by requiring that the product of the rows be < 2^31-1
# however, if any cases < 10, fisher.test can be more accurate

A_vs_B_fisher <- fisher.test(counts_tbl[c(1,2),])


### Augmented Dickey-Fuller test for time-series autoregressiveness
adf.test(dat1, alternative = "explosive")

# Comparing MFs given Mutation Depth, Total Depth and sample group
glm(Mut Depth ~ group, family = quasipoisson(), offset = log(Total Depth))



