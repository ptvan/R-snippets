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

### Kruskal-Wallis test
kruskal.test(dat1, dat2)

### Kendall's Tau
cor.test(dat1, dat2, method="kendall")

### ANOVA
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

### Augmented Dickey-Fuller test for time-series autoregressiveness
adf.test(dat1, alternative = "explosive")
