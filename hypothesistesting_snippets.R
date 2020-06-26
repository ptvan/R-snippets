library(MASS)

dat1 <- runif(1000)
dat2 <- runif(1000)
dat3 <- runif(1000)

# Shapiro-Wilk test for normality
shapiro.test(dat1)

# Student's t-test
t.test(dat1, dat2)

# ANOVA
df <- reshape2::melt(cbind(dat1, dat2, dat3))[,c("Var2","value")]
colnames(df) <- c("group","value")
anovaOut <- aov(value ~ group, data=df)
summary(anovaOut)

# Pearson's Chi-Squared
tbl <- table(survey$Smoke, survey$Exer) 
chisq.test(tbl) 

# Mann-Whitney U test
wilcox.test(dat1, dat2)
