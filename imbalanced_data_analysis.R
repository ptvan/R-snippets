library(smotefamily)
library(ROSE)
library(dplyr)

# data from the ROSE package
data(hacide)

# data imbalance: minority outcome only occurs in 2% of cases
hacide.train %>% group_by(cls) %>% count()

# SMOTE
smote_out <- SMOTE(hacide.train[,2:3],
             as.numeric(hacide.train$cls),
             K = 3, dup_size = 0)

# ADASYN
adasyn_out <- ADAS(hacide.train[,2:3],
      as.numeric(hacide.train$cls),
      K = 3)