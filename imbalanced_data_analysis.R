library(smotefamily)
library(ROSE)

data(hacide)

dat_plot = SMOTE(hacide.train[,2:3],
                 as.numeric(hacide.train$cls),
                 K = 3, dup_size = 0)