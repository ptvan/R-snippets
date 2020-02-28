##########################
# Naive Bayes
##########################
## using the naivebayes package
library(tidyverse)
library(naivebayes)

set.seed(1)

# abalone data from https://archive.ics.uci.edu/ml/datasets/Abalone
# NOTE: length/diameter are likely correlated, while naive Bayes assumes
# covars to be independent !

data <- read.csv("abalone.csv")
colnames(data) <- c("sex"
                  ,"length"
                  ,"diameter"
                  ,"height"
                  ,"weight_whole"
                  ,"weight_shucked"
                  ,"weight_viscera"
                  ,"weight_shell"
                  ,"n_rings")

n <- nrow(data)

train <- data[1:floor(n*0.9), ]
test <- data[(n-floor(n*0.1)):n,-1] # remember to exclude the labels !
nb <- naive_bayes(sex ~ ., train, usepoisson = TRUE)
summary(nb)
# classification
predict(nb, test, type = "class")
# probabilities
predict(nb, test, type = "prob")

## using the caret package
library(caret)

#####################################
# Linear Discriminant Analysis (LDA)
####################################
library(MASS)
r <- lda(formula = Species ~ .,
         data = iris,
         prior = c(1,1,1)/3)

r2 <- lda(formula = Species ~ .,
          data = iris,
          prior = c(1,1,1)/3,
          CV = TRUE)


train <- sample(1:150, 75)

plda = predict(object = r, # predictions
               newdata = iris[-train, ])
head(plda$class)
