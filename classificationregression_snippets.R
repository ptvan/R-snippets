##########################
# Naive Bayes
##########################
## using the naivebayes package
library(naivebayes)

n <- 100
set.seed(1)
data <- data.frame(class = sample(c("classA", "classB"), n, TRUE)
                   ,bern = sample(LETTERS[1:2], n, TRUE)
                   ,cat  = sample(letters[1:3], n, TRUE)
                   ,logical = sample(c(TRUE,FALSE), n, TRUE)
                   ,norm = rnorm(n),count = rpois(n, lambda = c(5,15)))

train <- data[1:95, ]
test <- data[96:100, -1]
nb <- naive_bayes(class ~ ., train, usepoisson = TRUE)
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
