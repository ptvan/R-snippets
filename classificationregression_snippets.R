##########################
# Naive Bayes
##########################
## using the naivebayes package
library(tidyverse)
library(ggplot2)

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

train <- sample(1:150, 75)

plda = predict(object = r, # predictions
               newdata = iris[-train, ])
head(plda$class)
ldahist(plda$x[,2], g = plda$class)

# ggplot is a little nicer...
irisProjection <- cbind(scale(as.matrix(iris[,-5]),scale=FALSE) %*% r$scaling,iris[,5,drop=FALSE])
ggplot(data=irisProjection,aes(x=LD1,y=LD2,col=Species)) +
 geom_point()  

##########################
# XGBOOST
##########################
library(xgboost)

data(agaricus.train, package='xgboost')
data(agaricus.test, package='xgboost')
train <- agaricus.train
test <- agaricus.test

# random forests
bst_rf <- xgboost(data = train$data
                  , label = train$label
                  , max_depth = 4
                  , num_parallel_tree = 1000
                  , subsample = 0.5
                  , colsample_bytree =0.5
                  , nrounds = 1
                  , objective = "binary:logistic")

bst <- xgboost(data = dtrain
               , booster = "gbtree"
               , max_depth = 2
               , eta = 1
               , nthread = 2
               , nrounds = 2
               , objective = "binary:logistic", verbose = 2)

# XGB's default is regression, so the predictions are probabilities
pred <- predict(bst, test$data)

prediction <- as.numeric(pred > 0.5)
err <- mean(prediction != test$label)
print(paste("test-error=", err))

xgb.importance(model = bst)
xgb.plot.shap(agaricus.test$data, model = bst, features = "odor=none")


#########################
# Support Vector Machines
#########################
library(e1071)

gauss2D <- mlbench.2dnormals(n=10000, cl=2, r=1)

train <- list()
train$x <- gauss2D$x[1:8000,]
train$labels <- gauss2D$classes[1:8000]

test <- list()
test$x <- gauss2D$x[8001:10000,]
test$labels <- gauss2D$classes[8001:10000]


classifier <- svm(train$x, train$labels, type="C-classification", kernel='linear')
y_pred <- predict(classifier, newdata = test$x) 

plot(test$x, col=test$labels)
cf <- coef(classifier)
abline(-cf[1]/cf[3], -cf[2]/cf[3], col = "red")