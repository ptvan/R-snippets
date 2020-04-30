##########################
# Naive Bayes
##########################
## using the naivebayes package
library(naivebayes)
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

train <- data.frame(cbind(gauss2D$x[1:8000,], gauss2D$classes[1:8000]))
colnames(train) <- c("X1","X2", "y")

test <- data.frame(cbind(gauss2D$x[8001:10000,], gauss2D$classes[8001:10000]))
colnames(test) <- c("X1","X2", "y")

### using linear kernel
cl_lin <- svm(y~., data=train,  type="C-classification", kernel='linear', scale=TRUE)
pred_lin <- predict(cl_lin, newdata = test) 

plot(cl_lin, train)
plot(cl_lin, test)

### using radial kernel
cl_rbf <- svm(y~., data=train, type="C-classification", kernel='radial')
pred_rbf <- predict(cl_rbf, newdata = test)

# confusion matrix
table(predicted = pred_rbf, actual = test$y)

# misclassification rate
mean(pred_rbf != test$y)*100

plot(cl_rbf, train)


#################### 
# Tree-based methods
#################### 

### classifying BreastCancer data using rpart
library(mlbench)
library(rpart)
data(BreastCancer)

BreastCancer$Cell.size <- as.numeric(BreastCancer$Cell.size)
BreastCancer$Cl.thickness <- as.numeric(BreastCancer$Cl.thickness)
BreastCancer$Marg.adhesion <- as.numeric(BreastCancer$Marg.adhesion)
BreastCancer$Cell.shape <- as.numeric(BreastCancer$Cell.shape)
BreastCancer$Mitoses <- as.numeric(BreastCancer$Mitoses)
BreastCancer$Epith.c.size <- as.numeric(BreastCancer$Epith.c.size)
BreastCancer$Bl.cromatin <- as.numeric(BreastCancer$Bl.cromatin)
BreastCancer$Normal.nucleoli <- as.numeric(BreastCancer$Normal.nucleoli)
BreastCancer$Bare.nuclei <- as.numeric(BreastCancer$Bare.nuclei)

# randomForest needs complete cases
BreastCancer <- BreastCancer[which(complete.cases(BreastCancer)),]

fit <- rpart(Class ~ Cell.size + Cl.thickness + Marg.adhesion + Cell.shape + Mitoses + Epith.c.size + Bl.cromatin + Normal.nucleoli + Bare.nuclei, method="class", data=BreastCancer)
printcp(fit)
plotcp(fit)
plot(fit, uniform=TRUE, main="breast cancer classification")
text(fit, use.n=TRUE, all=TRUE, cex=.8)

# prune the forest
pfit<- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
plot(pfit, uniform=TRUE,
     main="breast cancer classification, pruned")
text(pfit, use.n=TRUE, all=TRUE, cex=.8)

### classifying BreastCancer data using random forest
library(randomForest)
rffit <- randomForest(Class ~ Cell.size + Cl.thickness + Marg.adhesion + Cell.shape + Mitoses + Epith.c.size + Bl.cromatin + Normal.nucleoli + Bare.nuclei, data=BreastCancer)
importance(rffit)



#########################
# Partial Least Squares
#########################
library(mlbench)
library(pls)

### regression so response must be continuous, corrected median value from BostonHousing2 in this case
# when ncomp is not specified, the maximum number is used
plsfit <- plsr(cmedv ~ lon + lat + crim + indus + nox + rm + age + dis + ptratio , data=BostonHousing2, validation="LOO")

# plot Root Mean Squared Error of Prediction vs. no. of components
plot(RMSEP(plsfit), legendpos="topright")

# plot cross-validated prediction vs. observed values
plot(plsfit, ncomp = 2, asp = 1, line = TRUE)

# plot loadings
plot(plsfit, "loadings", comps=1:2, legendpos="topright")

# extract explained variances
explvar(plsfit)