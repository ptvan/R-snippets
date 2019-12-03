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