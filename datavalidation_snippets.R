library(validate)
data(women)
cf <- check_that(women, height > 0, weight > 0, height/weight > 0.5)
summary(cf)

rules <- validator(height == women_reference$height)
cf <- confront(women, rules, ref = list(women_reference = women1))
summary(cf)

rules <- validator( fruit %in% codelist )
fruits <-  c("apple", "banana", "orange")
dat <- data.frame(fruit = c("apple","broccoli","orange","banana"))
cf <- confront(dat, rules, ref = list(codelist = fruits))
summary(cf)