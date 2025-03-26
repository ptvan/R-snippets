library(validate)
library(simputation)
library(lubridate)
library(admiral)


####################
# Data validation
####################
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

## impute missing data
# read in steps data
steps <- read.csv("~/working/datasets/iphone_health/stepsData.csv") %>%
                    select(startDate, stepsWalked) %>%
                    mutate(startDate = as.Date(startDate)) %>%
                    group_by(startDate) %>%
                    summarize(dailySteps = sum(stepsWalked)) %>%
                    mutate(dayOfWeek = factor(weekdays(startDate),
                                              levels = c("Monday",
                                                         "Tuesday",
                                                         "Wednesday",
                                                         "Thursday",
                                                         "Friday",
                                                         "Saturday",
                                                         "Sunday"
                                                         )
                                              )) %>%
                    mutate(weekNumber = as.numeric(strftime(startDate, "%V")))

# insert new date with no step
steps_missing <- steps_df %>% 
            add_row(., 
                    startDate=as.Date("2019-10-10"), 
                    dailySteps=NA, 
                    dayOfWeek=factor("Thursday"), 
                    weekNumber=41)
tail(steps_missing)

# inpute by median weekNumber 
steps_imputed_weekNumber <- impute_median(steps_missing, dailySteps ~ weekNumber)
tail(steps_imputed_weekNumber)


