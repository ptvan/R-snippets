library(mongolite)
library(dplyr)
library(magrittr)

m <- mongo(collection = "test", db = "test", url = "mongodb://localhost",
           verbose = FALSE, options = ssl_options())

steps <- read.csv("~/working/datasets/iphone_health/stepsData.csv")

# day-level data
steps <- steps %>%
  mutate(startDate = as.Date(startDate)) %>%
  group_by(startDate) %>%
  summarize(stepsWalked = sum(stepsWalked))

iphone_steps <- mongo("steps")
iphone_steps$insert(steps)
iphone_steps$count()

print(iphone_steps$find('{}'))

iphone_steps$find('{"startDate": "2015-12-07"}')

iphone_steps$drop()