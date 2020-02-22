library(mongolite)
library(dplyr)
library(magrittr)

m <- mongo(collection = "test", db = "test", url = "mongodb://localhost",
           verbose = FALSE, options = ssl_options())

stepsData <- read.csv("~/working/datasets/iphone_health/stepsData.csv")

# day-level stepcount data
stepsData <- stepsData %>%
  mutate(startDate = as.Date(startDate)) %>%
  group_by(startDate) %>%
  summarize(stepsWalked = sum(stepsWalked))

# INSERT
stepsDB <- mongo("steps")
stepsDB$insert(stepsData)
stepsDB$count()

print(stepsDB$find('{}'))

# FIND
stepsDB$find('{"startDate": "2015-12-07"}')

# DROP
stepsDB$drop()

# run arbitrary commands
mongoCon <- mongo()
mongoCon$run('{"listCollections":1}')

# some commands can only be run against the `admin` database
admin <- mongo(db = "admin")
admin$run('{"listDatabases":1}')
