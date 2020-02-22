library(mongolite)
library(dplyr)
library(magrittr)

m <- mongo(collection = "test", db = "test", url = "mongodb://localhost",
           verbose = FALSE, options = ssl_options())

stepsData <- read.csv("~/working/datasets/iphone_health/stepsData.csv")

# day-level stepcount data
stepsData <- stepsData %>%
  mutate(startDate = strptime(startDate, "%Y-%m-%dT%H:%M:%SZ", 'UTC') %>%
  group_by(startDate) %>%
  summarize(stepsWalked = sum(stepsWalked))

# INSERT
stepsDB <- mongo("steps")
stepsDB$insert(stepsData)
stepsDB$count()

# FIND/SELECT
# show all
print(stepsDB$find('{}'))

# exact match
stepsDB$find(query = '{"startDate": "2015-12-07"}')

# sort, "_id" is included by default so we have to explicitly exclude it 
stepsDB$find(fields = '{"startDate":true, "stepsWalked":true, "_id":false}'
             , sort='{"stepsWalked": -1}')

# date range
stepsDB$find(
  query = '{"startDate": { "$gte" : { "$date" : "2017-01-01T00:00:00Z" }}}',
  fields = '{"startDate" : true, "stepsWalked":true, "_id": false}'
)

# UPDATE
stepsDB$update('{"startDate": "2015-12-07"}', '{"$set":{"stepsWalked":10000}}')

# DROP
stepsDB$drop()

# run arbitrary commands
mongoCon <- mongo()
mongoCon$run('{"listCollections":1}')

# some commands can only be run against the `admin` database
admin <- mongo(db = "admin")
admin$run('{"listDatabases":1}')
