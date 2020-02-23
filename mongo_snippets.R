library(mongolite)
library(dplyr)
library(magrittr)

m <- mongo(collection = "test", db = "test", url = "mongodb://localhost",
           verbose = FALSE, options = ssl_options())

# get iPhone health data
stepsData <- read.csv("~/working/datasets/iphone_health/stepsData.csv")
cyclingData <- read.csv("~/working/datasets/iphone_health/cyclingData.csv")

# aggregate into day-level data and convert into dates POSIXct for Mongo 

stepsData$startDate = substr(stepsData$startDate, 1, 10)
stepsData <- aggregate(stepsData$stepsWalked, by=list(stepsData$startDate), sum)
colnames(stepsData) <- c("startDate", "stepsWalked")
#stepsData$startDate <- as.POSIXct.Date(as.Date(stepsData$startDate))

cyclingData$startDate = substr(cyclingData$startDate, 1, 10)
cyclingData <- aggregate(cyclingData$milesCycled, by=list(cyclingData$startDate), sum)
colnames(cyclingData) <- c("startDate", "milesCycled")
# cyclingData$startDate <- as.POSIXct.Date(as.Date(cyclingData$startDate))

# connect 
stepsDB <- mongo(db="mydb", collection="steps")
cyclingDB <- mongo(db="mydb", collection="cycling")


# INSERT
stepsDB$insert(stepsData)
cyclingDB$insert(cyclingData)
stepsDB$count()
cyclingDB$count()


# FIND/SELECT
# show all
print(stepsDB$find('{}'))

stepsDB$find('{}', limit = 10)

# exact match
stepsDB$find(query = '{"startDate": "2016-10-07"}')
cyclingDB$find(query = '{"startDate": "2016-10-07"}')

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

# lookup, a LEFT OUTER JOIN for noSQL...
# NOTE: THIS DOESN'T APPEAR TO WORK AS OF FEBRUARY 2020:
# https://github.com/jeroen/mongolite/issues/183
out <- cyclingDB$aggregate(
  '[
    {
      "$lookup":
        {
          "from": "stepsDB",
          "localField": "startDate",
          "foreignField": "startDate",
          "as": "health"
        }
   }
]'
)



# RUN ARBITRARY COMMANDS
mongoCon <- mongo()
mongoCon$run('{"listCollections":1}')

# some commands can only be run against the `admin` database
admin <- mongo(db = "admin")
admin$run('{"listDatabases":1}')
