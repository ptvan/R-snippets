library(XML)
library(dplyr)
library(magrittr)
library(lubridate)

### parse the Apple Health XML
setwd("~/working/R-snippets/apple_health_export/")
xml <- xmlParse("export.xml")
summary(xml)
records <- XML:::xmlAttrsToDataFrame(xml["//Record"])
unique(records$type)

### extract Strava cycling data
cycling <- subset(records, sourceName == "Strava" & type == "HKQuantityTypeIdentifierDistanceCycling") %>%
           dplyr::select(creationDate, startDate, endDate, value) %>%
           mutate(value = as.numeric(as.character(value))) %>%
           rename(value = "milesCycled")
# merge so we get calories for each biking session
cycling <- cycling %>% 
           inner_join(subset(records, sourceName =="Strava" & type == "HKQuantityTypeIdentifierActiveEnergyBurned")) %>%
           dplyr::select(creationDate, startDate, endDate, value, milesCycled) %>%
           mutate(value = as.numeric(as.character(value))) %>%
           rename(value = "kcalBurned")

write.csv(cycling, "cyclingData.csv", row.names = FALSE, quote = FALSE)

#### extract RunKeeper running data  
running <- subset(records, sourceName == "Runkeeper" & type == "HKQuantityTypeIdentifierDistanceWalkingRunning") %>%
           dplyr::select(creationDate, startDate, endDate, value) %>%
           mutate(value = as.numeric(as.character(value))) %>%
           rename(value = "milesRun")

running <- running %>% 
          inner_join(subset(records, sourceName =="Runkeeper" & type == "HKQuantityTypeIdentifierActiveEnergyBurned")) %>%
          dplyr::select(creationDate, startDate, endDate, value, milesRun) %>%
          mutate(value = as.numeric(as.character(value))) %>%
          rename(value = "kcalBurned")

write.csv(running, "runningData.csv", row.names = FALSE, quote = FALSE)

#### extract step counts
steps <- subset(records, sourceName == "iPhu" & type == "HKQuantityTypeIdentifierStepCount") %>%
  dplyr::select(creationDate, startDate, endDate, value) %>%
  mutate(value = as.numeric(as.character(value))) %>%
  rename(value = "stepsWalked")

write.csv(steps, "stepsData.csv", row.names = FALSE, quote = FALSE)



  