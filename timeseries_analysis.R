library(lubridate)
library(tidyr)
library(dplyr)
library(zoo)
library(forecast)

# import data
steps <- read.csv("stepsData.csv")

# clean things up
# steps <- steps %>% 
#   mutate(creationDate = as_datetime(creationDate)) %>% 
#   mutate(startDate = as_datetime(startDate)) %>% 
#   mutate(endDate = as_datetime(endDate)) 

# day-level data
steps <- steps %>%
  mutate(creationDate = as.Date(creationDate)) %>%
  mutate(startDate = as.Date(startDate)) %>%
  mutate(endDate = as.Date(endDate)) %>%
  group_by(startDate) %>%
  summarize(stepsWalked = sum(stepsWalked)) %>%
  mutate(year = year(startDate)) %>%
  mutate(month = month(startDate)) %>%
  mutate(day = day(startDate)) 

as.ts(steps)

ma(steps, order=52, centre = T)