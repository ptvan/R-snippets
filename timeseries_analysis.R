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
  summarize(stepsWalked = sum(stepsWalked))

# create a ts
steps_ts <- ts(steps$stepsWalked, start = c(2015,yday(steps$startDate[1])), frequency = 365)

# find a weekly trend
steps_trend <- ma(steps_ts, order=7, centre = T)

# plot data with trend overlaid
plot(steps_ts, main="step counts")
lines(steps_trend, col="red")

# detrend
steps_detrended_ts <- steps_ts / steps_trend
plot(steps_detrended_ts, main="step counts, detrended")

# decompose into seasonal, trend and random
steps_stl <- stl(steps_ts, "periodic")

steps_stl_seasonal <- steps_stl$time.series[,1]
steps_stl_trend <- steps_stl$time.series[,2]
steps_stl_random <- steps_stl$time.series[,3]

par(mfrow=c(3,1))
plot(steps_stl_seasonal, main="seasonal component")
plot(steps_stl_trend, main="trend component")
plot(steps_stl_random, main="random component")