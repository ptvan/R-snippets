library(lubridate)
library(tidyr)
library(dplyr)
library(zoo)
library(forecast)
library(changepoint)

### STEPS DATA
steps <- read.csv("stepsData.csv")

# clean things up
# steps <- steps %>% 
#   mutate(creationDate = as_datetime(creationDate)) %>% 
#   mutate(startDate = as_datetime(startDate)) %>% 
#   mutate(endDate = as_datetime(endDate)) 

# day-level data
steps <- steps %>%
  mutate(startDate = as.Date(startDate)) %>%
  group_by(startDate) %>%
  summarize(stepsWalked = sum(stepsWalked))

# create a ts object
steps_ts <- ts(steps$stepsWalked, start = c(2015,yday(steps$startDate[1])), frequency = 365)

# basic plotting
monthplot(steps_ts)

# get weekly trends using moving average
steps_trend <- ma(steps_ts, order=52, centre = T)

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


# changepoint detection

# try to figure out our own penalty vector...
# penalties <- seq(0, 20,.2)
# 
# tryPenalty <- function(data, pen){
#   ans <- cpt.mean(data, test.stat="Normal", method = "PELT", penalty = "Manual", pen.value = penalties) 
#   return(length(cpts(ans)) +1)
# }
# 
# tryPenalty(steps_ts, penalties)
# 
# elbowplotData <- unlist(lapply(penalties, function(p) 
#   tryPenalty(data = steps_ts, pen = p)))

# ... or just let the function do it...
plot(cpt.mean(steps_ts, penalty = "Hannan-Quinn"), ylab="steps", main="change point detection by Hannan-Quinn")

### BIKING DISTANCE DATA
biking <- read.csv("cyclingData.csv", header=T)
biking <- biking %>%
  mutate(startDate = as.Date(startDate)) %>%
  group_by(startDate) %>%
  summarize(kcalBurned = sum(kcalBurned),milesCycled = sum(milesCycled)) %>%
  mutate(startDate = as.character(startDate)) 

# fill in data for missing dates & merge with real data
missing_dates <- as.character(seq(as.Date(biking$startDate[1]), as.Date(biking$startDate[nrow(biking)]), by=1)) %>%
  setdiff(.,  biking$startDate) %>%
  data.frame()
colnames(missing_dates) <- c("startDate")
missing_dates$kcalBurned = NA
missing_dates$milesCycled = NA  
biking <- rbind(biking, missing_dates) %>% 
  read.zoo(format="%Y-%m-%d")

# convert steps data.frame into a zoo object
steps <- read.zoo(steps,format="%Y-%m-%d")

#### COMBINED DATA

# merging steps and biking zoo objects, filling out NAs with zero
health <- merge(biking, steps) %>% 
  na.fill(0)

# calculate ACF & partial ACF to test for stationarity 
acf(health$milesCycled, na.action = na.pass, main="milesCycled ACF")
pacf(health$steps, na.action = na.pass, main="steps pACF")

