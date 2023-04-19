library(ggplot2)
library(ggrepel)
library(maps)
library(usmap)
library(tidycensus)
library(SpatialEpi)

census_api_key("MYAPIKEYGOESHERE", install = TRUE)
options(tigris_use_cache = TRUE)

# basic plotting of locations 
plot_usmap("states", exclude = c("AK","HI"), labels = T)

usa <- map_data("usa")

places <- data.frame(
  long = c(-79.995888, -122.306417, -122.306445),
  lat = c(40.440624, 47.644855, 47.644849),
  names = c("CMU\n(2009 - 2014)", "ISB\n(2006 - 2009)", "FredHutch\n(2015 - now)"),
  stringsAsFactors = FALSE
)  

continental_US <- ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group), fill = NA, color = "gray") +
  coord_fixed(1.3) +
  theme_void()

continental_US + 
  geom_point(data = places, aes(x = long, y = lat), color = "red", size = 3, alpha = 0.2) +
  geom_text_repel(data = places, aes(x = long, y = lat, label=names), color = "red", size = 3) +
  coord_fixed(1.3) +
  theme_void()

# load in nutria2007 sighting data and convert to 
nutria_sightings <- read.csv("~/working/nutria2007/nutria_obs.csv")
nutria_centroids <- latlong2grid(nutria_sightings[, 1:2])

# load in King County from US Census ACS 2020 data
king_county <- get_acs(
  state = "WA",
  county = "King",
  geography = "tract",
  variables = "B19013_001",
  geometry = TRUE,
  year = 2020
)

plot(king_county["estimate"])

king_county %>%
  ggplot(aes(fill = estimate)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c(option = "magma") 

