library(ggplot2)
library(ggrepel)
library(maps)
library(usmap)
library(tmap)
library(tidycensus)
library(SpatialEpi)
library(osmdata)
library(ggmap)
library(sf)

census_api_key("MYCENSUSAPIKEYGOESHERE", install = TRUE)
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

# load in nutria2007 sighting data and convert to Simple Features
# using WGS84 Coordinate Reference System
nutria_sightings <- read.csv("~/working/nutria2007/nutria_obs.csv")
nutria_centroids <- latlong2grid(nutria_sightings[, 1:2])
nutria_sf <- st_as_sf(nutria_sightings, 
                      coords = c("longitude", "latitude"),
                      crs = st_crs(4326)
)

tmap_mode("view")
tm_shape(nutria_sf) + tm_dots("type")

tmap_mode("plot")
nutria_map <- tm_shape(nutria_sf) + tm_dots("type")
tmap_save(nutria_map, "nutria_map.png", width=900, height=600)


WA_counties <- map_data('county', 'washington') %>% 
  select(lon = long, lat, group, id = subregion)

ggplot(WA_counties, aes(lon, lat)) + 
  geom_polygon(fill = "white", colour = "grey50") + 
  coord_quickmap()

ggplot(nutria_sf) + 
  geom_sf() + 
  coord_sf()

Seattle_bb <- getbb("Seattle") 
Shoreline_bb <- getbb("Shoreline") 
KC_bb <- getbb("King County")

Seattle_roads <- Seattle_bb %>%
  opq() %>%
  add_osm_feature(key = "highway", value = c("tertiary", "residential")) %>%
  osmdata_sf() 

Shoreline_roads <- Shoreline_bb %>%
  opq() %>%
  add_osm_feature(key = "highway", value = c("tertiary", "residential")) %>%
  osmdata_sf() 

Seattle_highways <- Seattle_bb %>%
  opq() %>%
  add_osm_feature(key = "highway", value = c("motorway", "primary", "secondary")) %>%
  osmdata_sf() 

KC_roads <- KC_bb %>%
  opq() %>%
  add_osm_feature(key = "highway", value = c("tertiary", "residential")) %>%
  osmdata_sf() 

KC_highways <- KC_bb %>%
  opq() %>%
  add_osm_feature(key = "highway", value = c("motorway", "primary", "secondary")) %>%
  osmdata_sf() 
  
ggplot() +
  geom_sf(data = Seattle_roads$osm_lines,
          inherit.aes = FALSE,
          color = "blue",
          size = 0.2) +
  geom_sf(data = Shoreline_roads$osm_lines,
          inherit.aes = FALSE,
          color = "green",  
          size = 0.1) +
  geom_sf(data = nutria_sf,
          color = "red",
          size = 0.5
          ) +
  theme_void()
