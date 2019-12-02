library(ggplot2)
library(ggrepel)
library(maps)
usa <- map_data("usa")


places <- data.frame(
  long = c(-79.995888, -122.306417, -122.306445),
  lat = c(40.440624, 47.644855, 47.644849),
  names = c("CMU\n(2009 - 2014)", "ISB\n(2006 - 2009)", "FredHutch\n(2015 - now)"),
  stringsAsFactors = FALSE
)  

continental_US <- ggplot() +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill=NA, color="gray") +
  coord_fixed(1.3) +
  theme_void()

continental_US + 
  geom_point(data = places, aes(x = long, y = lat), color = "red", size = 3, alpha = 0.2) +
  geom_text_repel(data = places, aes(x = long, y = lat, label=names), color = "red", size = 3
                  )

library(usmap)
plot_usmap("states", exclude=c("AK","HI"), labels = T)
