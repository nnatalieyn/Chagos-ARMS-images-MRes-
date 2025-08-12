library(ggplot2)     # plotting
library(sf)          # deal with spatial tables
library(raster)
library(sp)
library(terra)       # handle raster data
library(tidyterra)   # handeling terra
library(tidyverse)
library(ggspatial)   # some useful ggplot2 additions
library(ggrepel)     # handeling text on plots
library(patchwork)   # for assembling multiple ggplots
library(mapdata)     # higher resolution maps
library(rnaturalearth)
library(rnaturalearthdata)

# map of chagos in the world ----------------------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
chagos_epicenter <- st_point(c(72,-6.5))
chagos_epicenter <- st_sfc(chagos_epicenter, crs = 4326) #make into sf object with WGS84
chagos_epicenter <- st_transform(chagos_epicenter, 3857) #reproject to metric CRS like web mercator
chagos_circle <- st_buffer(chagos_epicenter, dist = 200000) #200km radius
chagos_circle <- st_transform(chagos_circle, 4356) #transform back to WGS84

ggplot() +
  geom_sf(data = world, fill = "black", color = "black") +
  geom_sf(data = chagos_circle, fill = NA, color = "red", linewidth = 1.5) +
  coord_sf(xlim = c(50, 85), ylim = c(-10, 25), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.3) +  # bottom-left corner
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

# map of chagos -----------------------------------------------------------

#chagos is too small to show up on standard map packages
my_region <- "Chagos Archipelago"
lca <- map_data("worldHires", region = my_region)
lca |> glimpse()

#sites
sites<- data.frame(
  name = c("Anglaise", "Moresby", "Coin", "Court"),
  lat = c(-5.3422, -5.237, -5.4506, -5.30883333),
  long = c(72.225, 71.8339, 71.7678, 72.25236667)
)

#check scale 
lca_sf <- st_as_sf(lca, coords = c("long", "lat"), crs = 4326) |>
  st_transform(crs = 32743) #south indian ocean/UTM zone 43S
ggplot() +
  geom_sf(data = lca_sf, fill = "black", color = "black") +
  coord_sf(xlim = c(71.75, 72.30), ylim = c(5.5, 5.2), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  coord_sf() +
  theme_minimal()

#plot
glp <- 
  map_data("worldHires", region = "Chagos Archipelago") |> 
  ggplot(aes(long, lat, group = group)) +
  geom_path()
glp

#add aspect ration (ratio of dimensions between x and y)
glp +
  geom_polygon(color = "black") +
  coord_quickmap(xlim = c(71.75, 72.35), ylim = c(-5.5, -5.2)) +
  geom_point(data = sites, aes(x = long, y = lat), color = "red", size = 3, inherit.aes = FALSE) +
  annotation_scale(
    location = "br", 
    width_hint = 0.2,    
    style = "bar"        
  ) +
  scale_x_continuous(
    name = "Longitude",
    breaks = seq(71.8, 72.3, 0.1),
    labels = function(x) paste0(x, "°E")
  ) +
  scale_y_continuous(
    name = "Latitude",
    breaks = seq(-5.5, -5.2, 0.1),
    labels = function(y) paste0(abs(y), "°S") 
  )+
  theme_bw()
