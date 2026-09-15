# Map of US with state outlines, north arrow, scale bar, and site markers

library(sf)
library(ggplot2)
library(ggspatial)
library(usmap)
library(ggrepel)
library(here)

# Get US map data with states (exclude PR)
us_map_data <- us_map()
us_map_data <- us_map_data[us_map_data$fips != "72", ]

# Define locations
sites <- data.frame(
  name = c("Ottumwa, IA", "Boise, ID", "Riverside, CA", "Boston, MA"),
  lon = c(-92.4186, -116.202, -117.396156, -71.0571),
  lat = c(41.0194, 43.615, 33.953350, 42.3611)
)

# Transform coordinates to usmap projection (returns an sf object)
sites_transformed <- usmap_transform(sites, input_names = c("lon", "lat"))

# Extract point coordinates for use with ggrepel (needs plain x/y, not an sf geom column)
pt_coords <- sf::st_coordinates(sites_transformed)
sites_transformed$x <- pt_coords[, "X"]
sites_transformed$y <- pt_coords[, "Y"]

# Create map
map <- ggplot() +
  geom_sf(data = us_map_data, fill = "white", col = "gray30", linewidth = 0.3) +
  geom_sf(data = sites_transformed, col = "darkred",
          size = 3, shape = 21, fill = "white", stroke = 1.2) +
  geom_label_repel(data = sites_transformed,
                    aes(x = x, y = y, label = name),
                    size = 3.2, col = "black", fill = "white",
                    box.padding = 0.6, min.segment.length = 0,
                    segment.col = "gray40", seed = 42) +
  annotation_north_arrow(
    location = "tl",
    which_north = "grid",
    pad_x = unit(0.15, "in"),
    pad_y = unit(0.15, "in"),
    style = north_arrow_fancy_orienteering()
  ) +
  annotation_scale(
    location = "br",
    pad_x = unit(0.3, "in"),
    pad_y = unit(0.3, "in"),
    unit_category = "metric"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  coord_sf()

print(map) # AK and HI are repositioned as insets and not to scale

ggsave(here("output/figures/Fig1_map.png"), width = 9.5, height = 3.5, dpi = 400)
