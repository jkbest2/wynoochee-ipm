## Packages --------------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(ggspatial)
library(sf)
library(USAboundaries)
library(patchwork)

source("R/geodata-functions.R")

## Load data -------------------------------------------------------------------
crs <- state_proj |>
  filter(state == "WA",
         zone == "north") |>
  pull(epsg) |>
  as.integer() |>
  st_crs()

## State boundaries and cities for context
wa_state <- us_states(state = "WA", resolution = "high") |>
  st_transform(crs)

## Load watershed boundary for Chehalis River as a convenient way to bound any
## other spatial data.
che_hu8 <- read_wbd("8") |>
  filter(grepl("Lower Chehalis", name)) |>
  st_transform(crs)

## Get the location of the Wynoochee Hydropower Dam
wyn_hydro <- read_sf("rawdata/geodata/Power_Plants/Power_Plants.shp") |>
  filter(Plant_Name == "Wynoochee") |>
  st_transform(crs)
## Get Wishkah Cr flowlines to filter out the "other" Big Creek
wis_fl <- read_nhd("Flowline") |>
  st_zm() |>
  st_transform(crs) |>
  filter(grepl("Wishkah", gnis_name)) |>
  select() |>
  st_union()

## Get NHD Flowlines
che_fl <- read_nhd("Flowline") |>
  st_zm() |>
  st_transform(crs) |>
  filter(grepl("Chehalis River|Wynoochee|Big Creek|^Satsop River|East Fork Satsop River|Bingham Creek", gnis_name),
  ## filter(grepl("Wynoochee|Bingham", gnis_name),
         ## This GNIS ID represents a separate Big Creek, which is a tributary
         ## of the Wishkah River
         gnis_id != "01516480") |>
  select(gnis_id, gnis_name) |>
  group_by(gnis_id, gnis_name) |>
  summarize(geometry = st_union(geometry),
            .groups = "drop")

che_fl2 <- read_nhd("Flowline") |>
  st_zm() |>
  st_transform(crs) |>
  filter(visibility >= 5000000) |>
  summarize(geometry = st_union(geometry),
            .by = c(gnis_id, gnis_name))

## Get waterbodies
che_wb <- read_nhd("Waterbody") |>
  st_zm() |>
  st_transform(crs) |>
  st_filter(che_fl, .predicate = st_intersects)

che_area <- read_nhd("Area") |>
  st_zm() |>
  st_transform(crs) |>
  st_filter(che_fl, .predicate = st_intersects)

che_nhd <- rbind(select(che_fl),
                 select(che_wb),
                 select(che_area))

## Get line features - the only one that intersects the flowline is the
## Wynoochee Dam!
che_line <- read_nhd("Line") |>
  st_zm() |>
  st_transform(crs) |>
  st_filter(che_fl, .predicate = st_intersects)

## Get cities in the vicinity of the map
wa_cities <- us_cities(states = "WA") |>
  st_transform(crs)|>
  st_crop(st_as_sfc(st_bbox(che_nhd)))

## Create a spatial point object for the Bingham Creek Fish Hatchery
bng_fh <- st_sf(geometry = st_sfc(st_point(c(-123.40086239525105, 47.14645959145278))),
                 crs = st_crs(4326)) |>
  st_transform(crs) |>
  mutate(name = "Bingham Creek Fish Hatchery")


## River map
che_labs <- che_fl |>
  select(gnis_name) |>
  filter(gnis_name %in% c("Chehalis River",
                          "Satsop River",
                          ## "Big Creek",
                          "Bingham Creek",
                          ## "East Fork Satsop River",
                          "Wynoochee River")) |>
  rbind(che_wb |>
          filter(gnis_name == "Wynoochee Lake") |>
          select(gnis_name)) |>
  rbind(wyn_hydro |>
          rename(gnis_name = Plant_Name) |>
          select(gnis_name) |>
          mutate(gnis_name = "Wynoochee Dam")) |>
  rename(name = gnis_name) |>
  rbind(mutate(bng_fh, name = "Bingham Hatchery"))
## Manual label nudges to keep things looking nice.
lab_nudge <- tribble(
  ~ name, ~ x, ~ y,
  "Chehalis River", -2300, -4500,
  "Satsop River", 7250, 2000,
  ## "Big Creek", -4000, -250,
  "Bingham Creek", 9500, 0,
  ## "East Fork Satsop River", 8500, 0,
  "Wynoochee River", 9500, 5000,
  "Wynoochee Lake", 10500, 0,
  "Wynoochee Dam", 10000, -1750,
  "Bingham Hatchery", 11000, -1250
)
bbox <- st_bbox(rbind(select(che_nhd), select(che_hu8)))
bbox_poly <- st_sf(geometry = st_as_sfc(bbox)) |>
  st_buffer(0)
che_map <- ggplot() +
  geom_sf(data = wa_state, alpha = 0.4, fill = NA) +
  ## geom_sf(data = che_hu8, fill = NA) +
  geom_sf(data = filter(che_fl2, !is.na(gnis_name)), color = "lightblue") +
  geom_sf(data = che_nhd, color = "skyblue", fill = "skyblue") +
  geom_sf(data = che_line, color = "black") +
  geom_sf(data = wyn_hydro) +
  geom_sf(data = bng_fh) +
  geom_label(data = che_labs,
                   aes(geometry = geometry,
                       label = name),
                   stat = "sf_coordinates",
                   nudge_x = lab_nudge$x,
                   nudge_y = lab_nudge$y) +
  annotation_scale(location = "br",
                   width_hint = 0.5) +
  coord_sf(xlim = bbox[c(1, 3)] + c(0, 5e3),
           ylim = bbox[c(2, 4)]) +
  labs(x = NULL, y = NULL)
che_map

wa_bbox <- st_bbox(wa_state)
wa_map <- ggplot() +
  geom_sf(data = wa_state) +#, fill = NA) +
  geom_sf(data = bbox_poly, fill = NA, size = 2) +
  # geom_sf(data = che_hu8) +
  geom_sf(data = che_nhd, color = "skyblue", fill = "skyblue") +
  geom_sf(data = che_line, color = "black") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(xlim = c(wa_bbox[1], mean(wa_bbox[c(1, 3)])),
           ylim = wa_bbox[c(2, 4)])
che_map +
  inset_element(wa_map, left = 0.025, bottom = 0.6, right = 0.4, top = 0.95)

ggsave("figs/che-map.png",
       width = 7, height = 7)

### Habitat map -----------------
wyn_hab <- read_rds("data/hab_df.rds") |>
  filter(name == "Upper Wynoochee")

che_fl3 <- read_nhd("Flowline") |>


bbox <- st_bbox(wyn_hab)
wdfw_darkblue <- "#003F6B"
ggplot() +
  geom_sf(data = che_wb, fill = wdfw_darkblue, color = wdfw_darkblue) +
  geom_sf(data = che_fl, color = wdfw_darkblue) +
  geom_sf(data = che_area, fill = wdfw_darkblue, color = wdfw_darkblue) +
  geom_sf(data = wyn_hab, linewidth = 1.5, color = "#92D4DA") +
  coord_sf(xlim = bbox[c(1, 3)] + c(-4e3, 2e3),
           ylim = bbox[c(2, 4)] + c(-1e3, 2e3)) +
  annotation_scale(location = "br")

ggsave("figs/upper_wyn_hab.png",
       width = 7, height = 7)
