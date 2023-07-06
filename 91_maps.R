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
  filter(grepl("Chehalis River|Wynoochee|Big Creek|Satsop|Bingham", gnis_name),
         ## This GNIS ID represents a separate Big Creek, which is a tributary
         ## of the Wishkah River
         gnis_id != "01516480") |>
  select(gnis_id, gnis_name) |>
  group_by(gnis_id, gnis_name) |>
  summarize(geometry = st_union(geometry),
            .groups = "drop")

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

che_labs <- che_fl |>
  select(gnis_name) |>
  filter(gnis_name %in% c("Chehalis River", "Satsop River", "Big Creek",
                          "Bingham Creek", "East Fork Satsop River",
                          "Wynoochee River")) |>
  rbind(che_wb |>
          filter(gnis_name == "Wynoochee Lake") |>
          select(gnis_name)) |>
  rbind(wyn_hydro |>
          rename(gnis_name = Plant_Name) |>
          select(gnis_name) |>
          mutate(gnis_name = "Wynoochee Dam")) |>
  rename(name = gnis_name) |>
  rbind(bng_fh)
## Manual label nudges to keep things looking nice.
lab_nudge <- tribble(
  ~ name, ~ x, ~ y,
  "Chehalis River", 0, -2500,
  "Satsop River", 4000, 2000,
  "Big Creek", -4000, -250,
  "Bingham Creek", 6000, 0,
  "East Fork Satsop River", 8500, 0,
  "Wynoochee River", -7000, 0,
  "Wynoochee Lake", 7000, 0,
  "Wynoochee Dam", 6000, -1250,
  "Bingham Fish Hatchery", 10000, -500
)

## Maps ------------------------------------------------------------------------
bbox <- st_bbox(rbind(select(che_nhd), select(che_hu8)))
bbox_poly <- st_sf(geometry = st_as_sfc(bbox)) |>
  st_buffer(5e3)
che_map <- ggplot() +
  geom_sf(data = wa_state, alpha = 0.4, fill = NA) +
  # geom_sf(data = che_hu8, fill = NA) +
  geom_sf(data = che_nhd, color = "skyblue", fill = "skyblue") +
  geom_sf(data = che_line, color = "black") +
  geom_sf(data = wyn_hydro) +
  geom_sf(data = bng_fh) +
  geom_sf(data = wa_cities, alpha = 0.4) +
  geom_text(data = wa_cities,
            aes(geometry = geometry,
                label = city),
            alpha = 0.4,
            vjust = "top",
            hjust = "left",
            stat = "sf_coordinates") +
  geom_label(data = che_labs,
                   aes(geometry = geometry,
                       label = name),
                   stat = "sf_coordinates",
                   nudge_x = lab_nudge$x,
                   nudge_y = lab_nudge$y) +
  annotation_scale(location = "br",
                   width_hint = 0.5) +
  coord_sf(xlim = bbox[c(1, 3)],
           ylim = bbox[c(2, 4)])

## Create a map of the Western half of the state to use as an inset
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

## Put them together and save them
che_map +
  inset_element(wa_map, left = 0.025, bottom = 0.6, right = 0.4, top = 0.975)
ggsave("figs/wyn-che-map.png",
       width = 12, height = 12)
