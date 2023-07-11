## Libraries -------------------------------------------------------------------
library(tidyverse)  # Data manipulation and plotting
library(sf)         # Spatial data handling
library(USAboundaries)

source("R/geodata-functions.R")

## Data import -----------------------------------------------------------------
## Get the coordinate reference system (CRS) for Washington state from the
## `USAboundaries` package.
crs <- state_proj |>
  filter(state == "WA",
         # Usage section includes list of counties
         zone == "south") |> 
  pull(epsg) |>
  as.integer() |>
  st_crs()

##  Load watershed boundary for Chehalis River as a
## convenient way to bound any other spatial data.
che_hu8 <- read_wbd("8") |>
  filter(grepl("Lower Chehalis", name)) |>
  st_transform(crs)

## Load WDFW Statewide Washington Integrated Fish Distribution layer, filter for
## Coho and for the HUC-8 of interest, then drop the `M` dimension (otherwise
## lots of errors down the line) and convert to the above CRS (in this case
## converting coordinates and distances from survey feet to meters).
swifd <- read_sf("rawdata/geodata/swifd/SWIFD.shp") |>
  filter(HUC_8 == che_hu8$huc8,
         SPECIES == "Coho Salmon") |>
  st_zm() |>
  st_transform(crs)

## Wynoochee River above the reservoir -----------------------------------------
## Get Wynoochee Reservoir
wyn_res <- read_nhd("Waterbody") |>
  st_zm() |>
  st_transform(crs) |>
  filter(gnis_name == "Wynoochee Lake")

## Find minimum rotated rectangle, extend northern boundary to create a buffer
## that can be used to extract features upstream of the reservoir.
wyn_res_mrr <- st_minimum_rotated_rectangle(wyn_res)
## Finds the northern edge of the minimum rotated rectangle
wyn_mrr_top <- wyn_res_mrr |>
  st_coordinates() |>
  as_tibble() |>
  select(easting = X, northing = Y) |>
  distinct() |>
  arrange(-northing) |>
  head(n = 2)
## Calculates the slope of the northern edge
wyn_mrr_s <- wyn_mrr_top |>
  summarize(s = diff(northing) / diff(easting)) |>
  pull(s)
## Define the distance to extend the line east-west
d <- 1e3
## Use this to create a new line feature
wyn_res_top <- wyn_mrr_top |>
  mutate(sgn = sign(easting - mean(easting)),
         new_easting = easting + sgn * d,
         new_northing = northing + sgn * d * wyn_mrr_s) |>
  select(new_easting, new_northing) |>
  as.matrix() |>
  st_linestring() |>
  st_sfc(crs = crs)
## And create a single-sided buffer
wyn_up_buff <- wyn_res_top |>
  st_buffer(dist = 6e3, singleSide = TRUE)
## Which is then used to extract the features above the reservoir
wynup <- swifd |>
  st_intersection(wyn_up_buff)

## Plots used to double check geometry selection and refine the buffer distances
# plot(wyn_up_buff, col = "gray")
# plot(wyn_res, col = "skyblue", add = TRUE)
# plot(wyn_up$geometry, col = "blue", add = TRUE)

## Get Bingham Creek and tributaries -------------------------------------------
bing <- swifd |>
  filter(GNIS_NAME == "Bingham Creek") |>
  iter_intersect(filter(swifd, !grepl("Satsop", GNIS_NAME)))

## Plots to double check operations.
# plot(bing$geometry)
# plot(swifd$geometry, add = TRUE)
# plot(bing$geometry, col = "skyblue", lwd = 3, add = TRUE)

## Get Big Creek for potential strays ------------------------------------------
bigcr <- swifd |>
  filter(GNIS_NAME == "Big Creek") |>
  st_filter(filter(swifd, GNIS_NAME == "Wynoochee River"),
            .predicate = st_intersects) |>
  iter_intersect(filter(swifd, GNIS_NAME != "Wynoochee River"))

## Check with plots
# plot(bigcr$geometry)
# plot(swifd$geometry, add = TRUE)
# plot(bigcr$geometry, col = "skyblue", lwd = 3, add = TRUE)

## Postprocess the river sections ----------------------------------------------
wynup_line <- wynup |>
  st_geometry() |>
  st_union() |>
  st_sf() |>
  st_set_geometry("geometry") |>
  mutate(name = "Upper Wynoochee")
bing_line <- bing |>
  st_geometry() |>
  st_union() |>
  st_sf() |>
  st_set_geometry("geometry") |>
  mutate(name = "Bingham Creek")
bigcr_line <- bigcr |>
  st_geometry() |>
  st_union() |>
  st_sf() |>
  st_set_geometry("geometry") |>
  mutate(name = "Big Creek")
hab_df <- bind_rows(wynup_line, bing_line, bigcr_line) |>
  mutate(habitat = set_units(st_length(geometry), "km"))

write_rds(hab_df, "data/hab_df.rds")

## Final plot for reference
hab_df_plot <- hab_df |>
  mutate(habkm = paste(signif(habitat, 3), "km"),
         System = paste(name, habkm, sep = ": "))
ggplot() +
  geom_sf(data = swifd, color = "gray50") +
  geom_sf(data = wyn_res, color = "gray50", fill = "gray70") +
  geom_sf(data = hab_df_plot, aes(color = System),
          linewidth = 1) +
  theme_bw()
ggsave("figs/habitat_map.png", width = 10, height = 10)
