library(tidyverse)
library(readxl)

if (!dir.exists("data")) {
  dir.create("data")
}

bing_coho_up <- read_xlsx("rawdata/2000-2019 Upstream CohoBingham.xlsx",
                     col_types = rep("numeric", 4),
                     range = "A1:D21",
                     # n_max = 20, # Skip the summary rows at the end
                     .name_repair = ~ str_split_i(tolower(.), pattern = "s ", i = 1)) |>
  pivot_longer(-year, names_to = "type", values_to = "count") |>
  mutate(year = factor(year),
         type = factor(type, levels = c("female", "male", "jack")))
write_rds(bing_coho_up, "data/bing_coho_up.rds")

ggplot(bing_coho_up, aes(x = year, y = count, fill = type)) +
  geom_col(position = position_dodge())
  
bing_coho_smolt_spawner <- read_xlsx("rawdata/SPAWNERSMOLT.xlsx", skip = 2,
                                col_types = c("numeric", "numeric", "text", "numeric", "numeric", "numeric"),
                                range = "A3:F27") |>
  slice(-1) |> # Drop first, empty row
  set_names("ocean_year", "smolts", "spawn_season", "males", "females", "jacks") |>
  mutate(spawn_year = as.numeric(str_split_i(spawn_season, "-", 1)))
write_rds(bing_coho_smolt_spawner, "data/bing_coho_smolt_spawner.rds")

ggplot(bing_coho_smolt_spawner, aes(x = females, y = smolts, color = ocean_year)) +
  geom_point() +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0, 0.05, 0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0, 0.05, 0))

bing_coho_smolt_spawner |>
  mutate(smolt_per_female = smolts / females) |>
  ggplot(aes(x = ocean_year, y = smolt_per_female)) +
  geom_point() +
  geom_line()

fix_wyn_trap_year <- function(year) {
  ## Split out the first of the two years in the season
  if (year == 1970) {
    ny <- 1969 # Not sure why 1970 is listed twice here, so push this one
                    # back a year?
  } else if (grepl("/", year)) {
    ny <- as.numeric(str_split_i(year, "/", 1))
  } else if (grepl("-", year)) {
    ny <- as.numeric(str_split_i(year, "-", 1))
  }
  
  ## Add the century where appropriate
  if (ny < 100) {
    if (ny >= 70) {
      ny <- ny + 1900
    } else if (ny < 23) {
      ny <- ny + 2000
    }
  }
  return(ny)
}


wyn_trap <- read_xls("rawdata/Wynoochee Trap Totals 1970-2023.xls",
                     range = "A2:E56",
                     .name_repair = tolower) |>
  mutate(year = map_dbl(year, fix_wyn_trap_year)) |>
  rename(steelhead_w = `w sthd`,
         steelhead_nor = `stlhd nor`) |>
  pivot_longer(-year, names_to = "species", values_to = "count",
               values_drop_na = TRUE) |>
  mutate(type = ifelse(species == "steelhead_nor", "NOR", "Wild"),
         species = fct_collapse(species, steelhead = c("steelhead_w", "steelhead_nor")))
write_rds(wyn_trap, "data/wyn_trap.rds")

ggplot(wyn_trap, aes(x = year, y = count, color = species, shape = type)) +
  geom_point() +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0, 0.1, 0))
