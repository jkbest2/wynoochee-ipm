# devtools::load_all("../../src/salmonIPM")
library(salmonIPM)
library(tidyverse)

bing <- read_rds("data/bing_coho_up.rds") |>
  mutate(sex = ifelse(type == "female", "female", "male"),
         age = ifelse(type == "jack", 2, 3)) |>
  select(-type) |>
  pivot_wider(id_cols = year,
              names_from = c(sex, age),
              values_from = count) |>
  mutate(age_2 = male_2,
         age_3 = female_3 + male_3,
         total = male_2 + male_3 + female_3)

bing_harvest <- read_rds("data/bing_harvest.rds") |>
  select(year = return_year,
         ocean_survival,
         harvest_total)

bing_df <- bing |>
  left_join(bing_harvest, by = join_by(year)) |>
  mutate(pop = factor("Bingham Creek"),
         A = 40,
         n_W_obs = 0,
         n_H_obs = 0,
         fit_p_HOS = FALSE,
         B_take_obs = 0,
         ) |>
  select(pop,
         year,
         A,
         S_obs = total,
         n_age2_obs = age_2,
         n_age3_obs = age_3,
         n_W_obs,
         n_H_obs,
         fit_p_HOS,
         F_rate = harvest_total,
         B_take_obs
) 

bing_fit <- salmonIPM(
  model = "IPM", life_cycle = "SS", pool_pops = FALSE,
  SR_fun = "BH",
  center = FALSE, scale = FALSE,
  fish_data = bing_df,
  age_F = c(1, 1),
  age_B = c(0, 0),
  # age_S_obs = c(TRUE, TRUE),
  # age_S_eff = c(TRUE, TRUE),
  prior = list(mu_p ~ dirichlet(c(0.1, 0.9)), # More informative prior on age distribution
               tau ~ gnormal(0, 0.01, 1)),    # Spawners are observed precisely due to trap-and-haul
  chains = 4, iter = 10000,
  control = list(adapt_delta = 0.99,
                 metric = "dense_e")
)

bing_fit
