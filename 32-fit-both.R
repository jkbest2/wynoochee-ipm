## devtools::load_all("../../src/salmonIPM")
library(salmonIPM)
options(mc.cores = parallel::detectCores(logical = FALSE))
rstan_options(auto_write = TRUE)
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
         year = as.integer(as.character(year)),
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

wyn <- read_rds("data/wyn_trap.rds") |>
  filter(species == "coho") |>
  left_join(select(bing_harvest, year, harvest_total),
            by = join_by(year)) |>
  mutate(pop = "Upper Wynoochee",
         A = 10,
         # n_age2_obs = 0.1,
         # n_age3_obs = 0.48,
         S_obs = count,
         n_W_obs = 0,
         n_H_obs = 0,
         fit_p_HOS = FALSE,
         B_take_obs = 0) |>
  left_join(select(bing, year, age_2, age_3), by = join_by(year)) |>
  mutate(n_age2_obs = age_2 / 10,
         n_age2_obs = replace_na(n_age2_obs, 0),
         n_age3_obs = age_3 / 10,
         n_age3_obs = replace_na(n_age3_obs, 0)) |>
  select(pop,
         year,
         A,
         S_obs,
         n_age2_obs,
         n_age3_obs,
         n_W_obs,
         n_H_obs,
         fit_p_HOS,
         F_rate = harvest_total,
         B_take_obs) |>
  filter(!is.na(F_rate))
         # year >= 2000)

pop_df <- bind_rows(bing_df, wyn)

wyn_fit <- salmonIPM(
  model = "IPM", life_cycle = "SS", pool_pops = FALSE,
  SR_fun = "BH",
  center = FALSE, scale = FALSE,
  fish_data = wyn,
  age_F = c(1, 1),
  age_B = c(0, 0),
  # age_S_obs = c(FALSE, TRUE),
  # age_S_eff = c(FALSE, TRUE),
  prior = list(mu_p ~ dirichlet(c(1, 9))),
  chains = 4, iter = 10000,
  control = list(adapt_delta = 0.999,
                 metric = "dense_e"))

wyn_fit

pop_fit <- salmonIPM(
  model = "IPM", life_cycle = "SS", pool_pops = TRUE,
  SR_fun = "BH",
  center = FALSE, scale = FALSE,
  fish_data = pop_df,
  age_F = c(1, 1),
  age_B = c(0, 0),
  age_S_obs = 
  prior = list(sigma_alpha ~ normal(0, 0.1),
               sigma_Rmax ~ normal(0, 0.1)),
  control = list(adapt_delta = 0.999,
                 metric = "dense_e")
)
pop_fit

shinystan::launch_shinystan(wyn_fit)
