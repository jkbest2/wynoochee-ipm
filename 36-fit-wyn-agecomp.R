library(salmonIPM)
options(mc.cores = parallel::detectCores(logical = FALSE))
rstan_options(auto_write = TRUE)
library(tidyverse)
library(posterior)
library(units)

if (!dir.exists("data/agecomp"))
  dir.create("data/agecomp")
data_dir <- "data/agecomp"

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
         survival_escapement, # Overall smolt -> adult survival
         harvest_total)

wyn_hab <- read_rds("data/hab_df.rds") |>
  sf::st_drop_geometry() |>
  filter(name == "Upper Wynoochee") |>
  pluck("habitat") |>
  units::drop_units()


## Use overall fraction jacks to back-calculate smolts
frac_jack <- sum(bing$age_2) / sum(bing$total)

## Get the data and prep for model
wyn <- read_rds("data/wyn_trap.rds") |>
  filter(species == "coho") |>
  left_join(bing_harvest, by = join_by(year)) |>
  filter(!is.na(harvest_total)) |>
  mutate(pop = "Upper Wynoochee",
         A = wyn_hab,
         n_W_obs = 0,
         n_H_obs = 0,
         fit_p_HOS = FALSE,
         B_take_obs = 0) |>
  rename(S_obs = count) |>
  left_join(select(bing, year, age_2, age_3), by = join_by(year)) |>
  ## Find age structure
  mutate(n_age2_obs = replace_na(age_2, 0),
         n_age3_obs = replace_na(age_3, 0))

## Use overall Bingham age composition to back-calculate smolt abundance
wyn_oac <- wyn |>
  mutate(
    frac_jack = sum(n_age2_obs) / sum(n_age2_obs + n_age3_obs),
    jack_smolts = lead(S_obs, 2) * lead(frac_jack, 2) /
      lead(survival_escapement, 3), # Best assumption of ocean survival is that it tracks the brood year adult ocean survival
    mature_smolts = lead(S_obs, 3) * (1 - lead(frac_jack, 3)) /
      lead(survival_escapement, 3),
    M_obs = jack_smolts + mature_smolts) |>
  select(pop,
         year,
         A,
         M_obs,
         S_obs,
         n_age2_obs,
         n_age3_obs,
         n_W_obs,
         n_H_obs,
         fit_p_HOS,
         F_rate = harvest_total,
         B_take_obs) |>
  filter(!is.na(M_obs),
         !is.na(F_rate),
         !is.na(S_obs))
write_rds(wyn_oac, file.path(data_dir, "wyn_oac.rds"))

## Use Bingham age composition where available to back-calculate smolt abundance
wyn_bac <- wyn |>
  mutate(
    frac_jack = age_2 / (age_2 + age_3),
    ## Fill in the fraction of jacks using the overall mean. This allows
    ## us to fill in many more years of data.
    frac_jack = replace_na(
      frac_jack,
      sum(age_2, na.rm = TRUE) / sum(age_2 + age_3, na.rm = TRUE)),
    ## Back-calculate the number of smolts we would have seen (if we were
    ## looking). This treats the overall smolt-to-adult survival as known, based
    ## on the Bingham data. It also assumes that jacks share the same mortality
    ## rate *as the mature fish they return with*.
    jack_smolts = lead(S_obs, 2) * lead(frac_jack, 2) /
      lead(survival_escapement, 3), # Best assumption of ocean survival is that it tracks the brood year adult ocean survival
    mature_smolts = lead(S_obs, 3) * (1 - lead(frac_jack, 3)) /
      lead(survival_escapement, 3),
    M_obs = jack_smolts + mature_smolts) |>
  select(pop,
         year,
         A,
         M_obs,
         S_obs,
         n_age2_obs,
         n_age3_obs,
         n_W_obs,
         n_H_obs,
         fit_p_HOS,
         F_rate = harvest_total,
         B_take_obs) |>
  filter(!is.na(M_obs),
         !is.na(F_rate),
         !is.na(S_obs))
write_rds(wyn_bac, file.path(data_dir, "wyn_bac.rds"))

## Compare smolt abundance between age compostion methods
ggplot() +
  geom_line(data = wyn_oac,
            aes(x = year, y = M_obs),
            linetype = "dashed") +
  geom_line(data = wyn_bac,
            aes(x = year, y = M_obs))
ggsave(file.path(data_dir, "smolts_agecomp.png"))

## Use a single observation of each age class for each year
wyn_oac1 <- wyn_oac |>
  mutate(n_age2_obs = 1,
         n_age3_obs = 9)

## Fit with Jeffreys prior on age composition
wyn_jeff_fit <- salmonIPM(
  model = "IPM", life_cycle = "SMS", pool_pops = FALSE,
  SR_fun = "BH",
  ages = list(M = 1),
  center = FALSE, scale = FALSE,
  fish_data = wyn_oac1,
  age_F = c(1, 1),
  age_B = c(0, 0),
  ## age_S_obs = c(FALSE, TRUE),
  ## age_S_eff = c(FALSE, TRUE),
  prior = list(
    ## Priors from the Barrowman et al. (2003) meta-analysis. alpha prior
    ## includes correction assuming 45% of spawners are female and 50% of
    ## smolts are female. Mmax does the same, doubling the capacity per unit
    ## habitat to account for male smolts.
    alpha ~ lognormal(4.27 + log(0.9), sqrt(0.18^2 + 0.43^2)),
    Mmax ~ lognormal(6.58 + log(2), sqrt(0.18^2 + 0.64^2)),
    ## mu_MS ~ beta(), # DEFAULT
    ## Age distribution is informed by Bingham data; this is a semi-informative
    ## prior to help keep things reasonable
    ## mu_p ~ dirichlet(c(1, 1)), # DEFAULT
    mu_p ~ dirichlet(c(0.5, 0.5)),
    ## Spawners are observed precisely during trap-and-haul upstream, but
    ## allowing variance to get too close to zero results in divergent
    ## transitions. This is why the default prior has an explicitly zero-avoiding
    ## prior.
    ## , tau_S ~ gnormal(1, 0.85, 30) # DEFAULT
    tau_S ~ gnormal(0.25, 0.075, 3) # Keeps values > 0.1, but smaller than default
    ## Smolts are very approximate because they are back-calculated from Bingham
    ## SAR rates, so
    ## , tau_M ~ gnormal(1, 0.85, 30) # DEFAULT
    ## , tau_M ~ gnormal(2, 1, 2) # Wider and larger than default, might need to
    ##                            # avoid small values more
  ),
  chains = 4, iter = 10000,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 10,
                 ## "diag_e" works better than "dense_e"
                 metric = "diag_e")
)
write_rds(wyn_jeff_fit, file.path(data_dir, "wyn_jeff_fit.rds"))

wyn_jeff_post <- as_draws_rvars(wyn_jeff_fit)
write_rds(wyn_jeff_post, file.path(data_dir, "wyn_jeff_post.rds"))

## Find the parameterization of the beta that will put 5% probability that jacks
## are more than 10% of the returners
b90 <- uniroot(\(b) qbeta(0.05, b, 1) - 0.9, interval = c(28, 29))$root

## Fit with informative prior on age composition
wyn_info_fit <- salmonIPM(
  model = "IPM", life_cycle = "SMS", pool_pops = FALSE,
  SR_fun = "BH",
  ages = list(M = 1),
  center = FALSE, scale = FALSE,
  fish_data = wyn_oac1,
  age_F = c(1, 1),
  age_B = c(0, 0),
  ## age_S_obs = c(FALSE, TRUE),
  ## age_S_eff = c(FALSE, TRUE),
  prior = list(
    ## Priors from the Barrowman et al. (2003) meta-analysis. alpha prior
    ## includes correction assuming 45% of spawners are female and 50% of
    ## smolts are female. Mmax does the same, doubling the capacity per unit
    ## habitat to account for male smolts.
    alpha ~ lognormal(4.27 + log(0.9), sqrt(0.18^2 + 0.43^2)),
    Mmax ~ lognormal(6.58 + log(2), sqrt(0.18^2 + 0.64^2)),
    ## mu_MS ~ beta(), # DEFAULT
    ## Age distribution is informed by Bingham data; this is a semi-informative
    ## prior to help keep things reasonable
    ## mu_p ~ dirichlet(c(1, 1)), # DEFAULT
    mu_p ~ dirichlet(c(1, b90)),
    ## Spawners are observed precisely during trap-and-haul upstream, but
    ## allowing variance to get too close to zero results in divergent
    ## transitions. This is why the default prior has an explicitly zero-avoiding
    ## prior.
    ## , tau_S ~ gnormal(1, 0.85, 30) # DEFAULT
    tau_S ~ gnormal(0.25, 0.075, 3) # Keeps values > 0.1, but smaller than default
    ## Smolts are very approximate because they are back-calculated from Bingham
    ## SAR rates, so
    ## , tau_M ~ gnormal(1, 0.85, 30) # DEFAULT
    ## , tau_M ~ gnormal(2, 1, 2) # Wider and larger than default, might need to
    ##                            # avoid small values more
  ),
  chains = 4, iter = 10000,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 10,
                 ## "diag_e" works better than "dense_e"
                 metric = "diag_e")
)
write_rds(wyn_info_fit, file.path(data_dir, "wyn_info_fit.rds"))

wyn_info_post <- as_draws_rvars(wyn_info_fit)
write_rds(wyn_info_post, file.path(data_dir, "wyn_info_post.rds"))

### Additional fits
## Reduce age composition observations to minimal, still back-calculating smolts from overall
wyn_oac2 <- wyn_oac |>
  mutate(n_age2_obs = 0.01,
         n_age3_obs = 0.01)

## Fit with Jeffreys prior on age composition
wyn_minobs_fit <- salmonIPM(
  model = "IPM", life_cycle = "SMS", pool_pops = FALSE,
  SR_fun = "BH",
  ages = list(M = 1),
  center = FALSE, scale = FALSE,
  fish_data = wyn_oac2,
  age_F = c(1, 1),
  age_B = c(0, 0),
  ## age_S_obs = c(FALSE, TRUE),
  ## age_S_eff = c(FALSE, TRUE),
  prior = list(
    ## Priors from the Barrowman et al. (2003) meta-analysis. alpha prior
    ## includes correction assuming 45% of spawners are female and 50% of
    ## smolts are female. Mmax does the same, doubling the capacity per unit
    ## habitat to account for male smolts.
    alpha ~ lognormal(4.27 + log(0.9), sqrt(0.18^2 + 0.43^2)),
    Mmax ~ lognormal(6.58 + log(2), sqrt(0.18^2 + 0.64^2)),
    ## mu_MS ~ beta(), # DEFAULT
    ## Age distribution is informed by Bingham data; this is a semi-informative
    ## prior to help keep things reasonable
    ## mu_p ~ dirichlet(c(1, 1)), # DEFAULT
    mu_p ~ dirichlet(c(0.5, 0.5)),
    ## Spawners are observed precisely during trap-and-haul upstream, but
    ## allowing variance to get too close to zero results in divergent
    ## transitions. This is why the default prior has an explicitly zero-avoiding
    ## prior.
    ## , tau_S ~ gnormal(1, 0.85, 30) # DEFAULT
    tau_S ~ gnormal(0.25, 0.075, 3) # Keeps values > 0.1, but smaller than default
    ## Smolts are very approximate because they are back-calculated from Bingham
    ## SAR rates, so
    ## , tau_M ~ gnormal(1, 0.85, 30) # DEFAULT
    ## , tau_M ~ gnormal(2, 1, 2) # Wider and larger than default, might need to
    ##                            # avoid small values more
  ),
  chains = 4, iter = 10000,
  control = list(adapt_delta = 0.99,
                 max_treedepth = 10,
                 ## "diag_e" works better than "dense_e"
                 metric = "diag_e")
)
write_rds(wyn_minobs_fit, file.path(data_dir, "wyn_minobs_fit.rds"))

wyn_minobs_post <- as_draws_rvars(wyn_minobs_fit)
write_rds(wyn_minobs_post, file.path(data_dir, "wyn_minobs_post.rds"))


## Back-calculate smolts from adults only, again use minimal observations


## age comp from wynoochee from Megan
library(readxl)
ac <- read_xlsx("rawdata/2016-2022WynoocheeCohoAdultNumbers.xlsx") |>
  select(year = Coho,
         wyn_adult = Above,
         between = Between,
         wyn_jack = Jacks)

wyn |>
  select(year = year,
         bing_adult = age_3,
         bing_jack = age_2) |>
  right_join(ac, by = join_by(year)) |>
  transmute(year = year,
            bing_pjack = bing_jack / (bing_jack + bing_adult),
            wyn_pjack = wyn_jack / (wyn_jack + wyn_adult))
