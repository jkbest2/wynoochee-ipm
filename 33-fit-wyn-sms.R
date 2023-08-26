library(salmonIPM)
options(mc.cores = parallel::detectCores(logical = FALSE))
rstan_options(auto_write = TRUE)
library(tidyverse)
library(posterior)

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
  filter(name == "Upper Wynoochee") |>
  pluck("habitat") |>
  drop_units()

## These are the values of dam survival that we are using to inflate the smolt
## counts.
ds <- seq(0.2, 1, by = 0.2)

## Currently need at least two age classes, so using the information from
## Bingham on age class observations to inform. We don't want *all* the
## information, so we're reducing the sample size by a factor of
## `age_obs_factor`. This means that if in a given year 10 jacks and 90 mature
## fish were observed and `age_obs_factor` is 10, then we will tell the model
## that 1 jack and 9 mature fish were observed at Wynoochee. This is also
## informed by the prior set on the age distribution, which is currently mildly
## informative as a Dirichlet(2, 18).
age_obs_factor <- 10

## Generate the data frames assuming different levels of smolt survival during
## downstream dam passage
dfs <- lapply(ds, \(dam_survival) {
  wyn <- read_rds("data/wyn_trap.rds") |>
    filter(species == "coho") |>
    left_join(bing_harvest, by = join_by(year)) |>
    filter(!is.na(harvest_total)) |>
    mutate(pop = "Upper Wynoochee",
           A = wyn_hab,
           S_obs = count,
           n_W_obs = 0,
           n_H_obs = 0,
           fit_p_HOS = FALSE,
           B_take_obs = 0) |>
    left_join(select(bing, year, age_2, age_3), by = join_by(year)) |>
    ## Find age structure
    mutate(n_age2_obs = age_2 / age_obs_factor,
           n_age2_obs = replace_na(n_age2_obs, 0),
           n_age3_obs = age_3 / age_obs_factor,
           n_age3_obs = replace_na(n_age3_obs, 0),
           frac_jack = age_2 / (age_2 + age_3),
           ## Fill in the fraction of jacks using the overall mean. This allows
           ## us to fill in many more years of data.
           frac_jack = replace_na(
             frac_jack,
             sum(age_2, na.rm = TRUE) / sum(age_2 + age_3, na.rm = TRUE))) |>
    ## Back-calculate the number of smolts we would have seen (if we were
    ## looking). This treats the overall smolt-to-adult survival as known, based
    ## on the Bingham data. It also assumes that jacks share the same mortality
    ## rate *as the mature fish they return with*.
    mutate(
      jack_smolts = lead(S_obs, 2) * lead(frac_jack, 2) /
        lead(survival_escapement, 3), # Best assumption of ocean survival is that it tracks the brood year adult ocean survival
      mature_smolts = lead(S_obs, 3) * (1 - lead(frac_jack, 3)) /
        lead(survival_escapement, 3),
      M_obs = (jack_smolts + mature_smolts) / dam_survival) |>
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
           B_take_obs)

  ## Check that fishing rates are between 0 and 1.
  ## FIXME Are the more checks that can be done here?
  stopifnot(all(wyn$F_rate >= 0,
                wyn$F_rate < 1))

  wyn
})
write_rds(dfs, "data/wyn_dfs.rds")

walk(seq_along(ds), \(idx) {
  ## Create the directory to save results if necessary
  data_dir <- paste0("data/ds", format(100 * ds[idx], trim = TRUE, digits = 2))
  if (!dir.exists(data_dir)) dir.create(data_dir)

  ## Extract the data set and double check that there are no missing
  ## observations that will crash R
  wyn_data <- dfs[[idx]] |>
    filter(!is.na(M_obs),
           !is.na(F_rate),
           !is.na(S_obs))
  write_rds(wyn_data, file.path(data_dir, "wyn_data.rds"))

  ## Use the `salmonIPM::stan_data` function to generate data. This isn't used,
  ## but we need to save it in order to extract the default prior on Mmax, which
  ## is derived from the data. There may be a better way to do this, either by
  ## extracting these values from the `stanfit` object at the end, or by
  ## treating this as the source of "truth" and extracting the necessary values
  ## below. This would prevent mismatched prior specifications in the final
  ## figures.
  wyn_stan_data <- stan_data(
    stan_model = "IPM_SS_np",
    SR_fun = "BH",
    ages = list(M = 1),
    center = FALSE, scale = FALSE,
    fish_data = wyn_data,
    age_F = c(1, 1),
    age_B = c(0, 0)
  )
  write_rds(wyn_stan_data, file.path(data_dir, "wyn_stan_data.rds"))

  wyn_fit <- salmonIPM(
    model = "IPM", life_cycle = "SMS", pool_pops = FALSE,
    SR_fun = "BH",
    ages = list(M = 1),
    center = FALSE, scale = FALSE,
    fish_data = wyn_data,
    age_F = c(1, 1),
    age_B = c(0, 0),
    ## age_S_obs = c(FALSE, TRUE),
    ## age_S_eff = c(FALSE, TRUE),
    prior = list(
      ## alpha ~ lognormal(2, 5), # DEFAULT
      ## Mmax ~ lognormal(wyn_stan_data$prior_Rmax[1], wyn_stan_data$prior_Rmax[2]) # DEFAULT (derived from data)

      ## Priors from the Barrowman et al. (2003) meta-analysis. alpha prior
      ## includes correction assuming 45% of spawners are female and 50% of
      ## smolts are female. Mmax does the same, doubling the capacity per unit
      ## habitat to account for male smolts.
      alpha ~ lognormal(4.27 + log(0.9), sqrt(0.18^2 + 0.43^2)),
      Mmax ~ lognormal(6.58 + log(2), sqrt(0.18^2 + 0.64^2)),
      ## , mu_MS ~ beta() # DEFAULT
      ## Age distribution is informed by Bingham data; this is a semi-informative
      ## prior to help keep things reasonable
      ## , mu_p ~ dirichlet(c(1, 1)) # DEFAULT
      mu_p ~ dirichlet(c(2, 18)),
      ## Spawners are observed precisely during trap-and-haul upstream, but
      ## allowing variance to get too close to zero results in divergent
      ## transitions. This is why the default prior has an explicitly zero-avoiding
      ## prior.
      ## , tau_S ~ gnormal(1, 0.85, 30) # DEFAULT
      tau_S ~ gnormal(0.25, 0.075, 3) # Keeps values < 0.1, but lower than default
      ## Smolts are very approximate because they are back-calculated from Bingham
      ## SAR rates, so
      ## , tau_M ~ gnormal(1, 0.85, 30) # DEFAULT
      ## , tau_M ~ gnormal(2, 1, 2) # Wider and larger than default, might need to avoid small values more
    ),
    chains = 4, iter = 10000,
    control = list(adapt_delta = 0.99,
                   metric = "dense_e")
  )
  attr(wyn_fit, "dam_survival") <- ds[idx]

  write_rds(wyn_fit, file.path(data_dir, "wyn_fit.rds"))

  wyn_post <- as_draws_rvars(wyn_fit)
  attr(wyn_post, "dam_survival") <- ds[idx]
  write_rds(wyn_post, file.path(data_dir, "wyn_post.rds"))
})
