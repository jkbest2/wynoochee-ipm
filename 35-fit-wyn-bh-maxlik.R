library(salmonIPM)
## options(mc.cores = parallel::detectCores(logical = FALSE))
## rstan_options(auto_write = TRUE)
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
         survival_escapement, # Overall smolt -> adult survival
         harvest_total)

## bing_osurv <- bing_harvest |>
##   mutate(smolt_year = year - 3) |>
##   select(smolt_year, ocean_survival)

bing_harvest2 <- read_rds("data/bing_harvest.rds")

## bing_df <- bing |>
##   left_join(bing_harvest, by = join_by(year)) |>
##   mutate(pop = factor("Bingham Creek"),
##          year = as.integer(as.character(year)),
##          A = 40,
##          n_W_obs = 0,
##          n_H_obs = 0,
##          fit_p_HOS = FALSE,
##          B_take_obs = 0,
##          ) |>
##   select(pop,
##          year,
##          A,
##          S_obs = total,
##          n_age2_obs = age_2,
##          n_age3_obs = age_3,
##          n_W_obs,
##          n_H_obs,
##          fit_p_HOS,
##          F_rate = harvest_total,
##          B_take_obs
## )

ds <- seq(0.2, 1, by = 0.2)
age_obs_factor <- 10

dfs <- lapply(ds, \(dam_survival) {

wyn <- read_rds("data/wyn_trap.rds") |>
  filter(species == "coho") |>
  left_join(bing_harvest, by = join_by(year)) |>
  filter(!is.na(harvest_total)) |>
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
  ## Find age structure
  mutate(n_age2_obs = age_2 / age_obs_factor,
         n_age2_obs = replace_na(n_age2_obs, 0),
         n_age3_obs = age_3 / age_obs_factor,
         n_age3_obs = replace_na(n_age3_obs, 0),
         frac_jack = age_2 / (age_2 + age_3),
         frac_jack = replace_na(frac_jack,
                                sum(age_2, na.rm = TRUE) / sum(age_2 + age_3, na.rm = TRUE))) |>
  ## Back-calculate the number of smolts we would have seen if we were looking
  mutate(jack_smolts = lead(S_obs, 2) * lead(frac_jack, 2) / lead(survival_escapement, 2),
         mature_smolts = lead(S_obs, 3) * (1 - lead(frac_jack, 3)) / lead(survival_escapement, 3),
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

## stopifnot(all(wyn$F_rate >= 0,
##               wyn$F_rate < 1))
})

bev_holt <- function(sp, alpha, rmax)
  alpha * sp / (1 + alpha * sp / rmax)

bh_lik <- function(pars, data) {
  m_hat <- bev_holt(data$S_obs, exp(pars[1]), exp(pars[2]))
  -sum(dlnorm(data$M_obs, log(m_hat), sdlog = exp(pars[3]), log = TRUE))
}

opts <- lapply(dfs, \(wyn) {
wyn_comp <- wyn |>
  select(year, M_obs, S_obs) |>
  filter(!is.na(M_obs),
         !is.na(S_obs))

pars <- c(log_alpha = log(200),
          log_Rmax = log(25000 / dam_survival),
          log_sd = log(1))
bh_lik(pars = pars, data = wyn_comp)
opt <- optim(pars, bh_lik, data = wyn_comp)
opt
})

plts <- lapply(seq_along(ds), \(idx) {

  opt <- opts[[idx]]
  wyn_comp <- dfs[[idx]] |>
    select(year, M_obs, S_obs) |>
    filter(!is.na(M_obs),
           !is.na(S_obs))

  bh_df <- tibble(S = seq(0, 5500, 25)) |>
    mutate(M_hat = SR(alpha = exp(opt$par[1]), Rmax = exp(opt$par[2]), S = S),
           M10 = qlnorm(0.1, log(M_hat), exp(opt$par[3])),
           M90 = qlnorm(0.9, log(M_hat), exp(opt$par[3])))

  ggplot() +
    geom_ribbon(data = bh_df,
                aes(x = S, ymin = M10, ymax = M90),
                alpha = 0.3) +
    geom_line(data = bh_df,
              aes(x = S, y = M_hat)) +
    geom_point(data = wyn_comp,
               aes(x = S_obs, y = M_obs, color = year)) +
    geom_path(data = wyn_comp,
              aes(x = S_obs, y = M_obs, color = year),
              alpha = 0.3) +
    geom_label(data = filter(wyn_comp, year %% 5 == 0),
               aes(x = S_obs, y = M_obs, label = year),
               nudge_x = 100) +
    coord_cartesian(xlim = c(0, NA),
                    ## Ensure that the label and point aren't *right* at the edge
                    ylim = c(0, 1.02 * max(wyn_comp$M_obs)),
                    expand = FALSE) +
    labs(x = "Spawners", y = "Smolts", color = "Year",
         title = paste0(ds[idx] * 100, "% Dam survival "))
})

pars <- lapply(seq_along(ds), \(idx) {
  opt <- opts[[idx]]
  p <- c(dam_survival = ds[idx],
         alpha = exp(opt$par[1]),
         Rmax = exp(opt$par[2]),
         sdlog = exp(opt$par[3]))
  names(p) <- c("dam_survival", "alpha", "Rmax", "sdlog")
  p
})
