library(tidyverse)
library(posterior)
library(ggdist)
library(distributional)

## Define a distributional-style interface for the generalized normal
## distribution in order to plot priors
source("dist_gen_normal.R")

data_dir <- "data/agecomp"
## Read in the "observations" for plotting against the posteriors
wyn_data <- read_rds(file.path(data_dir, "wyn_oac.rds"))
## Need the stan_data output to get the default prior specification for Mmax
## wyn_stan_data <- read_rds(file.path(data_dir, "wyn_stan_data.rds"))
## Read in the posteriors
wyn_info_post <- read_rds(file.path(data_dir, "wyn_info_post.rds"))
wyn_jeff_post <- read_rds(file.path(data_dir, "wyn_jeff_post.rds"))
wyn_minobs_post <- read_rds(file.path(data_dir, "wyn_minobs_post.rds"))

### Prior-posterior plots ======================================================
## Productivity parameters -----------------------------------------------------
prodpar_df <- tibble(
  parameter = c("M[max]", "alpha"),
  prior = c(dist_lognormal(6.58 + log(2), sqrt(0.18^2 + 0.64^2)),
            dist_lognormal(4.27 + log(0.9), sqrt(0.18^2 + 0.43^2))),
  Informative = c(wyn_info_post$Mmax, wyn_info_post$alpha),
  Jeffreys = c(wyn_jeff_post$Mmax, wyn_jeff_post$alpha)) |>
  pivot_longer(c(Informative, Jeffreys),
               names_to = "agecomp_prior",
               values_to = "posterior")

prodpar_df |>
  filter(parameter == "alpha") |>
  ggplot() +
  stat_slab(aes(xdist = prior), alpha = 0.5) +
  stat_slab(aes(xdist = posterior, fill = agecomp_prior),
            alpha = 0.5) +
  scale_thickness_shared() + # Put slab thickness on common scale!
  ## scale_x_log10(labels = scales::comma,
  ##               expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Smolt productivity at 100% dam survival",
       x = "Smolts per spawner at low density (alpha)",
       y = "Probability density",
       fill = "Age composition\nprior")
ggsave(file.path(data_dir, "01-alpha-post.png"))

prodpar_df |>
  filter(parameter == "M[max]") |>
  ggplot() +
  stat_slab(aes(xdist = prior), alpha = 0.5) +
  stat_slab(aes(xdist = posterior, fill = agecomp_prior),
            alpha = 0.5) +
  scale_thickness_shared() + # Put slab thickness on common scale!
  ## scale_x_log10(labels = scales::comma,
  ##               expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Smolt productivity at 100% dam survival",
       x = "Smolt capacity per km",
       y = "Probability density",
       fill = "Age composition\nprior")
ggsave(file.path(data_dir, "02-Mmax-post.png"))

## bev_holt <- function(sp, alpha, rmax, A = 20)
##   alpha * sp / (1 + alpha * sp / (A * rmax))

## bev_holt_rv <- rfun(bev_holt)

## prodpar_df <- tibble(
##   parameter = c("M[max]", "alpha"),
##   prior = c(dist_lognormal(6.58 + log(2), sqrt(0.18^2 + 0.64^2)),
##             dist_lognormal(4.27 + log(0.9), sqrt(0.18^2 + 0.43^2))),
##   Informative = c(wyn_info_post$Mmax, wyn_info_post$alpha),
##   Jeffreys = c(wyn_jeff_post$Mmax, wyn_jeff_post$alpha)) |>
## tibble(spawners = seq(0, 1.02 * max(wyn_data$S_obs, na.rm = TRUE), length.out = 256),
##        smolts = bev_holt(spawners, wyn_post$alpha, wyn_post$Mmax)) |>
##   ggplot() +
##   stat_lineribbon(aes(x = spawners, ydist = smolts),
##                   linewidth = 0.5,
##                   alpha = 0.7) +
##   geom_point(data = wyn_data,
##              aes(x = S_obs, y = M_obs)) +
##   scale_x_continuous(labels = scales::comma,
##                      expand = expansion(mult = c(0, 0.02))) +
##   scale_y_continuous(labels = scales::comma,
##                      expand = expansion(mult = c(0, 0.02))) +
##   scale_fill_brewer() +
##   guides(linewidth = "none") +
##   labs(x = "Spawners", y = "Smolts",
##        title = "Beverton-Holt spawner-smolt relationship",
##        subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
## ggsave(file.path(data_dir, "03-bevholt.png"), ...)

  ## Observation variances -------------------------------------------------------

## tau_S_prior <- dist_gen_normal(1, 0.85, 30) # DEFAULT
## tau_S_prior <- dist_gen_normal(0.25, 0.075, 3)
## tau_M_prior <- dist_gen_normal(1, 0.85, 30) # DEFAULT
## tau_M_prior <- dist_gen_normal(2, 1, 2)

## tibble(
##   parameter = c("tau[S]", "tau[M]"),
##   prior = c(tau_S_prior, tau_M_prior),
##   posterior = c(wyn_post$tau_S, wyn_post$tau_M)
## ) |>
##   ggplot() +
##   stat_slab(aes(xdist = prior), alpha = 0.5) +
##   stat_slab(aes(xdist = posterior),
##             density = "histogram",
##             color = "skyblue", fill = "skyblue", alpha = 0.7) +
##   geom_vline(xintercept = 0.1, linetype = "dashed", alpha = 0.5) +
##   facet_wrap(~ parameter, nrow = 1, scales = "free",
##              ## strip.position = "bottom",
##              labeller = label_parsed) +
##   scale_thickness_shared() + # Put slab thickness on common scale!
##   scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
##   expand_limits(x = 0) +
##   scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
##   ## coord_cartesian(xlim = c(0, NA), expand = FALSE) +
##   labs(x = "Value", y = "",
##        title = "Observation uncertainty parameters",
##        subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
## ggsave(file.path(data_dir, "04-obs-uncert-pars.png"), ...)

mod_post <- read_rds("data/ds100/wyn_post.rds")

b90 <- uniroot(\(b) qbeta(0.05, b, 1) - 0.9, interval = c(28, 29))$root

## Age distributions
agecomp_df <- tibble(
  year = wyn_data$year,
  prior = dist_beta(b90, 1),
  info_post = wyn_info_post$p[, 2, drop = TRUE],
  jeff_post = wyn_jeff_post$p[, 2, drop = TRUE],
  mod_post = mod_post$p[, 2, drop = TRUE]) |>
  pivot_longer(-year) |>
  mutate(name = factor(name, levels = c("prior", "mod_post", "jeff_post", "info_post")))

agecomp_df |>
  ggplot() +
  stat_pointinterval(aes(x = year, ydist = value, color = name), alpha = 0.5,
                     position = position_dodge()) +
  scale_thickness_shared() +
  coord_cartesian(ylim = c(0.9, 1.0)) +
  scale_y_continuous(labels = scales::percent)

  stat_slab(aes(x = year, ydist = posterior),
            density = "bounded",
            color = "skyblue", fill = "skyblue", alpha = 0.7) +
  geom_point(data = wyn_data,
             aes(x = year, y = n_age3_obs / (n_age3_obs + n_age2_obs))) +
  ## Not sharing thickness here or you can't see the prior because it is much
  ## more diffuse than the posterior
  ## scale_thickness_shared() + # Put slab thickness on common scale!
  coord_cartesian(ylim = c(0.75, 1), expand = FALSE) +
  scale_x_continuous(breaks = seq(1950, 2050, by = 5)) +
  labs(x = "Year", y = "Fraction Mature (vs. Jack)",
       title = "Spawner age distribution",
       subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
ggsave(file.path(data_dir, "05-age-dist.png"), ...)

## Smolts -----
smolt_df <- tibble(
  year = wyn_data$year,
  info_smolts = wyn_info_post$M,
  jeff_smolts = wyn_jeff_post$M,
  mod_smolts = mod_post$M
) |>
  pivot_longer(-year,
               names_to = "prior_type",
               values_to = "smolts")

ggplot(smolt_df,
       aes(x = year,
           ydist = smolts)) +
  stat_pointinterval(aes(fill = prior_type, color = prior_type),
                     position = position_dodge()) +
  geom_point(data = wyn_data,
             aes(x = year, y = M_obs),
             inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(1950, 2050, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), labels = scales::comma) +
  expand_limits(y = 0) +
  labs(x = "Year", y = "Smolts",
       title = "Smolt production")
## ggsave(file.path(data_dir, "06-smolts.png"), ...)

## Spawners -----
spawn_df <- tibble(
  year = wyn_data$year,
  info_smolts = wyn_info_post$S,
  jeff_smolts = wyn_jeff_post$S,
  mod_smolts = mod_post$S
) |>
  pivot_longer(-year,
               names_to = "prior_type",
               values_to = "spawners")

ggplot(spawn_df,
       aes(x = year,
             ydist = spawners)) +
  stat_pointinterval(aes(fill = prior_type, color = prior_type),
                     position = position_dodge()) +
  geom_point(data = wyn_data,
             aes(x = year, y = S_obs),
             inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(1950, 2050, by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), labels = scales::comma) +
  expand_limits(y = 0) +
  labs(x = "Year", y = "Spawners",
       title = "Spawners")
## ggsave(file.path(data_dir, "06-smolts.png"), ...)

wyn_data2 <- wyn_data |>
  select(year, smolt_obs = M_obs, spawn_obs = S_obs) |>


left_join(smolt_df, spawn_df, by = join_by(year, prior_type)) |>
  left_join(wyn_data2, join_by(year)) |>
  ggplot(aes(color = prior_type, fill = prior_type)) +
  stat_pointinterval(aes(xdist = spawners, y = median(smolts)), alpha = 0.5) +
  stat_pointinterval(aes(x = median(spawners), ydist = smolts), alpha = 0.5) +
  geom_point(aes(x = spawn_obs, y = smolt_obs), color = "black") +
  geom_segment(aes(x = spawn_obs, xend = median(spawners), y = smolt_obs, yend = median(smolts))) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(ylim = c(3e3, 62e3))

left_join(smolt_df, spawn_df, by = join_by(year, prior_type)) |>
  left_join(wyn_data2, join_by(year)) |>
  ggplot(aes(x = smolt_obs, ydist = smolts, color = prior_type)) +
  stat_pointinterval(position = position_dodge(width = 1000), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

left_join(smolt_df, spawn_df, by = join_by(year, prior_type)) |>
  left_join(wyn_data2, join_by(year)) |>
  filter(prior_type %in% c("jeff_smolts", "mod_smolts")) |>
  select(year, prior_type, smolts) |>
  pivot_wider(names_from = prior_type, values_from = smolts) |>
  ggplot() +
  stat_pointinterval(aes(x = median(jeff_smolts), ydist = mod_smolts), alpha = 0.5) +
  stat_pointinterval(aes(xdist = jeff_smolts, y = median(mod_smolts)), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(aes(x = median(jeff_smolts), y = median(mod_smolts)), se = FALSE, method = lm) +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Smolts: Jeffreys prior",
       y = "Smolts: Bingham age comp")

left_join(smolt_df, spawn_df, by = join_by(year, prior_type)) |>
  left_join(wyn_data2, join_by(year)) |>
  filter(prior_type %in% c("jeff_smolts", "mod_smolts")) |>
  select(year, prior_type, spawners) |>
  pivot_wider(names_from = prior_type, values_from = spawners) |>
  ggplot() +
  stat_pointinterval(aes(y = median(jeff_smolts), xdist = mod_smolts), alpha = 0.5) +
  stat_pointinterval(aes(ydist = jeff_smolts, x = median(mod_smolts)), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(aes(y = median(jeff_smolts), x = median(mod_smolts)), se = FALSE, method = lm) +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Spawners: Bingham age comp",
       y = "Spawners: Jeffreys prior")

## Alternative age comp formulations
tibble(year = wyn_oac$year,
       p = wyn_minobs_post$p[, 2, drop = TRUE]) |>
  ggplot(aes(xdist = p, y = year, fill = factor(year), color = factor(year))) +
  stat_halfeye() +
  scale_thickness_shared()

df2 <- tibble(year = wyn_oac$year,
       smolt_obs = wyn_oac$M_obs,
       smolt_minobs = wyn_minobs_post$M,
       smolt_mod = mod_post$M,
       spawn_obs = wyn_oac$S_obs,
       spawn_minobs = wyn_minobs_post$S,
       spawn_mod  = mod_post$S
       )


df2 |>
  select(year, starts_with("smolt")) |>
  pivot_longer(c(smolt_minobs, smolt_mod)) |>
  ggplot(aes(x = smolt_obs, ydist = value, color = name)) +
  stat_pointinterval(position = position_dodge(width = 0.01)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10()

df2 |>
  select(year, starts_with("spawn")) |>
  pivot_longer(c(spawn_minobs, spawn_mod)) |>
  ggplot(aes(x = spawn_obs, ydist = value, color = name)) +
  stat_pointinterval(position = position_dodge(width = 0.01)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10()

prod_df <- tibble(
mod = c("mod", "minobs"),
alpha = c(mod_post$alpha, wyn_minobs_post$alpha),
Mmax = c(mod_post$Mmax, wyn_minobs_post$Mmax))

ggplot(prod_df, aes(xdist = alpha, fill = mod)) +
  stat_halfeye(alpha = 0.4) +
  scale_thickness_shared()

ggplot(prod_df, aes(xdist = Mmax, fill = mod)) +
  stat_halfeye(alpha = 0.4) +
  scale_thickness_shared()
