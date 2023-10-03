library(tidyverse)
library(posterior)
library(ggdist)
library(distributional)

## Define a distributional-style interface for the generalized normal
## distribution in order to plot priors
source("dist_gen_normal.R")

plot_posteriors <- function(ds, ...) {
  data_dir <- paste0("data/ds", format(100 * ds, trim = TRUE, digits = 2))
  ## Read in the "observations" for plotting against the posteriors
  wyn_data <- read_rds(file.path(data_dir, "wyn_data.rds"))
  ## Need the stan_data output to get the default prior specification for Mmax
  wyn_stan_data <- read_rds(file.path(data_dir, "wyn_stan_data.rds"))
  ## Read in the posterior
  wyn_post <- read_rds(file.path(data_dir, "wyn_post.rds"))

### Prior-posterior plots ======================================================
  ## Productivity parameters -----------------------------------------------------
  tibble(
    parameter = c("M[max]", "alpha"),
    prior = c(dist_lognormal(6.58 + log(2), sqrt(0.18^2 + 0.64^2)),
              dist_lognormal(4.27 + log(0.9), sqrt(0.18^2 + 0.43^2))),
    posterior = c(wyn_post$Mmax, wyn_post$alpha)
  ) |>
    ggplot() +
    stat_slab(aes(xdist = prior), alpha = 0.5) +
    stat_slab(aes(xdist = posterior, group = parameter),
              density = "histogram",
              color = "skyblue", fill = "skyblue", alpha = 0.7) +
    facet_wrap(~ parameter, nrow = 1, scales = "free",
               labeller = label_parsed) +
    scale_thickness_shared() + # Put slab thickness on common scale!
    scale_x_log10(labels = scales::comma,
                  expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Value", y = "",
         title = "Marginal smolt production parameter posteriors",
         subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
  ggsave(file.path(data_dir, "01-prod-pars.png"), ...)

  tibble(
    Mmax = draws_of(wyn_post$Mmax),
    alpha = draws_of(wyn_post$alpha)
  ) |>
    ggplot(aes(x = alpha, y = Mmax)) +
    geom_density_2d_filled() +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    guides(fill = "none") +
    labs(x = expression(alpha), y = expression(M[max]),
         title = "Joint smolt production parameter posterior",
         subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
  ggsave(file.path(data_dir, "02-joint-prod-pars.png"), ...)

  bev_holt <- function(sp, alpha, rmax, A = 20)
    alpha * sp / (1 + alpha * sp / (A * rmax))

  bev_holt_rv <- rfun(bev_holt)

  tibble(spawners = seq(0, 1.02 * max(wyn_data$S_obs, na.rm = TRUE), length.out = 256),
         smolts = bev_holt(spawners, wyn_post$alpha, wyn_post$Mmax)) |>
    ggplot() +
    stat_lineribbon(aes(x = spawners, ydist = smolts),
                    linewidth = 0.5,
                    alpha = 0.7) +
    geom_point(data = wyn_data,
               aes(x = S_obs, y = M_obs)) +
    scale_x_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_brewer() +
    guides(linewidth = "none") +
    labs(x = "Spawners", y = "Smolts",
         title = "Beverton-Holt spawner-smolt relationship",
         subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
  ggsave(file.path(data_dir, "03-bevholt.png"), ...)

  ## Observation variances -------------------------------------------------------

  ## tau_S_prior <- dist_gen_normal(1, 0.85, 30) # DEFAULT
  tau_S_prior <- dist_gen_normal(0.25, 0.075, 3)
  tau_M_prior <- dist_gen_normal(1, 0.85, 30) # DEFAULT
  ## tau_M_prior <- dist_gen_normal(2, 1, 2)

  tibble(
    parameter = c("tau[S]", "tau[M]"),
    prior = c(tau_S_prior, tau_M_prior),
    posterior = c(wyn_post$tau_S, wyn_post$tau_M)
  ) |>
    ggplot() +
    stat_slab(aes(xdist = prior), alpha = 0.5) +
    stat_slab(aes(xdist = posterior),
              density = "histogram",
              color = "skyblue", fill = "skyblue", alpha = 0.7) +
    geom_vline(xintercept = 0.1, linetype = "dashed", alpha = 0.5) +
    facet_wrap(~ parameter, nrow = 1, scales = "free",
               ## strip.position = "bottom",
               labeller = label_parsed) +
    scale_thickness_shared() + # Put slab thickness on common scale!
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    expand_limits(x = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    ## coord_cartesian(xlim = c(0, NA), expand = FALSE) +
    labs(x = "Value", y = "",
         title = "Observation uncertainty parameters",
         subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
  ggsave(file.path(data_dir, "04-obs-uncert-pars.png"), ...)

  ## Setting the age composition prior using an informative Beta distribution.
  ## Parameterized such that Pr(prop_jacks < 0.1) == 0.95
  b90 <- uniroot(\(b) pbeta(0.1, 2, b) - 0.95, c(40, 50))
  age_comp_prior <- c(b90$root, 2)

  ## Age distributions
  tibble(
    year = wyn_data$year,
    prior = dist_beta(age_comp_prior[1] , age_comp_prior[2]),
    posterior = wyn_post$p[, 2, drop = TRUE]
  ) |>
    ggplot() +
    stat_slab(aes(x = year, ydist = prior), alpha = 0.5) +
    stat_slab(aes(x = year, ydist = posterior),
              density = "bounded",
              color = "skyblue", fill = "skyblue", alpha = 0.7) +
    geom_point(data = wyn_data,
               aes(x = year, y = n_age3_obs / (n_age3_obs + n_age2_obs))) +
    ## Not sharing thickness here or you can't see the prior because it is much
    ## more diffuse than the posterior
    ## scale_thickness_shared() + # Put slab thickness on common scale!
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(0.75, 1), expand = FALSE) +
    scale_x_continuous(breaks = seq(1950, 2050, by = 5)) +
    labs(x = "Year", y = "Fraction Mature (vs. Jack)",
         title = "Spawner age distribution",
         subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
  ggsave(file.path(data_dir, "05-age-dist.png"), ...)

  ## Smolts -----
  tibble(
    year = wyn_data$year,
    smolts = wyn_post$M
  ) |>
    ggplot(aes(x = year,
               ydist = smolts)) +
    stat_pointinterval(fill = "skyblue", color = "skyblue") +
    geom_point(data = wyn_data,
               aes(x = year, y = M_obs),
               inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(1950, 2050, by = 5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), labels = scales::comma) +
    expand_limits(y = 0) +
    labs(x = "Year", y = "Smolts",
         title = "Smolt production",
         subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
  ggsave(file.path(data_dir, "06-smolts.png"), ...)

  ## Spawners -----
  tibble(
    year = wyn_data$year,
    spawners = wyn_post$S
  ) |>
    ggplot(aes(x = year,
               ydist = spawners)) +
    stat_pointinterval(fill = "skyblue", color = "skyblue") +
    geom_point(data = wyn_data,
               aes(x = year, y = S_obs),
               inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(1950, 2050, by = 5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), labels = scales::comma) +
    expand_limits(y = 0) +
    labs(x = "Year", y = "Spawners",
         title = "Returning spawners",
         subtitle = paste0(ds * 100, "% Downstream Dam Passage Survival"))
  ggsave(file.path(data_dir, "07-spawners.png"), ...)

  invisible(data_dir)
}

walk(seq(0.2, 1, 0.2), plot_posteriors, width = 10)
