library(tidyverse)
library(RTMB)

dam_survival <- 1
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

wyn_comp <- wyn |>
  filter(!is.na(M_obs),
         !is.na(S_obs))

par_init <- list(
  log_m_hat = rep(0, nrow(wyn)),
  log_rec_dev = rep(0, nrow(wyn)),
  ocean_survival = rep(0.02, nrow(wyn))
)

## Density namespace in R; taken from `ar1_4D` example in RTMB
N01 <- function(x) sum(dnorm(x, log=TRUE))
AR1 <- function(phi, DMARG = N01) {
    function(x) {
        x <- as.array(x)
        d <- dim(x); ncolx <- tail(d, 1)
        xcol <- function(j) {
            dim(x) <- c(length(x) / ncolx, ncolx)
            ans <- x[,j]
            newd <- head(d, -1)
            if (length(newd)==0) newd <- NULL
            dim(ans) <- newd
            ans
        }
        ans <- DMARG(xcol(1))
        sd <- sqrt(1-phi*phi)
        for (i in 2:ncolx) {
            ans <- ans + DMARG( ( xcol(i) - phi * xcol(i-1) ) / sd )
        }
        ans <- ans - (length(x) / ncolx) * (ncolx - 1) *  log(sd)
        ans
    }
}

## Define the stock-recruit relationship. Here sp is number of spawners, alpha
## is the intrinsic productivity, and rmax is the maximum recruitment.
bev_holt <- function(sp, alpha, rmax)
  alpha * sp / (1 + alpha * sp / rmax)

bh_lik <- function(pars, data) {
  ## m_hat <- SR(SR_fun = "BH",
  ##             alpha = exp(pars[1]),
  ##             Rmax = exp(pars[2]),
  ##             S = data$S_obs,
  ##             A = 1,
  ##             R_per_S = FALSE)
  m_hat <- bev_holt(data$S_obs,
                    exp(pars[1]),
                    exp(pars[2]))
  -sum(dnorm(log(data$M_obs), log(m_hat), exp(pars[3]), log = TRUE))
}

wyn_nlp <- function(pars) {
  getAll(pars, wyn_comp, warn = FALSE)

  alpha <- exp(pars$log_alpha)
  rmax <- exp(pars$log_rmax)
  sdl <- exp(pars$log_sdlog)

  OBS(M_obs)

  m_hat <- bev_holt(wyn_comp$S_obs, alpha, rmax)

  ## -sum(dnorm(log(M_obs), log(m_hat), sdl, log = TRUE))
  log(M_obs) %~% dnorm(log(m_hat), sdl)
}

par_init <- list(log_alpha = log(70),
              log_rmax = log(28e3),
              log_sdlog = log(0.8))

wyn_nlp(par_init)

optim(unlist(par_init), \(p) wyn_nlp(list(log_alpha = p[1], log_rmax = p[2], log_sdlog = p[3])))

obj <- MakeADFun(wyn_nlp, par_init)

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt2 <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")
