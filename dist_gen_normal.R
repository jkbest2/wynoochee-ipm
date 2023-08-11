library(distributional)
library(vctrs)
library(rlang)

dist_gen_normal <- function(mu = 0, alpha = 1, beta = 2){
  mu <- vec_cast(mu, double())
  alpha <- vec_cast(alpha, double())
  beta <- vec_cast(beta, double())

  if (any(alpha[!is.na(alpha)] < 0))
    abort("alpha must be nonnegative")
  if (any(beta[!is.na(beta)] < 0))
    abort("beta must be nonnegative")

  new_dist(mu = mu, alpha = alpha, beta = beta, class = "dist_gen_normal")
}

format.dist_gen_normal <- function(x, digits = 2, ...) {
  sprintf(
    "gN(%s, %s, %s)",
    format(x[["mu"]], digits = digits, ...),
    format(x[["alpha"]], digits = digits, ...),
    format(x[["beta"]], digits = digits, ...)
  )
}

density.dist_gen_normal <- function(x, at, ...) {
  gnorm::dgnorm(at, x[["mu"]], x[["alpha"]], x[["beta"]])
}
#' @export
log_density.dist_gen_normal <- function(x, at, ...){
  gnorm::dgnorm(at, x[["mu"]], x[["alpha"]], log = TRUE)
}

#' @export
quantile.dist_gen_normal <- function(x, p, ...){
  gnorm::qgnorm(p, x[["mu"]], x[["alpha"]], x[["beta"]])
}
#' @export
log_quantile.dist_gen_normal <- function(x, p, ...){
  gnorm::qgnorm(p, x[["mu"]], x[["alpha"]], x[["beta"]], log.p = TRUE)
}

#' @export
cdf.dist_gen_normal <- function(x, q, ...){
  gnorm::pgnorm(q, x[["mu"]], x[["alpha"]], x[["beta"]])
}
#' @export
log_cdf.dist_gen_normal <- function(x, q, ...){
  gnorm::pgnorm(q, x[["mu"]], x[["alpha"]], x[["beta"]], log.p = TRUE)
}

#' @export
generate.dist_gen_normal <- function(x, times, ...){
  gnorm::rgnorm(times, x[["mu"]], x[["alpha"]], x[["beta"]])
}

#' @export
mean.dist_gen_normal <- function(x, ...){
  x[["mu"]]
}

#' @export
covariance.dist_gen_normal <- function(x, ...){
  x[["alpha"]] ^ 3 * gamma(3 / x[["beta"]]) /
    gamma(1 / x[["beta"]])
}

#' @export
skewness.dist_gen_normal <- function(x, ...) {
  0
}

#' @export
kurtosis.dist_gen_normal <- function(x, ...) {
  beta <- x[["beta"]]
  gamma(5 / beta) * gamma(1 / beta) /
    gamma(3 / beta) ^ 2 - 3
}

# make a normal distribution from a gen_normal distribution using the
# specified base
## normal_dist_with_base <- function(x, base = exp(1)) {
##   vec_data(dist_normal(x[["mu"]], x[["alpha"]]) / log(base))[[1]]
## }

## #' @method Math dist_gen_normal
## #' @export
## Math.dist_gen_normal <- function(x, ...) {
##   switch(.Generic,
##     # Shortcuts to get Normal distribution from log-normal.
##     log = normal_dist_with_base(x, ...),
##     log2 = normal_dist_with_base(x, 2),
##     log10 = normal_dist_with_base(x, 10),

##     NextMethod()
##   )
## }
