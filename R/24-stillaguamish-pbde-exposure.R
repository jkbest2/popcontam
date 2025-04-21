library(tidyverse)
library(posterior)

source("R/fit-lnorm-pop-exposure.R")

## This is the data frame containing only samples from the Puyallup estuary and
## nearshore. Currently it is processed and saved in R/exposure-map.R
contam <- read_rds(
  here::here("data", "stillaguamish", "pbde_exposure_obs.rds")
)

data <- lnorm_data(contam, bde_sum)

## Increased the numer of iterations in order to ensure Rhat < 1.01
fit <- fit_lnorm_exposure(data, iter = 4000)
post <- as_draws_rvars(fit)

rhat_check <- tibble(
  mod = "ww",
  par = c("pop_meanlog", "pop_sdlog")
) |>
  mutate(
    post = map_vec(par, \(p) pluck(post, p)),
    rhat = rhat(post),
    rhat_good = abs(rhat - 1) < 1.01
  )

if (any(!rhat_check$rhat_good)) {
  warning(
    "Rhat indicates MCMC may not have converged, consider more iterations"
  )
}

# par_check <- tibble(
#   mod = c("ww", "lw", "lw1"),
#   par = "pop_sdlog"
# ) |>
#   mutate(
#     post = map2_vec(mod, par, \(m, p) pluck(post, m, p)),
#     sdl_high = sum(post > 2.5),
#     sdl_max = max(post)
#   )

# if (any(par_check$sdl_high > 0)) {
#   warning("pop_sdlog > 2.5 in ", sum(par_check$sdl_high), "MCMC samples")
# }

write_rds(
  post,
  here::here("data", "stillaguamish", "pbde_exposure.rds")
)
