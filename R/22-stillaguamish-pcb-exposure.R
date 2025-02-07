library(tidyverse)
library(posterior)

source("R/fit-lnorm-pop-exposure.R")

## This is the data frame containing only samples from the Puyallup estuary and
## nearshore. Currently it is processed and saved in R/exposure-map.R
contam <- read_rds("data/stillaguamish/pcb_estns.rds")

data_ww <- lnorm_data(contam, pcb_ug_ww)
data_lw <- lnorm_data(contam, pcb_ug_lw)
data_lw1 <- lnorm_data(contam, pcb_ug_lw1)

fit_ww <- fit_lnorm_exposure(data_ww, iter = 4000)
post_ww <- as_draws_rvars(fit_ww)

fit_lw <- fit_lnorm_exposure(data_lw, iter = 4000)
post_lw <- as_draws_rvars(fit_lw)

## R hats were a bit high for some parameters with 2000 iterations
fit_lw1 <- fit_lnorm_exposure(data_lw1, iter = 4000)
post_lw1 <- as_draws_rvars(fit_lw1)

post <- list(ww = post_ww, lw = post_lw, lw1 = post_lw1)

rhat_check <- expand_grid(
  mod = c("ww", "lw", "lw1"),
  par = c("pop_meanlog", "pop_sdlog")
) |>
  mutate(
    post = map2_vec(mod, par, \(m, p) pluck(post, m, p)),
    rhat = rhat(post),
    rhat_good = abs(rhat - 1) < 1.01
  )

if (any(!rhat_check$rhat_good)) {
  warning(
    "Rhat indicates MCMC may not have converged, consider more iterations"
  )
}

par_check <- tibble(
  mod = c("ww", "lw", "lw1"),
  par = "pop_sdlog"
) |>
  mutate(
    post = map2_vec(mod, par, \(m, p) pluck(post, m, p)),
    sdl_high = sum(post > 2.5),
    sdl_max = max(post)
  )

if (any(par_check$sdl_high > 0)) {
  warning("pop_sdlog > 2.5 in ", sum(par_check$sdl_high), "MCMC samples")
}

write_rds(
  post,
  "data/stillaguamish/pcb_exposure.rds"
)
