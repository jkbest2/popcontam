library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/pcb-effects.R")
source("R/puyallup.R")
source("R/fit-lnorm-pop-exposure.R")

## This is the data frame containing only samples from the Puyallup estuary and
## nearshore. Currently it is processed and saved in R/exposure-map.R
contam <- read_rds("data/puy-pcb-estns.rds")

if (TRUE) {
  data_ww <- lnorm_data(contam, pcb_ug_ww)
  data_lw <- lnorm_data(contam, pcb_ug_lw)
  data_lw1 <- lnorm_data(contam, pcb_ug_lw1)

  fit_ww <- fit_lnorm_exposure(data_ww)
  post_ww <- as_draws_rvars(fit_ww)
  write_rds(post_ww, "data/pcb_ww_post.rds")

  fit_lw <- fit_lnorm_exposure(data_lw)
  post_lw <- as_draws_rvars(fit_lw)
  write_rds(post_lw, "data/pcb_lw_post.rds")

  fit_lw1 <- fit_lnorm_exposure(data_lw1)
  post_lw1 <- as_draws_rvars(fit_lw1)
  write_rds(post_lw1, "data/pcb_lw1_post.rds")
} else {
  post_ww <- read_rds("data/pcb_ww_post.rds")
  post_lw <- read_rds("data/pcb_lw_post.rds")
  post_lw1 <- read_rds("data/pcb_lw1_post.rds")
}

## Get baseline population distribution
puy0 <- eq_pop(
  puyallup_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5)
)
puy0_oa <- get_oceanadults(puy0)

pcb_effect <- function(exp_meanlog, exp_sdlog, wt_type) {
  exp_surv <- function(pcb) {
    dlnorm(pcb, exp_meanlog, exp_sdlog) *
      (1 - combo_mort(pcb, wt_type = wt_type))
  }
  integrate(
    exp_surv,
    lower = 0, upper = Inf,
    ## The default rel.tol occasionally gives incorrect results, including one
    ## set of population parameters that consistently shows an increase in
    ## survival due to PCB exposure!
    rel.tol = .Machine$double.eps^0.5
  )$value
}
rv_pcb_effect <- rfun(pcb_effect)

eff_ww <- rv_pcb_effect(post_ww$pop_meanlog, post_ww$pop_sdlog, "ww")
eff_lw <- rv_pcb_effect(post_lw$pop_meanlog, post_lw$pop_sdlog, "lw")
eff_lw1 <- rv_pcb_effect(post_lw1$pop_meanlog, post_lw1$pop_sdlog, "lw")

puyallup_exposed <- function(ns_surv, stage = NULL) {
  eqp <- eq_pop(puyallup_sim, nearshore_surv_adj = ns_surv, pop0 = rep(1000, 5))
  ## Don't always want the population age structure
  if (!is.null(stage)) {
    eqp <- attr(eqp, stage)
  }
  eqp
}
rv_puyallup_exposed <- rfun(puyallup_exposed)

puy_oa_ww <- rv_puyallup_exposed(eff_ww, "oceanadults")
puy_oa_lw <- rv_puyallup_exposed(eff_lw, "oceanadults")
puy_oa_lw1 <- rv_puyallup_exposed(eff_lw1, "oceanadults")

## Plot effects
puy_eff <- tibble(
  type = c("Wet", "Lipid", "1% Lipid"),
  eff = c(eff_ww, eff_lw, eff_lw1),
  oceanadults = c(puy_oa_ww, puy_oa_lw, puy_oa_lw1),
  oa0 = puy0_oa
) |>
  mutate(
    oa_change = (oceanadults - oa0) / oa0
  )
write_rds(puy_eff, "data/puy-eff.rds")

puy_eff |>
  ggplot(aes(xdist = eff, y = type)) +
  stat_histinterval() +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()

puy_eff |>
  ggplot(aes(xdist = oceanadults, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = puy0_oa, linetype = "dashed") +
  scale_x_continuous(
    name = "Number of 1+ Chinook in the ocean",
    labels = scales::comma
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()

puy_eff |>
  ggplot(aes(xdist = oa_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in adult Chinook ocean abundance",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
