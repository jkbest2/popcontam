library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/pcb-effects.R")
source("R/puyallup.R")

if (!file.exists("data/pcb_ww_ln_post.rds")) {
  stop("Need to run the exposure model first")
}
post_ln <- read_rds("data/pcb_ww_ln_post.rds")

puy0 <- eq_pop(
  puyallup_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5)
)
puy0_oa <- get_oceanadults(puy0)

pcb_effect <- function(exp_meanlog, exp_sdlog) {
  exp_surv <- function(pcb) {
    dlnorm(pcb, exp_meanlog, exp_sdlog) *
      (1 - combo_mort(pcb))
  }
  ns_surv <- integrate(
    exp_surv,
    lower = 0, upper = Inf,
    ## The default rel.tol occasionally gives incorrect results, including one
    ## set of population parameters that consistently shows an increase in
    ## survival due to PCB exposure!
    rel.tol = .Machine$double.eps^0.5
  )$value
}
rv_pcb_effect <- rfun(pcb_effect)
eff <- rv_pcb_effect(post_ln$pop_meanlog, post_ln$pop_sdlog)

as_draws_df(eff) |>
  ggplot(aes(x = x)) +
  geom_histogram()

## TODO Why are PCBs sometimes *increasing* survival?
puyallup_exposed <- function(exp_meanlog, exp_sdlog) {
  exp_surv <- function(pcb) {
    dlnorm(pcb, exp_meanlog, exp_sdlog) *
      (1 - combo_mort(pcb))
  }
  ns_surv <- integrate(
    exp_surv,
    lower = 0, upper = Inf,
    rel.tol = .Machine$double.eps^0.5
  )$value
  if (ns_surv > 1) {
    stop("PCBs increasing survival?")
  }
  eq_pop(puyallup_sim, nearshore_surv_adj = ns_surv, pop0 = rep(1000, 5))
}
rv_puyallup_exposed <- rfun(puyallup_exposed)

puy <- rv_puyallup_exposed(post_ln$pop_meanlog, post_ln$pop_sdlog)
puy_oa <- rvar_sum(puy[2:5])

as_draws_df(puy_oa) |>
  ggplot(aes(x = x)) +
  geom_density() +
  geom_vline(xintercept = get_oceanadults(puy0), linetype = "dashed")

oa_reduction <- (puy_oa - puy0_oa) / puy0_oa

as_draws_df(oa_reduction) |>
  ggplot(aes(x = x)) +
  geom_density() +
  scale_x_continuous(labels = scales::percent)
