library(tidyverse)
library(posterior)

source("R/utils.R")
source("R/pcb-effects.R")
source("R/white.R")

if (!file.exists("data/pcb_ww_ln_post.rds")) {
  stop("Need to run the exposure model first")
}
post_ln <- read_rds("data/pcb_ww_ln_post.rds")

white0 <- eq_pop(
  white_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 10)
)
white0_oa <- get_oceanadults(white0)

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
white_exposed <- function(exp_meanlog, exp_sdlog) {
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
  eq_pop(white_sim, nearshore_surv_adj = ns_surv, pop0 = rep(1000, 10))
}
rv_white_exposed <- rfun(white_exposed)

white <- rv_white_exposed(post_ln$pop_meanlog, post_ln$pop_sdlog)
white_oa <- rvar_sum(white[2:10])

write_rds(list(white0 = white0, white = white), "data/white_pop.rds")

as_draws_df(white_oa) |>
  ggplot(aes(x = x)) +
  geom_density() +
  geom_vline(xintercept = get_oceanadults(white0), linetype = "dashed")

oa_reduction <- (white_oa - white0_oa) / white0_oa

as_draws_df(oa_reduction) |>
  ggplot(aes(x = x)) +
  geom_density() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = scales::percent)
