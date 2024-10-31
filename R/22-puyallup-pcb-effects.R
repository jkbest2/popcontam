library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/pcb-effects.R")

post <- read_rds("data/puyallup/pcb_exposure.rds")

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

eff_ww <- rv_pcb_effect(post$ww$pop_meanlog, post$ww$pop_sdlog, "ww")
eff_lw <- rv_pcb_effect(post$lw$pop_meanlog, post$lw$pop_sdlog, "lw")
eff_lw1 <- rv_pcb_effect(post$lw1$pop_meanlog, post$lw1$pop_sdlog, "lw")

write_rds(
  list(
    ww = eff_ww,
    lw = eff_lw,
    lw1 = eff_lw1
  ),
  "data/puyallup/pcb_eff.rds"
)
