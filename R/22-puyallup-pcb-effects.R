library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/pcb-effects.R")

post <- read_rds("data/puyallup/pcb_exposure.rds")

pcb_effect <- function(
    pop_meanlog, pop_sdlog, wt_type, remove_pcbs = TRUE,
    rel.tol = .Machine$double.eps^0.5, subdivisions = 100) {
  expected_surv <- function(pcb) {
    dlnorm(pcb, pop_meanlog, pop_sdlog) *
      combo_surv(pcb, wt_type = wt_type, remove_pcbs = remove_pcbs)
  }
  integrate(
    expected_surv,
    lower = 0, upper = Inf,
    ## The default rel.tol occasionally gives incorrect results, including one
    ## set of population parameters that consistently shows an increase in
    ## survival due to PCB exposure!
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )$value
}
rv_pcb_effect <- rfun(pcb_effect)

eff_ww <- rv_pcb_effect(
  pop_meanlog = post$ww$pop_meanlog,
  pop_sdlog = post$ww$pop_sdlog,
  wt_type = "ww",
  remove_pcbs = TRUE
)
eff_lw <- rv_pcb_effect(
  pop_meanlog = post$lw$pop_meanlog,
  pop_sdlog = post$lw$pop_sdlog,
  wt_type = "lw",
  remove_pcbs = TRUE
)
eff_lw1 <- rv_pcb_effect(
  pop_meanlog = post$lw1$pop_meanlog,
  pop_sdlog = post$lw1$pop_sdlog,
  wt_type = "lw",
  remove_pcbs = TRUE
)

write_rds(
  list(
    ww = eff_ww,
    lw = eff_lw,
    lw1 = eff_lw1
  ),
  "data/puyallup/pcb_eff.rds"
)
