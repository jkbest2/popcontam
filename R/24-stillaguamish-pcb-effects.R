library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/pcb-effects.R")

post <- read_rds("data/stillaguamish/pcb_exposure.rds")

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
  post$ww$pop_meanlog,
  post$ww$pop_sdlog,
  "ww"
)
## Get roundoff errors using a rel.tol that is too low, so increased here.
eff_lw <- rv_pcb_effect(
  post$lw$pop_meanlog,
  post$lw$pop_sdlog,
  wt_type = "lw",
  rel.tol = 1e-6
)
eff_lw1 <- rv_pcb_effect(
  post$lw1$pop_meanlog,
  post$lw1$pop_sdlog,
  "lw",
  rel.tol = 1e-6
)

write_rds(
  list(
    ww = eff_ww,
    lw = eff_lw,
    lw1 = eff_lw1
  ),
  "data/stillaguamish/pcb_eff.rds"
)
