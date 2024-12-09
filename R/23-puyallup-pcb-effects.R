library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/pcb-effects.R")
source("R/puyallup.R")

## Get baseline population distribution
puy0 <- eq_pop(
  puyallup_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5),
  n_max = 500L
)
puy_sar <- attr(puy0, "spawners") /
  (attr(puy0, "parr_mig") + attr(puy0, "fry_mig"))

post <- read_rds("data/puyallup/pcb_exposure.rds")

eff_ww <- rv_pcb_effect(
  pop_meanlog = post$ww$pop_meanlog,
  pop_sdlog = post$ww$pop_sdlog,
  base_surv = puy_sar,
  wt_type = "ww",
  remove_pcbs = TRUE
)
eff_lw <- rv_pcb_effect(
  pop_meanlog = post$lw$pop_meanlog,
  pop_sdlog = post$lw$pop_sdlog,
  base_surv = puy_sar,
  wt_type = "lw",
  remove_pcbs = TRUE
)
eff_lw1 <- rv_pcb_effect(
  pop_meanlog = post$lw1$pop_meanlog,
  pop_sdlog = post$lw1$pop_sdlog,
  base_surv = puy_sar,
  wt_type = "lw",
  remove_pcbs = TRUE
)

write_rds(
  list(
    ww = eff_ww,
    lw = eff_lw,
    lw1 = eff_lw1
  ),
  "data/puyallup/puy_pcb_eff.rds"
)
