library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/pcb-effects.R")
source("R/white.R")

white0 <- eq_pop(
  white_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 10),
  n_max = 500L
)
white_sar <- attr(white0, "subyearspawners") /
  (attr(white0, "parr_mig") + attr(white0, "fry_mig") + attr(white0, "yr_mig"))

post <- read_rds("data/puyallup/pcb_exposure.rds")

eff_ww <- rv_pcb_effect(
  pop_meanlog = post$ww$pop_meanlog,
  pop_sdlog = post$ww$pop_sdlog,
  base_surv = white_sar,
  wt_type = "ww",
  remove_pcbs = TRUE
)
eff_lw <- rv_pcb_effect(
  pop_meanlog = post$lw$pop_meanlog,
  pop_sdlog = post$lw$pop_sdlog,
  base_surv = white_sar,
  wt_type = "lw",
  remove_pcbs = TRUE
)
eff_lw1 <- rv_pcb_effect(
  pop_meanlog = post$lw1$pop_meanlog,
  pop_sdlog = post$lw1$pop_sdlog,
  base_surv = white_sar,
  wt_type = "lw",
  remove_pcbs = TRUE
)

write_rds(
  list(
    ww = eff_ww,
    lw = eff_lw,
    lw1 = eff_lw1
  ),
  "data/puyallup/white_pcb_eff.rds"
)
