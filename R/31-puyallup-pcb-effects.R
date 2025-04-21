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

## Combination effect, removing PCBs, will be used in population model --------
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

## Direct mortality only, not used in population model ------------------------
eff_ww_dm <- rv_pcb_effect(
  pop_meanlog = post$ww$pop_meanlog,
  pop_sdlog = post$ww$pop_sdlog,
  base_surv = puy_sar,
  wt_type = "ww",
  eff_type = "dir_mort",
  remove_pcbs = FALSE
)
eff_lw_dm <- rv_pcb_effect(
  pop_meanlog = post$lw$pop_meanlog,
  pop_sdlog = post$lw$pop_sdlog,
  base_surv = puy_sar,
  wt_type = "lw",
  eff_type = "dir_mort",
  remove_pcbs = FALSE
)
eff_lw1_dm <- rv_pcb_effect(
  pop_meanlog = post$lw1$pop_meanlog,
  pop_sdlog = post$lw1$pop_sdlog,
  base_surv = puy_sar,
  wt_type = "lw",
  eff_type = "dir_mort",
  remove_pcbs = FALSE
)

write_rds(
  list(
    ww = eff_ww_dm,
    lw = eff_lw_dm,
    lw1 = eff_lw1_dm
  ),
  "data/puyallup/puy_pcb_dm_eff.rds"
)

## Growth-related mortality only, not used in population model ----------------
eff_ww_gr <- rv_pcb_effect(
  pop_meanlog = post$ww$pop_meanlog,
  pop_sdlog = post$ww$pop_sdlog,
  base_surv = puy_sar,
  wt_type = "ww",
  eff_type = "gr_mort",
  remove_pcbs = FALSE
)
eff_lw_gr <- rv_pcb_effect(
  pop_meanlog = post$lw$pop_meanlog,
  pop_sdlog = post$lw$pop_sdlog,
  base_surv = puy_sar,
  wt_type = "lw",
  eff_type = "gr_mort",
  remove_pcbs = FALSE
)
eff_lw1_gr <- rv_pcb_effect(
  pop_meanlog = post$lw1$pop_meanlog,
  pop_sdlog = post$lw1$pop_sdlog,
  base_surv = puy_sar,
  wt_type = "lw",
  eff_type = "gr_mort",
  remove_pcbs = FALSE
)

write_rds(
  list(
    ww = eff_ww_gr,
    lw = eff_lw_gr,
    lw1 = eff_lw1_gr
  ),
  "data/puyallup/puy_pcb_gr_eff.rds"
)
