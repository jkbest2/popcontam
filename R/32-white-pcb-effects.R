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

## Direct mortality only, not used in population model ------------------------
eff_ww_dm <- rv_pcb_effect(
  pop_meanlog = post$ww$pop_meanlog,
  pop_sdlog = post$ww$pop_sdlog,
  base_surv = white_sar,
  wt_type = "ww",
  eff_type = "dir_mort",
  remove_pcbs = FALSE
)
eff_lw_dm <- rv_pcb_effect(
  pop_meanlog = post$lw$pop_meanlog,
  pop_sdlog = post$lw$pop_sdlog,
  base_surv = white_sar,
  wt_type = "lw",
  eff_type = "dir_mort",
  remove_pcbs = FALSE
)
eff_lw1_dm <- rv_pcb_effect(
  pop_meanlog = post$lw1$pop_meanlog,
  pop_sdlog = post$lw1$pop_sdlog,
  base_surv = white_sar,
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
  "data/puyallup/white_pcb_dm_eff.rds"
)

## Growth-related mortality only, not used in population model ----------------
eff_ww_gr <- rv_pcb_effect(
  pop_meanlog = post$ww$pop_meanlog,
  pop_sdlog = post$ww$pop_sdlog,
  base_surv = white_sar,
  wt_type = "ww",
  eff_type = "gr_mort",
  remove_pcbs = FALSE
)
eff_lw_gr <- rv_pcb_effect(
  pop_meanlog = post$lw$pop_meanlog,
  pop_sdlog = post$lw$pop_sdlog,
  base_surv = white_sar,
  wt_type = "lw",
  eff_type = "gr_mort",
  remove_pcbs = FALSE
)
eff_lw1_gr <- rv_pcb_effect(
  pop_meanlog = post$lw1$pop_meanlog,
  pop_sdlog = post$lw1$pop_sdlog,
  base_surv = white_sar,
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
  "data/puyallup/white_pcb_gr_eff.rds"
)
