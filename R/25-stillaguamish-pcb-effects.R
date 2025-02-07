library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/pcb-effects.R")
source("R/stillaguamish.R")

## Get baseline population distribution, calculate smolt-to-adult return ratio
## in order to get expected July size based on Duffy and Beauchamp.
stilly0 <- eq_pop(
  stillaguamish_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5)
)
sar <- attr(stilly0, "spawners") /
  (attr(stilly0, "delta_fry") + attr(stilly0, "parr_mig"))

post <- read_rds("data/stillaguamish/pcb_exposure.rds")

## Integration success is sensitive to outliers in the estimated population
## parameters. In particular, values of `pop_sdlog` larger than ~2.5 can result
## in divergent integrals.
eff_ww <- rv_pcb_effect(
  pop_meanlog = post$ww$pop_meanlog,
  pop_sdlog = post$ww$pop_sdlog,
  wt_type = "ww",
  base_surv = sar,
  remove_pcbs = TRUE
)
eff_lw <- rv_pcb_effect(
  pop_meanlog = post$lw$pop_meanlog,
  pop_sdlog = post$lw$pop_sdlog,
  wt_type = "lw",
  base_surv = sar,
  remove_pcbs = TRUE,
  ## Get roundoff errors using a rel.tol that is too low, so increased here.
  rel.tol = 1e-6
)
eff_lw1 <- rv_pcb_effect(
  pop_meanlog = post$lw1$pop_meanlog,
  pop_sdlog = post$lw1$pop_sdlog,
  wt_type = "lw",
  base_surv = sar,
  remove_pcbs = TRUE
)

write_rds(
  list(
    ww = eff_ww,
    lw = eff_lw,
    lw1 = eff_lw1
  ),
  "data/stillaguamish/pcb_eff.rds"
)
