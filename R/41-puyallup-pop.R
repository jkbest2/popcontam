library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/puyallup.R")

## Get baseline population distribution
puy0 <- eq_pop(
  puyallup_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5),
  n_max = 500L
)
puy0_sp <- get_spawners(puy0)

puy_eff <- read_rds("data/puyallup/puy_pcb_eff.rds")

puyallup_exposed <- function(ns_surv, stage = NULL) {
  eqp <- eq_pop(
    puyallup_sim,
    nearshore_surv_adj = ns_surv,
    pop0 = rep(1000, 5),
    n_max = 500L
  )
  ## Don't always want the population age structure
  if (!is.null(stage)) {
    eqp <- attr(eqp, stage)
  }
  eqp
}
rv_puyallup_exposed <- rfun(puyallup_exposed)

puy_sp_ww <- rv_puyallup_exposed(puy_eff$ww, "spawners")
puy_sp_lw <- rv_puyallup_exposed(puy_eff$lw, "spawners")
puy_sp_lw1 <- rv_puyallup_exposed(puy_eff$lw1, "spawners")

puy_ww <- rv_puyallup_exposed(puy_eff$ww)
puy_lw <- rv_puyallup_exposed(puy_eff$lw)
puy_lw1 <- rv_puyallup_exposed(puy_eff$lw1)

puy_age00 <- puy0[1]
puy_age00_ww <- puy_ww[1]
puy_age00_lw <- puy_lw[1]
puy_age00_lw1 <- puy_lw1[1]

write_rds(
  list(
    sp = puy0_sp,
    ww = puy_sp_ww,
    lw = puy_sp_lw,
    lw1 = puy_sp_lw1
  ),
  "data/puyallup/puy_pcb_spawners.rds"
)

puy_eff <- tibble(
  type = c("Wet", "Lipid", "1% Lipid"),
  eff = c(eff$ww, eff$lw, eff$lw1),
  spawners = c(puy_sp_ww, puy_sp_lw, puy_sp_lw1),
  sp0 = puy0_sp
) |>
  mutate(
    sp_change = (spawners - sp0) / sp0
  )
write_rds(puy_eff, "data/puyallup/puy_pcb_eff_df.rds")
