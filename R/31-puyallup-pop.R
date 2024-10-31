library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/puyallup.R")

## Get baseline population distribution
puy0 <- eq_pop(
  puyallup_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5)
)
puy0_oa <- get_oceanadults(puy0)

eff <- read_rds("data/puyallup/pcb_eff.rds")

puyallup_exposed <- function(ns_surv, stage = NULL) {
  eqp <- eq_pop(puyallup_sim, nearshore_surv_adj = ns_surv, pop0 = rep(1000, 5))
  ## Don't always want the population age structure
  if (!is.null(stage)) {
    eqp <- attr(eqp, stage)
  }
  eqp
}
rv_puyallup_exposed <- rfun(puyallup_exposed)

puy_oa_ww <- rv_puyallup_exposed(eff$ww, "oceanadults")
puy_oa_lw <- rv_puyallup_exposed(eff$lw, "oceanadults")
puy_oa_lw1 <- rv_puyallup_exposed(eff$lw1, "oceanadults")

write_rds(
  list(
    oa0 = puy0_oa,
    ww = puy_oa_ww,
    lw = puy_oa_lw,
    lw1 = puy_oa_lw1
  ),
  "data/puyallup/puy_pcb_oceanadults.rds"
)

puy_eff <- tibble(
  type = c("Wet", "Lipid", "1% Lipid"),
  eff = c(eff$ww, eff$lw, eff$lw1),
  oceanadults = c(puy_oa_ww, puy_oa_lw, puy_oa_lw1),
  oa0 = puy0_oa
) |>
  mutate(
    oa_change = (oceanadults - oa0) / oa0
  )
write_rds(puy_eff, "data/puyallup/puy_pcb_eff_df.rds")
