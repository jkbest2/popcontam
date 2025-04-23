library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/stillaguamish.R")

## Get baseline population distribution
stilly0 <- eq_pop(
  stillaguamish_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5)
)
stilly0_sp <- get_spawners(stilly0)

eff <- read_rds("data/stillaguamish/pcb_eff.rds")

stillaguamish_exposed <- function(ns_surv, stage = NULL) {
  eqp <- eq_pop(
    stillaguamish_sim,
    nearshore_surv_adj = ns_surv,
    pop0 = rep(1000, 5)
  )
  ## Don't always want the population age structure
  if (!is.null(stage)) {
    eqp <- attr(eqp, stage)
  }
  eqp
}
rv_stillaguamish_exposed <- rfun(stillaguamish_exposed)

stilly_sp_ww <- rv_stillaguamish_exposed(eff$ww, "spawners")
stilly_sp_lw <- rv_stillaguamish_exposed(eff$lw, "spawners")
stilly_sp_lw1 <- rv_stillaguamish_exposed(eff$lw1, "spawners")

write_rds(
  list(
    sp0 = stilly0_sp,
    ww = stilly_sp_ww,
    lw = stilly_sp_lw,
    lw1 = stilly_sp_lw1
  ),
  "data/stillaguamish/stilly_pcb_spawners.rds"
)

stilly_eff <- tibble(
  type = c("Wet", "Lipid", "1% Lipid"),
  eff = c(eff$ww, eff$lw, eff$lw1),
  spawners = c(stilly_sp_ww, stilly_sp_lw, stilly_sp_lw1),
  sp0 = stilly0_sp
) |>
  mutate(
    sp_change = (spawners - sp0) / sp0
  )
write_rds(stilly_eff, "data/stillaguamish/stilly_pcb_eff_df.rds")
