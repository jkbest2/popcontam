library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/white.R")

white0 <- eq_pop(
  white_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 10),
  n_max = 500L
)
white0_sp <- get_spawners(white0)

white_eff <- read_rds("data/puyallup/white_pcb_eff.rds")

## TODO Why are PCBs sometimes *increasing* survival?
white_exposed <- function(ns_surv, stage = NULL) {
  eqp <- eq_pop(
    white_sim,
    nearshore_surv_adj = ns_surv,
    pop0 = rep(1000, 10),
    n_max = 500L
  )
  ## Don't always want population age structure
  if (!is.null(stage)) {
    eqp <- attr(eqp, stage)
  }
  eqp
}
rv_white_exposed <- rfun(white_exposed)

white_sp_ww <- rv_white_exposed(white_eff$ww, "spawners")
white_sp_lw <- rv_white_exposed(white_eff$lw, "spawners")
white_sp_lw1 <- rv_white_exposed(white_eff$lw1, "spawners")

white_ww <- rv_white_exposed(white_eff$ww)
white_lw <- rv_white_exposed(white_eff$lw)
white_lw1 <- rv_white_exposed(white_eff$lw1)

white_age00 <- white0[1]
white_age00_ww <- white_ww[1]
white_age00_lw <- white_lw[1]
white_age00_lw1 <- white_lw1[1]

tibble()

write_rds(
  list(
    sp0 = white0_sp,
    ww = white_sp_ww,
    lw = white_sp_lw,
    lw1 = white_sp_lw1
  ),
  "data/puyallup/white_pcb_spawners.rds"
)

white_eff <- tibble(
  type = c("Wet", "Lipid", "1% Lipid"),
  eff = c(eff$ww, eff$lw, eff$lw1),
  spawners = c(white_sp_ww, white_sp_lw, white_sp_lw1),
  sp0 = white0_sp
) |>
  mutate(
    sp_change = (spawners - sp0) / sp0
  )
write_rds(white_eff, "data/puyallup/white_pcb_eff_df.rds")
