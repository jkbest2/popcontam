library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/white.R")

white0 <- eq_pop(
  white_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 10)
)
white0_oa <- get_oceanadults(white0)

eff <- read_rds("data/puyallup/pcb_eff.rds")

## TODO Why are PCBs sometimes *increasing* survival?
white_exposed <- function(ns_surv, stage = NULL) {
  eqp <- eq_pop(white_sim, nearshore_surv_adj = ns_surv, pop0 = rep(1000, 10))
  ## Don't always want population age structure
  if (!is.null(stage)) {
    eqp <- attr(eqp, stage)
  }
  eqp
}
rv_white_exposed <- rfun(white_exposed)

white_oa_ww <- rv_white_exposed(eff$ww, "oceanadults")
white_oa_lw <- rv_white_exposed(eff$lw, "oceanadults")
white_oa_lw1 <- rv_white_exposed(eff$lw1, "oceanadults")

write_rds(
  list(
    oa0 = white0_oa,
    ww = white_oa_ww,
    lw = white_oa_lw,
    lw1 = white_oa_lw1
  ),
  "data/puyallup/white_pcb_oceanadults.rds"
)

white_eff <- tibble(
  type = c("Wet", "Lipid", "1% Lipid"),
  eff = c(eff$ww, eff$lw, eff$lw1),
  oceanadults = c(white_oa_ww, white_oa_lw, white_oa_lw1),
  oa0 = white0_oa
) |>
  mutate(
    oa_change = (oceanadults - oa0) / oa0
  )
write_rds(white_eff, "data/puyallup/white_pcb_eff_df.rds")
