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
stilly0_oa <- get_oceanadults(stilly0)

eff <- read_rds("data/stillaguamish/pcb_eff.rds")

stillaguamish_exposed <- function(ns_surv, stage = NULL) {
  eqp <- eq_pop(stillaguamish_sim, nearshore_surv_adj = ns_surv, pop0 = rep(1000, 5))
  ## Don't always want the population age structure
  if (!is.null(stage)) {
    eqp <- attr(eqp, stage)
  }
  eqp
}
rv_stillaguamish_exposed <- rfun(stillaguamish_exposed)

stilly_oa_ww <- rv_stillaguamish_exposed(eff$ww, "oceanadults")
stilly_oa_lw <- rv_stillaguamish_exposed(eff$lw, "oceanadults")
stilly_oa_lw1 <- rv_stillaguamish_exposed(eff$lw1, "oceanadults")

write_rds(
  list(
    oa0 = stilly0_oa,
    ww = stilly_oa_ww,
    lw = stilly_oa_lw,
    lw1 = stilly_oa_lw1
  ),
  "data/stillaguamish/stilly_pcb_oceanadults.rds"
)

stilly_eff <- tibble(
  type = c("Wet", "Lipid", "1% Lipid"),
  eff = c(eff$ww, eff$lw, eff$lw1),
  oceanadults = c(stilly_oa_ww, stilly_oa_lw, stilly_oa_lw1),
  oa0 = stilly0_oa
) |>
  mutate(
    oa_change = (oceanadults - oa0) / oa0
  )
write_rds(stilly_eff, "data/stillaguamish/stilly_pcb_eff_df.rds")
