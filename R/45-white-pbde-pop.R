library(tidyverse)
library(posterior)
library(ggdist)

source(here::here("R", "utils.R"))
source(here::here("R", "white.R"))

## Get baseline population distribution
white0 <- eq_pop(
  white_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 10),
  n_max = 500L
)
white0_sp <- get_spawners(white0)

## PBDE effects are assumed the same between Puyallup and White because there
## is no growth-related mortality here.
pbde_eff <- read_rds("data/puyallup/puy-pbde-eff.rds")

white_exposed <- function(ns_surv, stage = NULL) {
  eqp <- eq_pop(
    white_sim,
    nearshore_surv_adj = ns_surv,
    pop0 = rep(1000, 10),
    n_max = 500L
  )
  ## Don't always want the population age structure
  if (!is.null(stage)) {
    eqp <- attr(eqp, stage)
  }
  eqp
}
rv_white_exposed <- rfun(white_exposed)

white_sp <- rv_white_exposed(pbde_eff, "spawners")

write_rds(
  list(
    base = white0_sp,
    exposed = white_sp
  ),
  "data/puyallup/white_pbde_spawners.rds"
)

white_eff <- tibble(
  type = "PBDE",
  eff = pbde_eff,
  spawners = white_sp,
  sp0 = white0_sp
) |>
  mutate(
    sp_change = (spawners - sp0) / sp0
  )
write_rds(white_eff, "data/puyallup/white_pbde_eff_df.rds")
