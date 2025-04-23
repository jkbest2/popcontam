library(tidyverse)
library(posterior)
library(ggdist)

source(here::here("R", "utils.R"))
source(here::here("R", "puyallup.R"))

## Get baseline population distribution
puy0 <- eq_pop(
  puyallup_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5),
  n_max = 500L
)
puy0_sp <- get_spawners(puy0)

pbde_eff <- read_rds("data/puyallup/puy-pbde-eff.rds")

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

puy_sp <- rv_puyallup_exposed(pbde_eff, "spawners")

write_rds(
  list(
    base = puy0_sp,
    exposed = puy_sp
  ),
  "data/puyallup/puy_pbde_spawners.rds"
)

puy_eff <- tibble(
  type = "PBDE",
  eff = pbde_eff,
  spawners = puy_sp,
  sp0 = puy0_sp
) |>
  mutate(
    sp_change = (spawners - sp0) / sp0
  )
write_rds(puy_eff, "data/puyallup/puy_pbde_eff_df.rds")
