library(tidyverse)
library(posterior)
library(ggdist)

source(here::here("R", "utils.R"))
source(here::here("R", "stillaguamish.R"))

## Get baseline population distribution
stilly0 <- eq_pop(
  stillaguamish_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5),
  n_max = 500L
)
stilly0_sp <- get_spawners(stilly0)

pbde_eff <- read_rds("data/stillaguamish/stilly-pbde-eff.rds")

stillaguamish_exposed <- function(ns_surv, stage = NULL) {
  eqp <- eq_pop(
    stillaguamish_sim,
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
rv_stillaguamish_exposed <- rfun(stillaguamish_exposed)

stilly_sp <- rv_stillaguamish_exposed(pbde_eff, "spawners")

write_rds(
  list(
    base = stilly0_sp,
    exposed = stilly_sp
  ),
  "data/stillaguamish/pbde_spawners.rds"
)

stilly_eff <- tibble(
  type = "PBDE",
  eff = pbde_eff,
  spawners = stilly_sp,
  sp0 = stilly0_sp
) |>
  mutate(
    sp_change = (spawners - sp0) / sp0
  )
write_rds(stilly_eff, "data/stillaguamish/pbde_eff_df.rds")
