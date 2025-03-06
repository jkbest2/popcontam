library(tidyverse)
library(mgcv)
library(gratia)
library(posterior)
library(ggdist)

pbde_surv_full <- read_csv(
  here::here("data", "pbde_survival_raw.csv"),
  col_types = cols(
    strata = col_character(),
    time = col_double(),
    n.risk = col_double(),
    n.event = col_double(),
    surv = col_double(),
    std.err = col_double(),
    lower = col_double(),
    upper = col_double()
  )
)

pbde_conc <- read_csv(
  here::here("data", "pbde_survival.csv"),
  col_types = cols(
    strata = col_character(),
    time = col_skip(),
    n.risk = col_skip(),
    n.event = col_skip(),
    surv = col_skip(),
    std.err = col_skip(),
    lower = col_skip(),
    upper = col_skip(),
    concentration = col_double()
  )
)

## Need to account for individuals dropped from the study ------------
pbde_start <- pbde_surv_full |>
  slice_min(time, by = strata) |>
  select(strata, start_n_risk = n.risk)

pbde_end <- pbde_surv_full |>
  slice_max(time, by = strata) |>
  mutate(end_surv = n.risk - n.event) |>
  select(strata, end_surv)

pbde_surv <- pbde_surv_full |>
  summarize(n_event = sum(n.event), .by = strata) |>
  left_join(pbde_start, by = join_by(strata)) |>
  left_join(pbde_end, by = join_by(strata)) |>
  left_join(pbde_conc, by = join_by(strata)) |>
  mutate(
    vulnerable = end_surv + n_event,
    n_surv = end_surv,
    n_dead = n_event,
    p_surv = n_surv / vulnerable,
    q10 = qbeta(0.1, n_surv, n_dead),
    q90 = qbeta(0.9, n_surv, n_dead),
  ) |>
  select(strata, concentration, vulnerable, n_surv, n_dead, p_surv, q10, q90)

write_rds(
  pbde_surv,
  here::here("data", "pbde_surv.rds")
)
