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

## Fit model --------------------------------
mod <- gam(
  cbind(n_surv, n_dead) ~ s(sqrt(concentration), k = 6, bs = "bs", m = c(2, 2)),
  data = pbde_surv,
  family = binomial()
)

## Predict values with uncertainty
pred_df <- tibble(
  concentration = seq(0, sqrt(300), length.out = 257)^2
)

gam_pred_rv <- function(mod, newdata, n = 4e3, type = "response") {
  beta <- coef(mod)
  n_beta <- length(beta)

  v <- vcov(mod)
  vchol <- chol(v)

  x <- rvar_rng(rnorm, n_beta, ndraws = n)
  beta_sim <- beta + t(vchol) %*% x
  covar_sim <- predict(mod, newdata = newdata, type = "lpmatrix")
  pred <- covar_sim %**% beta_sim
  if (type == "response") {
    inv_link <- family(mod)$linkinv
    pred <- rfun(inv_link)(pred)
  }
  pred
}

# pred_df |>
#   mutate(
#     pred_rv = gam_pred_rv(mod, pred_df)
#   ) |>
#   point_interval(pred_rv, .width = 0.8) |>
#   ggplot(aes(x = concentration, y = pred_rv)) +
#   geom_ribbon(
#     aes(
#       ymin = .lower,
#       ymax = .upper
#     ),
#     alpha = 0.2
#   ) +
#   geom_line() +
#   geom_pointrange(
#     data = pbde_surv,
#     aes(
#       y = p_surv,
#       ymin = q10, ymax = q90,
#     )
#   ) +
#   scale_x_continuous(
#     name = "PBDE Concentration (ng/g ww)"
#   ) +
#   scale_y_continuous(
#     name = "Probability of survival",
#     labels = scales::percent
#   )
# ggsave(here::here("figs", "pbde-dose-response.png"))
