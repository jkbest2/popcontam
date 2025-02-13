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

pbde_surv <- read_csv(
  here::here("data", "pbde_survival.csv"),
  col_types = cols(
    strata = col_character(),
    time = col_double(),
    n.risk = col_double(),
    n.event = col_double(),
    surv = col_double(),
    std.err = col_double(),
    lower = col_double(),
    upper = col_double(),
    concentration = col_double()
  )
) |>
  left_join(pbde_surv_nrisk, by = join_by(strata)) |>
  mutate(
    n_surv = n.risk - n.event,
    n_dead = n_risk - n_surv,
    p_surv = n_surv / n_risk
  ) |>
  select(strata, concentration, time, n_risk, n_surv, n_dead, p_surv, surv)

ggplot(
  pbde_surv_full,
  aes(
    x = time,
    y = surv, ymin = lower, ymax = upper,
    color = strata, fill = strata
  )
) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_step()

# ggplot(pbde_surv, aes(x = concentration)) +
#   geom_errorbar(aes(ymin = lower, ymax = upper)) +
#   geom_point(aes(y = surv)) +
#   geom_line(aes(y = surv))

fake_threshold_data <- function(df, max_conc_quantile, n_supp, n_risk) {
  obs_conc <- df$concentration
  obs_risk <- df$n_risk
  obs_surv <- df$n_surv
  obs_dead <- df$n_dead

  threshold_p <- sum(obs_surv) / sum(obs_risk)

  conc_rate <- -log(1 - max_conc_quantile) / max(obs_conc)
  supp_conc <- rexp(n_supp, conc_rate)
  supp_surv <- rbinom(n_supp, n_risk, threshold_p)

  tibble(
    strata = "SupplementalPrior",
    concentration = supp_conc,
    time = NA,
    n_risk = n_risk,
    n_surv = supp_surv,
    n_dead = n_risk - supp_surv,
    p_surv = n_surv / n_risk,
    surv = NA
  )
}

supp_data <- fake_threshold_data(pbde_surv[1:2, ], 0.95, 0, 100)

dat <- rbind(pbde_surv, supp_data) |>
  mutate(
    obs_type = ifelse(
      strata == "SupplementalPrior",
      "Supplemental",
      "Observed"
    ),
    q025 = qbeta(0.025, n_surv, n_dead),
    q10 = qbeta(0.1, n_surv, n_dead),
    q90 = qbeta(0.9, n_surv, n_dead),
    q975 = qbeta(0.975, n_surv, n_dead)
  )

mod <- gam(
  cbind(n_surv, n_dead) ~ s(sqrt(concentration), k = 6, bs = "bs", m = c(2, 3)),
  data = dat,
  family = binomial()
)


# draw(mod)
# draw(smooth_estimates(mod))
# appraise(mod)

pred_df <- tibble(
  concentration = seq(0, sqrt(200), length.out = 257)^2
)
# pred_df$pred <- predict(mod, newdata = pred_df, type = "response")
# pred_df |>
#   ggplot(aes(x = concentration, y = pred)) +
#   geom_line() +
#   geom_point(data = dat, aes(y = p_surv, color = obs_type))

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

pred_df |>
  mutate(
    pred_rv = gam_pred_rv(mod, pred_df)
  ) |>
  point_interval(pred_rv, .width = 0.8) |>
  ggplot(aes(x = concentration, y = pred_rv)) +
  geom_ribbon(
    aes(
      ymin = .lower,
      ymax = .upper
    ),
    alpha = 0.2
  ) +
  geom_line() +
  geom_pointrange(
    data = dat,
    aes(
      y = p_surv,
      ymin = q10, ymax = q90,
      color = obs_type
    )
  ) +
  scale_x_continuous(
    name = "PBDE Concentration (ng/g ww)"
  ) +
  scale_y_continuous(
    name = "Probability of survival",
    labels = scales::percent
  )



# mod <- lm(
#   surv ~ bs(concentration, df = 4, degree = 2),
#   weights = 1 / (std.err),
#   data = pbde_surv
# )
# pred <- tibble(
#   concentration = seq(0, 180, length.out = 129)
# ) |>
#   mutate(
#     surv = predict(mod, newdata = tibble(concentration = concentration))
#   )

# plot(surv ~ concentration, data = pbde_surv, type = "b", ylim = c(0.5, 0.8))
# lines(pred$concentration, pred$surv, col = "blue")

# ggplot(pbde_surv, aes(x = concentration, y = surv)) +
#   geom_errorbar(aes(ymin = lower, ymax = upper)) +
#   geom_point() +
#   geom_line(data = pred) +
#   geom_vline(xintercept = 27) +
#   labs(y = "Survival", x = "Concentration") +
#   scale_y_continuous(labels = scales::percent)
