library(tidyverse)
library(posterior)
library(mgcv)
library(ggdist)

source(here::here("R", "pbde-effects.R"))

pbde_exp <- read_rds(here::here("data", "puyallup", "pbde_exposure.rds"))
pbde_surv <- read_rds(here::here("data", "pbde_surv.rds"))

pbde_mod <- read_rds(here::here("data", "pbde_model.rds"))
dr_mod <- pbde_mod$pbde_model
threshold <- pbde_mod$pbde_threshold

eff <- pbde_eff(pbde_exp, dr_mod, threshold, pbde_surv, thin = 10)

write_rds(eff, here::here("data", "puyallup", "puy-pbde-eff.rds"))

## Fit model --------------------------------
# dr_mod <- gam(
#   cbind(n_surv, n_dead) ~ s(sqrt(concentration), k = 6, bs = "bs", m = c(2, 2)),
#   data = pbde_surv,
#   family = binomial()
# )

# dr_coefs_rv <- function(mod = dr_mod, ndraws = 4e3) {
#   beta <- coef(mod)
#   n_beta <- length(beta)

#   v <- vcov(mod)
#   vchol <- chol(v)

#   x <- rvar_rng(rnorm, n_beta, ndraws = ndraws)
#   beta + t(vchol) %*% x
# }

# gen_base_surv <- function(threshold, surv_df = pbde_surv, ndraws = 4e3) {
#   thr_surv <- surv_df |>
#     filter(concentration < threshold) |>
#     summarize(n_surv = sum(n_surv), n_dead = sum(n_dead))
#   rvar_rng(rbeta, 1, thr_surv$n_surv + 1, thr_surv$n_dead + 1, ndraws = ndraws)
# }


## Geometric mean of second and third concentrations
# gm23 <- sqrt(prod(pbde_surv$concentration[2:3]))
# eff0 <- pbde_eff(thr = 0, thin = 10)
# eff3 <- pbde_eff(thr = 3, thin = 10)
# eff5.5 <- pbde_eff(thr = gm23, thin = 10)
# eff7 <- pbde_eff(thr = 7, thin = 10)

# thr_df <- tibble(
#   threshold = c("0 ng/g", "3 ng/g", "5.5 ng/g", "7 ng/g"),
#   rsurv = c(eff0, eff3, eff5.5, eff7)
# )

# ggplot(thr_df, aes(xdist = rsurv, y = threshold)) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   stat_halfeye(aes(fill = threshold)) +
#   scale_x_continuous(
#     name = "Relative Survival", labels = scales::percent,
#     breaks = seq(0, 1.2, 0.05), minor_breaks = seq(0, 1.2, 0.01)
#   ) +
#   labs(y = "Threshold") +
#   guides(color = "none", fill = "none")
# ggsave(here::here("figs", "pbde-effects-thresholds.png"), width = 8, height = 5)


# pred_df <- tibble(
#   concentration = seq(0, sqrt(250), length.out = 257)^2
# ) |>
#   mutate(
#     `0` = dr_pred_rv(concentration, n = 1e3, threshold = 0),
#     `3` = dr_pred_rv(concentration, n = 1e3, threshold = 3),
#     `5.5` = dr_pred_rv(concentration, n = 1000, threshold = gm23),
#     `7` = dr_pred_rv(concentration, n = 1e3, threshold = 7)
#   ) |>
#   pivot_longer(
#     c(`0`, `3`, `5.5`, `7`),
#     names_to = "threshold",
#     names_transform = as.numeric,
#     values_to = "survival"
#   ) |>
#   mutate(threshold_label = paste(threshold, "ng/g")) |>
#   point_interval(survival, .width = 0.8)

# pred_df |>
#   ggplot(aes(x = concentration, y = survival)) +
#   geom_ribbon(
#     aes(ymin = .lower, ymax = .upper, fill = threshold_label),
#     alpha = 0.2
#   ) +
#   geom_line(aes(color = threshold_label)) +
#   geom_point(
#     data = pbde_surv,
#     aes(y = p_surv),
#     alpha = 0.5
#   ) +
#   geom_errorbar(
#     data = filter(pred_df, concentration == 0),
#     aes(
#       x = threshold,
#       y = survival,
#       ymin = .lower, ymax = .upper,
#       color = threshold_label
#     )
#   ) +
#   scale_x_continuous(
#     name = "PBDE Concentration (ng/g ww)"
#   ) +
#   scale_y_continuous(
#     name = "Probability of survival",
#     labels = scales::percent
#   ) +
#   labs(color = "Threshold", fill = "Threshold")
# ggsave(here::here("figs", "pbde-dose-response.png"), width = 8, height = 5)
