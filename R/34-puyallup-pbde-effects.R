library(tidyverse)
library(posterior)
library(mgcv)
library(ggdist)

pbde_exp <- read_rds(here::here("data", "puyallup", "pbde_exposure.rds"))
pbde_surv <- read_rds(here::here("data", "pbde_surv.rds"))

## Fit model --------------------------------
dr_mod <- gam(
  cbind(n_surv, n_dead) ~ s(sqrt(concentration), k = 6, bs = "bs", m = c(2, 2)),
  data = pbde_surv,
  family = binomial()
)

dr_coefs_rv <- function(mod = dr_mod, ndraws = 4e3) {
  beta <- coef(mod)
  n_beta <- length(beta)

  v <- vcov(mod)
  vchol <- chol(v)

  x <- rvar_rng(rnorm, n_beta, ndraws = ndraws)
  beta + t(vchol) %*% x
}

gen_base_surv <- function(threshold, surv_df = pbde_surv, ndraws = 4e3) {
  thr_surv <- surv_df |>
    filter(concentration < threshold) |>
    summarize(n_surv = sum(n_surv), n_dead = sum(n_dead))
  rvar_rng(rbeta, 1, thr_surv$n_surv + 1, thr_surv$n_dead + 1, ndraws = ndraws)
}

exp_rel_surv <- function(conc, beta, pop_meanlog, pop_sdlog, base_surv, thr = 0) {
  map_dbl(
    conc,
    function(conc) {
      if (conc < thr) {
        surv <- base_surv
      } else {
        dm <- predict(
          dr_mod,
          newdata = data.frame(concentration = conc),
          type = "lpmatrix"
        )
        surv <- plogis(dm %*% beta)
      }
      surv / base_surv * dlnorm(conc, pop_meanlog, pop_sdlog)
    }
  )
}

pbde_eff <- function(
    exp_post = pbde_exp,
    mod = dr_mod,
    threshold = 0,
    surv = pbde_surv,
    thin = 10) {
  pop_meanlog <- exp_post$pop_meanlog |>
    merge_chains() |>
    thin_draws(thin) |>
    draws_of()
  pop_sdlog <- exp_post$pop_sdlog |>
    merge_chains() |>
    thin_draws(thin) |>
    draws_of()
  beta_draws <- dr_coefs_rv(mod, ndraws = length(pop_meanlog)) |>
    draws_of()
  beta <- map(
    seq_len(nrow(beta_draws)),
    ~ beta_draws[., , , drop = TRUE]
  )
  base_surv <- gen_base_surv(threshold, surv, ndraws = length(pop_meanlog)) |>
    draws_of()
  thr <- rep(threshold, length(pop_meanlog))

  pmap_dbl(
    list(
      beta = beta,
      pop_meanlog = pop_meanlog,
      pop_sdlog = pop_sdlog,
      base_surv = base_surv,
      thr = thr
    ),
    function(beta, pop_meanlog, pop_sdlog, base_surv, thr) {
      integrate(
        function(conc) {
          exp_rel_surv(
            conc, beta,
            pop_meanlog, pop_sdlog,
            base_surv, thr
          )
        },
        lower = 0, upper = Inf
      )$value
    }
  ) |>
    rvar()
}

## Geometric mean of second and third concentrations
gm23 <- sqrt(prod(pbde_surv$concentration[2:3]))
eff0 <- pbde_eff(thr = 0, thin = 10)
eff3 <- pbde_eff(thr = 3, thin = 10)
eff5.5 <- pbde_eff(thr = gm23, thin = 10)
eff7 <- pbde_eff(thr = 7, thin = 10)

thr_df <- tibble(
  threshold = c("0 ng/g", "3 ng/g", "5.5 ng/g", "7 ng/g"),
  rsurv = c(eff0, eff3, eff5.5, eff7)
)

ggplot(thr_df, aes(xdist = rsurv, y = threshold)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  stat_halfeye(aes(fill = threshold)) +
  scale_x_continuous(
    name = "Relative Survival", labels = scales::percent,
    breaks = seq(0, 1.2, 0.05), minor_breaks = seq(0, 1.2, 0.01)
  ) +
  labs(y = "Threshold") +
  guides(color = "none", fill = "none")
ggsave(here::here("figs", "pbde-effects-thresholds.png"), width = 8, height = 5)

dr_pred_rv <- function(conc, mod = dr_mod, n = 4e3, threshold = 0) {
  newdata <- tibble(
    concentration = conc
  )
  nd0 <- filter(newdata, concentration < threshold)
  nd1 <- filter(newdata, concentration >= threshold)

  pred0 <- rep(
    gen_base_surv(threshold, pbde_surv, ndraws = 1000),
    nrow(nd0)
  )

  beta <- coef(mod)
  n_beta <- length(beta)

  v <- vcov(mod)
  vchol <- chol(v)

  x <- rvar_rng(rnorm, n_beta, ndraws = n)
  beta_sim <- beta + t(vchol) %*% x
  covar_sim <- predict(mod, newdata = nd1, type = "lpmatrix")
  pred1 <- covar_sim %**% beta_sim
  inv_link <- family(mod)$linkinv
  pred1 <- rfun(inv_link)(pred1)
  c(pred0, pred1)
}

pred_df <- tibble(
  concentration = seq(0, sqrt(250), length.out = 257)^2
) |>
  mutate(
    `0` = dr_pred_rv(concentration, n = 1e3, threshold = 0),
    `3` = dr_pred_rv(concentration, n = 1e3, threshold = 3),
    `5.5` = dr_pred_rv(concentration, n = 1000, threshold = gm23),
    `7` = dr_pred_rv(concentration, n = 1e3, threshold = 7)
  ) |>
  pivot_longer(
    c(`0`, `3`, `5.5`, `7`),
    names_to = "threshold",
    names_transform = as.numeric,
    values_to = "survival"
  ) |>
  mutate(threshold_label = paste(threshold, "ng/g")) |>
  point_interval(survival, .width = 0.8)

pred_df |>
  ggplot(aes(x = concentration, y = survival)) +
  geom_ribbon(
    aes(ymin = .lower, ymax = .upper, fill = threshold_label),
    alpha = 0.2
  ) +
  geom_line(aes(color = threshold_label)) +
  geom_point(
    data = pbde_surv,
    aes(y = p_surv),
    alpha = 0.5
  ) +
  geom_errorbar(
    data = filter(pred_df, concentration == 0),
    aes(
      x = threshold,
      y = survival,
      ymin = .lower, ymax = .upper,
      color = threshold_label
    )
  ) +
  scale_x_continuous(
    name = "PBDE Concentration (ng/g ww)"
  ) +
  scale_y_continuous(
    name = "Probability of survival",
    labels = scales::percent
  ) +
  labs(color = "Threshold", fill = "Threshold")
ggsave(here::here("figs", "pbde-dose-response.png"), width = 8, height = 5)
