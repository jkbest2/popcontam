library(tidyverse)
library(posterior)
library(mgcv)
library(ggdist)

source(here::here("R", "pbde-effects.R"))

pbde_exp <- read_rds(here::here("data", "puyallup", "pbde_exposure.rds"))
pbde_surv <- read_rds(here::here("data", "pbde_surv.rds"))

## Fit model --------------------------------
dr_mod <- gam(
  cbind(n_surv, n_dead) ~ s(sqrt(concentration), k = 6, bs = "bs", m = c(2, 2)),
  data = pbde_surv,
  family = binomial()
)


beta <- dr_coefs_rv()
base_surv <- gen_base_surv(5)

benchresp <- base_surv - 0.05

dr_surv(1:10, beta)

surv_df <- tibble(
  concentration = seq(2.8, 7.8, length.out = 200)
) |>
  mutate(
    surv = dr_surv(concentration, beta),
    diff = base_surv - surv,
    diff_mean = mean(diff),
    thr = diff > 0.05,
    thr_mean = mean(thr)
  )

surv_df |>
  mutate(mean_diff = mean(diff)) |>
  ggplot(aes(x = concentration, ydist = diff)) +
  stat_lineribbon(alpha = 0.8) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = "dashed")

thr_idx <- detect_index(surv_df$diff_mean, \(d) d > 0.05)

surv_df[(thr_idx - 2):(thr_idx + 2), ]

pbde_threshold <- 6.6

write_rds(
  list(pbde_model = dr_mod, pbde_threshold = pbde_threshold),
  here::here("data", "pbde_model.rds")
)
