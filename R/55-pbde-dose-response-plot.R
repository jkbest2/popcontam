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


pred_mort <- function(conc, dr_mod, threshold, surv = pbde_surv, ndraws = 4e3) {
  base_mort <- 1 - gen_base_surv(threshold, surv, ndraws)
  beta <- dr_coefs_rv(dr_mod, ndraws)

  under_thr <- which(conc <= threshold)
  over_thr <- which(conc > threshold)

  mort <- rep(base_mort, length(conc))
  if (length(over_thr) > 0) {
    mort[over_thr] <- 1 - dr_surv(conc[over_thr], beta, threshold)
  }

  mort
}

mort_df <- tibble(
  conc = seq(0, 300, 1),
  mort = pred_mort(conc, dr_mod, threshold, pbde_surv, 4e2)
)

mort_df |>
  point_interval() |>
  ggplot(mort_df, aes(x = conc, y = mort, ymin = .lower, ymax = .upper)) +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  labs(
    x = "PBDE Concentration (ng/g wet weight)",
    y = "Predicted mortality"
  )
