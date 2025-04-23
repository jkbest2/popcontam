library(tidyverse)
library(posterior)
library(mgcv)
library(ggdist)

source(here::here("R", "pbde-effects.R"))

pbde_exp <- read_rds(here::here("data", "stillaguamish", "pbde_exposure.rds"))
pbde_surv <- read_rds(here::here("data", "pbde_surv.rds"))

pbde_mod <- read_rds(here::here("data", "pbde_model.rds"))
dr_mod <- pbde_mod$pbde_model
threshold <- pbde_mod$pbde_threshold

eff <- pbde_eff(pbde_exp, dr_mod, threshold, pbde_surv, thin = 10)

write_rds(eff, here::here("data", "stillaguamish", "stilly-pbde-eff.rds"))
