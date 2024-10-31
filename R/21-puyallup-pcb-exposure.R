library(tidyverse)
library(posterior)

source("R/fit-lnorm-pop-exposure.R")

## This is the data frame containing only samples from the Puyallup estuary and
## nearshore. Currently it is processed and saved in R/exposure-map.R
contam <- read_rds("data/puyallup/pcb_estns.rds")

data_ww <- lnorm_data(contam, pcb_ug_ww)
data_lw <- lnorm_data(contam, pcb_ug_lw)
data_lw1 <- lnorm_data(contam, pcb_ug_lw1)

fit_ww <- fit_lnorm_exposure(data_ww)
post_ww <- as_draws_rvars(fit_ww)

fit_lw <- fit_lnorm_exposure(data_lw)
post_lw <- as_draws_rvars(fit_lw)

fit_lw1 <- fit_lnorm_exposure(data_lw1)
post_lw1 <- as_draws_rvars(fit_lw1)

write_rds(
  list(
    ww = post_ww,
    lw = post_lw,
    lw1 = post_lw1
  ),
  "data/puyallup/pcb_exposure.rds"
)
