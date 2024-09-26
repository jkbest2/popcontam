library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(posterior)
library(bayesplot)

contam <- read_csv(
  "data/PCB_Puyallup_monitoring_data.csv",
  col_types =
    cols(
      SAMPLE_ID = col_character(),
      Study_ID = col_character(),
      Latitude = col_double(),
      Longitude = col_double(),
      Species = col_character(),
      RiverSystem = col_character(),
      COMPOUND = col_character(),
      UNIT = col_character(),
      `CONC_FOUND (ng/g ww)` = col_skip(),
      `Sampling Date` = col_date(format = "%m/%d/%Y"),
      `% Lipids` = col_double(),
      pcb_uggwet = col_double(),
      pcb_ugglw_sample = col_double(),
      pcb_ugglw_1 = col_double(),
      `effect-mort ww` = col_skip(),
      `effect-mort lipid-sample specific` = col_skip(),
      `effect-mort lipid - 1%` = col_skip(),
      Composite_Count = col_double()
    )) |>
  rename(
    id = SAMPLE_ID,
    study = Study_ID,
    latitude = Latitude,
    longitude = Longitude,
    species = Species,
    river = RiverSystem,
    compound = COMPOUND,
    unit = UNIT,
    date = `Sampling Date`,
    pct_lipid = `% Lipids`,
    pcb_ug_ww = pcb_uggwet,
    pcb_ug_lw = pcb_ugglw_sample,
    pcb_ug_lw1 = pcb_ugglw_1,
    n_composite = Composite_Count,
  )

data <- list(N = length(contam$n_composite),
             ncomp = contam$n_composite,
             conc = contam$pcb_ug_ww)

fit <- stan("gamma-exposure.stan",
            data = data,
            chains = 4,
            iter = 2000)

post <- as_draws_rvars(fit)

## This plot shows that the gamma model is not accounting for enough low values;
## this also means that it may be 
ppc_dens_overlay(data$conc, as_draws_matrix(post$conc_post)[1:50, ]) +
xlim(c(0,
 1))

## This seems to show that the PIT residuals don't fall outside the expected envelope
ppc_pit_ecdf(data$conc, as_draws_matrix(post$conc_post)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed")

## Mapping the observations shows (unsurprisingly) that the higher observations
## are in the estuary. So maybe considering distance from the estuary as a
## covariate would be helpful here? There aren't really any other options here.
## Though possibly the number of samples in the composite could be related to
## the size (and hence age), so time to accumulate exposure, and could be usefuL?
## Very noisy measure however.
library(sf)

st_as_sf(contam, coords = c("longitude", "latitude")) |>
ggplot(aes(color = pcb_ug_ww)) +
geom_sf(alpha = 0.7) +
scale_color_viridis_c()
