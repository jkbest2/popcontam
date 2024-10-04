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
    )
) |>
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

gamma_data <- function(contam, conc_col) {
  contam <- contam |>
    rename(conc = {{ conc_col }})

  list(
    N = nrow(contam),
    ncomp = contam$n_composite,
    conc = contam$conc
  )
}

data_g <- gamma_data(contam, pcb_ug_lw1)

fit_g <- stan("gamma-exposure.stan",
  data = data_g,
  chains = 4,
  iter = 2000
)

post_g <- as_draws_rvars(fit_g)

## This plot shows that the gamma model is not accounting for enough low values;
## this also means that it may be
ppc_dens_overlay(data_g$conc, as_draws_matrix(post_g$conc_post)[1:100, ]) +
  xlim(c(0, max(data_g$conc) * 1.1))

## This seems to show that the PIT residuals don't fall outside the expected
## envelope
ppc_pit_ecdf(data_g$conc, as_draws_matrix(post_g$conc_post)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

ci_g <- as_draws_df(fit_g) |>
  select(.draw, shape, scale) |>
  # slice_sample(n = 100) |>
  mutate(x = list(x)) |>
  unnest(x) |>
  mutate(y = dgamma(x, shape, scale)) |>
  group_by(x) |>
  curve_interval(y, .width = 0.5)

ci_g |>
  ggplot(aes(x = x, y = y)) +
  # geom_line(alpha = 0.1)
  geom_lineribbon(aes(ymin = .lower, ymax = .upper), fill = "skyblue") +
  coord_cartesian(xlim = c(0, 0.2))

ci_g |>
  ggplot(aes(x = x, y = y)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper), fill = "skyblue") +
  coord_cartesian(xlim = c(0.1, 0.3), ylim = c(0, 5))

## Mapping the observations shows (unsurprisingly) that the higher observations
## are in the estuary. So maybe considering distance from the estuary as a
## covariate would be helpful here? There aren't really any other options here.
## Though possibly the number of samples in the composite could be related to
## the size (and hence age), so time to accumulate exposure, and could be
## usefuL? Very noisy measure however.
library(sf)

st_as_sf(contam, coords = c("longitude", "latitude")) |>
  ggplot(aes(color = pcb_ug_ww)) +
  geom_sf(alpha = 0.7) +
  scale_color_viridis_c()

contam_sf <- st_as_sf(
  contam,
  coords = c("longitude", "latitude"),
  crs = st_crs(4326)
) |>
  st_transform(crs = st_crs(32149))

coord_svd <- st_coordinates(contam_sf) |>
  apply(2, scale, center = TRUE, scale = FALSE) |>
  svd()

### ---------------------------------------------------------------------------
### Simulated data fits
n <- 1000
sim_df <- tibble(
  n_composite = rpois(n, 5) + 1
) |>
  mutate(
    ## n_composite = 1,
    conc = map_dbl(n_composite, \(nc) mean(rgamma(nc, shape = 2, scale = 0.5)))
  )

sim_df |>
  ggplot(aes(x = conc)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 7, 0.1)) +
  facet_wrap(~n_composite, scales = "free")

data_sim <- gamma_data(sim_df, conc)

fit_sim <- stan(
  "gamma-exposure.stan",
  data = data_sim,
  ## pars = c("shape", "scale"),
  chains = 4,
  iter = 2000
)

post_sim <- as_draws_rvars(fit_sim)

post_sim$shape
post_sim$scale

ppc_dens_overlay(
  sim_df$conc,
  as_draws_matrix(post_sim$conc_post)[1:50, ]
)

ppc_pit_ecdf(
  sim_df$conc,
  as_draws_matrix(post_sim$conc_post)
) +
  geom_abline(linetype = "dashed") +
  coord_cartesian(expand = FALSE)
