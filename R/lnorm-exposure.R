library(tidyverse)

### ---------------------------------------------------------------------------
### Simulation
sumlnorm_pars <- function(n, pop_mu, pop_sig) {
  pop_sig2 <- pop_sig^2
  sum_sig2 <- log(1 / n * (exp(pop_sig2) - 1) + 1)
  sum_mu <- log(n) + pop_mu + pop_sig2 / 2 - sum_sig2 / 2
  c(meanlog = sum_mu, meansd = sqrt(sum_sig2))
}

sumlnorm_pars(2, 0, 1)

x <- matrix(rlnorm(10 * 1000), nrow = 1000)
## Need to apply the cumsum function over rows, but output of apply produces
## answers in columns, so need to transpose to get back to the original shape.
xs <- t(apply(x, 1, cumsum))
xm <- sapply(1:10, \(i) xs[, i] / i)

df <- tibble(n = rep(1:10, each = 1000)) |>
  mutate(
    xs = c(xs),
    xm = c(xm)
  )

pars <- map(1:10, sumlnorm_pars, pop_mu = 0, pop_sig = 1)
max_sum <- apply(xs, 2, max)
max_mean <- apply(xm, 2, max)

df_sum_approx <- tibble(n = 1:10) |>
  mutate(data = map(
    n,
    \(i) {
      tibble(
        xs = seq(0, max_sum[i], length.out = 1025),
        ds = dlnorm(xs, pars[[i]][1], pars[[i]][2])
      )
    }
  )) |>
  unnest(data)

df_mean_approx <- tibble(n = 1:10) |>
  mutate(data = map(
    n,
    \(i) {
      tibble(
        xm = seq(0, max_mean[i], length.out = 1025),
        dm = dlnorm(xm, pars[[i]][1] - log(i), pars[[i]][2])
      )
    }
  )) |>
  unnest(data)

## Sum plot
ggplot(df, aes(x = xs)) +
  geom_histogram(aes(y = after_stat(density))) +
  geom_line(data = df_sum_approx, aes(y = ds)) +
  facet_wrap(~n, scales = "free")

## Mean plot
ggplot(df, aes(x = xm)) +
  geom_histogram(
    aes(y = after_stat(density)),
    breaks = seq(0, 30, 0.1)
  ) +
  geom_line(data = df_mean_approx, aes(y = dm)) +
  coord_cartesian(xlim = c(0, 5), expand = FALSE) +
  facet_wrap(~n, ncol = 1)

### ---------------------------------------------------------------------------
### Stan model fit
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(posterior)
library(bayesplot)
library(ggdist)

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

contam |>
  ggplot(aes(x = pcb_ug_ww)) +
  geom_histogram() +
  facet_wrap(~n_composite) +
  scale_y_continuous(breaks = seq(0, 10, 2))

lnorm_data <- function(contam, conc_col) {
  ncomp <- sort(unique(contam$n_composite))
  N <- nrow(contam) # nolint: object_name_linter.
  contam <- contam |>
    mutate(
      ncomp_idx = match(contam$n_composite, ncomp),
      conc = {{ conc_col }}
    )

  list(
    C = length(ncomp),
    ncomp = as.array(ncomp), # In case there is only a single ncomp value
    N = N,
    ncomp_idx = contam$ncomp_idx,
    conc = contam$conc
  )
}

data_ln <- lnorm_data(contam, pcb_ug_ww)
write_rds(data_ln, "data/pcb_ww_ln_data.rds")

fit_ln <- stan("inst/lnorm-exposure.stan",
  data = data_ln,
  chains = 4,
  iter = 2000
)
write_rds(fit_ln, "data/pcb_ww_ln_fit.rds")

post_ln <- as_draws_rvars(fit_ln)
write_rds(post_ln, "data/pcb_ww_ln_post.rds")

post_ln$pop_meanlog
post_ln$pop_sdlog

## This plot shows that the gamma model is not accounting for enough low values;
## this also means that it may be
ppc_dens_overlay(data_ln$conc, as_draws_matrix(post_ln$conc_gen)[1:50, ]) +
  xlim(c(0, 0.25))
ppc_dens_overlay(data_ln$conc, as_draws_matrix(post_ln$conc_gen2)[1:50, ]) +
  xlim(c(0, 0.25))

## This seems to show that the PIT residuals don't fall outside the expected
## envelope
ppc_pit_ecdf(data_ln$conc, as_draws_matrix(post_ln$conc_gen)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
ppc_pit_ecdf(data_ln$conc, as_draws_matrix(post_ln$conc_gen2)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

x <- seq(0, 0.3, length = 1025)

ci_ln <- as_draws_df(fit_ln) |>
  select(.draw, pop_meanlog, pop_sdlog) |>
  # slice_sample(n = 100) |>
  mutate(x = list(x)) |>
  unnest(x) |>
  mutate(y = dlnorm(x, pop_meanlog, pop_sdlog)) |>
  group_by(x) |>
  curve_interval(y, .width = 0.8)

ci_ln |>
  ggplot(aes(x = x, y = y)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper), fill = "skyblue") +
  coord_cartesian(xlim = c(0, 0.2))

ci_ln |>
  ggplot(aes(x = x, y = y)) +
  geom_lineribbon(aes(ymin = .lower, ymax = .upper), fill = "skyblue") +
  xlim(0.1, 0.3) +
  ylim(0, 1)


### ---------------------------------------------------------------------------
### Simulated data model fits
n <- 1000
pop_meanlog <- 0
pop_sdlog <- 1
sim_df <- tibble(
  n_composite = rpois(n, 100) + 1
) |>
  mutate(
    conc = map_dbl(n_composite, \(nc) mean(rlnorm(nc, pop_meanlog, pop_sdlog)))
  )

sim_df |>
  ggplot(aes(x = conc)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 7, 0.1)) +
  facet_wrap(~n_composite, scales = "free")

data_sim <- lnorm_data(sim_df, conc)

fit_sim <- stan("lnorm-exposure.stan", data = data_sim, chains = 4, iter = 2000)

post_sim <- as_draws_rvars(fit_sim)

ppc_dens_overlay(
  sim_df$conc,
  as_draws_matrix(post_sim$conc_gen)[1:100, ]
)

ppc_pit_ecdf(
  sim_df$conc,
  as_draws_matrix(post_sim$conc_gen)
) +
  geom_abline(linetype = "dashed") +
  coord_cartesian(expand = FALSE)

library(patchwork)

(mcmc_hist(fit_sim, pars = c("pop_meanlog")) +
  vline_at(pop_meanlog, linetype = "dashed")) +
  (mcmc_hist(fit_sim, pars = c("pop_sdlog")) +
    vline_at(pop_sdlog, linetype = "dashed"))
