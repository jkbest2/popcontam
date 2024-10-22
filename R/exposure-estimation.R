library(tidyverse)
library(pracma)

## dsumlnorm <- function(x, n, meanlog, sdlog, log = FALSE, nsamp = 10e3) {
##   s <- replicate(nsamp, sum(rlnorm(n, meanlog, sdlog)))
##   ef <- ecdf(s)
##   d <- fderiv(e, x, 1)
## }

## d <-numdiff(e, seq(0,75, length = 1025))

## dens <- density(s)
## plot(dens)

## ss <- sort(s)
## which.min(abs(ss - 10))

## diff(e(ss[6096:6097])) / diff(ss[6096:6097])

## find_interval <- function(x, vec, n = 1, unique = TRUE) {
##   v_lower <- -sort(-vec[vec <= x], partial = n) |> head(n)
##   v_upper <- sort(vec[vec >= x], partial = n) |> head(n)
##   int <- c(v_lower, v_upper)
##   if (unique)
##     int <- unique(int)
##   return (int)
## }

## polygrad <- function(x, vec, n) {
##   vec0 <- vec - x
##   e <- ecdf(vec0)
##   vals <- find_interval(0, vec0, n, unique = TRUE)
##   val_ecdf <- e(vals)
##   m <- length(vals)

##   p <- polyfit(vals, val_ecdf, n = 3)
##   X <- Reduce()
## }

## n <- 6
## nsamp <- 1e7
## meanlog <- 0
## sdlog <- 1

## s <- replicate(nsamp, sum(rlnorm(n, meanlog, sdlog)))

## vec0 <- vec - x
## e <- ecdf(vec0)
## vals <- find_interval(0, vec0, 20, unique = TRUE)
## val_ecdf <- e(vals)
## m <- length(vals)

## p1 <- polyfit(vals, val_ecdf, n = 1)
## p2 <- polyfit(vals, val_ecdf, n = 2)
## p3 <- polyfit(vals, val_ecdf, n = 3)
## p4 <- polyfit(vals, val_ecdf, n = 4)

## curve(polyval(p1, x), from = min(vals), to = max(vals))
## curve(polyval(p2, x), from = min(vals), to = max(vals), add = TRUE)
## curve(polyval(p3, x), from = min(vals), to = max(vals), add = TRUE)
## curve(polyval(p4, x), from = min(vals), to = max(vals), add = TRUE)
## points(vals, val_ecdf)

dmeangamma <- function(x, n, shape, rate, log = FALSE) {
  dgamma(x, n * shape,  n * rate, log = log)
}
rmeangamma <- function(nsamp, n, shape, rate) {
  replicate(nsamp, mean(rgamma(n, shape, rate)))
}

curve(dgamma(x, 2, 1), from = 0, to = 8, n = 1025, ylim = c(0, 1))
curve(dmeangamma(x, 2, 2, 1), from = 0, to = 8, n = 1025, add = TRUE)
curve(dmeangamma(x, 3, 2, 1), from = 0, to = 8, n = 1025, add = TRUE)
curve(dmeangamma(x, 4, 2, 1), from = 0, to = 8, n = 1025, add = TRUE)
curve(dmeangamma(x, 5, 2, 1), from = 0, to = 8, n = 1025, add = TRUE)
curve(dmeangamma(x, 6, 2, 1), from = 0, to = 8, n = 1025, add = TRUE)
curve(dmeangamma(x, 7, 2, 1), from = 0, to = 8, n = 1025, add = TRUE)
curve(dmeangamma(x, 8, 2, 1), from = 0, to = 8, n = 1025, add = TRUE)
curve(dmeangamma(x, 9, 2, 1), from = 0, to = 8, n = 1025, add = TRUE)
curve(dmeangamma(x, 10, 2, 1), from = 0, to = 8, n = 1025, add = TRUE)

library(tidyverse)

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


contam |>
  select(species, pcb_ug_lw1, n_composite) |>
  ggplot(aes(x = factor(n_composite), y = pcb_ug_lw1)) +
  geom_violin()

nll <- function(pars, conc, ncomp) {
  shape <- exp(pars[1])
  rate <- exp(pars[2])
  ll <- map2_dbl(conc, ncomp, \(.x, .y) dmeangamma(.x, .y, shape, rate, log = TRUE))
  -sum(ll)
}

nll(c(0, 0), contam$pcb_ug_lw1, contam$n_composite)

### With lipid concentrations normalized to 1%
opt <- optim(c(0, 0), nll, conc = contam$pcb_ug_lw1, ncomp = contam$n_composite)

shape <- exp(opt$par[1])
rate <- exp(opt$par[2])

post_dens <- expand_grid(n_composite = contam$n_composite,
                         pcb_ug_lw1 = seq(7e-3, 8, length = 1025)) |>
  mutate(d = dmeangamma(pcb_ug_lw1, n_composite, shape, rate))

ggplot(post_dens, aes(x = pcb_ug_lw1, y = d)) +
  geom_histogram(aes(y = after_stat(density)),
                 data = contam,
                 bins = 30,
                 alpha = 0.5) +
  geom_line() +
  facet_wrap(~ n_composite) +
  scale_x_continuous(expand = expansion(c(0, 0.05))) +
  scale_y_continuous(limits = c(0, 2),
                     expand = expansion(c(0, 0.05)))

### Raw lipid concentrations
opt <- optim(c(0, 0), nll, conc = contam$pcb_ug_lw, ncomp = contam$n_composite)

shape <- exp(opt$par[1])
rate <- exp(opt$par[2])

post_dens <- expand_grid(n_composite = contam$n_composite,
                         pcb_ug_lw = seq(7e-3, 8, length = 1025)) |>
  mutate(d = dmeangamma(pcb_ug_lw, n_composite, shape, rate))

ggplot(post_dens, aes(x = pcb_ug_lw, y = d)) +
  geom_histogram(aes(y = after_stat(density)),
                 data = contam,
                 bins = 30,
                 alpha = 0.5) +
  geom_line() +
  facet_wrap(~ n_composite) +
  scale_x_continuous(expand = expansion(c(0, 0.05))) +
  scale_y_continuous(limits = c(0, 2),
                     expand = expansion(c(0, 0.05)))

### Water weight concentrations
opt <- optim(c(0, 0), nll, conc = contam$pcb_ug_ww, ncomp = contam$n_composite)

shape <- exp(opt$par[1])
rate <- exp(opt$par[2])

post_dens <- expand_grid(n_composite = contam$n_composite,
                         pcb_ug_ww = seq(1e-3, 0.25, length = 1025)) |>
  mutate(d = dmeangamma(pcb_ug_ww, n_composite, shape, rate))

ggplot(post_dens, aes(x = pcb_ug_ww, y = d)) +
  geom_histogram(aes(y = after_stat(density)),
                 data = contam,
                 bins = 30,
                 breaks = seq(0, 0.25, 0.01),
                 alpha = 0.5) +
  ## geom_vline(aes(xintercept = pcb_ug_ww), data = contam) +
  geom_line() +
  facet_wrap(~ n_composite) +
  scale_x_continuous(expand = expansion(c(0, 0.05))) +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(c(0, 0.05)))

library(RTMB)

contam2 <- contam |>
  select(n_composite,
         conc = pcb_ug_ww)


nll2 <- function(pars) {
  getAll(pars, contam2)
  shape <- exp(logshape)
  rate <- exp(lograte)
  conc <- OBS(conc)
  jnll <- 0
  for (i in seq_along(n_composite)) {
    n <- n_composite[i]
    jnll <- jnll - dgamma(conc[i], n * shape, n * rate, log = TRUE)
  }
  jnll
}

obj <- MakeADFun(nll2, list(logshape = log(shape), lograte = log(rate)))
