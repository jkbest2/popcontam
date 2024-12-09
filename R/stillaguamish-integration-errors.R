library(tidyverse)
library(posterior)
library(ggdist)

source("R/utils.R")
source("R/pcb-effects.R")
source("R/stillaguamish.R")

## Get baseline population distribution, calculate smolt-to-adult return ratio
## in order to get expected July size based on Duffy and Beauchamp.
stilly0 <- eq_pop(
  stillaguamish_sim,
  nearshore_surv_adj = 1,
  pop0 = rep(1000, 5)
)
sar <- attr(stilly0, "spawners") /
  (attr(stilly0, "delta_fry") + attr(stilly0, "parr_mig"))

post <- read_rds("data/stillaguamish/pcb_exposure.rds")

df <- as_draws_df(post$lw1) |>
  select(.draw, pop_meanlog, pop_sdlog) |>
  rowwise() |>
  mutate(
    eff = possibly(pcb_effect, NA)(
      pop_meanlog,
      pop_sdlog,
      base_surv = sar,
      wt_type = "lw",
      rel.tol = .Machine$double.eps^0.25),
    miss = is.na(eff)
  )

ggplot(df, aes(x = pop_meanlog, y = pop_sdlog)) +
  geom_hex() +
  geom_point(data = filter(df, miss), color = "red") +
  geom_point(data = filter(df, !miss), color = "green")

df2 <- df |>
  filter(miss) |>
  expand_grid(pcb = seq(0, 1, length.out = 1025)) |>
  mutate(eff2 = dlnorm(pcb, pop_meanlog, pop_sdlog) * combo_surv(pcb, base_surv = sar, wt_type = "lw"))

ggplot(df2, aes(x = pcb, y = eff2, group = .draw)) +
  geom_line()

expected_surv <- function(pcb, pop_meanlog, pop_sdlog) {
  dlnorm(pcb, pop_meanlog, pop_sdlog) *
    combo_surv(pcb, base_surv = sar, wt_type = "lw", remove_pcbs = TRUE)
}

df3 <- df |>
  filter(miss) |>
  select(pop_meanlog, pop_sdlog)

integrate(expected_surv,
  lower = 0, upper = Inf,
  pop_meanlog = df3$pop_meanlog[1],
  pop_sdlog = df3$pop_sdlog[1]
)
integrate(expected_surv,
  lower = 0, upper = Inf,
  pop_meanlog = df3$pop_meanlog[2],
  pop_sdlog = df3$pop_sdlog[2],
  # rel.tol = 1e-8,
  # subdivisions = 1000
)

eff_grid <- expand_grid(
  pop_meanlog = seq(-1, 0.75, length.out = 257),
  pop_sdlog = seq(0.4, 3.5, length.out = 257)
) |>
  rowwise() |>
  mutate(
    eff = possibly(integrate, list(value = NA))(
      expected_surv,
      lower = 0, upper = Inf,
      pop_meanlog = pop_meanlog,
      pop_sdlog = pop_sdlog
    )$value
  )

ggplot(eff_grid, aes(x = pop_meanlog, y = pop_sdlog, z = eff)) +
  geom_contour_filled() +
  geom_point(data = df, aes(x = pop_meanlog, y = pop_sdlog), color = "red") +
  geom_point(data = as_draws_df(post$lw) |> mutate(eff = 1.1), color = "green") +
  geom_abline(slope = -0.3, intercept = 2.45, linetype = "dashed") +
  geom_abline(slope = -0.15, intercept = 2.475, linetype = "dashed")
