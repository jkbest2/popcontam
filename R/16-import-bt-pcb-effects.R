library(tidyverse)
library(posterior)

bt_eff <- read_csv(
  here::here("data", "Berninger and Tillitt Raw Datapoints.csv"),
  col_types = cols(
    concentration = col_double(),
    effect = col_double(),
    endpoint = col_factor()
  )
) |>
  mutate(effect = effect / 100)

ggplot(bt_eff, aes(x = concentration, y = effect, color = endpoint)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_x_log10()

bt_mort <- bt_eff |>
  filter(endpoint == "Mortality")
bt_growth <- bt_eff |>
  filter(endpoint == "Growth")

mod_mort <- lm(effect ~ log10(concentration), data = bt_mort)
mod_growth <- lm(effect ~ log10(concentration), data = bt_growth)

write_rds(
  list(
    mod_mort = mod_mort,
    summ_mort = summary(mod_mort),
    mod_growth = mod_growth,
    summ_growth = summary(mod_growth)
  ),
  here::here("data", "pcb-effect-regressions.rds")
)

## TODO: Figure out how to convert these to lipid numbers
pred_eff <- function(conc, mod, min_val = 0) {
  dm <- cbind(rep(1, length(conc)), log10(conc))
  summ <- summary(mod)
  beta <- rvar_rng(rnorm, 2, summ$coefficients[, 1], summ$coefficients[, 2])
  eff <- rvar_rng(rnorm, length(conc), drop(dm %*% beta), summ$sigma)
  # Add a minimum value to keep e.g. effect sizes from going negative
  if (!is.null(min_val)) {
    eff <- draws_of(eff) |>
      pmax(eff_draws, min_val) |>
      rvar()
  }
  eff
}
