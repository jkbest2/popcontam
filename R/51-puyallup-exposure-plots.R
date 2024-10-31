library(tidyverse)
library(posterior)
library(ggdist)

# ## This plot shows that the gamma model is not accounting for enough low values;
# ## this also means that it may be
# ppc_dens_overlay(data_ln$conc, as_draws_matrix(post_ln$conc_gen)[1:50, ]) +
#   xlim(c(0, 0.25))
# ppc_dens_overlay(data_ln$conc, as_draws_matrix(post_ln$conc_gen2)[1:50, ]) +
#   xlim(c(0, 0.25))

# ## This seems to show that the PIT residuals don't fall outside the expected
# ## envelope
# ppc_pit_ecdf(data_ln$conc, as_draws_matrix(post_ln$conc_gen)) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed")
# ppc_pit_ecdf(data_ln$conc, as_draws_matrix(post_ln$conc_gen2)) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed")

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
