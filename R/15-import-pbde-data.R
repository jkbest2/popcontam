library(tidyverse)

pbde_surv <- read_csv(
  here::here("data", "pbde_survival.csv"),
  col_types = cols(
    strata = col_character(),
    time = col_double(),
    n.risk = col_double(),
    n.event = col_double(),
    surv = col_double(),
    std.err = col_double(),
    lower = col_double(),
    upper = col_double(),
    concentration = col_double()
  )
)

ggplot(pbde_surv, aes(x = concentration)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point(aes(y = surv)) +
  geom_line(aes(y = surv))

library(splines)

mod <- lm(
  surv ~ bs(concentration, df = 4, degree = 2),
  weights = 1 / (std.err),
  data = pbde_surv
)
pred <- tibble(
  concentration = seq(0, 180, length.out = 129)
) |>
  mutate(
    surv = predict(mod, newdata = tibble(concentration = concentration))
  )

plot(surv ~ concentration, data = pbde_surv, type = "b", ylim = c(0.5, 0.8))
lines(pred$concentration, pred$surv, col = "blue")

ggplot(pbde_surv, aes(x = concentration, y = surv)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point() +
  geom_line(data = pred) +
  geom_vline(xintercept = 27) +
  labs(y = "Survival", x = "Concentration") +
  scale_y_continuous(labels = scales::percent)
