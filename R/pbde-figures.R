library(tidyverse)

### PBDE dose-response in terms of mortality ==================================
mort_df <- read_rds("pbde_mort_df.rds")

mort_df |>
  ggplot(aes(x = conc, y = mort, ymin = .lower, ymax = .upper)) +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  scale_x_continuous(
    breaks = seq(0, 250, 50),
    expand = expansion(c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.1),
    labels = scales::percent,
    limits = c(0, 1), expand = expansion(c(0, 0))
  ) +
  labs(
    x = "PBDE Concentration (ng/g wet weight)",
    y = "Predicted mortality"
  ) +
  theme_minimal()

ggsave("pbde_mort.png", width = 12, height = 8)

### PBDE Change in spawners ===================================================
## Puyallup -------------------------------------------------------------------
puy_eff <- read_rds("puy_pbde_eff_df.rds")

puy_eff |>
  ggplot(aes(xdist = sp_change)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
ggsave("puy_pbde_relpop.png", width = 12, height = 8)

## White ----------------------------------------------------------------------
white_eff <- read_rds("white_pbde_eff_df.rds")

white_eff |>
  ggplot(aes(xdist = sp_change)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
ggsave("white_pbde_relpop.png", width = 12, height = 8)

## Combined Puyallup/White ----------------------------------------------------
combo_eff <- tibble(
  type = puy_eff$type,
  eff = puy_eff$eff,
  spawners = puy_eff$spawners + white_eff$spawners,
  sp0 = puy_eff$sp0 + white_eff$sp0,
  sp_change = (spawners - sp0) / sp0
)

combo_eff |>
  ggplot(aes(xdist = sp_change)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
ggsave("combo_pbde_relpop.png", width = 12, height = 8)

## Stillaguamish --------------------------------------------------------------
stilly_eff <- read_rds("pbde_eff_df.rds")

stilly_eff |>
  ggplot(aes(xdist = sp_change)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  theme_minimal() +
  theme(axis.title.y = element_blank())
ggsave("stilly_pbde_relpop.png", width = 12, height = 8)
