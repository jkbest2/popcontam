library(tidyverse)
library(posterior)
library(ggdist)
library(patchwork)

if (!dir.exists("figs/puyallup")) dir.create("figs/puyallup")

### Puyallup exposures ----------------------------------------------------
puy_exp <- read_rds(here::here("data", "puyallup", "pbde_exposure.rds")) |>
  map(thin_draws, 1)

puy_pbde <- tibble(pbde = seq(0, 40, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pbde,
    puy_exp$pop_meanlog,
    puy_exp$pop_sdlog)) |>
  curve_interval(popdens)
puy_pbde_plt <- puy_pbde |>
  ggplot(aes(x = pbde)) +
  geom_vline(xintercept = 6.6, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(
    title = "Puyallup PBDE Exposure",
    x = "PBDE concentration (ng/g wet weight)"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
puy_pbde_plt
ggsave(
  here::here("figs", "puyallup", "pbde_exposure_pdf.png"),
  puy_pbde_plt
)

puy_ccdf <- tibble(
  pbde = seq(0, 40, length.out = 1025)
) |>
  mutate(
    popdens = rfun(plnorm)(
      pbde,
      puy_exp$pop_meanlog,
      puy_exp$pop_sdlog,
      lower.tail = FALSE),
  ) |>
  curve_interval(popdens)
puy_ccdf_plt <- puy_ccdf |>
  ggplot(aes(x = pbde)) +
  geom_vline(xintercept = 6.6, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(
    expand = expansion(c(0, 0.02)),
    minor_breaks = seq(0, 1, 0.05),
    label = scales::percent
  ) +
  labs(
    title = "Puyallup PBDE Exposure",
    x = "PBDE concentration (ng/g wet weight)",
    y = "Proportion with maximum exposure"
  ) +
  theme_minimal() +
  theme(
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank()
  )
puy_ccdf_plt
ggsave(
  here::here("figs", "puyallup", "pbde_exposure_ccdf.png"),
  puy_ccdf_plt
)

### Puyallup effects ----------------------------------------------------------
puy_eff <- read_rds("data/puyallup/puy_pbde_eff_df.rds")

puy_eff |>
  ggplot(aes(xdist = eff)) +
  stat_histinterval() +
  geom_vline(xintercept = 1.00, linetype = "dashed") +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
# ggsave("figs/puyallup/puy_pbde_surv.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/puy_pbde_surv.png", width = 11, height = 8.5)

puy_eff |>
  ggplot(aes(xdist = spawners)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = sp0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of spawners",
    labels = scales::comma
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
# ggsave("figs/puyallup/puy_pbde_pop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/puy_pbde_pop.png", width = 11, height = 8.5)

puy_eff |>
  ggplot(aes(xdist = sp_change)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
# ggsave("figs/puyallup/puy_pbde_relpop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/puy_pbde_relpop.png", width = 11, height = 8.5)

### White effects -------------------------------------------------------------
white_eff <- read_rds("data/puyallup/white_pbde_eff_df.rds")

white_eff |>
  ggplot(aes(xdist = eff)) +
  stat_histinterval() +
  geom_vline(xintercept = 1.00, linetype = "dashed") +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
# ggsave("figs/puyallup/white_pbde_surv.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/white_pbde_surv.png", width = 11, height = 8.5)

white_eff |>
  ggplot(aes(xdist = spawners)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = sp0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of spawners",
    labels = scales::comma
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
# ggsave("figs/puyallup/white_pbde_pop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/white_pbde_pop.png", width = 11, height = 8.5)

white_eff |>
  ggplot(aes(xdist = sp_change)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
# ggsave("figs/puyallup/white_pbde_relpop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/white_pbde_relpop.png", width = 11, height = 8.5)

### Combined effects ----------------------------------------------------------
combo_eff <- tibble(
  type = puy_eff$type,
  eff = puy_eff$eff,
  spawners = puy_eff$spawners + white_eff$spawners,
  sp0 = puy_eff$sp0 + white_eff$sp0,
  sp_change = (spawners - sp0) / sp0
)

combo_eff |>
  ggplot(aes(xdist = eff)) +
  stat_histinterval() +
  geom_vline(xintercept = 1.00, linetype = "dashed") +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
# ggsave("figs/puyallup/combo_pbde_surv.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/combo_pbde_surv.png", width = 11, height = 8.5)

combo_eff |>
  ggplot(aes(xdist = spawners)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = sp0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of spawners",
    labels = scales::comma
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
# ggsave("figs/puyallup/combo_pbde_pop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/combo_pbde_pop.png", width = 11, height = 8.5)

combo_eff |>
  ggplot(aes(xdist = sp_change)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank()
  )
# ggsave("figs/puyallup/combo_pbde_relpop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/combo_pbde_relpop.png", width = 11, height = 8.5)
