library(tidyverse)
library(posterior)
library(ggdist)
library(patchwork)

if (!dir.exists("figs/stillaguamish")) dir.create("figs/stillaguamish")

### Puyallup exposures ----------------------------------------------------
stilly_exp <- read_rds(here::here("data", "stillaguamish", "pbde_exposure.rds")) |>
  map(thin_draws, 1)

stilly_pbde <- tibble(pbde = seq(0, 40, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pbde,
    stilly_exp$pop_meanlog,
    stilly_exp$pop_sdlog)) |>
  curve_interval(popdens)
stilly_pdf_plt <- stilly_pbde |>
  ggplot(aes(x = pbde)) +
  geom_vline(xintercept = 6.6, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(
    title = "Stillaguamish PBDE Exposure",
    x = "PBDE concentration (ng/g wet weight)"
  ) +
  # theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
stilly_pdf_plt
ggsave(
  here::here("figs", "stillaguamish", "pbde_exposure_pdf.png"),
  stilly_pdf_plt
)

stilly_cdf <- tibble(
  pbde = seq(0, 40, length.out = 1025)
) |>
  mutate(
    popdens = rfun(plnorm)(
      pbde,
      stilly_exp$pop_meanlog,
      stilly_exp$pop_sdlog,
      lower.tail = FALSE),
  ) |>
  curve_interval(popdens)
stilly_ccdf_plt <- stilly_cdf |>
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
    title = "Stillaguamish PBDE Exposure",
    x = "PBDE concentration (ng/g wet weight)",
    y = "Proportion with maximum exposure"
  ) +
  # theme_minimal() +
  theme(
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank()
  )
stilly_ccdf_plt
ggsave(
  here::here("figs", "stillaguamish", "pbde_exposure_ccdf.png"),
  stilly_ccdf_plt
)

### Stillaguamish effects ------------------------------------------------------
stilly_eff <- read_rds("data/stillaguamish/pbde_eff_df.rds")

stilly_eff |>
  ggplot(aes(xdist = eff)) +
  stat_histinterval() +
  geom_vline(xintercept = 1.00, linetype = "dashed") +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(axis.title.y = element_blank())
# ggsave("figs/stillaguamish/stilly_pbde_surv.pdf", width = 11, height = 8.5)
ggsave("figs/stillaguamish/stilly_pbde_surv.png", width = 11, height = 8.5)

stilly_eff |>
  ggplot(aes(xdist = spawners)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = sp0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of spawners",
    labels = scales::comma
  ) +
  # ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(axis.title.y = element_blank())
# ggsave("figs/stillaguamish/stilly_pbde_pop.pdf", width = 11, height = 8.5)
ggsave("figs/stillaguamish/stilly_pbde_pop.png", width = 11, height = 8.5)

stilly_eff |>
  ggplot(aes(xdist = sp_change)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal() +
  theme(axis.title.y = element_blank())
# ggsave("figs/stillaguamish/stilly_pbde_relpop.pdf", width = 11, height = 8.5)
ggsave("figs/stillaguamish/stilly_pbde_relpop.png", width = 11, height = 8.5)
