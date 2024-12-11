library(tidyverse)
library(posterior)
library(ggdist)
library(patchwork)

if (!dir.exists("figs/stillaguamish")) dir.create("figs/stillaguamish")

### Stillaguamish exposures ----------------------------------------------------
stilly_exp <- read_rds("data/stillaguamish/pcb_exposure.rds") |>
  map(thin_draws, 1)

stilly_ww <- tibble(pcb = seq(0, 0.3, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pcb,
    stilly_exp$ww$pop_meanlog,
    stilly_exp$ww$pop_sdlog)) |>
  curve_interval(popdens)
stilly_ww_plt <- stilly_ww |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g wet weight)") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

stilly_lw <- tibble(pcb = seq(0, 6.6, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pcb,
    stilly_exp$lw$pop_meanlog,
    stilly_exp$lw$pop_sdlog)) |>
  curve_interval(popdens)
stilly_lw_plt <- stilly_lw |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 2.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g observed lipid weight)") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

stilly_lw1 <- tibble(pcb = seq(0, 6.6, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pcb,
    stilly_exp$lw1$pop_meanlog,
    stilly_exp$lw1$pop_sdlog)) |>
  curve_interval(popdens)
stilly_lw1_plt <- stilly_lw1 |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 2.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g 1% lipid weight)") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

stilly_ww_plt / stilly_lw_plt / stilly_lw1_plt
ggsave("figs/stillaguamish/stilly_pop_exposure.pdf", width = 11, height = 8.5)
ggsave("figs/stillaguamish/stilly_pop_exposure.png", width = 11, height = 8.5)

## Stillaguamish proportion affected ------------------------------------------
stilly_aff <- map2(
  stilly_exp, c(0.1, 2.2, 2.2),
  ~ rfun(plnorm)(.y, .x$pop_meanlog, .x$pop_sdlog, lower.tail = FALSE)
)

norm_levels <- c("Wet weight", "Observed lipids", "1% lipids")

stilly_aff_df <- tibble(
  type = factor(norm_levels, levels = rev(norm_levels)),
  prop_over = c(stilly_aff$ww, stilly_aff$lw, stilly_aff$lw1)
)

stilly_aff_plt <- ggplot(stilly_aff_df, aes(xdist = prop_over, y = type)) +
  stat_histinterval() +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Proportion of population over threshold", y = "Normalization") +
  theme_minimal()
ggsave("figs/stillaguamish/stilly_prop_affected.pdf", width = 11, height = 8.5)
ggsave("figs/stillaguamish/stilly_prop_affected.png", width = 11, height = 8.5)

### Stillaguamish effects ------------------------------------------------------
stilly_eff <- read_rds("data/stillaguamish/stilly_pcb_eff_df.rds")

stilly_eff |>
  ggplot(aes(xdist = eff, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 1.00, linetype = "dashed") +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/stillaguamish/stilly_pcb_surv.pdf", width = 11, height = 8.5)
ggsave("figs/stillaguamish/stilly_pcb_surv.png", width = 11, height = 8.5)

stilly_eff |>
  ggplot(aes(xdist = spawners, y = type)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = sp0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of spawners",
    labels = scales::comma
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/stillaguamish/stilly_pcb_pop.pdf", width = 11, height = 8.5)
ggsave("figs/stillaguamish/stilly_pcb_pop.png", width = 11, height = 8.5)

stilly_eff |>
  ggplot(aes(xdist = sp_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/stillaguamish/stilly_pcb_relpop.pdf", width = 11, height = 8.5)
ggsave("figs/stillaguamish/stilly_pcb_relpop.png", width = 11, height = 8.5)
