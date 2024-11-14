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
  labs(x = "PCB concentration (µg/g 1% lipid wet weight)") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

stilly_aff <- map2(
  stilly_exp, c(0.1, 2.2, 2.2),
  ~ rfun(plnorm)(.y, .x$pop_meanlog, .x$pop_sdlog, lower.tail = FALSE)
)

map(
  stilly_aff,
  function(aff) {
    as_draws_df(aff) |>
      ggplot(aes(x = x)) +
      geom_density() +
      scale_x_continuous(
        limits = c(0, 0.3),
        labels = scales::percent,
        expand = expansion(c(0, 0.02))) +
      labs(x = "Proportion over threshold") +
        theme_minimal() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank()
        )
      )
  }
)

norm_levels <- c("1% Lipid", "Observed Lipid", "Wet Weight")

map2(stilly_aff, rev(norm_levels), ~ as_draws_df(.x) |> mutate(norm = .y)) |>
  reduce(bind_rows) |>
  mutate(norm = factor(norm, levels = rev(norm_levels))) |>
  ggplot(aes(x = x, fill = norm, color = norm)) +
  geom_density(alpha = 0.3) +
  # facet_wrap(~norm, ncol = 1, scales = "free") +
  scale_x_continuous(
    labels = scales::percent,
    expand = expansion(c(0, 0.02))
  ) +
  labs(
    x = "Proportion over threshold",
    fill = "Normalization",
    color = "Normalization"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

stilly_aff_df <- tibble(
  type = factor(rev(norm_levels), levels = norm_levels),
  prop_over = c(stilly_aff$ww, stilly_aff$lw, stilly_aff$lw1)
)

stilly_aff_plt <- ggplot(stilly_aff_df, aes(xdist = prop_over, y = type)) +
  stat_halfeye() +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Proportion of population over threshold", y = "Normalization")

stilly_ww_plt / stilly_lw_plt / stilly_lw1_plt

stilly_ww |>
  filter(pcb >= 0.1) |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g wet weight") +
  theme_minimal()

### Stillaguamish effects ------------------------------------------------------
stilly_eff <- read_rds("data/stillaguamish/stilly_pcb_eff_df.rds")

stilly_eff |>
  ggplot(aes(xdist = eff, y = type)) +
  stat_histinterval() +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/stillaguamish/stilly_pcb_surv.pdf", width = 11, height = 8.5)

stilly_eff |>
  ggplot(aes(xdist = oceanadults, y = type)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = oa0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of 1+ Chinook in the ocean",
    labels = scales::comma
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/stillaguamish/stilly_pcb_pop.pdf", width = 11, height = 8.5)

stilly_eff |>
  ggplot(aes(xdist = oa_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in adult Chinook ocean abundance",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/stillaguamish/stilly_pcb_relpop.pdf", width = 11, height = 8.5)
