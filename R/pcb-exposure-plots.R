library(tidyverse)
library(posterior)
library(ggdist)
library(patchwork)

### Puyallup exposures ----------------------------------------------------
puy_exp <- read_rds("data/puyallup/pcb_exposure.rds") |>
  map(thin_draws, 1)

puy_ww <- tibble(pcb = seq(0, 0.3, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pcb,
    puy_exp$ww$pop_meanlog,
    puy_exp$ww$pop_sdlog)) |>
  curve_interval(popdens) |>
  mutate(river = "Puyallup/White")

puy_ww_plt <- puy_ww |>
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


## Stillaguamish exposures ----------------------------------------------------
stilly_exp <- read_rds("data/stillaguamish/pcb_exposure.rds") |>
  map(thin_draws, 1)

stilly_ww <- tibble(pcb = seq(0, 0.3, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pcb,
    stilly_exp$ww$pop_meanlog,
    stilly_exp$ww$pop_sdlog)) |>
  curve_interval(popdens) |>
  mutate(river = "Stillaguamish")

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

stilly_ww_plt / puy_ww_plt

bind_rows(stilly_ww, puy_ww) |>
  mutate(river = factor(river, levels = c("Stillaguamish", "Puyallup/White"))) |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3, fill = "#619CCF") +
  geom_line(aes(y = popdens), color = "#619CCF") +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g wet weight)") +
  facet_wrap(~river, ncol = 1) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
ggsave(here::here("figs", "pcb-exposure-ug.png"), height = 5.5, width = 3)

bind_rows(stilly_ww, puy_ww) |>
  mutate(
    river = factor(river, levels = c("Stillaguamish", "Puyallup/White")),
    pcb = pcb * 1e3,
  ) |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 100, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3, fill = "#619CCF") +
  geom_line(aes(y = popdens), color = "#619CCF") +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (ng/g wet weight)") +
  facet_wrap(~river, ncol = 1) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
ggsave(here::here("figs", "pcb-exposure-ng.png"), height = 5.5, width = 3)
