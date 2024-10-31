library(tidyverse)
library(posterior)
library(ggdist)

if (!dir.exists("figs/puyallup")) dir.create("figs/puyallup")

### Puyallup effects
puy_eff <- read_rds("data/puyallup/puy_pcb_eff_df.rds")

puy_eff |>
  ggplot(aes(xdist = eff, y = type)) +
  stat_histinterval() +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/puy_pcb_surv.pdf", width = 11, height = 8.5)

puy_eff |>
  ggplot(aes(xdist = oceanadults, y = type)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = oa0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of 1+ Chinook in the ocean",
    labels = scales::comma
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/puy_pcb_pop.pdf", width = 11, height = 8.5)

puy_eff |>
  ggplot(aes(xdist = oa_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in adult Chinook ocean abundance",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/puy_pcb_relpop.pdf", width = 11, height = 8.5)

### White effects
white_eff <- read_rds("data/puyallup/white_pcb_eff_df.rds")

white_eff |>
  ggplot(aes(xdist = eff, y = type)) +
  stat_histinterval() +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/white_pcb_surv.pdf", width = 11, height = 8.5)

white_eff |>
  ggplot(aes(xdist = oceanadults, y = type)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = oa0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of 1+ Chinook in the ocean",
    labels = scales::comma
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/white_pcb_pop.pdf", width = 11, height = 8.5)

white_eff |>
  ggplot(aes(xdist = oa_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in adult Chinook ocean abundance",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/white_pcb_relpop.pdf", width = 11, height = 8.5)

### Combined effects
combo_eff <- tibble(
  type = puy_eff$type,
  eff = puy_eff$eff,
  oceanadults = puy_eff$oceanadults + white_eff$oceanadults,
  oa0 = puy_eff$oa0 + white_eff$oa0,
  oa_change = (oceanadults - oa0) / oa0
)

combo_eff |>
  ggplot(aes(xdist = eff, y = type)) +
  stat_histinterval() +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/combo_pcb_surv.pdf", width = 11, height = 8.5)

combo_eff |>
  ggplot(aes(xdist = oceanadults, y = type)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = oa0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of 1+ Chinook in the ocean",
    labels = scales::comma
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/combo_pcb_pop.pdf", width = 11, height = 8.5)

combo_eff |>
  ggplot(aes(xdist = oa_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in adult Chinook ocean abundance",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/combo_pcb_relpop.pdf", width = 11, height = 8.5)
