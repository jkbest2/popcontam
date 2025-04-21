library(tidyverse)
library(posterior)
library(ggdist)
library(patchwork)

if (!dir.exists("figs/puyallup")) dir.create("figs/puyallup")

### Puyallup exposures ----------------------------------------------------
puy_exp <- read_rds("data/puyallup/pcb_exposure.rds") |>
  map(thin_draws, 1)

puy_ww <- tibble(pcb = seq(0, 0.3, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pcb,
    puy_exp$ww$pop_meanlog,
    puy_exp$ww$pop_sdlog)) 
    
    |>
  # curve_interval(popdens)
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

puy_lw <- tibble(pcb = seq(0, 6.6, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pcb,
    puy_exp$lw$pop_meanlog,
    puy_exp$lw$pop_sdlog)) |>
  curve_interval(popdens)
puy_lw_plt <- puy_lw |>
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

puy_lw1 <- tibble(pcb = seq(0, 6.6, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pcb,
    puy_exp$lw1$pop_meanlog,
    puy_exp$lw1$pop_sdlog)) |>
  curve_interval(popdens)
puy_lw1_plt <- puy_lw1 |>
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

puy_ww_plt / puy_lw_plt / puy_lw1_plt
ggsave("figs/puyallup/puy_pop_exposure.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/puy_pop_exposure.png", width = 11, height = 8.5)

## Stillaguamish proportion affected ------------------------------------------
puy_aff <- map2(
  puy_exp, c(0.1, 2.2, 2.2),
  ~ rfun(plnorm)(.y, .x$pop_meanlog, .x$pop_sdlog, lower.tail = FALSE)
)

norm_levels <- c("Wet weight", "Observed lipids", "1% lipids")

puy_aff_df <- tibble(
  type = factor(norm_levels, levels = rev(norm_levels)),
  prop_over = c(puy_aff$ww, puy_aff$lw, puy_aff$lw1)
)

puy_aff_plt <- ggplot(puy_aff_df, aes(xdist = prop_over, y = type)) +
  stat_histinterval() +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Proportion of population over threshold", y = "Normalization") +
  theme_minimal()
ggsave("figs/puyallup/puy_prop_affected.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/puy_prop_affected.png", width = 11, height = 8.5)


### Puyallup effects ----------------------------------------------------------
puy_eff <- read_rds("data/puyallup/puy_pcb_eff_df.rds")

puy_eff |>
  ggplot(aes(xdist = eff, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 1.00, linetype = "dashed") +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/puy_pcb_surv.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/puy_pcb_surv.png", width = 11, height = 8.5)

puy_eff |>
  ggplot(aes(xdist = spawners, y = type)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = sp0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of spawners",
    labels = scales::comma
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/puy_pcb_pop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/puy_pcb_pop.png", width = 11, height = 8.5)

puy_eff |>
  ggplot(aes(xdist = sp_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/puy_pcb_relpop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/puy_pcb_relpop.png", width = 11, height = 8.5)

### White effects -------------------------------------------------------------
white_eff <- read_rds("data/puyallup/white_pcb_eff_df.rds")

white_eff |>
  ggplot(aes(xdist = eff, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 1.00, linetype = "dashed") +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/white_pcb_surv.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/white_pcb_surv.png", width = 11, height = 8.5)

white_eff |>
  ggplot(aes(xdist = spawners, y = type)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = sp0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of spawners",
    labels = scales::comma
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/white_pcb_pop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/white_pcb_pop.png", width = 11, height = 8.5)

white_eff |>
  ggplot(aes(xdist = sp_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/white_pcb_relpop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/white_pcb_relpop.png", width = 11, height = 8.5)

### Combined effects ----------------------------------------------------------
combo_eff <- tibble(
  type = puy_eff$type,
  eff = puy_eff$eff,
  spawners = puy_eff$spawners + white_eff$spawners,
  sp0 = puy_eff$sp0 + white_eff$sp0,
  sp_change = (spawners - sp0) / sp0
)

combo_eff |>
  ggplot(aes(xdist = eff, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 1.00, linetype = "dashed") +
  scale_x_continuous(
    name = "Nearshore survival adjustment",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/combo_pcb_surv.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/combo_pcb_surv.png", width = 11, height = 8.5)

combo_eff |>
  ggplot(aes(xdist = spawners, y = type)) +
  stat_histinterval() +
  geom_vline(aes(xintercept = sp0), linetype = "dashed") +
  scale_x_continuous(
    name = "Number of spawners",
    labels = scales::comma
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/combo_pcb_pop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/combo_pcb_pop.png", width = 11, height = 8.5)

combo_eff |>
  ggplot(aes(xdist = sp_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_minimal()
ggsave("figs/puyallup/combo_pcb_relpop.pdf", width = 11, height = 8.5)
ggsave("figs/puyallup/combo_pcb_relpop.png", width = 11, height = 8.5)

## Mortality contributions ----------------------------------------------------
puy_dm_eff <- read_rds(
  here::here("data", "puyallup", "puy_pcb_dm_eff.rds")) |>
  map_vec(~ pluck(.))
puy_gr_eff <- read_rds(
  here::here("data", "puyallup", "puy_pcb_gr_eff.rds")) |>
  map_vec(~ pluck(.))
white_dm_eff <- read_rds(
  here::here("data", "puyallup", "white_pcb_dm_eff.rds")) |>
  map_vec(~ pluck(.))
white_gr_eff <- read_rds(
  here::here("data", "puyallup", "white_pcb_gr_eff.rds")) |>
  map_vec(~ pluck(.))
stilly_dm_eff <- read_rds(
  here::here("data", "stillaguamish", "pcb_dm_eff.rds")) |>
  map_vec(~ pluck(.))
stilly_gr_eff <- read_rds(
  here::here("data", "stillaguamish", "pcb_gr_eff.rds")) |>
  map_vec(~ pluck(.))


eff_df <- expand_grid(
  river = factor(c("Puyallup", "White", "Stillaguamish"),
    levels = c("Puyallup", "White", "Stillaguamish")
  ),
  eff_type = factor(
    c("Direct", "Growth"),
    levels = rev(c("Direct", "Growth", "Combined"))),
  wt_type = factor(
    c("Wet Weight", "Lipid Weight", "1% Lipid Weight"),
    levels = c("Wet Weight", "Lipid Weight", "1% Lipid Weight"))
) |>
  mutate(
    eff = rvar(c(
      puy_dm_eff, puy_gr_eff,
      white_dm_eff, white_gr_eff,
      stilly_dm_eff, stilly_gr_eff
    ))
  )

eff_df2 <- eff_df |>
  summarize(
    eff = rvar_sum(eff),
    .by = c(river, wt_type)
  ) |>
  mutate(
    eff_type = factor(
      "Combined",
      levels = rev(c("Direct", "Growth", "Combined"))
    )
  )

bind_rows(eff_df, eff_df2) |>
  ggplot(aes(xdist = eff, y = eff_type, color = eff_type)) +
  stat_slabinterval() +
  facet_grid(wt_type ~ river) +
  scale_x_continuous(
    name = "Mortality rate",
    labels = scales::percent
  ) +
  labs(
    y = "Mortality source",
  ) +
  guides(color = "none")
ggsave(here::here("figs", "mort_sources.png"), width = 7.5, height = 5)


# eff_df |>
#   left_join(
#     select(eff_df2, river, wt_type, combo_eff = eff),
#     by = join_by(river, wt_type)
#   ) |>
#   mutate(
#     prop_eff = eff / combo_eff
#   ) |>
#   ggplot(aes(xdist = prop_eff, y = eff_type, color = eff_type)) +
#   stat_slabinterval() +
#   facet_grid(wt_type ~ river) +
#   scale_x_continuous(
#     name = "Proportion of mortality",
#     labels = scales::percent
#   )


eff_df3 <- eff_df |>
  point_interval(eff) |>
  select(-.point, -.interval, -.width)
write_csv(data.frame(eff_df3), here::here("data", "pcb_effect_table.csv"))
  
eff_df3 |>
  mutate(
    across(
      c(eff, .lower, .upper),
      scales::percent
    )
  )