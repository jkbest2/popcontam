library(tidyverse)
library(rstan)
library(posterior)
library(patchwork)
library(ggdist)

chk_pcb <- read_xlsx("data/CB_Stilly_Puyallup_monitoring_data.xlsx") |>
  separate(Species, c("lifestage", "species")) |>
  mutate(
    river = str_to_title(RiverSystem),
    unit = "ng/g (wet)"
  ) |>
  select(river, species, lifestage, pcbs = `Max of CONC_FOUND`, unit) |>
  mutate(gt_lla = pcbs >= 100) |>
  filter(lifestage == "juvenile")

## PCB exposure observations
chk_pcb <- bind_rows(
  read_rds(here::here("data/puyallup/pcb_estns.rds")),
  read_rds(here::here("data/stillaguamish/pcb_estns.rds"))
) |>
  mutate(
    river = str_to_title(river)
  )

ggplot(chk_pcb, aes(x = pcb_ug_ww)) +
  geom_histogram(binwidth = 0.01, boundary = 0) +
  facet_wrap(~river) +
  labs(x = "PCB concentration (μg/g ww)", y = "Count") +
  theme_bw()

ggsave(
  here::here("figs/presentation/pcb-obs.png"),
  width = 2500,
  height = 1000,
  units = "px"
)

## Population exposures
puy_exp <- read_rds("data/puyallup/pcb_exposure.rds") |>
  map(thin_draws, 1)

puy_ww <- tibble(pcb = seq(0, 0.3, length.out = 1025)) |>
  mutate(
    popdens = rfun(dlnorm)(
      pcb,
      puy_exp$ww$pop_meanlog,
      puy_exp$ww$pop_sdlog
    )
  ) |>
  curve_interval(popdens)
puy_ww_plt <- puy_ww |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g wet weight)") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

puy_lw <- tibble(pcb = seq(0, 6.6, length.out = 1025)) |>
  mutate(
    popdens = rfun(dlnorm)(
      pcb,
      puy_exp$lw$pop_meanlog,
      puy_exp$lw$pop_sdlog
    )
  ) |>
  curve_interval(popdens)
puy_lw_plt <- puy_lw |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 2.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g observed lipid weight)") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

puy_lw1 <- tibble(pcb = seq(0, 6.6, length.out = 1025)) |>
  mutate(
    popdens = rfun(dlnorm)(
      pcb,
      puy_exp$lw1$pop_meanlog,
      puy_exp$lw1$pop_sdlog
    )
  ) |>
  curve_interval(popdens)
puy_lw1_plt <- puy_lw1 |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 2.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g 1% lipid weight)") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

puy_ww_plt / puy_lw_plt / puy_lw1_plt
ggsave(
  here::here("figs/presentation/puy-pop-exposure.png"),
  width = 2500,
  height = 1000,
  units = "px"
)

stilly_exp <- read_rds("data/stillaguamish/pcb_exposure.rds") |>
  map(thin_draws, 1)

stilly_ww <- tibble(pcb = seq(0, 0.3, length.out = 1025)) |>
  mutate(
    popdens = rfun(dlnorm)(
      pcb,
      stilly_exp$ww$pop_meanlog,
      stilly_exp$ww$pop_sdlog
    )
  ) |>
  curve_interval(popdens)
stilly_ww_plt <- stilly_ww |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g wet weight)") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

stilly_lw <- tibble(pcb = seq(0, 6.6, length.out = 1025)) |>
  mutate(
    popdens = rfun(dlnorm)(
      pcb,
      stilly_exp$lw$pop_meanlog,
      stilly_exp$lw$pop_sdlog
    )
  ) |>
  curve_interval(popdens)
stilly_lw_plt <- stilly_lw |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 2.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g observed lipid weight)") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

stilly_lw1 <- tibble(pcb = seq(0, 6.6, length.out = 1025)) |>
  mutate(
    popdens = rfun(dlnorm)(
      pcb,
      stilly_exp$lw1$pop_meanlog,
      stilly_exp$lw1$pop_sdlog
    )
  ) |>
  curve_interval(popdens)
stilly_lw1_plt <- stilly_lw1 |>
  ggplot(aes(x = pcb)) +
  geom_vline(xintercept = 2.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(x = "PCB concentration (µg/g 1% lipid weight)") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

stilly_ww_plt / stilly_lw_plt / stilly_lw1_plt

ggsave(
  here::here("figs/presentation/still-pop-exposure.png"),
  width = 2500,
  height = 1000,
  units = "px"
)


## Demographic effects
source(here::here("R/fish-size.R")) # Get data from D&B2011

direct_mortality <- function(pcb_ug) {
  pmax(0.1894 + 0.2115 * log10(pcb_ug), 0)
}
## FIXME This is *only* the growth effect, *does not include* the translation to survival yet!
growth_restriction <- function(pcb_ug) {
  pmax(0.1676 + 0.07580 * log10(pcb_ug), 0)
}
db2011_survival <- function(mass) {
  l <- -3.071 + 0.041 * mass
  10^l
}
growth_to_mort <- function(
  pcb,
  base_size = db_size$july_mass,
  noexp_surv = db_size$pred_surv
) {
  ## Calculate the expected reduction in size given exposure
  gred <- growth_restriction(pcb)
  exp_mass <- (1 - gred) * base_size
  1 - db2011_survival(exp_mass) / noexp_surv
}

## Example combined effects plot
combo_eff_df <- tibble(
  pcb_ng = seq(0, 1.025e3, length.out = 2049),
  ## pcb_ng = 10^seq(1, 6, length.out = 1025),
  pcb_ug = pcb_ng / 1000,
  `Direct Mortality` = direct_mortality(pcb_ug),
  `Growth Restriction` = growth_to_mort(pcb_ug),
  `Combined Naive` = `Direct Mortality` + `Growth Restriction`,
  `Combined Conditional` = `Direct Mortality` +
    (1 - `Direct Mortality`) * `Growth Restriction`
)

combo_ex <- tibble(
  pcb_ng = 500,
  pcb_ug = pcb_ng / 1000,
  `Direct Mortality` = direct_mortality(pcb_ug),
  `Growth Restriction` = growth_to_mort(pcb_ug),
  `Combined Naive` = `Direct Mortality` + `Growth Restriction`,
  `Combined Conditional` = `Direct Mortality` +
    (1 - `Direct Mortality`) * `Growth Restriction`
)

combo_eff_df |>
  pivot_longer(
    cols = c(`Direct Mortality`, `Combined Naive`, `Combined Conditional`),
    names_to = "aop",
    values_to = "effect"
  ) |>
  mutate(naive = ifelse(grepl("Naive", aop), TRUE, FALSE)) |>
  ggplot(aes(x = pcb_ng, y = effect, color = aop)) +
  annotate(
    "ribbon",
    x = c(0, 100),
    ymin = 0,
    ymax = Inf,
    fill = "gray50",
    alpha = 0.5
  ) +
  geom_line() +
  annotate(
    "errorbar",
    x = combo_ex$pcb_ng - 5,
    ymin = 0,
    ymax = combo_ex$`Direct Mortality`,
    width = 10
  ) +
  annotate(
    "text",
    x = combo_ex$pcb_ng - 15,
    y = combo_ex$`Direct Mortality` / 2,
    label = paste0(
      "Direct Mortality\n",
      scales::percent(combo_ex$`Direct Mortality`, 0.1)
    ),
    hjust = "right"
  ) +
  annotate(
    "errorbar",
    x = combo_ex$pcb_ng - 5,
    ymin = combo_ex$`Direct Mortality`,
    ymax = combo_ex$`Combined Naive`,
    width = 10
  ) +
  annotate(
    "text",
    x = combo_ex$pcb_ng - 15,
    y = mean(c(combo_ex$`Direct Mortality`, combo_ex$`Combined Naive`)),
    label = paste0(
      "Growth Restriction\n",
      scales::percent(combo_ex$`Growth Restriction`, 0.1)
    ),
    hjust = "right"
  ) +
  annotate(
    "errorbar",
    x = combo_ex$pcb_ng + 5,
    ymin = 0,
    ymax = combo_ex$`Combined Conditional`,
    width = 10
  ) +
  annotate(
    "text",
    x = combo_ex$pcb_ng + 15,
    y = combo_ex$`Combined Conditional` / 2 + 0.035,
    label = paste0(
      "Combined\n",
      scales::percent(combo_ex$`Combined Conditional`, 0.1),
      " = ",
      scales::percent(combo_ex$`Direct Mortality`, 0.1),
      " + (1 - ",
      scales::percent(combo_ex$`Direct Mortality`, 0.1),
      ") x ",
      scales::percent(combo_ex$`Growth Restriction`, 0.1)
    ),
    hjust = "left"
  ) +
  labs(x = "Tissue PCB Concentration (ng/g ww)", y = "Demographic effect") +
  #  title = "Growth-selective mortality occurs after direct mortality") +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = seq(0, 1000, 200),
    labels = scales::comma,
    expand = expansion()
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::percent,
    expand = expansion()
  ) +
  scale_color_discrete(guide = FALSE) +
  theme_bw()
ggsave(
  "figs/presentation/conditional-effects.png",
  width = 2500,
  height = 1000,
  units = "px"
)

### Results
## Puyallup-White effects
puy_eff <- read_rds("data/puyallup/puy_pcb_eff_df.rds")
white_eff <- read_rds("data/puyallup/white_pcb_eff_df.rds")
combo_eff <- tibble(
  type = puy_eff$type,
  eff = puy_eff$eff,
  spawners = puy_eff$spawners + white_eff$spawners,
  sp0 = puy_eff$sp0 + white_eff$sp0,
  sp_change = (spawners - sp0) / sp0
)

combo_eff |>
  ggplot(aes(xdist = sp_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  ylab("PCB concentration normalization") +
  theme_bw()
ggsave(
  "figs/presentation/combo-pcb-relpop.png",
  width = 2500,
  height = 1000,
  units = "px"
)

## Stillaguamish effects
stilly_eff <- read_rds("data/stillaguamish/stilly_pcb_eff_df.rds")

stilly_eff |>
  ggplot(aes(xdist = sp_change, y = type)) +
  stat_histinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Change in spawners",
    labels = scales::percent
  ) +
  coord_cartesian(xlim = c(0, 0.1)) +
  ylab("PCB concentration normalization") +
  theme_bw()
ggsave(
  "figs/presentation/stilly_pcb_relpop.png",
  width = 2500,
  height = 1000,
  units = "px"
)

## Mortality contributions ----------------------------------------------------
puy_dm_eff <- read_rds(
  here::here("data", "puyallup", "puy_pcb_dm_eff.rds")
) |>
  map_vec(~ pluck(.))
puy_gr_eff <- read_rds(
  here::here("data", "puyallup", "puy_pcb_gr_eff.rds")
) |>
  map_vec(~ pluck(.))
white_dm_eff <- read_rds(
  here::here("data", "puyallup", "white_pcb_dm_eff.rds")
) |>
  map_vec(~ pluck(.))
white_gr_eff <- read_rds(
  here::here("data", "puyallup", "white_pcb_gr_eff.rds")
) |>
  map_vec(~ pluck(.))
stilly_dm_eff <- read_rds(
  here::here("data", "stillaguamish", "pcb_dm_eff.rds")
) |>
  map_vec(~ pluck(.))
stilly_gr_eff <- read_rds(
  here::here("data", "stillaguamish", "pcb_gr_eff.rds")
) |>
  map_vec(~ pluck(.))


eff_df <- expand_grid(
  river = factor(
    c("Puyallup", "White", "Stillaguamish"),
    levels = c("Puyallup", "White", "Stillaguamish")
  ),
  eff_type = factor(
    c("Direct", "Growth"),
    levels = rev(c("Direct", "Growth", "Combined"))
  ),
  wt_type = factor(
    c("Wet Weight", "Lipid Weight", "1% Lipid Weight"),
    levels = c("Wet Weight", "Lipid Weight", "1% Lipid Weight")
  )
) |>
  mutate(
    eff = rvar(c(
      puy_dm_eff,
      puy_gr_eff,
      white_dm_eff,
      white_gr_eff,
      stilly_dm_eff,
      stilly_gr_eff
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
  guides(color = "none") +
  theme_bw()
ggsave(
  here::here("figs", "mort_sources.png"),
  width = 2500,
  height = 1000,
  units = "px"
)
