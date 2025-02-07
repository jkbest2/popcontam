library(tidyverse)
library(posterior)
library(ggdist)
library(patchwork)

if (!dir.exists("figs/puyallup")) dir.create("figs/puyallup")

### Puyallup exposures ----------------------------------------------------
puy_exp <- read_rds(here::here("data", "puyallup", "pbde_exposure.rds")) |>
  map(thin_draws, 1)

puy_ww <- tibble(pcb = seq(0, 40, length.out = 1025)) |>
  mutate(popdens = rfun(dlnorm)(
    pcb,
    puy_exp$pop_meanlog,
    puy_exp$pop_sdlog)) |>
  curve_interval(popdens)
puy_ww_plt <- puy_ww |>
  ggplot(aes(x = pcb)) +
  # geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line(aes(y = popdens)) +
  scale_x_continuous(expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.02))) +
  labs(
    title = "Puyallup PBDE Exposure",
    x = "PBDE concentration (ng/g wet weight)"
  ) +
  # theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
puy_ww_plt
ggsave(
  here::here("figs", "puyallup", "pbde_exposure_pdf.png"),
  puy_ww_plt
)

puy_ccdf <- tibble(
  pcb = seq(0, 40, length.out = 1025)
) |>
  mutate(
    popdens = rfun(plnorm)(
      pcb,
      puy_exp$pop_meanlog,
      puy_exp$pop_sdlog,
      lower.tail = FALSE),
  ) |>
  curve_interval(popdens)
puy_ccdf_plt <- puy_ccdf |>
  ggplot(aes(x = pcb)) +
  # geom_vline(xintercept = 0.1, linetype = "dashed") +
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
  # theme_minimal() +
  theme(
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank()
  )
puy_ccdf_plt
ggsave(
  here::here("figs", "puyallup", "pbde_exposure_ccdf.png"),
  puy_ccdf_plt
)
