library(tidyverse)
library(posterior)
library(distributional)
library(ggdist)
library(patchwork)

pop_mu <- 0
pop_sigma <- 1
nsamp <- 10e3 # Number of samples in Monte Carlo approximation

meanlnorm_pars <- function(n, pop_mu, pop_sig) {
  pop_sig2 <- pop_sig^2
  mean_sig <- sqrt(log(1 / n * (exp(pop_sig2) - 1) + 1))
  mean_mu <- pop_mu + pop_sig / 2 - mean_sig / 2
  c(meanlog = mean_mu, meansd = mean_sig)
}

rng_meanlnorm <- function(n, pop_mu, pop_sigma, nsamp = 4000) {
  replicate(nsamp, mean(rlnorm(n, pop_mu, pop_sigma))) |>
    rvar()
}

comp_df <- tibble(
  n = 1:25,
  approx_mu = map_dbl(n, \(n) meanlnorm_pars(n, pop_mu, pop_sigma)[1]),
  approx_sigma = map_dbl(n, \(n) meanlnorm_pars(n, pop_mu, pop_sigma)[2]),
  approx_dist = dist_lognormal(approx_mu, approx_sigma),
  mc_dist = map_vec(n, \(n) rng_meanlnorm(n, pop_mu, pop_sigma, nsamp))
)

ggplot(comp_df) +
  stat_slab(aes(xdist = approx_dist, color = n), fill = NA) +
  stat_slab(aes(xdist = mc_dist, color = n), linetype = "dashed", fill = NA) +
  facet_wrap(~n) +
  # scale_x_continuous(limits = c(0, 8), expand = FALSE) +
  scale_y_continuous(limits = c(0, NA), expand = FALSE) +
  coord_cartesian(xlim = c(0, 8), expand = FALSE) +
  scale_thickness_shared()

puy_pcb_obs <- bind_rows(
  read_rds(here::here("data", "puyallup", "pcb_estns.rds")) |>
    mutate(river = str_to_title(river)),
  read_rds(here::here("data", "stillaguamish", "pcb_estns.rds"))
)

ggplot(
  puy_pcb_obs,
  aes(x = pcb_ug_ww, fill = factor(n_composite))
) +
  stat_dots(position = position_dodge(), layout = "bin")

ggplot(
  puy_pcb_obs,
  aes(x = pcb_ug_ww)
) +
  geom_histogram(
    position = position_stack(reverse = TRUE),
    binwidth = 0.01,
    boundary = 0
  ) +
  scale_x_continuous(
    name = "μg PCB per g wet weight",
    limits = c(0, NA),
    expand = expansion(c(0, 0.025))
  ) +
  scale_y_continuous(
    name = "Number of observations",
    limits = c(0, NA),
    expand = expansion(c(0, 0.025))
  ) +
  facet_wrap(~river, nrow = 1) +
  theme_bw()
ggsave(
  here::here("docs", "figs", "ww_exposure_obs.png"),
  width = 12,
  height = 5
)

ggplot(
  puy_pcb_obs,
  aes(x = pcb_ug_ww, fill = factor(n_composite))
) +
  geom_histogram(
    position = position_stack(reverse = TRUE),
    binwidth = 0.01,
    boundary = 0
  ) +
  scale_x_continuous(
    name = "μg PCB per g wet weight",
    limits = c(0, NA),
    expand = expansion(c(0, 0.025))
  ) +
  scale_y_continuous(
    name = "Number of observations",
    limits = c(0, NA),
    expand = expansion(c(0, 0.025))
  ) +
  scale_fill_viridis_d(name = "# Composite", option = "E") +
  facet_wrap(~river, nrow = 1) +
  theme_bw()
ggsave(
  here::here("docs", "figs", "ww_exposure_comps.png"),
  width = 12,
  height = 5
)
