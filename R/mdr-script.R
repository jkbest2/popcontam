library(tidyverse)

## growth_morts <- 0.12

### Stillaguamish no effects
source("R/utils.R")
source("R/fish-size.R")
source("R/pcb-effects.R")

stilly_0 <- eq_pop(stilly_sim, pop0 = rep(1000, 5))
oa_0 <- get_oceanadults(stilly_0)

## Calculate smolt-to-adult ratio
smolts <- attr(stilly_0, "parr_mig") + attr(stilly_0, "fry_mig")
sar_0 <- get_spawners(stilly_0) / smolts
## Back-calculate July size from Duffy & Beauchamp 2011
inv_db2011_survival(sar_0)

## Now look at just parr migrants. Assuming that smolts are counted/estimated
## above the estuary, so the total number is parr_mig + fry_mig. All parr
## migrant stages after this are density-independent, so to look at parr to
## adult survival we can look at the proportion of parr migrants out of parr
## migrants plus delta fry
## prop_parrmig <- attr(stilly_0, "parr_mig") / (attr(stilly_0, "parr_mig") + attr(stilly_0, "delta_fry"))
## parr_spawners <- prop_parrmig * get_spawners(stilly_0)
## parr_sar <- parr_spawners / attr(stilly_0, "parr_mig")
## parr_size <- inv_db2011_survival(parr_sar)

## prop_frymig <- 1 - prop_parrmig
## fry_spawners <- prop_frymig * get_spawners(stilly_0)
## fry_sar <- fry_spawners / attr(stilly_0, "fry_mig")
## fry_size <- inv_db2011_survival(fry_sar)


dexposure_fun <- function(means, sdlogs, props) {
  function(xs) {
    vapply(
      xs, \(x) sum(props * dlnorm(x, means, sdlogs)),
      0.0
    )
  }
}

## Effect plots
tibble(
  ## pcb_ng = seq(0, 1e6, length.out = 2049),
  pcb_ng = 10^seq(1, 6, length.out = 1025),
  pcb_ug = pcb_ng / 1000,
  `Direct Mortality` = 0.1702 + 0.221 * log10(pcb_ug),
  `Growth Restriction` = 0.15 + 0.0938 * log10(pcb_ug),
  noeff = ifelse(pcb_ug >= 1, 1, 0.3)
) |>
  pivot_longer(cols = c(`Growth Restriction`, `Direct Mortality`), names_to = "aop", values_to = "effect") |>
  ggplot(aes(x = pcb_ng, y = effect)) +
  annotate("ribbon", x = c(0, 100), ymin = 0, ymax = Inf, fill = "gray50", alpha = 0.5) +
  geom_line() +
  labs(x = "Tissue PCB Concentration (ng/g ww)", y = "Effect") +
  scale_x_log10(limits = c(10, NA), labels = scales::comma, expand = expansion()) +
  scale_y_continuous(limits = c(0, NA), labels = scales::percent, expand = expansion()) +
  facet_wrap(~aop, ncol = 1) +
  theme_bw()
ggsave("figs/db2011-quantile-regressions.png", width = 6.5, height = 4.5)

dexposure <- dexposure_fun(c(log(0.100), log(0.350)), c(0.75, 0.5), c(0.4, 0.6))
curve(dexposure(x), from = 0.000, to = 1.000, n = 1025)
abline(v = 0.100, lty = "dashed")

exp_df <- tibble(
  pcb = seq(0, 1025, length.out = 1025),
  dens = dexposure(pcb / 1000)
)
exp_df |>
  ggplot(aes(x = pcb, y = dens)) +
  annotate("ribbon", x = c(0, 100), ymin = 0, ymax = Inf, fill = "gray50", alpha = 0.5) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 1000, 200), expand = expansion(c(0, 0.01))) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(x = "PCB Concentration ng/g (ww)", y = "Population exposure density") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave("figs/example-exposure-plot.png",
  width = 6.5, height = 4.5
)

## Calculate direct mortality rate based on the exposure distribution and the
## Berninger and Tillitt direct mortality relationship
morts <- integrate(\(pcb) dexposure(pcb) * mort_qreg(pcb), lower = 0, upper = Inf)$value
## Apply direct mortality to the to-ocean 0.0 transition (decreasing survival)
stilly_mort <- eq_pop(stilly_sim, nearshore_surv_adj = 1 - morts, pop0 = rep(1000, 5))
spawner_mort <- get_spawners(stilly_mort)
oa_mort <- get_oceanadults(stilly_mort)
red_mort <- 1 - (oa_0 - oa_mort) / oa_0

## Alternatively, *remove* direct mortalities from the system (increasing survival)
stilly_imort <- eq_pop(stilly_sim, nearshore_surv_adj = 1 / (1 - morts), pop0 = rep(1000, 5))
oa_imort <- get_oceanadults(stilly_imort)
red_imort <- 1 - (oa_imort - oa_0) / oa_imort

## Growth reduction
growth_reduction <- integrate(\(pcb) dexposure(pcb) * growth_qreg(pcb), lower = 0, upper = Inf)$value
base_size <- db_size$july_mass
noexp_surv <- db_size$pred_surv
## exp_surv <- integrate(\(pcb) dexposure(pcb) * db2011_survival(gro))

growth_eff <- integrate(\(pcb) dexposure(pcb) * growth_to_surv(pcb), lower = 0, upper = Inf)$value
stilly_growth <- eq_pop(stilly_sim, nearshore_surv_adj = growth_eff, pop0 = rep(1000, 5))
oa_growth <- get_oceanadults(stilly_growth)
red_growth <- 1 - (oa_0 - oa_growth) / oa_0
