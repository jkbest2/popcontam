## ## From Beechie et al. 2023 (Stillaguamish and Snohomish HARP model description), page 39
## beverton_holt <- function(n, prod, cap) {
##    n * prod / (1 + n * prod / cap)
## }

## get_spawners <- function(pop) {
##   attr(pop, "spawners")
## }

stilly_sim <- function(pop, nearshore_surv_adj = 1) {
  ## Contributions to prespawners from each age class
  prespawners <- 0.05 * pop[2] + 0.34 * pop[3] + 0.92 * pop[4] + pop[5]
  ## Prespawn mortality
  spawners <- 0.696937 * prespawners
  ## Egg production. A fecundity of 5400 is per female, and we assume (as in
  ## HARP) that returners have an even sex ratio.
  eggs <- beverton_holt(spawners / 2, 5400, 152859370)
  ## Survival to fry
  fry <- eggs * 0.548651
  ## Colonization is density-dependent. Pre-parr are fish that will choose to
  ## remain in the river to rear. Capacity is doubled to allow for two cohorts
  ## of fish, non-simultaneously occupying habitats.
  pre_parr <- beverton_holt(fry, 0.302312^(1/9), 1682111.1 * 2)
  ## Parr migrants rear for about nine weeks before moving downstream based on
  ## Snohomish screw trap data. HARP adjust this to impose a late June
  ## temperature penalty
  parr_mig <- beverton_holt(pre_parr, 0.302312^(8/9), 1682111.1 * 2)
  ## Migration "choice" occurs after the first week, so is based on the number
  ## of pre-parr. See also line 50 of lcm-chinook.R, where natal_fry are
  ## subtracted from surviving pre_fry. Migration survival (0.128, estimated) is
  ## applied at the same time here. This cannot be negative because the fry
  ## survival here is equal to the maximum possible pre_parr survival above.
  fry_mig <- (fry * 0.302312^(1/9) - pre_parr) * 0.128
  ## Transition to first marine year
  delta_fry <- beverton_holt(fry_mig, 0.35, 45335)
  ## Note that this age structure shouldn't be taken literally - we apply all
  ## mortality effects in the nearshore, when we know that part of it is
  ## attributable to harvest.
  age00 <- (delta_fry + parr_mig) * 0.01645436 * nearshore_surv_adj
  age01 <- pop[1] * 0.6
  age02 <- pop[2] * (1 - 0.05) * 0.7
  age03 <- pop[3] * (1 - 0.34) * 0.8
  age04 <- pop[4] * (1 - 0.92) * 0.9

  structure(c(age00, age01, age02, age03, age04, use.names = FALSE),
            names = paste0("age0", 0:4),
            oceanadults = unname(age01 + age02 + age03 + age04),
            prespawners = unname(prespawners),
            spawners = unname(spawners),
            eggs = unname(eggs),
            fry = unname(fry),
            parr_mig = unname(parr_mig),
            fry_mig = unname(fry_mig),
            delta_fry = unname(delta_fry))
}

## growth_morts <- 0.12

### Stillaguamish no effects
source("R/utils.R")
source("R/fish-size.R")


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

db2011_survival <- function(mass) {
  ## The Duffy and Beauchamp survival regression predicts percent survival. I
  ## prefer to work with survival rates, which are the percent divided by 100.
  ## Because log10(0.01) = -2, I included an extra -2 in the intercept of the
  ## regression. This is why the Duffy and Beauchamp intercept is -1.071 but
  ## we're using -3.071 here. Given that we're using relative changes in this
  ## survival it won't actually make a difference.
  l <- -3.071 + 0.041 * mass
  10^l
}

inv_db2011_survival <- function(survival) {
  ls <- log(survival, 10)
  (ls + 3.071) / 0.041
}

dexposure_fun <- function(means, sdlogs, props) {
  function(xs) {
    vapply(xs, \(x) sum(props * dlnorm(x, means, sdlogs)),
           0.0)
  }
}

mort_qreg <- function(pcb) {
  ## Effect threshold is 100 ng/g (ww), so return zero effect below this
  ifelse(
    pcb < 0.100,
    0,
    pmax(0.1702 + 0.221 * log10(pcb), 0))
}

growth_qreg <- function(pcb) {
  ifelse(
    pcb < 0.100,
    0,
    pmax(0.15 + 0.0938 * log10(pcb), 0)
  )
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
  facet_wrap(~ aop, ncol = 1) +
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
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("figs/example-exposure-plot.png",
       width = 6.5, height = 4.5)

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

growth_to_surv <- function(pcb, base_size = db_size$july_mass, noexp_surv = db_size$pred_surv) {
  ## Calculate the expected reduction in size given exposure
  gred <- growth_qreg(pcb)
  exp_mass <- (1 - gred) * base_size
  db2011_survival(exp_mass) / noexp_surv
}

growth_eff <- integrate(\(pcb) dexposure(pcb) * growth_to_surv(pcb), lower = 0, upper = Inf)$value
stilly_growth <- eq_pop(stilly_sim, nearshore_surv_adj = growth_eff, pop0 = rep(1000, 5))
oa_growth <- get_oceanadults(stilly_growth)
red_growth <- 1 - (oa_0 - oa_growth) / oa_0
