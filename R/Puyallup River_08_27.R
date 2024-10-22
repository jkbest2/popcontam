beverton_holt <- function(n, prod, cap) {
  n * prod / (1 + n * prod / cap)
}




library(tidyverse)

source("R/utils.R")
source("R/fish-size.R")

puyallup_0 <- eq_pop(puyallup_sim, pop0 = rep(1000, 5))

puyallup_0

spawners <- get_spawners(puyallup_0)

## Smolt to adult returns
smolts <- attr(puyallup_0, "parr_mig") + attr(puyallup_0, "fry_mig")
sar_0 <- get_spawners(puyallup_0) / smolts
sar_0
smolts
get_eggs <- function(pop) {
  attr(pop, "eggs")
}
## Calculate freshwater survival
freshsurv <- smolts / get_eggs(puyallup_0)
freshsurv

## Ocean adults
oa_0 <- get_oceanadults(puyallup_0)

## Berninger and Tillitt regressions

mort_qreg <- function(pcb) {
  ## Effect threshold is 100 ng/g (ww), so return zero effect below this
  ifelse(
    pcb < 0.100,
    0,
    pmax(0.1702 + 0.221 * log10(pcb), 0)
  )
}

growth_qreg <- function(pcb) {
  ifelse(
    pcb < 0.100,
    0,
    pmax(0.15 + 0.0938 * log10(pcb), 0)
  )
}

## added in adjusted lipid regression

mort_qreg_lipid <- function(pcb) {
  ## Effect threshold is 100 ng/g (ww) (2200 ng/g lw), so return zero effect below this
  ifelse(
    pcb < 3.7, # need to confirm the limit of applicability for lipid weight
    0,
    pmax(-0.1253 + 0.221 * log10(pcb), 0)
  )
}


pcb_data <- read.csv("data/PCB_Puyallup_monitoring_data.csv")

dexposure_fun_2 <- function(concentrations) { ## not sure i need to do this if I have a pre-generated exposure
  means <- log(concentrations)
  sdlogs <- rep(0.5, length(concentrations)) # Set a default standard deviation value
  props <- rep(1 / length(concentrations), length(concentrations)) # Equal proportions for each concentration

  function(xs) {
    vapply(
      xs, \(x) sum(props * dlnorm(x, means, sdlogs)),
      0.0
    )
  }
}

# Wet weight PCB concentrations
concentrations <- pcb_data$pcb_uggwet

dexposure_ww <- dexposure_fun_2(concentrations)

# Wet weight PCB concentrations with sample-specific lipid corrections
concentrations_lw_sample <- pcb_data$pcb_ugglw_sample

dexposure_lw <- dexposure_fun_2(concentrations_lw_sample)

# # Wet weight PCB concentrations using 1% lipid corrections
concentrations_lw_1 <- pcb_data$pcb_ugglw_1

dexposure_lw_1 <- dexposure_fun_2(concentrations_lw_1)

## Calculate direct mortality rate based on the exposure distribution and the
## Berninger and Tillitt direct mortality relationship
morts_ww <- integrate(\(pcb) dexposure_ww(pcb) * mort_qreg(pcb), lower = 0, upper = Inf)$value
morts_ww
## Apply direct mortality to the to-ocean 0.0 transition (decreasing survival)
## ug/g wet weight
puyallup_mort <- eq_pop(puyallup_sim, nearshore_surv_adj = 1 - morts_ww, pop0 = rep(1000, 5))
spawner_mort <- get_spawners(puyallup_mort)
oa_mort <- get_oceanadults(puyallup_mort)
red_mort <- 1 - (oa_0 - oa_mort) / oa_0
red_mort
Percent_spawner_reduction <- (spawner_mort - spawners) / spawners * 100
Percent_spawner_reduction

## ug/g lipid weight (sample specific lipids)
morts_lw <- integrate(\(pcb) dexposure_lw(pcb) * mort_qreg_lipid(pcb), lower = 0, upper = Inf)$value
morts_lw
## Apply direct mortality to the to-ocean 0.0 transition (decreasing survival)
## ug/g wet weight
puyallup_mort <- eq_pop(puyallup_sim, nearshore_surv_adj = 1 - morts_lw, pop0 = rep(1000, 5))
spawner_mort <- get_spawners(puyallup_mort)
oa_mort <- get_oceanadults(puyallup_mort)
red_mort <- 1 - (oa_0 - oa_mort) / oa_0
red_mort
Percent_spawner_reduction <- (spawner_mort - spawners) / spawners * 100
Percent_spawner_reduction

## ug/g lipid weight (1% lipids)
morts_lw_1 <- integrate(\(pcb) dexposure_lw_1(pcb) * mort_qreg_lipid(pcb), lower = 0, upper = Inf)$value
morts_lw_1
## Apply direct mortality to the to-ocean 0.0 transition (decreasing survival)
## ug/g wet weight
puyallup_mort <- eq_pop(puyallup_sim, nearshore_surv_adj = 1 - morts_lw_1, pop0 = rep(1000, 5))
spawner_mort <- get_spawners(puyallup_mort)
oa_mort <- get_oceanadults(puyallup_mort)
red_mort <- 1 - (oa_0 - oa_mort) / oa_0
red_mort
Percent_spawner_reduction <- (spawner_mort - spawners) / spawners * 100
Percent_spawner_reduction
