puyallup_sim <- function(pop, nearshore_surv_adj = 1) {
  ## Get number of spawners
  prespawners <- c(0.093 * pop[3], 0.647 * pop[4], pop[5])
  prespawner_prop <- proportions(prespawners)
  ## Migrant and holding prespawner steps added to account for associated
  ## capacities
  migrantprespawn <- beverton_holt(sum(prespawners), 0.722, 12657408.21) *
    prespawner_prop
  holdingprespawn <- beverton_holt(sum(migrantprespawn), 0.872, 729726.40) *
    prespawner_prop
  spawners <- beverton_holt(holdingprespawn, 0.971, 1139372348.62) *
    prespawner_prop

  ## Fecundity is age-dependent and assumes 50/50 sex ratio
  fecundity <- c(
    age02 = spawners[1] * 4600 / 2,
    age03 = spawners[2] * 5700 / 2,
    age04 = spawners[3] * 6600 / 2
  )
  egg_prod <- sum(fecundity)
  ## Account for capacity; production occurs in the fecundity step above
  eggs <- beverton_holt(egg_prod, 1, 623792859.31)

  ## Allow for different life histories
  fry <- beverton_holt(eggs, 0.391, 9629410.93)
  ## Assume 31% fry migrants
  fry_mig <- beverton_holt(fry, 0.626 * 0.31 * 0.924 * 0.310, 1558678.02)
  ## Assume 69% parr migrants
  parr_mig <- beverton_holt(fry, 0.626 * 0.69 * 0.35 * 0.551, 308326.34)

  ## Calculate number of fish at each age
  age00 <- beverton_holt(
    fry_mig + parr_mig, 0.08, 67139801.46
  ) * nearshore_surv_adj
  age01 <- beverton_holt(pop[1], 0.591, 963877008.33)
  age02 <- beverton_holt((1 - 0.093) * pop[2], 0.710, 20077276.13)
  age03 <- beverton_holt((1 - 0.647) * pop[3], 0.812, 8809670.23)
  age04 <- beverton_holt(pop[4], 0.881, 1889294.85)

  structure(c(age00, age01, age02, age03, age04),
    names = paste0("age0", 0:4),
    oceanadults = unname(age01 + age02 + age03 + age04),
    prespawners = unname(prespawners),
    holdingprespawn = unname(holdingprespawn),
    migrantprespawn = unname(migrantprespawn),
    spawners = unname(spawners[1] + spawners[2] + spawners[3]),
    egg_prod = unname(egg_prod),
    eggs = unname(eggs),
    fry = unname(fry),
    parr_mig = unname(parr_mig),
    fry_mig = unname(fry_mig)
  )
}
