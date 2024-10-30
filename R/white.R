white_sim <- function(pop, nearshore_surv_adj = 1) {
  ## Find number of spawners
  prespawners <- c(
    ## Age 0.* contributions
    0.093 * pop[3], 0.647 * pop[4], pop[5],
    ## Age 1.* contributions
    0.093 * pop[8], 0.647 * pop[9], pop[10]
  )
  ## Added migrant and holding prespawners so that the relevant capacities can
  ## be included.
  migrantprespawn <- beverton_holt(prespawners, 0.837, 22972996.51)
  holdingprespawn <- beverton_holt(migrantprespawn, 0.867, 680441.18)
  spawners <- beverton_holt(migrantprespawn, 0.971, 522439767.58)

  ## Calculate number of eggs; Fecundity is age dependent and assumes 50/50 sex
  ## ratio
  fecundity <- c(
    age02 = spawners[1] * 4600 / 2,
    age03 = spawners[2] * 5700 / 2,
    age04 = spawners[3] * 6600 / 2,
    age12 = spawners[4] * 4600 / 2,
    age13 = spawners[5] * 5700 / 2,
    age14 = spawners[6] * 6600 / 2
  )
  egg_prod <- sum(fecundity)
  ## Account for density dependence in egg production
  eggs <- beverton_holt(egg_prod, 1, 411853741.56)

  ## Find number in each (fry, parr, yearling) life history
  fry <- beverton_holt(eggs, 0.441, 10212032.14)
  ## Assume 88% subyearlings, 16% of which are fry migrants
  fry_mig <- beverton_holt(fry, 0.599 * 0.88 * 0.16 * 0.803 * 0.356, 920795.40)
  ## Assume 88% subyearlings, 72% of which are parr migrant
  parr_mig <- beverton_holt(fry, 0.599 * 0.88 * 0.72 * 0.240 * 0.563, 153767.53)
  # Assume 12% yearling migrants. Freshwater rearing stages are split because
  # they have separate capacities
  yr_0agerear <- beverton_holt(fry, 0.599 * 0.12 * 0.596, 768048.17)
  yr_0agemig <- beverton_holt(yr_0agerear, 0.900, 184719306.08)
  yr_0ageinac <- beverton_holt(yr_0agemig, 0.604, 970597.87)
  yr_1agerear <- beverton_holt(yr_0ageinac, 0.959, 6436546.88)
  yr_mig <- beverton_holt(yr_1agerear, 0.757, 197943007.58)

  ## Calculate the number at each age
  ## Subyearling migrants
  age00 <- beverton_holt(fry_mig + parr_mig, 0.08, 60660147.66) *
    nearshore_surv_adj
  age01 <- beverton_holt(pop[1], 0.590, 687697841.82)
  age02 <- beverton_holt((1 - 0.093) * pop[2], 0.731, 14835699.80)
  age03 <- beverton_holt((1 - 0.647) * pop[3], 0.849, 6838144.42)
  ## TODO Check on this mortality if all are expected to return
  age04 <- beverton_holt(pop[4], 0.882, 1890200.90)
  ## Yearling migrants
  age10 <- beverton_holt(yr_mig, 0.0978, 30203352.87)
  age11 <- beverton_holt(pop[6], 0.590, 660837621.60)
  age12 <- beverton_holt((1 - 0.093) * pop[7], 0.731, 14092487.28)
  age13 <- beverton_holt((1 - 0.647) * pop[8], 0.849, 6940007.42)
  ## TODO Check on this mortality if all are expected to return
  age14 <- beverton_holt(pop[9], 0.882, 1648312.88)

  structure(
    c(
      age00, age01, age02, age03, age04,
      age10, age11, age12, age13, age14
    ),
    names = c(paste0("age0", 0:4), paste0("age1", 0:4)),
    ## TODO Check on how to count "adults" for this model
    oceanadults = unname(
      age01 + age02 + age03 + age04 +
        age10 + age11 + age12 + age13 + age14
    ),
    prespawners = unname(prespawners),
    holdingprespawn = unname(holdingprespawn),
    migrantprespawn = unname(migrantprespawn),
    spawners = unname(sum(spawners)),
    subyearspawners = unname(sum(spawners[1:3])),
    yrlspawners = unname(spawners[4:6]),
    egg_prod = unname(egg_prod),
    eggs = unname(eggs),
    eggssubyr = unname(sum(fecundity[1:3])),
    fry = unname(fry),
    parr_mig = unname(parr_mig),
    yr_mig = unname(yr_mig),
    fry_mig = unname(fry_mig)
  )
}
