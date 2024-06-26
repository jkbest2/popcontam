## From Beechie et al. 2023 (Stillaguamish and Snohomish HARP model description), page 39
beverton_holt <- function(n, prod, cap) {
   n * prod / (1 + n * prod / cap)
}

get_spawners <- function(pop) {
  attr(pop, "spawners")
}

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
            oceanadults = age01 + age02 + age03 + age04,
            prespawners = unname(prespawners),
            spawners = unname(spawners),
            eggs = unname(eggs),
            fry = unname(fry),
            parr_mig = unname(parr_mig),
            fry_mig = unname(fry_mig))
}

## Set an initial population
pop0 <- rep(1e3, 5)
stilly_sim(pop0)

## Run it for a while to get to equilibrium
pop <- pop0
for (t in 2:300) {
  pop <- stilly_sim(pop)
}

pop
get_spawners(pop)

## Run the simulation function iteratively until the population vector reaches
## equilibrium
eq_pop <- function(sim_fun, ..., pop0, n_max = 300, abstol = 1e-10) {
  pop1 <- sim_fun(pop0, ...)
  n <- 1
  while (max(abs(pop1 - pop0)) > abstol) {
    pop0 <- pop1
    pop1 <- sim_fun(pop1, ...)
    n <- n + 1
    if (n > n_max) break
  }
  message(n)
  pop1
}

eq_pop(stilly_sim, pop0 = rep(3e3, 5))

growth_morts <- 0.12
