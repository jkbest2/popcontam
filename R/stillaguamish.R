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
  ## Egg production
  eggs <- beverton_holt(spawners / 2, 5400, 152859370)
  ## Survival to fry
  fry <- eggs * 0.548651
  # colonization is density-dependent
  # capacity doubled to allow for two cohorts of fish, non-simultaneously
  # occupying habitats
  pre_parr <- beverton_holt(fry, 0.302312^(1/9), 1682111.1 * 2)
  ## Split out fry and and parr migrants
  # about nine weeks of rearing for parr migrants - based on Snohomish screw
  # trap data. HARP adjust this to impose a late June temperature penalty
  parr_mig <- beverton_holt(pre_parr, 0.302312^(8/9), 1682111.1 * 2)
  ## Survival to subyearling rearing
  subyr <- fry * 0.302312^(1/9) # only 1 out of 9 weeks represented
  ## Using overall parr_mig as the number who stay indicates a constant-ish migration downstream?
  fry_mig <- (subyr - parr_mig) * 0.128 # migration mortality is estimated
  ## Subtracting off pre_parr rather than parr_mig seems more consistent with the lcm_chinook model
  ## fry_mig <- (subyr - pre_parr) * 0.128 # migration mortality is estimated
  stopifnot(fry_mig >= 0) # Shouldn't be a problem, but the check makes me feel better
  ## Transition to first marine year
  delta_fry <- beverton_holt(fry_mig, 0.35, 45335)
  # note that this age structure shouldn't be taken literally - we apply all
  # mortality effects in the nearshore, when we know that part of it is
  # attributable to harvest.
  age00 <- (delta_fry + parr_mig) * 0.01645436 * nearshore_surv_adj
  age01 <- pop[1] * 0.6
  age02 <- pop[2] * (1 - 0.05) * 0.7
  age03 <- pop[3] * (1 - 0.34) * 0.8
  age04 <- pop[4] * (1 - 0.92) * 0.9

  structure(c(age00, age01, age02, age03, age04, use.names = FALSE),
            names = paste0("age0", 0:4),
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
