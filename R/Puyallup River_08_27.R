beverton_holt <- function(n, prod, cap) {
n * prod / (1 + n * prod / cap)
}

puyallup_sim <- function(pop) {
  prespawners <- c(0.093 * pop[3], 0.647 * pop[4], pop[5])
  migrantprespawn <- beverton_holt(prespawners, 0.722, 12657408.21) #added in this step as there was capacity associated with holding prespawner
  holdingprespawn <- beverton_holt(migrantprespawn, 0.872, 729726.40)#added in this step as there was capacity associated with migrant prespawner
  spawners <- beverton_holt(migrantprespawn, 0.971, 1139372348.62)
  fecundity_age02 <- spawners[1] *4600/2 # fecundity specific for returning two year old females, assumed sex ratio is 50/50 female/male
  fecundity_age03 <- spawners[2] *5700/2 # fecundity specific for returning three year old females, assumed sex ratio is 50/50 female/male
  fecundity_age04 <- spawners[3] *6600/2 # fecundity specific for returning four year old females, assumed sex ratio is 50/50 female/male
  egg_prod <- sum(fecundity_age02 + fecundity_age03 + fecundity_age04) # summing eggs across spawners
  eggs <- beverton_holt(egg_prod, 1, 623792859.31) # applying capacity to egg production, using productivity of 1 since fecundity was captured in previous three equations
  fry <- beverton_holt(eggs, 0.391, 9629410.93)
  fry_mig <- beverton_holt(fry,0.626 * 0.31 * 0.924 * 0.310, 1558678.02) # capturing fry/parr split (31%/69%)
  parr_mig <- beverton_holt( fry, 0.626 * 0.69 * 0.35 * 0.551, 308326.34) # capturing fry/parr split (31%/69%)
  age00 <- beverton_holt(fry_mig + parr_mig, 0.08, 67139801.46) 
  age01 <- beverton_holt(pop[1], 0.591, 963877008.33)
  age02 <- beverton_holt((1-0.093) * pop[2], 0.710, 20077276.13)
  age03 <- beverton_holt((1 - 0.647) * pop[3], 0.812, 8809670.23)
  age04 <- beverton_holt(pop[4], 0.881, 1889294.85)
  

  structure(c(age00, age01, age02, age03, age04),
            names = paste0("age0", 0:4),
            oceanadults = unname(age01 + age02 + age03 + age04),
            prespawners = unname(prespawners),
            holdingprespawn = unname(holdingprespawn),
            migrantprespawn = unname(migrantprespawn),
            spawners = unname(spawners[1] + spawners[2] + spawners[3]), #spawners coming through as age02, age03, and age04 so summing for model output
            egg_prod = unname(egg_prod),
            eggs = unname(eggs),
            fry = unname(fry),
            parr_mig = unname(parr_mig),
            fry_mig = unname(fry_mig))
}




puyallup_0 <- eq_pop(puyallup_sim, pop0 = rep(1000, 5))

puyallup_0
smolts <- attr(puyallup_0, "parr_mig") + attr(puyallup_0, "fry_mig")
sar_0 <- get_spawners(puyallup_0) / smolts
sar_0
smolts
get_eggs <- function(pop) {
    attr(pop, "eggs")
  }

freshsurv <- smolts/get_eggs(puyallup_0)
freshsurv
