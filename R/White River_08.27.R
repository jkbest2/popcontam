beverton_holt <- function(n, prod, cap) {
  n * prod / (1 + n * prod / cap)
}

white_sim <- function(pop) {
  prespawners <- c(0.093 * pop[3], 0.647 * pop[4], pop[5], 0.093 * pop[8], 0.647 * pop[9], pop[10]) # additional pops related to yearlings
  migrantprespawn <- beverton_holt(prespawners, 0.837, 22972996.51) #added in this step as there was capacity associated with migrant prespawner
  holdingprespawn <- beverton_holt(migrantprespawn, 0.867, 680441.18)#added in this step as there was capacity associated with holding prespawner
  spawners <- beverton_holt(migrantprespawn, 0.971, 522439767.58)
  fecundity_age02 <- spawners[1] *4600/2 # fecundity specific for returning two year old females (0.2), assumed sex ratio is 50/50 female/male
  fecundity_age03 <- spawners[2] *5700/2 # fecundity specific for returning three year old females (0.3), assumed sex ratio is 50/50 female/male
  fecundity_age04 <- spawners[3] *6600/2 # fecundity specific for returning four year old females (0.4), assumed sex ratio is 50/50 female/male
  fecundity_age12 <- spawners[4] *4600/2 # fecundity specific for returning two year old females (1.2), assumed sex ratio is 50/50 female/male
  fecundity_age13 <- spawners[5] *5700/2 # fecundity specific for returning three year old females (1.3), assumed sex ratio is 50/50 female/male
  fecundity_age14 <- spawners[6] *6600/2 # fecundity specific for returning four year old females (1.4), assumed sex ratio is 50/50 female/male
  egg_prod <- sum(fecundity_age02 + fecundity_age03 + fecundity_age04 + fecundity_age12 + fecundity_age13 + fecundity_age14) # summing eggs across spawners
  eggs <- beverton_holt(egg_prod, 1, 411853741.56) # applying capacity to egg production, using productivity of 1 since fecundity was captured in previous three equations
  fry <- beverton_holt(eggs, 0.441, 10212032.14)
  fry_mig <- beverton_holt(fry, 0.599 * 0.88 *0.16* 0.803 * 0.356, 920795.40) #88% subyearlings, 16% fry migrant
  parr_mig <- beverton_holt( fry, 0.599 * 0.88* 0.72 * 0.240 * 0.563, 153767.53) #88% subyearlings, 72% parr migrant
  yr_0agerear <- beverton_holt(fry, 0.599 * 0.12 * 0.596, 768048.17) #12% yearling migrant. I am splitting out freshwater rearing stages (below) for yearlings as each stage has a separate carrying capacity
  yr_0agemig <- beverton_holt(yr_0agerear, 0.900, 184719306.08) ## this survival rate was missing from the original conceptual model, now included
  yr_0ageinac <- beverton_holt(yr_0agemig, 0.604, 970597.87)
  yr_1agerear <- beverton_holt(yr_0ageinac, 0.959, 6436546.88)
  yr_mig <- beverton_holt(yr_1agerear, 0.757, 197943007.58)
  age00 <- beverton_holt(fry_mig + parr_mig, 0.08, 60660147.66) 
  age01 <- beverton_holt(pop[1], 0.590, 687697841.82)
  age02 <- beverton_holt((1-0.093) * pop[2], 0.731, 14835699.80)
  age03 <- beverton_holt((1 - 0.647) * pop[3], 0.849, 6838144.42)
  age04 <- beverton_holt(pop[4], 0.882, 1890200.90) ## not sure how to handle this mortality if all are expected to return)
  age10 <- beverton_holt(yr_mig, 0.0978, 30203352.87) ## this survival rate was missing from the original conceptual model for age 1.0 yearlings
  age11 <- beverton_holt(pop[6], 0.590, 660837621.60) # Age 1.1 Yearlings
  age12 <- beverton_holt((1 - 0.093) * pop[7], 0.731, 14092487.28) # Age 2.1 Yearlings
  age13 <- beverton_holt((1 - 0.647) * pop[8], 0.849, 6940007.42) #Age 3.1 Yearlings
  age14 <- beverton_holt(pop[9], 0.882, 1648312.88) ## Age 4.1 Yearlings...not sure how to handle this mortality if all are expected to return)
  
  structure(c(age00, age01, age02, age03, age04, age10, age11, age12, age13, age14),
            names = paste0("age0", 0:4),
            oceanadults = unname(age01 + age02 + age03 + age04 + age10 + age11 + age12 + age13 + age14),
            prespawners = unname(prespawners),
            holdingprespawn = unname(holdingprespawn),
            migrantprespawn = unname(migrantprespawn),
            spawners = unname(spawners[1] + spawners[2] + spawners[3] + spawners[4] + spawners[5] + spawners[6]),
            subyearspawners = unname(spawners[1] + spawners[2] + spawners[3]),
            yrlspawners = unname(spawners[4] + spawners[5] + spawners[6]),
            egg_prod = unname(egg_prod),
            eggs = unname(eggs),
            eggssubyr = unname(fecundity_age02 + fecundity_age03 + fecundity_age04),
            fry = unname(fry),
            parr_mig = unname(parr_mig),
            yr_mig = unname(yr_mig),
            fry_mig = unname(fry_mig))
}

library(tidyverse)

source("R/utils.R")
source("R/fish-size.R")

white_0 <- eq_pop(white_sim, pop0 = rep(1000, 10))

white_0




##subyearling SAR
smolts <- attr (white_0, "parr_mig")+ attr(white_0, "fry_mig")
smolts
subyrspawn <- attr(white_0, "subyearspawners")
subyrspawn
sar_0 <- subyrspawn/ smolts
sar_0

##Yearling SAR
yrspawn <- attr(white_0, "yrlspawners")
sar_yrl <- yrspawn / attr(white_0, "yr_mig")
sar_yrl

get_eggs <- function(pop) {
    attr(pop, "eggs")
  }

##combined freshwater survival (yearlins and subyearlings combined)
freshsurv <- (smolts + attr(white_0, "yr_mig"))/get_eggs(white_0)
freshsurv


##freshwater survival of subyearlings
fresurvsub <- smolts / attr(white_0, "eggssubyr")
fresurvsub

#total outmigrants
total_outmigrants <- attr(white_0, "parr_mig") + attr(white_0, "fry_mig") + attr(white_0, "yr_mig")
total_outmigrants
