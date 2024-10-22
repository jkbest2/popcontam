library(tidyverse)

source("R/utils.R")
source("R/fish-size.R")

white_0 <- eq_pop(white_sim, pop0 = rep(1000, 10))

white_0




## subyearling SAR
smolts <- attr(white_0, "parr_mig") + attr(white_0, "fry_mig")
smolts
subyrspawn <- attr(white_0, "subyearspawners")
subyrspawn
sar_0 <- subyrspawn / smolts
sar_0

## Yearling SAR
yrspawn <- attr(white_0, "yrlspawners")
sar_yrl <- yrspawn / attr(white_0, "yr_mig")
sar_yrl

get_eggs <- function(pop) {
  attr(pop, "eggs")
}

## combined freshwater survival (yearlins and subyearlings combined)
freshsurv <- (smolts + attr(white_0, "yr_mig")) / get_eggs(white_0)
freshsurv


## freshwater survival of subyearlings
fresurvsub <- smolts / attr(white_0, "eggssubyr")
fresurvsub

# total outmigrants
total_outmigrants <- attr(white_0, "parr_mig") + attr(white_0, "fry_mig") + attr(white_0, "yr_mig")
total_outmigrants
