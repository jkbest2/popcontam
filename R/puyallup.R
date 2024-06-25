puyallup_sim <- function(pop) {
  prespawners <- c(0.093 * pop[3], 0.647 * pop[4], pop[5])
  spawners <- 0.722 * 0.873 * 0.971 * prespawners
  fecundity <- c(4600, 5700, 6600)
  eggs <- sum(fecundity * spawners / 2)
  fry <- 0.391 * eggs
  fry_mig <- 0.626 * 0.31 * 0.924 * fry
  parr_mig <- 0.626 * 0.69 * 0.35 * fry
  age00 <- 0.311 * fry_mig + 0.551 * parr_mig
  age01 <- 0.080 * pop[1]
  age02 <- 0.591 * pop[2]
  age03 <- 0.710 * (1 - 0.093) * pop[3]
  age04 <- 0.812 * (1 - 0.647) * pop[4]

  structure(c(age00, age01, age02, age03, age04),
            names = paste0("age0", 0:4),
            prespawners = prespawners,
            spawners = spawners,
            fry = fry,
            parr = parr)
}

pop <- rep(1e3, 6)
for (t in 2:300) {
  pop0 <- pop
  pop <- white_sim(pop)
}


return <- c(0, 0, 0.093, 0.647, 1)
spawners <- 0.722 * 0.873 * 0.971 * return
fecund <- c(0, 0, 4600 / 2, 5700 / 2, 6600 / 2)
eggs <- spawners * fecund
fry <- 0.391 * eggs
fry_mig <- 0.626 * 0.31 * 0.924 * fry
parr_mig <- 0.626 * 0.69 * 0.35 * fry
age00 <- 0.311 * fry_mig + 0.551 * parr_mig
age01 <- c(0.08, 0    , 0                  , 0                  , 0)
age02 <- c(0   , 0.591, 0                  , 0                  , 0)
age03 <- c(0   , 0    , 0.710 * (1 - 0.093), 0                  , 0)
age04 <- c(0   , 0    , 0                  , 0.812 * (1 - 0.647), 0)

puyallup_mat <- rbind(age00, age01, age02, age03, age04)
colnames(puyallup_mat) <- rownames(puyallup_mat)

puyallup_mat %*% rep(1e3, 6)
