white_sim <- function(pop) {
  prespawners <- c(0.093 * pop[4], 0.647 * pop[5], pop[6])
  spawners <- 0.838 * 0.867 * 0.971 * prespawners
  fecundity <- c(4600, 5700, 6600)
  eggs <- sum(fecundity * spawners / 2)
  fry <- eggs * 0.441
  fry_mig <-  0.16 * 0.599 * 0.803 * 0.356 * fry
  parr_mig <-  0.72 * 0.599 * 0.240 * 0.563 * fry
  yr_mig <- 0.12 * 0.599 * 0.596 * 0.604 * 0.959 * fry
  age00 <- 0.356 * fry_mig + 0.563 * parr_mig
  age01 <- 0.757 * pop[1] + 0.082 * pop[2]
  age02 <- 0.590 * pop[3]
  age03 <- 0.710 * (1 - 0.093) * pop[4]
  age04 <- 0.812 * (1 - 0.647) * pop[5]

  structure(c(yr_mig, age00, age01, age02, age03, age04),
            names = c("yr_mig", paste0("age0", 0:4)),
            oceanadults = age01 + age02 + age03 + age04,
            prespawners = prespawners,
            spawners = spawners,
            eggs = eggs,
            fry = fry,
            fry_migr = fry_mig,
            parr_migr = parr_mig)
}

pop <- rep(1e3, 6)
for (t in 2:300) {
  pop0 <- pop
  pop <- white_sim(pop)
}


return <- c(0, 0, 0, 0.093, 0.647, 1)
spawners <- 0.838 * 0.867 * 0.971 * return
fecund <- c(0, 0, 0, 4600 / 2, 5700 / 2, 6600 / 2)
eggs <- spawners * fecund
fry <- 0.441 * eggs
fry_mig <- 0.88 * 0.16 * 0.599 * 0.803 * 0.356 * fry
parr_mig <- 0.88 * 0.72 * 0.599 * 0.240 * 0.563 * fry

yr_mig <- 0.12 * 0.599 * 0.596 * 0.604 * 0.959 * fry
age00 <- 0.356 * fry_mig + 0.563 * parr_mig
age01 <- c(0.757, 0.082, 0, 0, 0, 0)
age02 <- c(0, 0, 0.590, 0, 0, 0)
age03 <- c(0, 0, 0, 0.710 * (1 - 0.093), 0, 0)
age04 <- c(0, 0, 0, 0, 0.812 * (1 - 0.647), 0)

white_mat <- rbind(yr_mig, age00, age01, age02, age03, age04)
colnames(white_mat) <- rownames(white_mat)

white_mat %*% rep(1e3, 6)
