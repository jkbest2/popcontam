## Individual-based model

ocean_entry_pr <- c(1, 0)
return_pr <- c(0, 0, 0.08, 0.36, 0.52, 0.04)

## States to consider
## 1. egg - by origin tributary
## 2. fry - by origin tributary
## 3.
##

state <- c(eggs = 10e3,
           fry_migrants = 0,
           parr_migrants = 0,
           nearshore_fry = 0,
           nearshore_parr = 0,
           estuary_fry = 0,
           estuary_parr = 0,
           age_01 = 0,
           age_02 = 0,
           age_03 = 0,
           age_04 = 0,
           age_05 = 0,
           adults = 0,
           spawners = 0)

tmat <- matrix(0, nrow = length(state), ncol = length(state),
               dimnames = list(from = names(state), to = names(state)))

settr <- function(tmat, ...) {
  trs <- list(...)
  ## frto <- stringr::str_split(names(trs), "\\.")
  frto <- strsplit(names(trs), "\\.")
  for (idx in seq_along(trs)) {
    tmat[frto[[idx]][1], frto[[idx]][2]] <- trs[[idx]]
  }
  tmat
}

mar_prod <- c(age_1.age_2 = 0.6,
              age_2.age_3 = 0.7,
              age_3.age_4 = 0.8,
              age_4.age_5 = 0.9,
              age_5.age_5 = 0.9) # Meaning of each of these is not clear
mat_rate <- c(bsy_2 = 0.050,
              bsy_3 = 0.34,
              bsy_4 = 0.92,
              bsy_5 = 1.0)

## Need to think about how to structure within-year steps correctly. For example,
settr(tmat,
      spawners.eggs = 5400 / 2,
      eggs.fry_migrants = NA,
      eggs.parr_migrants = NA,
      fry_migrants.nearshore_fry = NA,
      parr_migrants.nearshore_parr = NA,
      fry_migrants.estuary_fry = NA,
      parr_migrants.estuary_parr = NA,
      nearshore_fry.age_01 = 0, # based on 0.2% survival of Skagit nearshore fry
      nearshore_parr.age_01 = NA,
      estuary_fry.age_01 = NA,
      estuary_parr.age_01 = NA,
      age_01.age_02 = 0.6,
      age_02.adults = 0.050, # bsy_2
      age_02.age_03 = 0.7 (1 - ),
      age_03.adults = NA,
      age_03.age_04 = NA,
      age_04.adults = NA,
      age_04.age_05 = NA,
      age_05.adults = 1,
      adults.spawners = 0.95
      )
