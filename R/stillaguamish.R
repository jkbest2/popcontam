library(tidyverse)

## stilly_states <-c(
##   egg = NA,
##   fry = NA,
##   fry_migrant = NA,
##   parr_migrant = NA,
##   fry_migrant_nearshore_rearing = NA,
##   parr_migrant_nearshore_rearing = NA,
##   fry_migrant_estuary_rearing = NA,
##   parr_migrant_estuary_rearing = NA,
##   )
##

vl <- function(v) {
  as.vector(v, mode = "list")
}


stilly <- tribble(
  ~ state    , ~ tr_from             , ~ tr_to                    , ~ form,
  "egg"      , c("spawner")          , c("fry")                   , NA    ,
  "fry"      , c("egg")              , c("fry_migr", "parr_migr") , NA    ,
  "fry_migr" , c("fry")              , c("fm_near", "fm_est")     , NA    ,
  "parr_migr", c("fry")              , c("pm_near", "pm_est")     , NA    ,
  "fm_near"  , c("fry_migr")         , c("age_01")                , NA    ,
  "fm_est"   , c("fry_migr")         , c("age_01")                , NA    ,
  "pm_near"  , c("parr_migr")        , c("age_01")                , NA    ,
  "pm_est"   , c("parr_migr")        , c("age_01")                , NA    ,
  "age_01"   , c("fm_near", "fm_est",
                 "pm_near", "pm_est"), c("age_02")                , NA    ,
  "age_02"   , c("age_01")           , c("age_03", "return_02")   , NA    ,
  "age_03"   , c("age_02")           , c("age_04", "return_03")   , NA    ,
  "age_04"   , c("age_03")           , c("age_05", "return_04")   , NA    ,
  "age_05"   , c("age_04")           , c("return_05")             , NA    ,
  "return_02", c("age_02")           , c("spawn_02")              , NA    ,
  "return_03", c("age_03")           , c("spawn_03")              , NA    ,
  "return_04", c("age_04")           , c("spawn_04")              , NA    ,
  "return_05", c("age_05")           , c("spawn_05")              , NA    ,
  "spawn_02" , c("return_02")        , c("egg")                   , NA    ,
  "spawn_03" , c("return_03")        , c("egg")                   , NA    ,
  "spawn_04" , c("return_04")        , c("egg")                   , NA    ,
  "spawn_05" , c("return_05")        , c("egg")                   , NA    ,
)

stilly_edges <- stilly |>
  select(state, tr_to) |>
  unnest(tr_to) |>
  rename(from = state, to = tr_to)

stilly_edges |>
  ## filter(to != "egg") |>
ggplot(aes(from_id = from, to_id = to)) +
  geom_net(directed = TRUE,
           layout.alg = "target")

stilly_net <- network.initialize(n = nrow(stilly))
add.vertices(stilly_net, )
add.edges(stilly_net, stilly_edges$from, stilly_edges$to)

beverton_holt <- function(n, prod, cap) {
   n * prod / (1 + n * prod / cap)
}

density_independent <- function(n, prod) {
  n * prod
}

surv_bh <- function(n, prod, cap) {
  beverton_holt(n / cap, prod, 1)
}

surv_di <- function(n, prod) {
  prod
}

stg_pars <- tribble(
~ pop         , ~ lifestage     , ~ type        , ~ v_curr   ,
"fall_chinook", "Prespawn (N-O)", "Productivity", 0.696937384,
"fall_chinook", "Prespawn (H-O)", "Productivity", 0.521801472,
"fall_chinook", "Winter"        , "Productivity", 0.353555703,
"fall_chinook", "Winter"        , "Capacity"    , 166235.6105,
"fall_chinook", "Colonization"  , "Productivity", 0.302311982,
"fall_chinook", "Spawning"      , "Capacity"    , 152859371  ,
"fall_chinook", "Incubation"    , "Productivity", 0.548651115,
"fall_chinook", "Subyearling"   , "Capacity"    , 1682111.117,
"fall_chinook", "Subyearling"   , "Productivity", 0.30231197 ,
"fall_chinook", "Summer"        , "Capacity"    , 278121.5919,
"fall_chinook", "Summer"        , "Productivity", 0.645911874,
)
