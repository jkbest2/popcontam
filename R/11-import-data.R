library(tidyverse)

## This data set includes the PCB concentrations normalized by wet weight,
## lipid  weight, and 1% lipid weight. It also includes spatial locations of
## the samples.
contam <- read_csv(
  "data/PCB_Puyallup_monitoring_data.csv",
  col_types =
    cols(
      SAMPLE_ID = col_character(),
      Study_ID = col_character(),
      Latitude = col_double(),
      Longitude = col_double(),
      Species = col_character(),
      RiverSystem = col_character(),
      COMPOUND = col_character(),
      UNIT = col_character(),
      `CONC_FOUND (ng/g ww)` = col_skip(),
      `Sampling Date` = col_date(format = "%m/%d/%Y"),
      `% Lipids` = col_double(),
      pcb_uggwet = col_double(),
      pcb_ugglw_sample = col_double(),
      pcb_ugglw_1 = col_double(),
      `effect-mort ww` = col_skip(),
      `effect-mort lipid-sample specific` = col_skip(),
      `effect-mort lipid - 1%` = col_skip(),
      Composite_Count = col_double()
    )
) |>
  rename(
    id = SAMPLE_ID,
    study = Study_ID,
    latitude = Latitude,
    longitude = Longitude,
    species = Species,
    river = RiverSystem,
    compound = COMPOUND,
    unit = UNIT,
    date = `Sampling Date`,
    pct_lipid = `% Lipids`,
    pcb_ug_ww = pcb_uggwet,
    pcb_ug_lw = pcb_ugglw_sample,
    pcb_ug_lw1 = pcb_ugglw_1,
    n_composite = Composite_Count,
  )

## This data set includes labels for the sample locations such as hatchery,
## river, nearshore, and estuary. We are primarily interested in the latter two
## so we use the sample ID to associate each sample with a location.
contam2 <- read_csv(
  "data/PCB_Stilly_Puyallup_monitoring_data.csv",
  col_types = cols(
    SAMPLE_ID = col_character(),
    Study_ID = col_character(),
    Species = col_character(),
    RiverSystem = col_character(),
    COMPOUND = col_character(),
    UNIT = col_character(),
    `Max of CONC_FOUND` = col_double(),
    PERCENT_LIPID_CONTENT = col_character(),
    SAMPLE_LOCATION = col_character()
  )
) |>
  rename(
    id = SAMPLE_ID,
    study = Study_ID,
    species = Species,
    river = RiverSystem,
    compound = COMPOUND,
    unit = UNIT,
    pcb_ng_ww = `Max of CONC_FOUND`,
    pct_lipid = `PERCENT_LIPID_CONTENT`,
    location = SAMPLE_LOCATION
  )

## Filter thse to remove any that were not gathered in the estuary or nearshore
contam3 <- contam |>
  left_join(
    select(contam2, id, location),
    by = join_by(id)
  ) |>
  filter(grepl("Estuary|Nearshore", location))

write_rds(contam3, "data/puy-pcb-estns.rds")
