library(tidyverse)
library(readxl)

still <- read_xlsx(
  "data/PCB_juvenile_Stilly_monitoring_data.xlsx",
  col_types = c(
    "text", # SAMPLE_ID
    "text", # Study_ID
    "text", # Species
    "text", # RiverSystem
    "skip", # COMPOUND
    "skip", # UNIT
    "numeric", # CONC_FOUND (ng/g ww)
    "date", # Sampling Date
    "numeric", # % Lipids
    "numeric", # Composite
    "text" # Sample Location
  )
) |>
  rename(
    id = SAMPLE_ID,
    study = Study_ID,
    species = Species,
    river = RiverSystem,
    pcb_ng_ww = `CONC_FOUND (ng/g ww)`,
    date = `Sampling Date`,
    pct_lipid = `% Lipids`,
    n_composite = Composite,
    location = `Sample Location`
  ) |>
  mutate(
    pcb_ug_ww = 1e-3 * pcb_ng_ww,
    pcb_ug_lw = pcb_ug_ww / pct_lipid / 0.01,
    pcb_ug_lw1 = pcb_ug_ww / 0.01,
    date = as.Date(date)
  ) |>
  filter(
    location %in% c("Estuary", "Nearshore")
  ) |>
  select(
    id,
    study,
    species,
    river,
    date,
    pcb_ug_ww,
    pcb_ug_lw,
    pcb_ug_lw1,
    n_composite,
    location
  )

if (!dir.exists("data/stillaguamish")) dir.create("data/stillaguamish")
write_rds(still, "data/stillaguamish/pcb_estns.rds")
