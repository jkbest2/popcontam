library(tidyverse)
library(readxl)

chk_pcb <- read_xlsx("data/CB_Stilly_Puyallup_monitoring_data.xlsx") |>
  separate(Species, c("lifestage", "species")) |>
  mutate(river = str_to_title(RiverSystem),
         unit = "ng/g (wet)") |>
  select(river, species, lifestage,
         pcbs = `Max of CONC_FOUND`,
         unit) |>
  mutate(gt_lla = pcbs >= 100)
