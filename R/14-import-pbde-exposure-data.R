library(tidyverse)
library(readxl)

pbde <- read_xlsx(
  here::here("data", "PBDE Sum 47 and 99_Stillaguamish and Puyallup.xlsx")
) |>
  rename(
    id = SampleID,
    river = RiverSystem,
    n_composite = CompositeN,
    bde47_ww = `Conc_Found BDE_47 (ng/g wet)`,
    bde99_ww = `Conc_Found BDE_99 (ng/g wet)`,
    bd99_loq = `BDE_99 LOQ`,
    bde_sum = `Conc_BDE47+BDE99 (ng/g ww)`,
    ratio = Ratio,
    location = `Sample Location`,
    pct_lipid = `% lipids`,
    matrix = Matrix
  ) |>
  separate(
    location,
    into = c("river", "location"),
    sep = "_"
  ) |>
  mutate(
    river = str_to_title(river)
  )

pbde |>
  filter(
    river == "Puyallup",
    location %in% c("Estuary", "Nearshore")
  ) |>
  write_rds(
    here::here("data", "puyallup", "pbde_exposure_obs.rds")
  )

pbde |>
  filter(
    river == "Stillaguamish",
    location %in% c("Estuary", "Nearshore")
  ) |>
  write_rds(
    here::here("data", "stillaguamish", "pbde_exposure_obs.rds")
  )




pbde |>
  ggplot(aes(x = bde_sum, color = location, fill = location)) +
  geom_density(color = NA, alpha = 0.5, position = position_stack()) +
  facet_wrap(~river, nrow = 1)
