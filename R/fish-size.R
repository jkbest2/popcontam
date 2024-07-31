db_size <- read_csv("data/DuffyBeauchamp2011-SizeTable.csv",
                    col_types =
                      cols(
                        region = col_character(),
                        hatchery = col_character(),
                        survival_pct = col_double(),
                        release_date = col_character(),
                        release_mass = col_double(),
                        july_n = col_double(),
                        july_fl = col_double(),
                        july_fl_se = col_double(),
                        july_mass = col_double(),
                        september_n = col_double(),
                        september_fl = col_double(),
                        september_fl_se = col_double(),
                        september_mass = col_double()
                      )) |>
   summarize(july_mass = weighted.mean(july_mass, july_n),
             n = sum(july_n)) |>
             ## .by = region) |>
  mutate(pred_surv = 10 ^ (-3.071 + 0.041 * july_mass))
