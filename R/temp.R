base_df <- expand_grid(
  stage = c("Age 00", "Spawners"),
  river = c("Puyallup", "White"),
) |>
  mutate(
    baseline = c(puy_age00, white_age00, puy0_sp, white0_sp)
  )


df <- expand_grid(
  stage = c("Age 00", "Spawners"),
  river = c("Puyallup", "White"),
  norm = c("ww", "lw", "lw1"),
) |>
  mutate(
    n = c(
      puy_age00_ww, puy_age00_lw, puy_age00_lw1,
      white_age00_ww, white_age00_lw, white_age00_lw1,
      puy_sp_ww, puy_sp_lw, puy_sp_lw1,
      white_sp_ww, white_sp_lw, white_sp_lw1
    )
  ) |>
  left_join(
    base_df,
    by = join_by(stage, river)
  ) |>
  mutate(
    change = (n - baseline) / baseline
  )

df |>
  mutate(
    n = mean(n),
    change = mean(change)
  ) |>
  write_csv("white_puyallup_table.csv")

combo_df <- df |>
  filter(river != "Stillaguamish") |>
  summarize(
    n = rvar_sum(n),
    baseline = sum(baseline),
    .by = c(stage, norm)
  ) |>
  mutate(
    change = (n - baseline) / baseline,
    river = "Puyallup/White"
  )

bind_rows(df, combo_df) |>
  mutate(
    river = factor(river, levels = c("Puyallup", "White", "Puyallup/White")),
    n = mean(n),
    change = mean(change),
    change_pct = scales::percent(change)
  ) |>
  arrange(stage, river) |>
  write_csv("combined_table.csv")

puy_eff <- read_rds("data/puyallup/puy_pcb_eff.rds")
white_eff <- read_rds("data/puyallup/white_pcb_eff.rds")
stilly_eff <- read_rds("data/stillaguamish/pcb_eff.rds")

eff_df <- expand_grid(
  river = c("Puyallup", "White", "Stillaguamish"),
  norm = c("ww", "lw", "lw1")
) |>
  mutate(
    effect = c(
      puy_eff$ww, puy_eff$lw, puy_eff$lw1,
      white_eff$ww, white_eff$lw, white_eff$lw1,
      stilly_eff$ww, stilly_eff$lw, stilly_eff$lw1
    )
  ) |>
  mutate(
    addl_mort = mean(1 - 1 / effect),
    addl_mort_pct = scales::percent(addl_mort),
    effect = mean(effect - 1),
    eff_pct = scales::percent(effect),
  ) |>
  select(
    river, norm, effect, eff_pct, addl_mort, addl_mort_pct
  )

write_csv(eff_df, "effect_table.csv")
