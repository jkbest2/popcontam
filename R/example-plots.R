library(tidyverse)
source("R/utils.R")

prod <- 0.5
cap <- 100

xrange <- c(0, 400)
xseq <- seq(xrange[1], xrange[2])

bh_df <- tibble(
  x = xseq,
  ymin = 0,
  ymax = beverton_holt(xseq, prod, cap),
  type = "Remain"
)
di_df <- tibble(
  x = xseq,
  ymin = bh_df$ymax,
  ymax = xseq * prod,
  type = "Migrate"
)

ex_loc <- 175
bind_rows(bh_df, di_df) |>
  ggplot(aes(x = x, ymin = ymin, ymax = ymax, fill = type, color = type)) +
  geom_ribbon(alpha = 0.8) +
  annotate("errorbar", x = ex_loc, ymin = 0, ymax = beverton_holt(ex_loc, prod, cap),
           linetype = "solid", width = 4) +
  annotate("label", x = ex_loc, y = beverton_holt(ex_loc, prod, cap) / 2,
           label = paste(round(beverton_holt(ex_loc, prod, cap)), "fish remain")) +
  annotate("errorbar", x = ex_loc, ymin = beverton_holt(ex_loc, prod, cap), ymax = ex_loc * prod,
           width = 4, linetype = "solid") +
  annotate("label", x = ex_loc, y = (beverton_holt(ex_loc, prod, cap) + ex_loc * prod) / 2,
           label = paste(round(ex_loc * prod - beverton_holt(ex_loc, prod, cap)), "fish migrate")) +
  labs(x = "Fish before density-dependent migration",
       y = "Fish after density-dependent migration",
       fill = "Action", color = "Action") +
  coord_cartesian(expand = FALSE)
ggsave("figs/dd_mig_example.png")
