##' @title Plot an illustration of density-dependent migration
##'
##' This plot follows the form of density-dependent migration presented in the
##' HARP model.
##'
##' @param prod Shared productivity
##' @param cap Density dependent capacity
##' @param xrange Range of populations to plot
##' @param ex_pop Example population with labeled percentages
##' @export
plot_dd_mig <- function(prod = 1, cap = 100, xrange = c(0, 250), ex_pop = 150) {
  ## xrange <- c(0, 250)
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

  ## ex_loc <- 150
  ex_rem <- beverton_holt(ex_loc, prod, cap)
  ex_mig <- ex_loc * prod - ex_rem
  ex_props <- proportions(c(ex_rem, ex_mig))

  bind_rows(bh_df, di_df) |>
    ggplot(aes(x = x, ymin = ymin, ymax = ymax, fill = type, color = type)) +
    geom_ribbon(alpha = 0.8) +
    geom_hline(yintercept = cap, linetype = "dashed", alpha = 0.8) +
    annotate("label", x = 25, cap, label = "Capacity") +
    annotate("errorbar", x = ex_loc, ymin = 0, ymax = ex_rem,
             linetype = "solid", width = 4) +
    annotate("label", x = ex_loc, y = ex_rem / 2,
             label = paste0(round(ex_props[1] * 100), "% remain")) +
    annotate("errorbar", x = ex_loc, ymin = ex_rem, ymax = ex_mig + ex_rem,
             width = 4, linetype = "solid") +
    annotate("label", x = ex_loc, y = (2 * ex_rem + ex_mig) / 2,
             label = paste0(round(ex_props[2] * 100), "% migrate")) +
    labs(x = "Fish before density-dependent migration",
         y = "Fish after density-dependent migration",
         fill = "Action", color = "Action") +
    coord_cartesian(expand = FALSE)
}

## ggsave("figs/dd_mig_example.png")

## curve(beverton_holt(x, prod, cap), from = 0, to = 10000, ylim = c(0, 105))
## abline(h = 100, lty = "dashed")

## T <- 15
## p <- matrix(NA, nrow = T, ncol = 9)
## p[1, ] <- c(1, 25, 50, 75, 100, 200, 300, 400, 500)
## for (t in 2:T) {
##   p[t, ] <- beverton_holt(p[t - 1, ], 2, cap)
## }
## matplot(p, type = 'l', ylim = c(0, 500))
## abline(h = cap, lty = "dashed")
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Plot the Beverton-Holt production function
##' @param prod Productivity parameter
##' @param cap Capacity parameter
##' @param xrange Range of populations to plot
##' @return A figure (ggplot) illustrating the B-H production function
##' @export
plot_bh <- function(prod = 0.5, cap = 100, xrange = c(0, 1025)) {
  xseq <- seq(xrange[1], xrange[2])

  bh_df <- tibble(
    x = xseq,
    y = beverton_holt(xseq, prod, cap),
  )

  ## ## ex_loc <- 150
  ## ex_rem <- beverton_holt(ex_loc, prod, cap)
  ## ex_mig <- ex_loc * prod - ex_rem
  ## ex_props <- proportions(c(ex_rem, ex_mig))

  bh_df |>
    ggplot(aes(x = x, y = y)) +
    geom_line() +
    ## geom_ribbon(alpha = 0.8) +
    geom_hline(yintercept = cap, linetype = "dashed", alpha = 0.8) +
    annotate("label", x = diff(xrange) / 4, cap, label = "Capacity") +
    ## annotate("errorbar", x = ex_loc, ymin = 0, ymax = ex_rem,
    ##          linetype = "solid", width = 4) +
    ## annotate("label", x = ex_loc, y = ex_rem / 2,
    ##          label = paste0(round(ex_props[1] * 100), "% remain")) +
    ## annotate("errorbar", x = ex_loc, ymin = ex_rem, ymax = ex_mig + ex_rem,
    ##          width = 4, linetype = "solid") +
    ## annotate("label", x = ex_loc, y = (2 * ex_rem + ex_mig) / 2,
    ##          label = paste0(round(ex_props[2] * 100), "% migrate")) +
    ## labs(x = "Fish before density-dependent migration",
    ##      y = "Fish after density-dependent migration",
    ##      fill = "Action", color = "Action") +
    scale_x_continuous(expand = expansion(0, 0)) +
    scale_y_continuous(expand = expansion(c(0, 0.05), 0)) +
    theme_bw()

    ## coord_cartesian(expand = FALSE)
}

cap_df <- tibble(
  Transition = "Density-dependent",
  capacity = 100,
  label = "Capacity",
  lab_x = 500
)

tibble(
  n = 0:1020,
  `Density-independent` = 0.5 * n,
  `Density-dependent` = beverton_holt(n, 0.5, 100)
) |>
  pivot_longer(cols = c(`Density-independent`, `Density-dependent`),
               names_to = "Transition", values_to = "n1") |>
  ggplot(aes(x = n, y = n1, color = Transition)) +
  geom_line() +
  geom_hline(data = cap_df, aes(yintercept = capacity, color = Transition),
             linetype = "dashed", alpha = 0.8) +
  geom_label(data = cap_df, aes(x = lab_x, y = capacity, label = label),
             alpha = 0.8, show.legend = FALSE) +
  labs(x = "Fish before transition", y = "Fish after transition") +
  scale_x_continuous(breaks = seq(0, 1000, 200), expand = expansion()) +
  scale_y_continuous(expand = expansion()) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("figs/transition-plot.png", width = 6.5, height = 4.5)


source("R/fish-size.R") # Get data from D&B2011

direct_mortality <- function(pcb_ug) {
  pmax(0.1702 + 0.221 * log10(pcb_ug), 0)
}
## FIXME This is *only* the growth effect, *does not include* the translation to survival yet!
growth_restriction <- function(pcb_ug) {
  pmax(0.15 + 0.0938 * log10(pcb_ug), 0)
}
db2011_survival <- function(mass) {
  l <- -3.071 + 0.041 * mass
  10^l
}
growth_to_mort <- function(pcb, base_size = db_size$july_mass, noexp_surv = db_size$pred_surv) {
  ## Calculate the expected reduction in size given exposure
  gred <- growth_restriction(pcb)
  exp_mass <- (1 - gred) * base_size
  1 - db2011_survival(exp_mass) / noexp_surv
}

## Example combined effects plot
combo_eff_df <- tibble(
  pcb_ng = seq(0, 1.025e3, length.out = 2049),
  ## pcb_ng = 10^seq(1, 6, length.out = 1025),
  pcb_ug = pcb_ng / 1000,
  `Direct Mortality` = direct_mortality(pcb_ug),
  `Growth Restriction` = growth_to_mort(pcb_ug),
  `Combined Naive` = `Direct Mortality` + `Growth Restriction`,
  `Combined Conditional` = `Direct Mortality` + (1 - `Direct Mortality`) * `Growth Restriction`
)

combo_ex <- tibble(
  pcb_ng = 500,
  pcb_ug = pcb_ng / 1000,
  `Direct Mortality` = direct_mortality(pcb_ug),
  `Growth Restriction` = growth_to_mort(pcb_ug),
  `Combined Naive` = `Direct Mortality` + `Growth Restriction`,
  `Combined Conditional` = `Direct Mortality` + (1 - `Direct Mortality`) * `Growth Restriction`
)

combo_eff_df |>
  pivot_longer(cols = c(`Direct Mortality`, `Combined Naive`, `Combined Conditional`),
               names_to = "aop", values_to = "effect") |>
  mutate(naive = ifelse(grepl("Naive", aop), TRUE, FALSE))  |>
  ggplot(aes(x = pcb_ng, y = effect, color = aop)) +
  annotate("ribbon", x = c(0, 100), ymin = 0, ymax = Inf, fill = "gray50", alpha = 0.5) +
  geom_line() +
  annotate("errorbar",
           x = combo_ex$pcb_ng - 5,
           ymin = 0, ymax = combo_ex$`Direct Mortality`, width = 10) +
  annotate("text",
           x = combo_ex$pcb_ng - 15, y = combo_ex$`Direct Mortality` / 2,
           label = paste0("Direct Mortality\n", scales::percent(combo_ex$`Direct Mortality`, 0.1)),
           hjust = "right")+
  annotate("errorbar",
           x = combo_ex$pcb_ng - 5,
           ymin = combo_ex$`Direct Mortality`, ymax = combo_ex$`Combined Naive`, width = 10) +
  annotate("text",
           x = combo_ex$pcb_ng - 15, y = mean(c(combo_ex$`Direct Mortality`, combo_ex$`Combined Naive`)),
           label = paste0("Growth Restriction\n", scales::percent(combo_ex$`Growth Restriction`, 0.1)),
           hjust = "right")+
  annotate("errorbar",
           x = combo_ex$pcb_ng + 5,
           ymin = 0, ymax = combo_ex$`Combined Conditional`, width = 10) +
  annotate("text",
           x = combo_ex$pcb_ng + 15,
           y = combo_ex$`Combined Conditional` / 2 + 0.035,
           label = paste0("Combined\n", scales::percent(combo_ex$`Combined Conditional`, 0.1), " = ",
                          scales::percent(combo_ex$`Direct Mortality`, 0.1), " + (1 - ",
                          scales::percent(combo_ex$`Direct Mortality`, 0.1), ") x ",
                          scales::percent(combo_ex$`Growth Restriction`, 0.1)),
           hjust = "left") +
  labs(x = "Tissue PCB Concentration (ng/g ww)", y = "Effect",
       title = "Growth-selective mortality occurs after direct mortality") +
  scale_x_continuous(limits = c(0, NA), breaks = seq(0, 1000, 200), labels = scales::comma, expand = expansion()) +
  scale_y_continuous(limits = c(0, NA), labels = scales::percent, expand = expansion()) +
  scale_color_discrete(guide = FALSE) +
  ## facet_wrap(~ aop, ncol = 1) +
  theme_bw()
ggsave("figs/conditional-effects.png", width = 10, height = 6)
