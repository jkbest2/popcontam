source("R/fish-size.R")

db2011_survival <- function(mass) {
  ## The Duffy and Beauchamp survival regression predicts percent survival. I
  ## prefer to work with survival rates, which are the percent divided by 100.
  ## Because log10(0.01) = -2, I included an extra -2 in the intercept of the
  ## regression. This is why the Duffy and Beauchamp intercept is -1.071 but
  ## we're using -3.071 here. Given that we're using relative changes in this
  ## survival it won't actually make a difference.
  l <- -3.071 + 0.041 * mass
  10^l
}

## Calculates expected mass (g) given the survival rate.
inv_db2011_survival <- function(survival) {
  ls <- log(survival, 10)
  (ls + 3.071) / 0.041
}

mort_qreg <- function(pcb, wt_type = "ww") {
  ## Effect threshold is 100 ng/g (ww), so return zero effect below this
  if (wt_type == "ww") {
    eff <- ifelse(
      pcb < 0.100,
      0,
      pmax(0.1702 + 0.221 * log10(pcb), 0)
    )
  } else if (wt_type == "lw") {
    eff <- ifelse(
      pcb < 2.2,
      0,
      pmax(-0.1253 + 0.221 * log10(pcb), 0)
    )
  } else {
    stop("wt_type must be \"ww\" for wet weight or \"lw\" for lipid weight")
  }
  eff
}

growth_qreg <- function(pcb, wt_type = "ww") {
  eff <- if (wt_type == "ww") {
    ifelse(
      pcb < 0.100,
      0,
      pmax(0.15 + 0.0938 * log10(pcb), 0)
    )
  } else if (wt_type == "lw") {
    eff <- ifelse(
      pcb < 2.2,
      0,
      pmax(0.0246 + 0.0938 * log10(pcb), 0)
    )
  } else {
    stop("wt_type must be \"ww\" for wet weight or \"lw\" for lipid weight")
  }
  eff
}

growth_mort <- function(
    pcb,
    base_size = NULL,
    base_surv = db_size$pred_surv,
    wt_type = "ww",
    remove_pcbs = TRUE) {
  ## If base size is not provided, calculate its value from the provided
  ## survival rate
  base_size <- base_size %||% inv_db2011_survival(base_surv)

  ## Calculate the reduction in growth due to PCB exposure.
  gr_red <- 1 - growth_qreg(pcb, wt_type)

  if (remove_pcbs) {
    exp_size <- base_size
    noexp_size <- base_size / gr_red
    exp_surv <- base_surv
    noexp_surv <- db2011_survival(noexp_size)
  } else {
    exp_size <- base_size * gr_red
    noexp_size <- base_size
    exp_surv <- db2011_survival(exp_size)
    noexp_surv <- base_surv
  }
  ## Calculate mortality due to PCB exposure. The proportion here gives the
  ## fraction of fish that that survive after exposure relative to those that
  ## would have survived with no exposure. The complement is then taken to get
  ## the relative change in mortality.
  1 - exp_surv / noexp_surv
}

combo_surv <- function(
    pcb,
    base_size = db_size$july_mass,
    base_surv = db_size$pred_surv,
    wt_type = "ww",
    remove_pcbs = TRUE) {
  dir_mort <- mort_qreg(pcb, wt_type)
  gr_mort <- growth_mort(pcb, base_size, base_surv, wt_type, remove_pcbs)

  mort <- dir_mort + (1 - dir_mort) * gr_mort
  surv <- 1 - mort
  if (remove_pcbs) {
    surv <- 1 / surv
  }
  surv
}
