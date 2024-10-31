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

growth_to_surv <- function(
    pcb,
    base_size = db_size$july_mass,
    noexp_surv = db_size$pred_surv,
    wt_type = "ww") {
  ## Calculate the expected reduction in size given exposure
  gred <- growth_qreg(pcb, wt_type)
  exp_mass <- (1 - gred) * base_size
  db2011_survival(exp_mass) / noexp_surv
}

combo_mort <- function(pcb, base_size = db_size$july_mass, noexp_surv = db_size$pred_surv, wt_type = "ww") {
  dir_mort <- mort_qreg(pcb, wt_type)
  grwth_mort <- 1 - growth_to_surv(pcb, base_size, noexp_surv, wt_type)
  ## Only individuals who are not affected by direct mortality are vulnerable
  ## to growth-restriction mortality
  dir_mort + (1 - dir_mort) * grwth_mort
}
