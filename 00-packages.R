# install.packages("pak")

pak::pkg_install(
  c(
    "tidyverse",
    "rstan",
    # "BH", # Boost for Stan
    # "RcppEigen", # Eigen also not installed for Stan when using pak
    "posterior",
    "bayesplot",
    "ggdist",
    "readxl",
    "splines",
    "patchwork",
    # "RTMB",
    "sf",
    "logKDE",
    "writexl",
    "here",
    "knitr",
    "scales"
  ),
  dependencies = TRUE # Does this fix the BH and RcppEigen problems?
)
