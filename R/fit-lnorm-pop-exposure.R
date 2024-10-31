library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

lnorm_data <- function(contam, conc_col) {
  ncomp <- sort(unique(contam$n_composite))
  N <- nrow(contam) # nolint: object_name_linter.
  contam <- contam |>
    mutate(
      ncomp_idx = match(contam$n_composite, ncomp),
      conc = {{ conc_col }}
    )

  list(
    C = length(ncomp),
    ncomp = as.array(ncomp), # In case there is only a single ncomp value
    N = N,
    ncomp_idx = contam$ncomp_idx,
    conc = contam$conc
  )
}

fit_lnorm_exposure <- function(data, ...) {
  stan("inst/lnorm-exposure.stan",
    data = data,
    ...
  )
}
