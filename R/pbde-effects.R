# pbde_exp <- read_rds(here::here("data", "puyallup", "pbde_exposure.rds"))
# pbde_surv <- read_rds(here::here("data", "pbde_surv.rds"))

# pbde_mod <- read_rds(here::here("data", "pbde_model.rds"))
# dr_mod <- pbde_mod$pbde_model
# threshold <- pbde_mod$pbde_threshold

dr_coefs_rv <- function(mod = dr_mod, ndraws = 4e3) {
  beta <- coef(mod)
  n_beta <- length(beta)

  v <- vcov(mod)
  vchol <- chol(v)

  x <- rvar_rng(rnorm, n_beta, ndraws = ndraws)
  beta + t(vchol) %*% x
}

gen_base_surv <- function(threshold, surv_df = pbde_surv, ndraws = 4e3) {
  thr_surv <- surv_df |>
    filter(concentration < threshold) |>
    summarize(n_surv = sum(n_surv), n_dead = sum(n_dead))
  rvar_rng(rbeta, 1, thr_surv$n_surv + 1, thr_surv$n_dead + 1, ndraws = ndraws)
}

dr_surv <- function(conc, beta, thr = 0) {
  # map_vec(
  #   conc,
  #   function(conc) {
  dm <- predict(
    dr_mod,
    newdata = data.frame(concentration = conc),
    type = "lpmatrix"
  )
  rfun(plogis)(dm %**% beta)[, 1, drop = TRUE]
  #   }
  # )
}

exp_rel_surv <- function(conc, beta, pop_meanlog, pop_sdlog, base_surv, thr = threshold) {
  map_dbl(
    conc,
    function(conc) {
      if (conc < thr) {
        surv <- base_surv
      } else {
        dm <- predict(
          dr_mod,
          newdata = data.frame(concentration = conc),
          type = "lpmatrix"
        )
        surv <- plogis(dm %*% beta)
      }
      surv / base_surv * dlnorm(conc, pop_meanlog, pop_sdlog)
    }
  )
}

pbde_eff <- function(
    exp_post = pbde_exp,
    mod = dr_mod,
    threshold = 0,
    surv = pbde_surv,
    thin = 10) {
  pop_meanlog <- exp_post$pop_meanlog |>
    merge_chains() |>
    thin_draws(thin) |>
    draws_of()
  pop_sdlog <- exp_post$pop_sdlog |>
    merge_chains() |>
    thin_draws(thin) |>
    draws_of()
  beta_draws <- dr_coefs_rv(mod, ndraws = length(pop_meanlog)) |>
    draws_of()
  beta <- map(
    seq_len(nrow(beta_draws)),
    ~ beta_draws[., , , drop = TRUE]
  )
  base_surv <- gen_base_surv(threshold, surv, ndraws = length(pop_meanlog)) |>
    draws_of()
  thr <- rep(threshold, length(pop_meanlog))

  pmap_dbl(
    list(
      beta = beta,
      pop_meanlog = pop_meanlog,
      pop_sdlog = pop_sdlog,
      base_surv = base_surv,
      thr = thr
    ),
    function(beta, pop_meanlog, pop_sdlog, base_surv, thr) {
      integrate(
        function(conc) {
          exp_rel_surv(
            conc, beta,
            pop_meanlog, pop_sdlog,
            base_surv, thr
          )
        },
        lower = 0, upper = Inf
      )$value
    }
  ) |>
    rvar()
}

dr_pred_rv <- function(conc, mod = dr_mod, n = 4e3, threshold = 0) {
  newdata <- tibble(
    concentration = conc
  )
  nd0 <- filter(newdata, concentration < threshold)
  nd1 <- filter(newdata, concentration >= threshold)

  pred0 <- rep(
    gen_base_surv(threshold, pbde_surv, ndraws = 1000),
    nrow(nd0)
  )

  beta <- coef(mod)
  n_beta <- length(beta)

  v <- vcov(mod)
  vchol <- chol(v)

  x <- rvar_rng(rnorm, n_beta, ndraws = n)
  beta_sim <- beta + t(vchol) %*% x
  covar_sim <- predict(mod, newdata = nd1, type = "lpmatrix")
  pred1 <- covar_sim %**% beta_sim
  inv_link <- family(mod)$linkinv
  pred1 <- rfun(inv_link)(pred1)
  c(pred0, pred1)
}
