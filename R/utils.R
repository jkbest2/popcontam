##' Calculates a density-dependent transition using the Beverton-Holt
##' productivity curve.
##'
##' Uses the paramterization from Beechie et al. 2023 (Stillaguamish and
##' Snohomish HARP model description), page 39
##' @title Beverton-Holt suvival function
##' @param n Number of individuals in
##' @param prod Intrinsic productivity
##' @param cap Capacity
##' @return Number of individuals out
##' @author John Best
##' @export
beverton_holt <- function(n, prod, cap) {
   n * prod / (1 + n * prod / cap)
}

##' Convenience function to access the \code{"spawners"} attribute
##'
##' @title Get number of spawners
##' @param pop Population vector with \code{"spawners"} attribute
##' @return Number of spawners
##' @author John Best
##' @export
get_spawners <- function(pop) {
  attr(pop, "spawners")
}

##' Convenience function to access \code{"oceanadults"} attribute, which
##' represents available food for SRKWs
##'
##' @title Get number of ocean adults
##' @param pop Population vector with \code{"oceanadults"} attribute
##' @author John Best
##' @export
get_oceanadults <- function(pop) {
  attr(pop, "oceanadults")
}

##' Iterates lifecycle models to find an equilibrium population. Currently works
##' with \link{\code{stilly_sim}}. Should also work with \link{\code{white_sim}}
##' and \link{\code{puyallup_sim}} once those have density-dependent stages.
##'
##' @title Determine equilibrium population state through iteration
##' @param sim_fun Lifecycle model simulation function
##' @param ...Extra parameters to pass to \code{sim_fun}
##' @param pop0 Initial population
##' @param n_max Maximum number of iterations
##' @param abstol Absolute tolerance for maximum difference between two
##'   population state iterations
##' @return Population state at equilibrium, including relevant
##' @author John Best
##' @export
eq_pop <- function(sim_fun, ..., pop0, n_max = 300, abstol = 1e-10) {
  pop1 <- sim_fun(pop0, ...)
  n <- 1
  while (max(abs(pop1 - pop0)) > abstol) {
    pop0 <- pop1
    pop1 <- sim_fun(pop1, ...)
    n <- n + 1
    if (n >= n_max) {
      warning("abstol not met after ", n,
              " iterations, max abs difference is ",
              max(abs(pop1 - pop0)))
      break
      }
  }
  pop1
}
