##' @describeIn pc_statevec Create the state vector
##' @export
new_pc_statevec <- function(state, statenames = names(state), t = 0) {
  if (!is.null(statenames) && (length(state) != length(statenames))) {
    stop("If provided, names must be same length as state")
  }
  structure(state,
            names = statenames,
            class = "statevec")
}

##' @describeIn pc_statevec Validate the state vector
##' @export
validate_pc_statevec <- function(statevec) {
  if (any(is.na(statevec) | (statevec < 0))) {
    stop("All states must be non-missing and non-negative",
         call. = FALSE)
  }
  if (!is.numeric(attr(statevec, "t"))) {
    stop("Time must be numeric")
  }
  statevec
}

##' Create a named vector representing the state of the population at a given
##' point in time
##'
##' @title Create a state vector for the population model
##' @param state Number of individuals in each state
##' @param statenames Name of each state (e.g. lifestages)
##' @param t Time
##' @return A named vector with class \code{statevec}
##' @export
pc_statevec <- function(state, statenames = names(state), t = 0) {
  new_statevec(state, statenames, t)
}

pc_time.statevec <- function(state) {
  attr(state, "t")
}

new_pc_transition <- function(fn, timestep) {


}

validate_pc_transition <- function(tr) {


}

pc_transition <- function() {

}
