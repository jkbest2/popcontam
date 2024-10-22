kldiv <- function(modfun, modpars, truefun, truepars, xdom = c(0, Inf)) {
  fn <- function(x) {
    truefun(x, truepars[1], truepars[2]) *
      (truefun(x, truepars[1], truepars[2], log = TRUE) -
        modfun(x, modpars[1], modpars[2], log = TRUE))
  }
  integrate(fn, lower = xdom[1], upper = xdom[2])
}

kldiv(dlnorm, c(2, 1), dlnorm, c(0, 1), c(0, Inf))

fn <- function(x) {
  dlnorm(x, 0, 1, log = FALSE) *
    (dlnorm(x, 0, 1, log = TRUE) - dlnorm(x, 1, 1, log = TRUE))
}

curve(fn(x), from = 0, to = 5)
abline(h = 0, lty = "dashed")

integrate(fn, lower = 0, upper = Inf)

totvardist <- function(modfun, modpars, truefun, truepars, xdom = c(0, Inf)) {
  fn <- function(x) {
    abs(truefun(x, truepars[1], truepars[2]) -
      modfun(x, truepars[1], truepars[2]))
  }
  td <- integrate(fn, lower = xdom[1], upper = xdom[2])
  td$value <- td$value / 2
  td$abs.error <- td$abs.error / 2
  td
}
