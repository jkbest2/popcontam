pcb <- seq(0, 1000)

pop

x <- seq(-4, 4, length.out = 1025)
d1 <- dnorm(x)
d2 <- dnorm(x)
y2 <- convolve(d1, d2, type = "filter")

plot(x, d1, type = "l")
lines(x, y2)

plot(seq(-8, 8, length.out = 2049), y2, type = "l")
lines(x, d1)

s2 <- function(t) {
  integrate(\(tau) dnorm(tau) * dnorm(t - tau), lower = -Inf, upper = Inf)$value
}

cv <- sapply(seq(-8, 8, length.out = 2049), s2)

plot(seq(-8, 8, length.out = 2049), cv, type = "l")
lines(x, d1)
lines(x, 2 * dnorm(2 * x, 0, sqrt(2)))

var(replicate(10000, mean(rnorm(3))))

rnlnorm <- function(n, m, meanlog, sdlog) {
  replicate(n, sum(rlnorm(m, meanlog, sdlog)), simplify = TRUE)
}

sig_n <- function(n, meanlog, sdlog) {
  sig2 <- sdlog^2
  log((exp(sig2) - 1) * n * exp(2 * meanlog) / (n * exp(meanlog))^2 + 1)
}

mu_n <- function(n, meanlog, sdlog) {
  sig2 <- sdlog^2
  log(n * exp(meanlog)) + sig2 / 2 - sig_n(n, meanlog, sdlog) / 2
}

n <- 10
ml <- 10
sdl <- 1
y10 <- rnlnorm(10e3, n, ml, sdl)

hist(y10, probability = TRUE)
curve(dlnorm(x, mu_n(n, ml, sdl), sig_n(n, ml, sdl)),
  from = 0, to = 100e5, n = 1025, add = TRUE
)

rnlnorm <- function(n, m, meanlog, sdlog) {
  replicate(n, sum(rlnorm(m, meanlog, sdlog)), simplify = TRUE)
}

sig_n <- function(n, meanlog, sdlog) {
  sig2 <- sdlog^2
  log((exp(sig2) - 1) * n * exp(2 * meanlog) / (n * exp(meanlog))^2 + 1)
}

mu_n <- function(n, meanlog, sdlog) {
  sig2 <- sdlog^2
  log(n * exp(meanlog)) + sig2 / 2 - sig_n(n, meanlog, sdlog) / 2
}

n <- 10
ml <- 10
sdl <- 1
y10 <- rnlnorm(10e3, n, ml, sdl)

hist(y10, probability = TRUE)
curve(dlnorm(x, mu_n(n, ml, sdl), sig_n(n, ml, sdl)), from = 0, to = 100e5, n = 1025, add = TRUE)
