#' State space model with AR(1) structure
#' Y_t = X_t + e_t,
#' X_t = \phi X_{t-1} + w_t,
#' X_0 ~ N[m_0, C_0]
#' where e_t ~ N[0, \sigma^2] and w_t ~ N[0, \tau^2] and \phi \in (0, 1)
#' The parameters are \theta = (\sigma^2, \tau^2, \phi). For now we keep the parameters
#' fixed.

rm(list = ls())
gc()
Rcpp::sourceCpp(file = "./R/ar1/ar1.cpp")
source("./R/wquantile.R")
source("./R/ar1/BFfilterAR1.R")


# Simulate data ----------------------------------------------------------------

set.seed(1212)
n <- 100L
phi <- 0.5
sigma <- .4
tau <- 1
m0 <- 0.0
C0 <- sqrt(tau^2 / (1 - phi^2)) # stationary variance 1.0
y <- numeric(n)
x <- numeric(n)
x[1] <- stats::rnorm(1L, mean = m0, sd = C0)
y[1] <- x[1] + stats::rnorm(1L, sd = sigma)
for (t in seq_len(n)[-1]) {
  x[t] <- phi * x[t - 1] + stats::rnorm(1L, sd = tau)
  y[t] <- x[t] + stats::rnorm(1L, sd = sigma)
}
# par(mfrow = c(1, 2))
plot(y, type = "l")
lines(x, col = "blue")


# Run BF and KF ----------------------------------------------------------------

parameters <- c(phi, sigma, tau)
initial_ss <- c(m0, C0)

# Initialise the module
mod <- new(StateSpaceAR1, y, initial_ss, 2)
mod$var_state
mod$mean_state
mod$RunKalmanFilter(parameters)

graphics.off()
plot(x, type = "l")
lines(mod$mean_ss, col = "red")
lines(mod$mean_ss + 2 * sqrt(mod$var_ss), col = "red", lty = 3)
lines(mod$mean_ss - 2 * sqrt(mod$var_ss), col = "red", lty = 3)

# Bootstrap filter implemented in R
n_particles <- 1000L
out_bf <- BFfilterAR1(y = y, parms = parameters, initial_ss = initial_ss,
                      n_particles = n_particles)

# Bootstrap filter implemented in C++
#Rcpp::sourceCpp(file = "./R/ar1.cpp")
mod2 <- new(StateSpaceAR1, y, initial_ss, n_particles)
mod2$RunBootstrapFilter(parameters)
mod2$lml
# mod2$particles
# mod2$ancestors
# mod2$weights

particles <- matrix(mod2$particles, nrow = n, ncol = n_particles, byrow = TRUE)
ancestors <- matrix(mod2$ancestors, nrow = n, ncol = n_particles, byrow = TRUE)
weights <- matrix(mod2$weights, nrow = n, ncol = n_particles, byrow = TRUE)
rowSums(weights)

# Get the quantiles of the filtering distribution for each time
qs <- matrix(nrow = n, ncol = 3)
for (i in 1:n) {
  qs[i,] <- wquantile(particles[i, ], weights[i, ], prob = c(0.05, 0.5, 0.95))
}

# Comparison -------------------------------------------------------------------


# Estimated log-likelihood
mod$lml # KF
mod2$lml # BF in C++
out_bf$lml # BF in R


# Filtering distribution under the KF and BF
png(filename = "./figures/comparison_kf_bf_ar1.png", width = 1000, height = 700)

par(mfrow = c(3, 1))
plot(x, main = "KF")
lines(mod$mean_ss, col = "red")
lines(mod$mean_ss + 2 * sqrt(mod$var_ss), col = "red", lty = 3)
lines(mod$mean_ss - 2 * sqrt(mod$var_ss), col = "red", lty = 3)

plot(x, type = "p", main = "BF (C++)")
lines(qs[, 1], col = "blue", lty = 2)
lines(qs[, 2], col = "blue")
lines(qs[, 3], col = "blue", lty = 2)

plot(x, type = "p", main = "BF (R)")
lines(out_bf$mean, col = "blue")
lines(out_bf$mean - 2*sqrt(out_bf$var), col = "blue", lty = 2)
lines(out_bf$mean + 2*sqrt(out_bf$var), col = "blue", lty = 2)


graphics.off()
