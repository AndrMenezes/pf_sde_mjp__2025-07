#' State space model with AR(1) structure
#' Y_t = X_t + e_t,
#' X_t = \phi X_{t-1} + w_t,
#' X_0 ~ N[m_0, C_0]
#' where e_t ~ N[0, \sigma^2] and w_t ~ N[0, \tau^2] and \phi \in (0, 1)
#' The parameters are \theta = (\sigma^2, \tau^2, \phi). For now we keep the parameters
#' fixed.

rm(list = ls())
Rcpp::sourceCpp(file = "./R/ar1.cpp")

# Simulate some data
set.seed(1212)
n <- 100L
phi <- 0.5
sigma <- 1
tau <- 0.5
m0 <- 0.0
C0 <- 0.2
y <- numeric(n)
x <- numeric(n)
x[1] <- stats::rnorm(1L, mean = m0, sd = sqrt(C0))
y[1] <- x[1] + stats::rnorm(1L, sd = sigma)
for (t in seq_len(n)[-1]) {
  x[t] <- phi * x[t - 1] + stats::rnorm(1L, sd = tau)
  y[t] <- x[t] + stats::rnorm(1L, sd = sigma)
}
par(mfrow = c(1, 2))
plot(x, type = "l")
plot(y, type = "l")

# Initialise the module
ml <- Rcpp::Module(module = "StateSpaceAR1")
parameters <- c(phi, sigma*sigma, tau*tau)
initial_ss <- c(m0, C0)
mod <- new(StateSpaceAR1, y, parameters, initial_ss)
mod$var_state
mod$mean_state
mod$RunKalmanFilter()

graphics.off()
plot(x, type = "l", ylim = c(-4, 4))
lines(mod$mean_ss, col = "red")
lines(mod$mean_ss + 2 * sqrt(mod$var_ss), col = "red", lty = 3)
lines(mod$mean_ss - 2 * sqrt(mod$var_ss), col = "red", lty = 3)




