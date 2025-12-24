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
sigma <- .4
tau <- 1
m0 <- 0.0
C0 <- tau^2 / (1 - phi^2) # stationary variance 1.0
y <- numeric(n)
x <- numeric(n)
x[1] <- stats::rnorm(1L, mean = m0, sd = sqrt(C0))
y[1] <- x[1] + stats::rnorm(1L, sd = sigma)
for (t in seq_len(n)[-1]) {
  x[t] <- phi * x[t - 1] + stats::rnorm(1L, sd = tau)
  y[t] <- x[t] + stats::rnorm(1L, sd = sigma)
}
# par(mfrow = c(1, 2))
plot(y, type = "l")
lines(x, col = "blue")

# Initialise the module
parameters <- c(phi, sigma*sigma, tau*tau)
initial_ss <- c(m0, 1)
mod <- new(StateSpaceAR1, y, parameters, initial_ss, 2)
mod$var_state
mod$mean_state
mod$RunKalmanFilter()

graphics.off()
plot(x, type = "l")
lines(mod$mean_ss, col = "red")
lines(mod$mean_ss + 2 * sqrt(mod$var_ss), col = "red", lty = 3)
lines(mod$mean_ss - 2 * sqrt(mod$var_ss), col = "red", lty = 3)

# Try the bootstrap filter now

Rcpp::sourceCpp(file = "./R/ar1.cpp")
n_particles <- 100L
mod2 <- new(StateSpaceAR1, y, parameters, initial_ss, n_particles)
mod2$RunBootstrapFilter()
mod2$particles
mod2$ancestors
mod2$weights

sum(mod2$weights_curr)

particles <- matrix(mod2$particles, nrow = n, ncol = n_particles, byrow = TRUE)
ancestors <- matrix(mod2$ancestors, nrow = n, ncol = n_particles, byrow = TRUE)
weights <- matrix(mod2$weights, nrow = n, ncol = n_particles, byrow = TRUE)
rowSums(weights)

wquantile<-function(xs,ws,prob) {
  ws = ws[order(xs)]
  xs = xs[order(xs)]
  sws = sum(ws); cws=cumsum(ws)
  np = length(prob); qs=rep(0,np)
  for (i in 1:np) {
    qs[i] <- xs[which.max(cws / sws >= prob[i])] ## which.max gives 1st instance
  }
  return(qs)
}

qs = matrix(nrow = n, ncol = 3)
for (i in 1:n) {
  qs[i,] = wquantile(particles[i, ], weights[i, ], prob = c(0.05, 0.5, 0.95))
}

qlo=min(qs); qhi=max(qs)
plot(x, type = "l", ylim = c(qlo, qhi))
lines(qs[, 1], col = "red", lty = 2)
lines(qs[, 2], col = "red")
lines(qs[, 3], col = "red", lty = 2)

