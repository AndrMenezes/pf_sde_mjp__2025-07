#' State space model with AR(1) structure
#' Y_t = X_t + e_t,
#' X_t = \phi X_{t-1} + w_t,
#' X_0 ~ N[m_0, C_0]
#' where e_t ~ N[0, \sigma^2] and w_t ~ N[0, \tau^2] and \phi \in (0, 1)
#' Here, we assumed that \phi is unknown and we want to estimate it

rm(list = ls())
gc()
Rcpp::sourceCpp(file = "./R/ar1.cpp")
# source("./R/wquantile")
# source("./R/BFfilterAR1.R")


# Simulate data ----------------------------------------------------------------

set.seed(1212)
n <- 100L
true_phi <- 0.5
true_sigma <- .4
true_tau <- 1
m0 <- 0.0
C0 <- sqrt(true_tau^2 / (1 - true_phi^2)) # stationary variance
y <- numeric(n)
x <- numeric(n)
x[1] <- stats::rnorm(1L, mean = m0, sd = C0)
y[1] <- x[1] + stats::rnorm(1L, sd = true_sigma)
for (t in seq_len(n)[-1]) {
  x[t] <- true_phi * x[t - 1] + stats::rnorm(1L, sd = true_tau)
  y[t] <- x[t] + stats::rnorm(1L, sd = true_sigma)
}
# par(mfrow = c(1, 2))
plot(y, type = "l")
lines(x, col = "blue")


# Run particle MCMC with KF ----------------------------------------------------


# Initialise the module
initial_ss <- c(m0, C0)
n_particles <- 500L
mod <- new(StateSpaceAR1, y, initial_ss, n_particles)
# Run KF fixing the parameters at true values
# mod$RunKalmanFilter(parameters)


# Number of MCMC iterations
n_mcmc <- 1000L

# Prior for log(phi/(1-phi))
sd_prior_logit_phi <- 1.0

# initial value for phi
phi_cur <- 0.4
parameters <- c(phi_cur, true_sigma, true_tau)

# transformation
logit_phi_cur <- log(phi_cur / (1 - phi_cur))

# Get estimate of the likelihood
mod$RunBootstrapFilter(parameters)
# mod$RunKalmanFilter(parameters)
lml_cur <- mod$lml + logit_phi_cur

# Proposal is random-walk: logit(phi)' ~ N[phi_cur, sd_prop_phi]
sd_prop_logit_phi <- 0.5

# Container for the draws of phi
phi_draws <- numeric(n_mcmc)
accept_rate <- 0L
for (k in seq_len(n_mcmc)) {

  # Sample from proposal at logit scale and transform back
  logit_phi_prop <- stats::rnorm(n = 1L, mean = logit_phi_cur, sd = sd_prop_logit_phi)
  phi_prop <- 1 / (1 + exp(-logit_phi_prop))

  # Compute log marginal likelihood using KF
  parameters[1L] <- phi_prop
  # mod$RunKalmanFilter(parameters)
  mod$RunBootstrapFilter(parameters)
  lml_prop <- mod$lml + logit_phi_prop # plus jacobian

  # log-prior
  lprior_prop <- stats::dnorm(x = logit_phi_prop, mean = 0, sd = sd_prior_logit_phi)
  lprior_cur <- stats::dnorm(x = logit_phi_cur, mean = 0,
                             sd = sd_prior_logit_phi)
  # log-proposals
  lq_prop <- stats::dnorm(x = logit_phi_prop, mean = logit_phi_cur,
                          sd = sd_prop_logit_phi)
  lq_cur <- stats::dnorm(x = logit_phi_cur, mean = logit_phi_prop,
                         sd = sd_prop_logit_phi)

  # MH ratio
  lr <- lml_prop - lml_cur + lprior_prop - lprior_cur + lq_cur - lq_prop
  if (log(stats::runif(1L)) < lr) {
    # Accept
    phi_cur <- phi_prop
    logit_phi_cur <- log(phi_cur / (1 - phi_cur))
    lml_cur <- lml_prop
    accept_rate <- accept_rate + 1L
  }
  phi_draws[k] <- phi_cur
}
accept_rate / n_mcmc

phi_draws_mat <- matrix(phi_draws, ncol = 1)
colnames(phi_draws_mat) <- "phi"
coda::effectiveSize(coda::as.mcmc(phi_draws))
bayesplot::mcmc_acf(x = phi_draws_mat)
bayesplot::mcmc_trace(x = phi_draws_mat)
bayesplot::mcmc_recover_hist(x = phi_draws_mat, true = true_phi)







