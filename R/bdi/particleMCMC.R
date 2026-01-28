#' Run particle MCMC to perform inference on the parameters of BDI system.
#' @param y observed data.
#' @param x0 initial vale for the state.
#' @param n_particles number of particles.
#' @param dt discretised time interval, where 1/dt give number of points to discretise
#' the interval (t, t+1].
#' @param ndpost,nskip number of posterior and skip draws.
#' @param sd_prior standard deviation for the normal prior on the log-scale of the
#' parameters.
#' @param sd_prop standard deviation of the normal proposal.
#' @param printevery integer to print number of iteration.
#' @description
#' This function run the particle MCMC algorithm to perform inference on the
#' Birth-Death-Immigration process. So far, it works when estimating only when
#' parameter and fixing the remaining.
#' To run the Bootstrap filter a C++ class named BDI is called.
#'
particleMCMC_bdi <- function(y, x0, n_particles, dt, ndpost, nskip,
                             theta_fix = c(0.8, 40.0),
                             wh_fix = c(2L, 3L),
                             sd_prior = rep(1.0, 3L - length(theta_fix)),
                             S_prop = rep(0.25, 3L - length(theta_fix)),
                             printevery = 10L) {
  if (is.matrix(S_prop)) {
    S_prop <- t(chol(S_prop))
    rprop <- function(np, mu, chol_cov_prop) drop(mu + chol_cov_prop %*% stats::rnorm(np))
  } else {
    rprop <- function(np, mu, sd_prop) stats::rnorm(np, mu, sd_prop)
  }

  # u=replicate(1e4,rprop(2, c(0,0), chol_cov_prop ))
  # cov(t(u))
  # cov_prop

  # Number of parameters
  nparms <- 3L - length(theta_fix)
  # Draw from the prior
  theta_cur <- exp(stats::rnorm(nparms, 0.0, sd = sd_prior[-wh_fix]))
  log_theta_cur <- log(theta_cur)

  # Initialise the parameter vector, some are fix other we infer
  parms <- rep(0.0, 3L)
  parms[wh_fix] <- theta_fix
  parms[-wh_fix] <- theta_cur

  # Initialise BID class
  mod <- new(BDI, y, x0, n_particles, dt)

  # Container for the draws of the parameters
  niter <- ndpost + nskip
  draws <- matrix(0.0, nrow = niter, ncol = nparms)
  accept_rate <- 0L

  # Get estimate of the likelihood
  mod$RunBootstrapFilter(parms)
  lml_cur <- mod$lml

  # Start MCMC
  ini <- proc.time()
  for (k in seq_len(niter)) {
    if (k %% printevery == 0L) cat(k, "\n")
    # Sample from proposal at log scale
    log_theta_prop <- rprop(nparms, log_theta_cur, S_prop)
    theta_prop <- exp(log_theta_prop)

    # Compute log marginal likelihood using BF
    parms[-wh_fix] <- theta_prop
    mod$RunBootstrapFilter(parms)
    lml_prop <- mod$lml

    # log-prior
    # lprior_prop <- sum(stats::dnorm(x = log_theta_prop, mean = 0.0, sd = sd_prior,
    #                                 log = TRUE))
    # lprior_cur <- sum(stats::dnorm(x = log_theta_cur, mean = 0.0, sd = sd_prior,
    #                                log = TRUE))

    lprior_prop <- sum(-0.5 * (log_theta_prop/sd_prior)^2)
    lprior_cur <- sum(-0.5 * (log_theta_cur/sd_prior)^2)

    # MH ratio
    lr <- lml_prop - lml_cur# + lprior_prop - lprior_cur
    if (log(stats::runif(1L)) < lr) {
      # Accept
      theta_cur <- theta_prop
      log_theta_cur <- log(theta_cur)
      lml_cur <- lml_prop
      accept_rate <- accept_rate + 1L
    }
    draws[k, ] <- theta_cur
  }
  end <- proc.time() - ini
  attr(draws, "AcceptRate") <- accept_rate / niter
  attr(draws, "time") <- end
  return(draws)
}
