#' Run particle MCMC to perform inference on the parameters of pure Death system.
#' @param x observed exact observations.
#' @param x0 initial vale for the state.
#' @param n_particles number of particles.
#' @param dt discretised time interval, where 1/dt give number of points to discretise
#' the interval (t, t+1].
#' @param ndpost,nskip number of posterior and skip draws.
#' @param sd_prior standard deviation for the normal prior on the log-scale of the
#' parameters.
#' @param sd_prop standard deviation of the normal proposal.
#' @param printevery integer to print number of iteration.
#'
particleMCMC_death <- function(x, x0,
                               n_particles = 150L,
                               max_trials = 1e4L,
                               target_sucess = 51L,
                               dt, ndpost, nskip,
                               filter = c("bootstrap", "alive", "franken"),
                               shape_prior = 10L, rate_prior = 1000L,
                               sd_prop = 1.0,
                               printevery = 10L) {

  filter <- match.arg(filter)
  # Draw from the prior
  theta_cur <- stats::rgamma(1L, shape = shape_prior, rate = rate_prior)
  log_theta_cur <- log(theta_cur)
  # cat(theta_cur, "\n")

  # Initialise Death class
  mod <- new(Death, x, x0, dt, n_particles, target_sucess, max_trials)

  # Container for the draws of the parameters
  niter <- ndpost + nskip
  draws <- numeric(niter)
  accept_rate <- 0L

  # Get estimate of the likelihood

  lml_cur <- switch(filter,
    "bootstrap" = mod$RunBootstrapFilter(theta_cur),
    "alive" = mod$RunAliveFilter(theta_cur),
    "franken" = mod$RunFrankenFilter(theta_cur)
  )

  # Start MCMC
  ini <- proc.time()
  for (k in seq_len(niter)) {
    if (k %% printevery == 0L) cat(k, "\n")
    # Sample from proposal at log scale
    log_theta_prop <- stats::rnorm(1L, mean = log_theta_cur, sd = sd_prop)
    theta_prop <- exp(log_theta_prop)

    # Compute log marginal likelihood using BF
    # mod$RunBootstrapFilter(theta_prop, n_particles)
    # cat(theta_prop, "\n")
    lml_prop <- switch(filter,
                      "bootstrap" = mod$RunBootstrapFilter(theta_prop),
                      "alive" = mod$RunAliveFilter(theta_prop),
                      "franken" = mod$RunFrankenFilter(theta_prop)
    )

    # log-prior
    # lprior_prop <- -0.5 * (log_theta_prop/sd_prior)^2
    # lprior_cur <- -0.5 * (log_theta_cur/sd_prior)^2
    lprior_prop <- stats::dgamma(x = theta_prop, shape = shape_prior,
                                 rate = rate_prior, log = TRUE)
    lprior_cur <- stats::dgamma(x = theta_cur, shape = shape_prior,
                                rate = rate_prior, log = TRUE)

    # MH ratio
    lr <- lml_prop - lml_cur + lprior_prop - lprior_cur
    if (log(stats::runif(1L)) < lr) {
      # Accept
      theta_cur <- theta_prop
      log_theta_cur <- log(theta_cur)
      lml_cur <- lml_prop
      accept_rate <- accept_rate + 1L
    }
    draws[k] <- theta_cur
  }
  end <- proc.time() - ini
  attr(draws, "AcceptRate") <- accept_rate / niter
  attr(draws, "time") <- end
  return(draws)
}

###########
#' Exact log-likelihood for pure death process
.lml_death <- function(x, theta) {
  tmax <- length(x)
  sum(stats::dbinom(x = x[-1], size = x[-tmax], prob = exp(-theta), log = TRUE))
}
exactMCMC_death <- function(y, x0, niter = 50000L, shape_prior = 10L,
                            rate_prior = 1000L, sd_prop = 1.0, printevery = 10L) {

  # Draw from the prior
  theta_cur <- stats::rgamma(1L, shape = shape_prior, rate = rate_prior)
  log_theta_cur <- log(theta_cur)
  # Container for the draws of the parameters
  draws <- numeric(niter)
  accept_rate <- 0L

  # Get estimate of the likelihood
  lml_cur <- .lml_death(x = y, theta = theta_cur)

  # Start MCMC
  ini <- proc.time()
  for (k in seq_len(niter)) {
    if (k %% printevery == 0L) cat(k, "\n")
    # Sample from proposal at log scale
    log_theta_prop <- stats::rnorm(1L, mean = log_theta_cur, sd = sd_prop)
    theta_prop <- exp(log_theta_prop)

    # Compute log likelihood
    lml_prop <- .lml_death(x = y, theta = theta_prop)
    # log-prior
    lprior_prop <- stats::dgamma(x = theta_prop, shape = shape_prior,
                                 rate = rate_prior, log = TRUE)
    lprior_cur <- stats::dgamma(x = theta_cur, shape = shape_prior,
                                rate = rate_prior, log = TRUE)

    # MH ratio
    lr <- lml_prop - lml_cur + lprior_prop - lprior_cur
    if (log(stats::runif(1L)) < lr) {
      # Accept
      theta_cur <- theta_prop
      log_theta_cur <- log(theta_cur)
      lml_cur <- lml_prop
      accept_rate <- accept_rate + 1L
    }
    draws[k] <- theta_cur
  }
  end <- proc.time() - ini
  attr(draws, "AcceptRate") <- accept_rate / niter
  attr(draws, "time") <- end
  return(draws)
}


