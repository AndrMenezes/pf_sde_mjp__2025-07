#' Run particle MCMC to perform inference on the parameters of BDI system with noise observations
#' generated from a gamma distribution.
#' @param y observed data.
#' @param x0 initial vale for the state.
#' @param dt discretised time interval, where 1/dt give number of points to discretise
#' the interval (t, t+1]. This is used to simulate from the hiden process using the tau-leap scheme.
#' @param filter character. Indicates which particle filter use. The options are: "\code{bootstrap}"
#' for Bootstrap filter and "\code{franken}" for the Franken filter.
#' @param n_particles number of particles use in the Bootstrap filter.
#' @param target_success,max_trials target success and the maximum number of trials
#' for the Franken Filter.
#' @param shape The shape parameter \eqn{(\alpha)} of gamma distribution for the
#' noise observations, i.e., \eqn{y_t | x_t \sim Gamma[\alpha, \alpha/x_t]}.
#' @param ndpost,nskip number of posterior and skip draws.
#' @param theta_fix,wh_fix provides the value of the rate parameters of the BDI
#' model, and which \eqn{(\mu, \lambda, \gamma)} parameters are considered fixed.
#' We recall that the model is non-identifiable if you try to estimate all the three parameters.
#' @param sd_prior standard deviation for the normal prior on the log-scale of the
#' parameters.
#' @param sd_prop standard deviation of the normal proposal.
#' @param printevery integer to print number of iteration.
#' @description
#' This function run the particle MCMC algorithm to perform inference on the
#' Birth-Death-Immigration process.
#'
particleMCMC_bdi <- function(y, x0, dt,
                             filter = c("bootstrap", "franken"),
                             n_particles = 50L,
                             target_success = 10L,
                             max_trials = 50L,
                             shape = 1.0,
                             ndpost = 1000L, nskip = 1000L, nadapt = floor(nskip/2),
                             theta_fix = c(0.8, 40.0),
                             wh_fix = c(2L, 3L),
                             sd_prior = rep(1.0, 3L - length(theta_fix)),
                             S_prop = rep(0.25, 3L - length(theta_fix)),
                             printevery = 10L,
                             theta_init = NULL) {
  filter <- match.arg(filter)

  if (is.matrix(S_prop)) {
    S_prop <- t(chol(S_prop))
    rprop <- function(np, mu, chol_cov_prop) drop(mu + chol_cov_prop %*% stats::rnorm(np))
  } else {
    rprop <- function(np, mu, sd_prop) stats::rnorm(np, mu, sd_prop)
  }

  # Number of parameters
  nparms <- 3L - length(theta_fix)
  # Draw from the prior
  if (is.null(theta_init)) theta_cur <- exp(stats::rnorm(nparms, 0.0, sd = sd_prior[-wh_fix]))
  else theta_cur <- theta_init
  log_theta_cur <- log(theta_cur)

  # Initialise the parameter vector, some are fix other we infer
  parms <- rep(0.0, 3L)
  parms[wh_fix] <- theta_fix
  parms[-wh_fix] <- theta_cur

  # Initialise BID class
  mod <- new(BDI, y, x0, dt, shape, n_particles, target_success, max_trials)
  # Container for the draws of the parameters
  niter <- ndpost + nskip
  draws <- matrix(0.0, nrow = ndpost, ncol = nparms)
  accept_rate <- 0L

  ini <- proc.time()
  # Get estimate of the likelihood
  if (filter == "bootstrap") mod$RunBootstrapFilter(parms)
  else mod$RunFrankenFilter(parms)
  lml_cur <- mod$lml

  # Run the burn-in/adapt phase to estimate the covariance of the proposal
  if (nskip > 0) {
    draws_adapt <- matrix(0.0, nrow = nadapt, ncol = nparms)
    for (k in seq_len(nskip)) {
      if (k %% printevery == 0L) cat("Adapt phase:", k, "of", nskip, "\n")
      # Sample from proposal at log scale
      log_theta_prop <- rprop(nparms, log_theta_cur, S_prop)
      theta_prop <- exp(log_theta_prop)
      # Compute log marginal likelihood using BF
      parms[-wh_fix] <- theta_prop
      if (filter == "bootstrap") mod$RunBootstrapFilter(parms)
      else mod$RunFrankenFilter(parms)
      lml_prop <- mod$lml
      # log-prior
      lprior_prop <- sum(-0.5 * (log_theta_prop/sd_prior)^2)
      lprior_cur <- sum(-0.5 * (log_theta_cur/sd_prior)^2)
      # MH ratio
      lr <- lml_prop - lml_cur + lprior_prop - lprior_cur
      if (log(stats::runif(1L)) < lr) {
        # Accept
        theta_cur <- theta_prop
        log_theta_cur <- log(theta_cur)
        lml_cur <- lml_prop
        # accept_rate <- accept_rate + 1L
      }
      if (k > nadapt) draws_adapt[k - nadapt, ] <- log_theta_cur
    }
    # Estimate the proposal variance
    if (is.matrix(S_prop)) {
      S_prop <- t(chol(2.5*2.5*cov(draws_adapt) / nparms))
    } else {
      S_prop <- 1.5*sd(draws_adapt)
    }
  }

  # Run the posterior sampling
  number_trials <- 0.0
  # lml <- numeric(ndpost)
  for (k in seq_len(ndpost)) {
    if (k %% printevery == 0L) cat("Posterior sampling:", k,  "of", ndpost, "\n")
    # Sample from proposal at log scale
    log_theta_prop <- rprop(nparms, log_theta_cur, S_prop)
    theta_prop <- exp(log_theta_prop)
    # Compute log marginal likelihood using BF
    parms[-wh_fix] <- theta_prop
    if (filter == "bootstrap") mod$RunBootstrapFilter(parms)
    else mod$RunFrankenFilter(parms)
    lml_prop <- mod$lml


    # log-prior
    lprior_prop <- sum(-0.5 * (log_theta_prop/sd_prior)^2)
    lprior_cur <- sum(-0.5 * (log_theta_cur/sd_prior)^2)
    # MH ratio
    lr <- lml_prop - lml_cur + lprior_prop - lprior_cur

    # cat(theta_cur, theta_prop, lml_prop, lml_cur, lr, "\n")

    if (log(stats::runif(1L)) < lr) {
      # Accept
      theta_cur <- theta_prop
      log_theta_cur <- log(theta_cur)
      lml_cur <- lml_prop
      accept_rate <- accept_rate + 1L
    }
    draws[k, ] <- theta_cur
    number_trials <- number_trials + mean(mod$number_trials)
    # lml[k] <- lml_cur
  }
  end <- proc.time() - ini
  attr(draws, "acceptance_rate") <- accept_rate / ndpost
  attr(draws, "time") <- end
  attr(draws, "avg_trials") <- number_trials / ndpost
  # attr(draws, "lml") <- lml_cur
  return(draws)
}
