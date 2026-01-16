
#' Exact algorithm due to Gillespie to simulate coupled chemical reaction systems.
#' @param theta vector of parameters that control the system.
#' @param hazard a function-vector for the hazard function associated to the system.
#' @param x0 initial condition for the state.
#' @param S stoichiometry matrix.
#' @param tmax maximum time.
#' @description Right now, this function only works for unidimensional states.
sim_gillespie <- function(theta, x0, hazard, S, tmax) {
  x <- numeric()
  tt <- numeric()
  x[1L] <- x0
  tt[1L] <- 0L
  # state dimension and number of reactions
  d <- length(x0)
  u <- length(theta)
  # Compute hazard and combined hazard for initial state
  h <- hazard(x0, theta)
  h0 <- sum(h)
  i <- 2L
  k <- 0.0
  while (k < tmax) {
    t_ <- stats::rexp(n = 1, rate = h0)
    j <- sample.int(n = u, size = 1L, prob = h)
    tt[i] <- tt[i - 1] + t_
    x[i] <- x[i - 1] + S[, j]
    k <- tt[i]
    h <- hazard(x[i], theta)
    h0 <- sum(h0)
    i <- i + 1L
  }
  list(x = x, tt = tt)
}

#' First reaction method to simulate coupled chemical reaction systems.
#' @param theta vector of parameters that control the system.
#' @param hazard a function-vector for the hazard function associated to the system.
#' @param x0 initial condition for the state.
#' @param S stoichiometry matrix.
#' @param tmax maximum time.
#' @description Right now, this function only works for unidimensional states.
sim_frm <- function(theta, x0, hazard, S, tmax) {
  x <- numeric()
  tt <- numeric()
  x[1L] <- x0
  tt[1L] <- 0L
  # state dimension and number of reactions
  d <- length(x0)
  u <- length(theta)
  # Compute hazard and generate putative times
  h <- hazard(x0, theta)
  # Auxliary indices
  i <- 2L
  k <- 0.0
  while (k < tmax) {
    t_ <- stats::rexp(n = u, rate = h)
    j <- which.min(t_)
    x[i] <- x[i - 1] + S[, j]
    tt[i] <- tt[i - 1] + t_[j]
    k <- tt[i]
    h <- hazard(x[i], theta)
    i <- i + 1L
  }
  list(x = x, tt = tt)
}

sim_poisson_leap <- function(theta, x0, hazard, S, tmax, dt = 0.1) {
  d <- length(x0)
  u <- length(theta)
  n <- tmax %/% dt
  x <- matrix(0.0, nrow = n, ncol = d)
  x_cur <- x0
  x[1L, ] <- x0
  for (i in seq_len(n)) {
    h <- hazard(x_cur, theta)
    r <- stats::rpois(n = u, lambda = h * dt)
    x_cur <- x_cur + drop(S %*% r)
    x[i, ] <- x_cur
  }
  list(x = x, tt = seq(0, tmax, length.out = n))
}

#' ODE solution for the Birth-Death-Immigration model
#' @param theta vector of parameters that control the model behaviour.
#' @param tt time expand.
#' @param x0 initial value for the population.
ode_solution_bdi <- function(theta, tt, x0) {
  mu <- theta[1]
  lambda <- theta[2]
  gama <- theta[3]
  gama / (lambda - mu) + (x0 - gama / (lambda - mu)) * exp((mu - lambda) * tt)
}





# n_obs <- 10L
# n_particles <- 1L
# n_dt <- 10L
# x_state <- matrix(data = 0.0, nrow = n_obs * n_dt, ncol = n_particles)
# for (i in seq_len(n_obs)) {
#   for (j in seq_len(n_particles)) {
#     # Call Poisson-leap for each particles, then appending
#     for (k in seq_len(n_dt)) {
#       row_index <- (i - 1) * n_dt + k
#       cat(row_index, "\n")
#       # x_state[i * k, j]
#     }
#   }
# }
