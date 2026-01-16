# Model:
#   x_0 ~ N[0, tau^2 / (1 - phi^2)]
#   x_t ~ N[phi x_{t-1}, tau^2]
#   y_t ~ N[x_{t}, sigma^2]

BFfilterAR1 <- function(y, parms, initial_ss, n_particles) {

  n <- length(y)
  phi <- parms[1]
  sigma <- parms[2]
  tau <- parms[3]

  # Initialize particles
  particles <- stats::rnorm(n_particles, initial_ss[1], sd = initial_ss[2])

  ## Initial weights
  w <- stats::dnorm(y[1], particles, sd = sigma)
  logl <- log(mean(w))
  w <- w / sum(w)

  ## Storage
  M.st <- V.st <- numeric(n) #matrix(0, nrow = n, ncol = d)
  M.st[1] <- sum(w * particles)
  V.st[1] <- sum(w * particles^2) - M.st[1]^2

  ## Filtering loop
  for (i in 2:n) {

    ## (1) Resample
    idx <- sample(1:n_particles, size = n_particles, replace = TRUE, prob = w)
    particles <- particles[idx]

    ## (2) Propagate
    particles <- stats::rnorm(n_particles, phi * particles, sd = tau)

    ## (3) Reweight
    w <- dnorm(y[i], particles, sd = sigma)
    logl <- logl + log(mean(w))
    w <- w / sum(w)

    ## Store moments
    M.st[i] <- sum(w * particles)
    V.st[i] <- sum(w * particles^2) - M.st[i]^2
  }

  list(mean = M.st, var = V.st, lml = logl)
}
