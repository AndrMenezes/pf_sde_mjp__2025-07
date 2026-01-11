############################################################
# SIR Filter for SV model
#
# Model:
#   x_t = phi x_{t-1} + N(0, sigma^2)
#   y_t ~ N(0, exp{2(gamma + x_t)})
#
# Inputs:
#   N      : number of particles
#   y      : observations
#   par    : (phi, sigma, gamma)
#   prior  : prior mean and sd for x_1
#   truth  : true latent state (optional, for plotting)
#   PLOT   : logical, produce plots?
#
# Output:
#   Filtered mean, variance, and log-likelihood
############################################################

SIRfilterSV <- function(N, y, par, prior, truth = NULL, PLOT = TRUE) {

  n <- length(y)
  d <- 1

  ## (1) Initialize particles
  particles <- matrix(
    rnorm(N, prior[1], prior[2]),
    ncol = d
  )

  ## (2) Initial weights
  w <- dnorm(y[1], 0, exp(par[3] + particles[, 1]))
  logl <- log(mean(w))
  w <- w / sum(w)

  ## Storage
  M.st <- V.st <- matrix(0, nrow = n, ncol = d)
  M.st[1, 1] <- sum(w * particles[, 1])
  V.st[1, 1] <- sum(w * particles[, 1]^2) - M.st[1, 1]^2

  ## Filtering loop
  for (i in 2:n) {

    ## (1) Resample
    idx <- sample(1:N, size = N, replace = TRUE, prob = w)
    particles <- matrix(particles[idx, ], ncol = d)

    ## (2) Propagate
    particles[, 1] <- par[1] * particles[, 1] +
      rnorm(N, 0, par[2])

    ## (3) Reweight
    w <- dnorm(y[i], 0, exp(par[3] + particles[, 1]))
    logl <- logl + log(mean(w))
    w <- w / sum(w)

    ## Store moments
    M.st[i, 1] <- sum(w * particles[, 1])
    V.st[i, 1] <- sum(w * particles[, 1]^2) - M.st[i, 1]^2

  }

  list(mean = M.st, var = V.st, l = logl)
}

n <- 100L
true_parms <- c(0.9, sqrt(1 - 0.9^2), 1)
x <- rep(0, n)
x[1] <- rnorm(1,0,true_parms[2]/sqrt(1-true_parms[1]^2))
for (i in 2:n) x[i] = x[i - 1]*true_parms[1] + rnorm(1,0,true_parms[2])
y <- rnorm(n, 0, exp(true_parms[3] + x))
plot(y, type = "l")
lines(x, type = "l", col = "blue")

out <- SIRfilterSV(N = 100, y = y, par = true_parms,
                   prior = c(0, 1))

out$mean

plot(x, type = "l")
lines(out$mean, type = "l", col = "red")


