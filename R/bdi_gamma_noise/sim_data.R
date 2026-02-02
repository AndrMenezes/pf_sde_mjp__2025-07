rm(list = ls()); gc()
library(parallel)

path_local <- "./R/bdi_gamma_noise/"
# path_fig <- file.path(path_local, "figures")
source("./R/utils.R")
# source(file.path(path_local, "particleMCMC.R"))


# set.seed(1212)
set.seed(11)

# Example: Birth-Death-Immigration (BDI) model
# X --> 2X (birth: \mu)
# X --> 0 (death: \lambda)
# 0 --> X (immigration: \lambda)
# theta = (mu[birth], lambda[death], gamma[immigration])

# Settings ---------------------------------------------------------------------
tmax <- 30L
dt <- 1 / 30 # Use 1/dt observations for each time period
x0 <- 500L
true_theta <- c(0.5, 0.8, 40)
hazard_bdi <- function(x, theta) c(x * theta[1], x * theta[2], theta[3])
S_bdi <- matrix(c(1, -1, 1), nrow = 1)

# Simulating -------------------------------------------------------------------
sim_xt <- sim_poisson_leap(tmax = tmax, x0 = x0, theta = true_theta,
                           hazard = hazard_bdi, S = S_bdi, dt = dt)

# Adding a observational layer
xt <- sim_xt$x[, 1L]
t_continuous <- sim_xt$tt
n <- length(xt)
retain <- c(1, seq.int(11, length(t_continuous), length.out = tmax - 1))
t_retain <- t_continuous[retain]
xt_retain <- xt[retain]

par(mfrow = c(3, 2))
alpha <- c(0.5, 1.0, 3.0, 5.0, 10, 100)
for (k in seq_along(alpha)) {
  yt <- stats::rgamma(n = tmax, shape = alpha[k], rate = alpha[k] / xt_retain)

  plot(0:29, yt, col = "red", ylim = c(0, max(max(yt), max(xt))), pch = 19,
       ylab = "", xlab = "t", main = paste0("alpha = ", alpha[k]))
  lines(t_continuous, xt)
  lines(0:29, xt_retain, col = "blue", pch = 18)
}

