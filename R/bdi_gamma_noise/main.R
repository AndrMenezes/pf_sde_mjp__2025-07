rm(list = ls()); gc()
library(parallel)

path_local <- "./R/bdi_gamma_noise/"
# path_fig <- file.path(path_local, "figures")
source("./R/utils.R")
Rcpp::sourceCpp(file = file.path(path_local, "bdi.cpp"))
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
graphics.off()
#
true_alpha <- 50.0
yt <- stats::rgamma(n = tmax, shape = true_alpha, rate = true_alpha / xt_retain)
plot(0:29, yt, col = "red", ylim = c(0, max(max(yt), max(xt))), pch = 19,
     ylab = "", xlab = "t", main = paste0("alpha = ", true_alpha))
lines(t_continuous, xt)
lines(0:29, xt_retain, col = "blue", pch = 18)


# Simple check
# Rcpp::sourceCpp(file = file.path(path_local, "bdi.cpp"))
n_particles <- 600L
max_trials <- 4000L
target_success <- 200L
mod <- new(BDI, yt, x0, 1/10, true_alpha, n_particles, target_success, max_trials)
mod$RunBootstrapFilter(true_theta)
mod$lml
mod$RunFrankenFilter(true_theta)
mod$lml
mod$number_trials

# Estimate Var(log(p)) for different number of particles -----------------------

mc <- 2000L

# BF
N_particles <- c(10L, 20L, 30L, 35, 40L, 50L, 100L, 200L, 300L)
list_lmls <- vector(mode = "list", length = length(N_particles))
var_lml <- numeric(length = length(N_particles))
for (j in seq_along(N_particles)) {
  cat(j, "\n")
  mod <- new(BDI, yt, x0, 1/10, true_alpha, N_particles[j], target_success, max_trials)
  lmls <- mclapply(seq_len(mc), function(i) {
    mod$RunBootstrapFilter(true_theta)
    return(mod$lml)
  }, mc.cores = 20L)
  lmls <- unlist(lmls)
  list_lmls[[j]] <- lmls
  var_lml[j] <- var(lmls)
}
tab_bf <-cbind(N_particles, var_lml, mean_lml = sapply(list_lmls, function(x) mean(x)))
plot(N_particles, var_lml)

# Franken filtering
TARGET_SUCCESS <- c(21L, 51L, 101L, 201L)
MAX_TRIALS <- c(100L, 500L, 1000L, 2000L)
grid <- expand.grid(s = TARGET_SUCCESS, mp = MAX_TRIALS)
nr <- nrow(grid)
var_lml_franken <- numeric(length = nr)
mean_lml_franken <- numeric(length = nr)
for (j in seq_len(nr)) {
  cat(j, "\n")
  s <-  grid$s[j]; M <-  grid$m[j]
  lmls <- mclapply(seq_len(mc), function(i) {
    mod <- new(BDI, yt, x0, 1/10, true_alpha, 10, s, M)
    mod$RunFrankenFilter(true_theta)
    return(mod$lml)
  }, mc.cores = 20L)
  lmls <- unlist(lmls)
  mean_lml_franken[j] <- mean(lmls)
  var_lml_franken[j] <- var(lmls)
}
tab_ff <- cbind(grid, var_lml_franken, exp_lml = mean_lml_franken)

out <- replicate(1000, {
  mod <- new(BDI, yt, x0, 1/10, true_alpha, 10, 21, 500)
  mod$RunFrankenFilter(true_theta)
  c(mod$lml, mean(mod$number_trials), max(mod$number_trials))
})
# do.call(rbind, out)
t(out)

ids <- which(is.infinite(out[1, ]))
t(out[, ids])
ids <- which(out[3, ] == 500)
t(out[, ids])


# Run particle MCMC using BF ---------------------------------------------------
n_particles <- 600L
max_trials <- 4000L
target_success <- 200L

# Particle MCMC using BF
mod1 <- new(BDI, yt, x0, 1/10, true_alpha, n_particles, 1, 1)
system.time(
  mod1$RunParticleMCMC_mu(true_theta[2L], true_theta[3L], 1.0, 0.04, 5000L))
plot(mod1$draws_mu[-(1:1000L)], type = "l"); abline(h = true_theta[1], col = "red")
sd_opt <- sd(log(mod1$draws_mu[-(1:1000L)]))
coda::effectiveSize(coda::as.mcmc(mod1$draws_mu[-(1:1000L)]))/60

# Particle MCMC using FF
mod2 <- new(BDI, yt, x0, 1/10, true_alpha, n_particles, target_success, max_trials)
system.time(
  mod2$RunParticleMCMCFranken_mu(true_theta[2L], true_theta[3L], 1.0, 0.04, 5000L))
plot(mod2$draws_mu[-(1:1000L)], type = "l"); abline(h = true_theta[1], col = "red")
sd_opt <- sd(log(mod2$draws_mu[-(1:1000L)]))
coda::effectiveSize(coda::as.mcmc(mod1$draws_mu[-(1:1000L)]))/60






