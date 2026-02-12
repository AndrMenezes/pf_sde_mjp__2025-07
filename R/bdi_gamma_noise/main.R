rm(list = ls()); gc()
library(parallel)

path_local <- "./R/bdi_gamma_noise/"
# path_fig <- file.path(path_local, "figures")
source("./R/utils.R")
Rcpp::sourceCpp(file = file.path(path_local, "bdi.cpp"))
source(file.path(path_local, "particleMCMC.R"))


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
alpha <- c(3.0, 5.0, 10, 50, 80, 100)
y <- matrix(0.0, nrow = tmax, ncol = length(alpha))
for (k in seq_along(alpha)) {
  y[, k] <- stats::rgamma(n = tmax, shape = alpha[k], rate = alpha[k] / xt_retain)
  plot(0:29, y[, k], col = "red", ylim = c(0, max(max(y[, k]), max(xt))), pch = 19,
       ylab = "", xlab = "t", main = paste0("alpha = ", alpha[k]))
  lines(t_continuous, xt)
  lines(0:29, xt_retain, col = "blue", pch = 18)
}
graphics.off()

id_ <- which(alpha == 10.0)
true_alpha <- alpha[id_]
yt <- y[, id_]
plot(ts(yt))
# true_alpha <- 50.0
# yt <- stats::rgamma(n = tmax, shape = true_alpha, rate = true_alpha / xt_retain)
# plot(0:29, yt, col = "red", ylim = c(0, max(max(yt), max(xt))), pch = 19,
#      ylab = "", xlab = "t", main = paste0("alpha = ", true_alpha))
# lines(t_continuous, xt)
# lines(0:29, xt_retain, col = "blue", pch = 18)


# Simple check
# Rcpp::sourceCpp(file = file.path(path_local, "bdi.cpp"))
n_particles <- 30L
max_trials <- 100L
target_success <- tmax
mod <- new(BDI, yt, x0, 1/10, true_alpha, n_particles, target_success, max_trials)
mod$RunFrankenFilter(true_theta)
mod$lml
mod$number_trials
mod$pct_m_max_reached
mod$RunBootstrapFilter(true_theta)#c(0.09, .8, 0.1)
mod$lml

# Estimate Var(log(p)) for different number of particles -----------------------

mc <- 2000L

# BF
N_particles <- c(10L, 20L, 50L, 100L, 200L, 300L)
var_lml <- mean_lml <- numeric(length = length(N_particles))
for (j in seq_along(N_particles)) {
  cat(j, "\n")
  lmls <- mclapply(seq_len(mc), function(i) {
    mod <- new(BDI, yt, x0, 1/10, true_alpha, N_particles[j], 1, 1)
    mod$RunBootstrapFilter(true_theta)
    return(mod$lml)
  }, mc.cores = 20L)
  lmls <- unlist(lmls)
  mean_lml[j] <- mean(lmls)
  var_lml[j] <- var(lmls)
}
tab_bf <- cbind(N_particles, var_lml, mean_lml)
tab_bf
plot(N_particles, var_lml)

# Franken filtering
TARGET_SUCCESS <- c(20L, 30L, 50L, 100L, 200L)
MAX_TRIALS <- c(1000L, 2000L)
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
tab_ff

# out <- replicate(1000, {
#   mod <- new(BDI, yt, x0, 1/10, true_alpha, 10, 21, 500)
#   mod$RunFrankenFilter(true_theta)
#   c(mod$lml, mean(mod$number_trials), max(mod$number_trials))
# })
# # do.call(rbind, out)
# t(out)
#
# ids <- which(is.infinite(out[1, ]))
# t(out[, ids])
# ids <- which(out[3, ] == 500)
# t(out[, ids])




# Particle MCMC to infer only \mu ----------------------------------------------

res_bf <- particleMCMC_bdi(y = yt, x0 = x0, dt = 1/10,
                        filter = "bootstrap",
                        n_particles = 30L,
                        shape = true_alpha,
                        ndpost = 10000L, nskip = 10000L,
                        theta_fix = true_theta[-1L],
                        wh_fix = c(2L, 3L),
                        sd_prior = 1.0,
                        S_prop = 0.5,
                        printevery = 100L, theta_init = true_theta[1L])

res_ff <- particleMCMC_bdi(y = yt, x0 = x0, dt = 1/10,
                           filter = "franken",
                           target_success = 30L, max_trials = 100L,
                           shape = true_alpha,
                           ndpost = 4000L, nskip = 5000L,
                           theta_fix = true_theta[-1L],
                           wh_fix = c(2L, 3L),
                           sd_prior = 1.0,
                           S_prop = 0.5,
                           printevery = 100L, theta_init = true_theta[1L])
plot(res_ff[, 1], type = "l"); abline(h = true_theta[1L], col = "red")
plot(res_bf[, 1], type = "l"); abline(h = true_theta[1L], col = "red")
ess_bf <- coda::effectiveSize(coda::as.mcmc(res_bf))
ess_ff <- coda::effectiveSize(coda::as.mcmc(res_ff))

attr(res_ff, "avg_trials")
ess_bf/attr(res_bf, "time")[3L]
ess_ff/attr(res_ff, "time")[3L]



# Particle MCMC to infer \mu and \gamma ----------------------------------------

res_bf <- particleMCMC_bdi(y = yt, x0 = x0, dt = 1/10,
                           filter = "bootstrap",
                           n_particles = 40L,
                           shape = true_alpha,
                           ndpost = 20000L, nskip = 20000L,
                           theta_fix = true_theta[2L],
                           wh_fix = 2L,
                           sd_prior = c(1.0, 1.0),
                           S_prop = diag(1, 2),
                           printevery = 100L, theta_init = true_theta[-2L])
par(mfrow = c(1, 2))
plot(res_bf[, 1], type = "l"); abline(h = true_theta[1L], col = 2)
plot(res_bf[, 2], type = "l"); abline(h = true_theta[3L], col = 2)

res_ff <- particleMCMC_bdi(y = yt, x0 = x0, dt = 1/10,
                           filter = "franken",
                           target_success = 30L, max_trials = 100L,
                           shape = true_alpha,
                           ndpost = 20000L, nskip = 20000L,
                           theta_fix = true_theta[2L],
                           wh_fix = 2L,
                           sd_prior = c(1.0, 1.0),
                           S_prop = diag(1, 2),
                           printevery = 100L, theta_init = true_theta[-2L])
# cov(log(res_ff))
# 2.5^2*cov(log(res_ff))/2

par(mfrow = c(1, 2))
plot(res_ff[, 1], type = "l"); abline(h = true_theta[1L], col = 2)
plot(res_ff[, 2], type = "l"); abline(h = true_theta[3L], col = 2)
plot(density(res_ff[, 1]))
plot(density(res_ff[, 2]))

head(cbind(res_ff, attr(res_ff, "lml")))
head(cbind(res_bf, attr(res_bf, "lml")))


par(mfrow = c(2, 2))
plot(res_bf[, 1], type = "l"); abline(h = true_theta[1L], col = 2)
plot(res_bf[, 2], type = "l"); abline(h = true_theta[3L], col = 2)
plot(res_ff[, 1], type = "l"); abline(h = true_theta[1L], col = 2)
plot(res_ff[, 2], type = "l"); abline(h = true_theta[3L], col = 2)
ess_bf <- coda::effectiveSize(coda::as.mcmc(res_bf))
ess_ff <- coda::effectiveSize(coda::as.mcmc(res_ff))

ess_bf/attr(res_bf, "time")[3L]
ess_ff/attr(res_ff, "time")[3L]



# Particle MCMC to perform inference on \mu ------------------------------------

n_particles <- 600L
max_trials <- 4000L
target_success <- 200L

# Particle MCMC using BF
mod1 <- new(BDI, yt, x0, 1/10, true_alpha, n_particles, 1, 1)
system.time(
  mod1$RunParticleMCMC_mu(true_theta[2L], true_theta[3L], 1.0, 0.04, 5000L))
plot(mod1$draws_mu[-(1:1000L)], type = "l"); abline(h = true_theta[1], col = "red")
plot(mod1$draws_mu[-(1:1000L)], type = "l"); abline(h = true_theta[1], col = "red")
sd_opt <- sd(log(mod1$draws_mu[-(1:1000L)]))
coda::effectiveSize(coda::as.mcmc(mod1$draws_mu[-(1:1000L)]))

# Particle MCMC using FF
mod2 <- new(BDI, yt, x0, 1/10, true_alpha, n_particles, target_success, max_trials)
system.time(
  mod2$RunParticleMCMCFranken_mu(true_theta[2L], true_theta[3L], 1.0, 0.04, 5000L))
plot(mod2$draws_mu[-(1:1000L)], type = "l"); abline(h = true_theta[1], col = "red")
plot(mod2$draws_mu, type = "l"); abline(h = true_theta[1], col = "red")
sd_opt <- sd(log(mod2$draws_mu[-(1:1000L)]))

coda::effectiveSize(coda::as.mcmc(mod1$draws_mu))
coda::effectiveSize(coda::as.mcmc(mod2$draws_mu))






