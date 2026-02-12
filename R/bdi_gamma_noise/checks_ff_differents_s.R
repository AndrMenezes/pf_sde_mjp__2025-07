rm(list = ls()); gc()
library(parallel)

path_local <- "./R/bdi_gamma_noise/"
# path_fig <- file.path(path_local, "figures")
source("./R/utils.R")
Rcpp::sourceCpp(file = file.path(path_local, "bdi.cpp"))
source(file.path(path_local, "particleMCMC.R"))


set.seed(696)

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

par(mfrow = c(1, 3))
true_alpha <- c(5.0, 10, 50)
y <- matrix(0.0, nrow = tmax, ncol = length(true_alpha))
for (k in seq_along(true_alpha)) {
  y[, k] <- stats::rgamma(n = tmax, shape = true_alpha[k],
                          rate = true_alpha[k] / xt_retain)
  plot(0:29, y[, k], col = "red", ylim = c(0, max(max(y[, k]), max(xt))),
       pch = 19, ylab = "", xlab = "t", main = paste0("alpha = ", true_alpha[k]))
  lines(t_continuous, xt)
  lines(0:29, xt_retain, col = "blue", pch = 18)
}

# Settings
MAX_TRIALS <- 400L
NDPOST <- 10000L
NSKIP <- 10000L


# FF for \alpha = 5.0 ----------------------------------------------------------

res_10_5 <- particleMCMC_bdi(y = y[, 1L], x0 = x0, dt = 1/10,
                              filter = "franken",
                              target_success = 10L,
                              max_trials = MAX_TRIALS,
                              shape = true_alpha[1L],
                              ndpost = NDPOST, nskip = NSKIP,
                              theta_fix = true_theta[-1L],
                              wh_fix = c(2L, 3L),
                              sd_prior = 1.0,
                              S_prop = 0.5,
                              printevery = 100L, theta_init = true_theta[1L])
attr(res_10_5, "avg_trials")
plot(ts(res_10_5[, 1L])); abline(h = true_theta[1L], col = "red")

res_30_5 <- particleMCMC_bdi(y = y[, 1L], x0 = x0, dt = 1/10,
                              filter = "franken",
                              target_success = 30L,
                              max_trials = MAX_TRIALS,
                              shape = true_alpha[1L],
                              ndpost = NDPOST, nskip = NSKIP,
                              theta_fix = true_theta[-1L],
                              wh_fix = c(2L, 3L),
                              sd_prior = 1.0,
                              S_prop = 0.5,
                              printevery = 100L, theta_init = true_theta[1L])
attr(res_30_5, "avg_trials")
plot(ts(res_30_5[, 1L])); abline(h = true_theta[1L], col = "red")

res_50_5 <- particleMCMC_bdi(y = y[, 1L], x0 = x0, dt = 1/10,
                              filter = "franken",
                              target_success = 50L,
                              max_trials = MAX_TRIALS,
                              shape = true_alpha[1L],
                              ndpost = NDPOST, nskip = NSKIP,
                              theta_fix = true_theta[-1L],
                              wh_fix = c(2L, 3L),
                              sd_prior = 1.0,
                              S_prop = 0.5,
                              printevery = 100L, theta_init = true_theta[1L])

attr(res_50_5, "avg_trials")
plot(ts(res_50_5[, 1L])); abline(h = true_theta[1L], col = "red")

res_70_5 <- particleMCMC_bdi(y = y[, 1L], x0 = x0, dt = 1/10,
                              filter = "franken",
                              target_success = 70L,
                              max_trials = MAX_TRIALS,
                              shape = true_alpha[1L],
                              ndpost = NDPOST, nskip = NSKIP,
                              theta_fix = true_theta[-1L],
                              wh_fix = c(2L, 3L),
                              sd_prior = 1.0,
                              S_prop = 0.5,
                              printevery = 100L, theta_init = true_theta[1L])

attr(res_70_5, "avg_trials")
plot(ts(res_70_5[, 1L])); abline(h = true_theta[1L], col = "red")

res_100_5 <- particleMCMC_bdi(y = y[, 1L], x0 = x0, dt = 1/10,
                               filter = "franken",
                               target_success = 100L,
                               max_trials = MAX_TRIALS,
                               shape = true_alpha[1L],
                               ndpost = NDPOST, nskip = NSKIP,
                               theta_fix = true_theta[-1L],
                               wh_fix = c(2L, 3L),
                               sd_prior = 1.0,
                               S_prop = 0.5,
                               printevery = 100L, theta_init = true_theta[1L])

res_bf_5 <- particleMCMC_bdi(y = y[, 1L], x0 = x0, dt = 1/10,
                             filter = "bootstrap",n_particles = 100L,
                             shape = true_alpha[1L],
                             ndpost = NDPOST, nskip = NSKIP,
                             theta_fix = true_theta[-1L],
                             wh_fix = c(2L, 3L),
                             sd_prior = 1.0,
                             S_prop = 0.5,
                             printevery = 100L, theta_init = true_theta[1L])
plot(ts(res_bf_5[, 1L])); abline(h = true_theta[1L], col = "red")

coda::effectiveSize(coda::as.mcmc(res_10_5))
coda::effectiveSize(coda::as.mcmc(res_30_5))
coda::effectiveSize(coda::as.mcmc(res_50_5))
coda::effectiveSize(coda::as.mcmc(res_70_5))
coda::effectiveSize(coda::as.mcmc(res_100_5))

e <- c(
  coda::effectiveSize(coda::as.mcmc(res_10_5)) / attr(res_10_5, "time")[3L],
  coda::effectiveSize(coda::as.mcmc(res_30_5)) / attr(res_30_5, "time")[3L],
  coda::effectiveSize(coda::as.mcmc(res_50_5)) / attr(res_50_5, "time")[3L],
  coda::effectiveSize(coda::as.mcmc(res_70_5)) / attr(res_70_5, "time")[3L],
  coda::effectiveSize(coda::as.mcmc(res_100_5)) / attr(res_100_5, "time")[3L]
)
s <- c(10L, 30L, 50L, 70L, 100L)
tab_5 <- cbind(s, e)
plot(s, e, type = "o", pch = 19)

par(mfrow = c(2, 3))
plot(ts(res_10_5[, 1L]), main = "s = 10"); abline(h = true_theta[1L], col = "red")
plot(ts(res_30_5[, 1L]), main = "s = 30"); abline(h = true_theta[1L], col = "red")
plot(ts(res_50_5[, 1L]), main = "s = 50"); abline(h = true_theta[1L], col = "red")
plot(ts(res_70_5[, 1L]), main = "s = 70"); abline(h = true_theta[1L], col = "red")
plot(ts(res_100_5[, 1L]), main = "s = 100"); abline(h = true_theta[1L], col = "red")
plot(ts(res_bf_5[, 1L]), main = "BF"); abline(h = true_theta[1L], col = "red")
graphics.off()

plot(density(res_10_5))
plot(density(res_bf_5))


# FF for \alpha = 50 -----------------------------------------------------------

res_10_50 <- particleMCMC_bdi(y = y[, 3L], x0 = x0, dt = 1/10,
                              filter = "franken",
                              target_success = 10L,
                              max_trials = MAX_TRIALS,
                              shape = true_alpha[3L],
                              ndpost = NDPOST, nskip = NSKIP,
                              theta_fix = true_theta[-1L],
                              wh_fix = c(2L, 3L),
                              sd_prior = 1.0,
                              S_prop = 0.5,
                              printevery = 100L, theta_init = true_theta[1L])
attr(res_10_50, "avg_trials")
plot(ts(res_10_50[, 1L])); abline(h = true_theta[1L], col = "red")

res_30_50 <- particleMCMC_bdi(y = y[, 3L], x0 = x0, dt = 1/10,
                              filter = "franken",
                              target_success = 30L,
                              max_trials = MAX_TRIALS,
                              shape = true_alpha[3L],
                              ndpost = NDPOST, nskip = NSKIP,
                              theta_fix = true_theta[-1L],
                              wh_fix = c(2L, 3L),
                              sd_prior = 1.0,
                              S_prop = 0.5,
                              printevery = 100L, theta_init = true_theta[1L])
attr(res_30_50, "avg_trials")
plot(ts(res_30_50[, 1L])); abline(h = true_theta[1L], col = "red")

res_50_50 <- particleMCMC_bdi(y = y[, 3L], x0 = x0, dt = 1/10,
                              filter = "franken",
                              target_success = 50L,
                              max_trials = MAX_TRIALS,
                              shape = true_alpha[3L],
                              ndpost = NDPOST, nskip = NSKIP,
                              theta_fix = true_theta[-1L],
                              wh_fix = c(2L, 3L),
                              sd_prior = 1.0,
                              S_prop = 0.5,
                              printevery = 100L, theta_init = true_theta[1L])

attr(res_50_50, "avg_trials")
plot(ts(res_50_50[, 1L])); abline(h = true_theta[1L], col = "red")

res_70_50 <- particleMCMC_bdi(y = y[, 3L], x0 = x0, dt = 1/10,
                              filter = "franken",
                              target_success = 70L,
                              max_trials = MAX_TRIALS,
                              shape = true_alpha[3L],
                              ndpost = NDPOST, nskip = NSKIP,
                              theta_fix = true_theta[-1L],
                              wh_fix = c(2L, 3L),
                              sd_prior = 1.0,
                              S_prop = 0.5,
                              printevery = 100L, theta_init = true_theta[1L])

attr(res_70_50, "avg_trials")
plot(ts(res_70_50[, 1L])); abline(h = true_theta[1L], col = "red")

res_100_50 <- particleMCMC_bdi(y = y[, 3L], x0 = x0, dt = 1/10,
                              filter = "franken",
                              target_success = 100L,
                              max_trials = MAX_TRIALS,
                              shape = true_alpha[3L],
                              ndpost = NDPOST, nskip = NSKIP,
                              theta_fix = true_theta[-1L],
                              wh_fix = c(2L, 3L),
                              sd_prior = 1.0,
                              S_prop = 0.5,
                              printevery = 100L, theta_init = true_theta[1L])

attr(res_100_50, "avg_trials")
plot(ts(res_100_50[, 1L])); abline(h = true_theta[1L], col = "red")

coda::effectiveSize(coda::as.mcmc(res_10_50))
coda::effectiveSize(coda::as.mcmc(res_30_50))
coda::effectiveSize(coda::as.mcmc(res_50_50))
coda::effectiveSize(coda::as.mcmc(res_70_50))
coda::effectiveSize(coda::as.mcmc(res_100_50))

e <- c(
  coda::effectiveSize(coda::as.mcmc(res_10_50)) / attr(res_10_50, "time")[3L],
  coda::effectiveSize(coda::as.mcmc(res_30_50)) / attr(res_30_50, "time")[3L],
  coda::effectiveSize(coda::as.mcmc(res_50_50)) / attr(res_50_50, "time")[3L],
  coda::effectiveSize(coda::as.mcmc(res_70_50)) / attr(res_70_50, "time")[3L],
  coda::effectiveSize(coda::as.mcmc(res_100_50)) / attr(res_100_50, "time")[3L]
)
s <- c(10L, 30L, 50L, 70L, 100L)
tab_50 <- cbind(s, e)
plot(s, e, type = "o", pch = 19)

par(mfrow = c(1, 2))
plot(tab_5[, 1L], tab_5[, 2L], type = "o", pch = 19)
plot(tab_50[, 1L], tab_5[, 2L], type = "o", pch = 19)


par(mfrow = c(2, 3))
plot(ts(res_10_50[, 1L]), main = "s = 10"); abline(h = true_theta[1L], col = "red")
plot(ts(res_30_50[, 1L]), main = "s = 30"); abline(h = true_theta[1L], col = "red")
plot(ts(res_50_50[, 1L]), main = "s = 50"); abline(h = true_theta[1L], col = "red")
plot(ts(res_70_50[, 1L]), main = "s = 70"); abline(h = true_theta[1L], col = "red")
plot(ts(res_100_50[, 1L]), main = "s = 100"); abline(h = true_theta[1L], col = "red")
# plot(ts(res_bf_5[, 1L]), main = "BF"); abline(h = true_theta[1L], col = "red")
graphics.off()
