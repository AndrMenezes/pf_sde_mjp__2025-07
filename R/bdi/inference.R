rm(list = ls()); gc()
library(parallel)

path_local <- "./R/bdi/"
path_fig <- file.path(path_local, "figures")
source("./R/utils.R")
source(file.path(path_local, "particleMCMC.R"))
Rcpp::sourceCpp(file = file.path(path_local, "bdi.cpp"))


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
yt <- stats::rpois(n = tmax, lambda = xt_retain)
# yt <- stats::rnorm(n = tmax, mean = xt_retain, sd = 1.0)
# x11()
# png(filename = file.path(path_fig, "data_sim1.png"), units = "cm",
#     width = 16, height = 14, res = 1000)
# par(mar = c(2.2, 2.2, 1.0, 0.5), cex = 1.1)
plot(0:29, yt, col = "red", ylim = c(0, max(max(yt), max(xt))), pch = 19, ylab = "", xlab = "t")
lines(t_continuous, xt)
lines(0:29, xt_retain, type = "o", col = "blue", pch = 18)
# graphics.off()


# Get an estimate of the variance of the log-likelihood  -----------------------
mc <- 2000L
N_particles <- c(10L, 20L, 30L, 35, 40L, 50L, 100L, 200L, 300L, 400L, 500L)
list_lmls <- vector(mode = "list", length = length(N_particles))
var_lml <- numeric(length = length(N_particles))
for (j in seq_along(N_particles)) {
  cat(j, "\n")
  mod <- new(BDI, yt, x0, N_particles[j], 1/10)
  lmls <- mclapply(seq_len(mc), function(i) {
    mod$RunBootstrapFilter(true_theta)
    return(mod$lml)
  }, mc.cores = 20L)
  lmls <- unlist(lmls)
  list_lmls[[j]] <- lmls
  var_lml[j] <- var(lmls)
}
cbind(N_particles, var_lml)
plot(N_particles, var_lml)

# Set the number of particles according to Var(log(lml)) \approx 1
n_particles_opt <- 400L

# Estimate the log-likelihood varying mu and fixing the other parameter --------
mod <- new(BDI, yt, x0, n_particles_opt, dt)
n_grid <- 500L
mu_grid <- seq(0.1, 1.0, length.out = n_grid)

parameters <- true_theta
lmls <- sapply(seq_len(n_grid), function(i) {
  cat(i, "\n")
  parameters[1L] <- mu_grid[i]
  mod$RunBootstrapFilter(parameters)
  return(mod$lml)
})

rmv <- which(mu_grid > 0.8)
plot(mu_grid[-rmv], lmls[-rmv])
abline(v = true_theta[1])

plot(log(mu_grid[-rmv]), lmls[-rmv])
abline(v = true_theta[1])

plot(log(mu_grid), lmls)
abline(v = true_theta[1])


# Estimate log-likelihood for varying \mu and \lambda --------------------------
n_particles <- 100L
mod <- new(BDI, yt, x0, n_particles, dt)
n_grid <- 50L
parms_grid <- expand.grid(mu = seq(0.1, 1.0, length.out = n_grid),
                          lambda = seq(0.1, 1.5, length.out = n_grid))
# \mu < \lambda
parms_grid <- parms_grid[which(parms_grid$mu < parms_grid$lambda), ]
nrow(parms_grid)

parameters <- true_theta
lmls <- sapply(seq_len(nrow(parms_grid)), function(i) {
  parameters[1L] <- parms_grid[i, ]$mu
  parameters[2L] <- parms_grid[i, ]$lambda
  mod$RunBootstrapFilter(parameters)
  cat(i, mod$lml, "\n")
  return(mod$lml)
})

parms_grid$lml <- lmls

library(ggplot2)
ggplot(parms_grid) +
  geom_raster(aes(mu, lambda, fill = lml), alpha = 0.8) +
  # geom_contour(aes(mu, lambda, z = lml), alpha = 0.8) +
  scale_fill_viridis_c(option = "C", limits = range(parms_grid$lml))

ggplot(parms_grid, aes(x = mu, y = lambda, z = lml)) +
  geom_raster(aes(fill = lml)) +
  geom_contour(col = "black") +
  geom_point(data = data.frame(mu = true_theta[1], lambda = true_theta[2],
                               lml = NA_real_),
             aes(x = mu, y = lambda), col = "black", fill = "firebrick2",
             shape = 21, size = 3.0) +
  scale_fill_gradient2(low ="dodgerblue2", high ="firebrick2", mid = "white",
                       midpoint = quantile(parms_grid$lml, 0.25),
                       limits = range(parms_grid$lml)) +
  labs(x = expression(mu), y = expression(lambda)) +
  theme_bw() +
  theme(text             = element_text(family = "Palatino", colour = "black",
                                        size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# particleMCMC to estimate only \mu --------------------------------------------
Rcpp::sourceCpp(file = file.path(path_local, "bdi.cpp"))
mod_1 <- new(BDI, yt, x0, 500L, 1/10)

system.time(
  mod_1$RunParticleMCMC_mu(true_theta[2L], true_theta[3L], 1.0, 0.25, 5000L))
plot(mod_1$draws_mu[-(1:1000L)], type = "l"); abline(h = true_theta[1], col = "red")
sd_opt <- sd(log(mod_1$draws_mu[-(1:1000L)]))

system.time(
  mod_1$RunParticleMCMC_mu(true_theta[2L], true_theta[3L], 10, sd_opt*2.5, 10000L))
plot(mod_1$draws_mu[-(1:2000L)], type = "l"); abline(h = true_theta[1], col = "red")

rmv <- 1:2000L
png(filename = file.path(path_fig, "trace_density_mu_bdi__sim2.png"), units = "cm",
    width = 20, height = 10, res = 1000)
par(mfrow = c(1, 2))
plot(mod_1$draws_mu[-rmv], type = "l", ylab = expression(mu[t]), xlab = "t")
abline(h = true_theta[1L], col = "red")
plot(density(mod_1$draws_mu[-rmv]), main = "", xlab = expression(mu))
points(x = true_theta[1L], y = 0)
graphics.off()


# plot(mod_1$draws_mu, type = "l"); abline(h = true_theta[1], col = "red")

# Pilot run to estimate the posterior variance
pilot_res <- particleMCMC_bdi(y = yt, x0 = x0, n_particles = 500L,
                              dt = 1/10,
                              ndpost = 5000L, nskip = 0L,
                              theta_fix = true_theta[-1L], wh_fix = c(2L, 3L),
                              sd_prior = 1.0, S_prop = 0.25, printevery = 10L)
attr(pilot_res, "AcceptRate")
attr(pilot_res, "time")

# Estimate posterior variance
rmv <- 1:1000L
var(log(pilot_res[-rmv, 1L]))
(sd_opt <- sd(log(pilot_res[-rmv, 1L])))
hist(pilot_res[-rmv, 1L]); abline(v=true_theta[1L], col = "red")
plot(pilot_res, type = "l"); abline(h=true_theta[1L], col = "red")

# Final run with the optimal proposal variance
res <- particleMCMC_bdi(y = yt, x0 = x0, n_particles = 500L, dt = 1/10,
                        ndpost = 50000L, nskip = 0L,
                        theta_fix = true_theta[-1L], wh_fix = c(2L, 3L),
                        sd_prior = 0.10, S_prop = sd_opt*2.5, printevery = 100L)
attr(res, "AcceptRate")
attr(res, "time")[3L] / 60

rmv <- 1:10000L
coda::effectiveSize(coda::as.mcmc(res[-rmv, 1L]))

png(filename = file.path(path_fig, "trace_density_mu_bdi__sim1.png"), units = "cm",
    width = 20, height = 10, res = 1000)
par(mfrow = c(1, 2))
plot(res[-rmv, 1L], type = "l", ylab = expression(mu[t]), xlab = "t")
abline(h = true_theta[1L], col = "red")
plot(density(res[-rmv, 1L]), main = "", xlab = expression(mu))
points(x = true_theta[1L], y = 0)
graphics.off()


# particleMCMC to estimate \mu and \gamma -------------------------------------

# Pilot run to estimate the posterior variance of the log(\mu) and log(\gamma)
pilot_res <- particleMCMC_bdi(y = yt, x0 = x0, n_particles = 50L, dt = 1/10,
                              ndpost = 50000L, nskip = 0L,
                              theta_fix = true_theta[2L], wh_fix = c(2L),
                              sd_prior = rep(1.0, 2L),
                              S_prop = rep(1.0, 2L),
                              printevery = 10L)
attr(pilot_res, "AcceptRate")
attr(pilot_res, "time")

# Estimate posterior variance
rmv <- 1:30000L
var(log(pilot_res[-rmv, 2L]))
cov_opt <- cov(log(pilot_res[-rmv, ]))*2.5^2 / 2

hist(pilot_res[-rmv, 1L])
hist(pilot_res[-rmv, 2L])
matplot(pilot_res)
plot(pilot_res[-rmv, 1L], type = "l");abline(h = true_theta[1L], col = "red")
plot(pilot_res[-rmv, 2L], type = "l");abline(h = true_theta[3L], col = "red")


# Final run with the optimal proposal variance
# cov_opt <- matrix(c(0.005450674, -0.00263535,-0.00263535,0.53284964), nrow=2L)
res <- particleMCMC_bdi(y = yt, x0 = x0, n_particles = 50L, dt = 1/10,
                        ndpost = 50000L, nskip = 0L,
                        theta_fix = true_theta[2L], wh_fix = c(2L),
                        sd_prior = rep(1.0, 2L), S_prop = cov_opt,
                        printevery = 10L)
attr(res, "AcceptRate")
attr(res, "time")[3L] / 60

plot(res[, 1L], type = "l", ylab = expression(mu[t]), xlab = "t")
abline(h = true_theta[1L], col = "red")
plot(res[, 2L], type = "l", ylab = expression(gamma[t]), xlab = "t")
abline(h = true_theta[3L], col = "red")

rmv <- 1:20000L
plot(res[-rmv, 1L], type = "l", ylab = expression(mu[t]), xlab = "t")
abline(h = true_theta[1L], col = "red")
plot(res[-rmv, 2L], type = "l", ylab = expression(gamma[t]), xlab = "t")
abline(h = true_theta[3L], col = "red")


png(filename = file.path(path_fig, "trace_density_mu_both__sim2.png"), units = "cm",
    width = 20, height = 10, res = 1000)
par(mfrow = c(1, 2))
plot(res[-rmv, 1L], type = "l", ylab = expression(mu[t]), xlab = "t")
abline(h = true_theta[1L], col = "red")
plot(density(res[-rmv, 1L]), main = "", xlab = expression(mu))
points(x = true_theta[1L], y = 0)
graphics.off()

png(filename = file.path(path_fig, "trace_density_gamma_bdi_both__sim2.png"), units = "cm",
    width = 20, height = 10, res = 1000)
par(mfrow = c(1, 2))
plot(res[-rmv, 2L], type = "l", ylab = expression(gamma[t]), xlab = "t")
abline(h = true_theta[3L], col = "red")
plot(density(res[-rmv, 2L]), main = "", xlab = expression(gamma))
points(x = true_theta[3L], y = 0)
graphics.off()


png(filename = "./figures/trace_mu_gamma_bdi.png", units = "cm",
    width = 16, height = 14, res = 1000)

par(mfrow = c(2, 1))
plot(res[-rmv, 1L], type = "l", ylab = expression(mu[t]), xlab = "t")
abline(h = true_theta[1L], col = "red")
# plot(density(res[-rmv, 1L]), main = "", xlab = expression(mu))
# points(x = true_theta[1L], y = 0)

plot(res[-rmv, 2L], type = "l", ylab = expression(gamma[t]), xlab = "t")
abline(h = true_theta[3L], col = "red")
# plot(density(res[-rmv, 2L]), main = "", xlab = expression(lambda))
# points(x = true_theta[2L], y = 0)

graphics.off()


#

# particleMCMC to estimate \mu and \lambda -------------------------------------

# Pilot run to estimate the posterior variance
pilot_res <- particleMCMC_bdi(y = yt, x0 = x0, n_particles = 200L, dt = 1/10,
                              ndpost = 10000L, nskip = 0L,
                              theta_fix = true_theta[3], wh_fix = c(3L),
                              sd_prior = rep(0.10, 2L), S_prop = rep(0.10, 2L),
                              printevery = 10L)
attr(pilot_res, "AcceptRate")
attr(pilot_res, "time")

# Estimate posterior variance
rmv <- 1:1000L
var(pilot_res[-rmv, 1L])
(cov_opt <- cov(pilot_res[-rmv, ]))

hist(pilot_res[-rmv, 1L])
hist(pilot_res[-rmv, 2L])
matplot(pilot_res)
plot(pilot_res[, 1L], type = "l")
lines(pilot_res[, 2L], type = "l", col = "red")


# Final run with the optimal proposal variance
res <- particleMCMC_bdi(y = yt, x0 = x0, n_particles = 200L, dt = 1/10,
                        ndpost = 50000L, nskip = 0L,
                        theta_fix = true_theta[3], wh_fix = c(3L),
                        sd_prior = 0.10, S_prop = cov_opt*2.5^2 / 2,
                        printevery = 100L)
attr(res, "AcceptRate")
attr(res, "time")[3L] / 60

rmv <- 1:10000L
png(filename = "./figures/trace_mu_lambda_bdi.png", units = "cm",
    width = 16, height = 14, res = 1000)

par(mfrow = c(2, 1))
plot(res[-rmv, 1L], type = "l", ylab = expression(mu[t]), xlab = "t")
abline(h = true_theta[1L], col = "red")
# plot(density(res[-rmv, 1L]), main = "", xlab = expression(mu))
# points(x = true_theta[1L], y = 0)

plot(res[-rmv, 2L], type = "l", ylab = expression(lambda[t]), xlab = "t")
abline(h = true_theta[2L], col = "red")
# plot(density(res[-rmv, 2L]), main = "", xlab = expression(lambda))
# points(x = true_theta[2L], y = 0)

graphics.off()

