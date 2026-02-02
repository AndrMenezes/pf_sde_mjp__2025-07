rm(list = ls())
source("./R/death/particleMCMC.R")
source("./R/utils.R")
Rcpp::sourceCpp(file = "./R/death/death.cpp")
library(parallel)
path_fig <- "./R/death/figures"
# Example: pure Death model
# X --> 0 (death: \lambda)
# theta = (mu[birth], lambda[death], gamma[immigration])

set.seed(12)

# Settings ---------------------------------------------------------------------
tmax <- 50L
dt <- 1 / 50 # Use 10 observations for each time period
x0 <- 100L
true_theta <- 0.01
hazard_bdi <- function(x, theta) x * theta[1L]
S_d <- matrix(-1, nrow = 1)

# Simulating -------------------------------------------------------------------
sim_xt <- sim_poisson_leap(tmax = tmax, x0 = x0, theta = true_theta,
                           hazard = hazard_bdi, S = S_d, dt = dt)
plot(sim_xt$x, type = "l")
t_continuous <- sim_xt$tt
retain <- c(1, seq.int(11, length(t_continuous), length.out = tmax - 1))
xt <- sim_xt$x[retain]
plot(xt)


# Check filtering distribution -------------------------------------------------

n_particles <- 100L
target_sucess <- 51L
max_trials <- 1e4L

mod <- new(Death, xt, x0, 1/10, 100, target_sucess, max_trials)
lml_boots <- mod$RunBootstrapFilter(true_theta)

# Check the alive PF
lml_alive <- mod$RunAliveFilter(true_theta)

# Check the Franken Filter
lml_franken <- mod$RunFrankenFilter(true_theta)

# Exact likelihood p(x_t | x_{t-1})
lml_exact <- sum(stats::dbinom(x = xt[-1], size = xt[-tmax],
                               prob = exp(-true_theta),
                               log = TRUE))

cbind(exact = lml_exact, boot = lml_boots, alive = lml_alive)

# Check log-likelihood estimation ----------------------------------------------

mc <- 2000L

# BF
N_particles <- c(10L, 50L, 100L, 400L, 500L)
list_lmls <- vector(mode = "list", length = length(N_particles))
var_lml <- numeric(length = length(N_particles))
for (j in seq_along(N_particles)) {
  cat(j, "\n")
  mod <- new(Death, xt, x0, 1/10, N_particles[j], target_sucess, max_trials)
  lmls <- mclapply(seq_len(mc), function(i) {
    return(mod$RunBootstrapFilter(true_theta))
  }, mc.cores = 20L)
  lmls <- unlist(lmls)
  list_lmls[[j]] <- lmls
  var_lml[j] <- var(lmls)
}
tab_bf <- cbind(N_particles, var_lml, exp_lml = sapply(list_lmls, mean),
                exact = lml_exact)

# Alive filtering
TARGET_SUCCESS <- c(21L, 51L, 101L, 201L)
MAX_TRIALS <- c(1e3L, 1e4L)
grid <- expand.grid(s = TARGET_SUCCESS, mp = MAX_TRIALS)
nr <- nrow(grid)
list_lmls_alive <- vector(mode = "list", length = nr)
var_lml_alive <- numeric(length = nr)
for (j in seq_len(nr)) {
  cat(j, "\n")
  lmls <- mclapply(seq_len(mc), function(i) {
    mod <- new(Death, xt, x0, 1/10, 1, grid$s[j], grid$mp[j])
    return(mod$RunAliveFilter(true_theta))
  }, mc.cores = 20L)
  lmls <- unlist(lmls)
  list_lmls_alive[[j]] <- lmls
  var_lml_alive[j] <- var(lmls)
}
tab_af <- cbind(grid, var_lml_alive, exp_lml = sapply(list_lmls_alive, mean),
                exact = lml_exact)

# Franken filtering
list_lmls_franken <- vector(mode = "list", length = nr)
var_lml_franken <- numeric(length = nr)
for (j in seq_len(nr)) {
  cat(j, "\n")
  lmls <- mclapply(seq_len(mc), function(i) {
    mod <- new(Death, xt, x0, 1/10, 1, grid$s[j], grid$mp[j])
    return(mod$RunFrankenFilter(true_theta))
  }, mc.cores = 20L)
  lmls <- unlist(lmls)
  list_lmls_franken[[j]] <- lmls
  var_lml_franken[j] <- var(lmls)
}

tab_ff <- cbind(grid, var_lml_franken, exp_lml = sapply(list_lmls_franken, mean),
                exact = lml_exact)
tab_bf
tab_af
tab_ff


# particleMCMC with Bootstrap filter -------------------------------------------

# pilot_res <- particleMCMC_death(x = xt, x0 = x0, n_particles = 100L, dt = 1/10,
#                                 filter = "bootstrap", ndpost = 5000L, nskip = 0L,
#                                 sd_prop = 1, printevery = 100L)
# plot(pilot_res[-(1:1000)])
# plot(density(pilot_res[-(1:1000)]))
# sd_opt <- 2.5 * sd(log(pilot_res[-(1:2000L)]))

res_bf <- particleMCMC_death(x = xt, x0 = x0, n_particles = 100L, dt = 1/40,
                             ndpost = 50000L, nskip = 0L, filter = "bootstrap",
                             sd_prop = 0.4, printevery = 100L)
attr(res_bf, "AcceptRate")
coda::effectiveSize(coda::as.mcmc(res_bf))
coda::effectiveSize(coda::as.mcmc(res_bf[-(1:20000L)]))

rmv <- 1:20000L
mean(res_bf[-rmv]/true_theta)
plot(log(res_bf), type = "l")
plot(res_bf, type = "l"); abline(h = true_theta, col = "red")
plot(density(res_bf[-rmv]))

# Save plot
png(filename = file.path(path_fig, "trace_density_bf__sim1.png"), units = "cm",
    width = 20, height = 10, res = 1000)
par(mfrow = c(1, 2))
plot(res_bf[-rmv], type = "l", ylab = expression(theta[t]), xlab = "t")
abline(h = true_theta[1L], col = "red")
plot(density(res_bf[-rmv]), main = "", xlab = expression(theta))
points(x = true_theta[1L], y = 0)
graphics.off()


# particle MCMC with alive filter
res_alive <- particleMCMC_death(x = xt, x0 = x0, n_particles = 1L,
                                max_trials = 1000L,
                                target_sucess = 50L, dt = 1/10,
                                ndpost = 50000L, nskip = 0L, filter = "alive",
                                sd_prop = 2, printevery = 100L)
attr(res_alive, "AcceptRate")
coda::effectiveSize(coda::as.mcmc(res_alive))
coda::effectiveSize(coda::as.mcmc(res_alive[-(1:20000L)]))
plot(log(res_alive), type = "l")
plot(res_alive, type = "l");  abline(h = true_theta, col = "red")
mean(res_alive[-(1:1000L)]/true_theta)
plot(density(res_alive[-(1:1000L)]))

png(filename = file.path(path_fig, "trace_density_alive__sim1.png"), units = "cm",
    width = 20, height = 10, res = 1000)
par(mfrow = c(1, 2))
plot(res_alive[-rmv], type = "l", ylab = expression(theta[t]), xlab = "t")
abline(h = true_theta[1L], col = "red")
plot(density(res_alive[-rmv]), main = "", xlab = expression(theta))
points(x = true_theta[1L], y = 0)
graphics.off()


# particle MCMC with franken filter
res_franken <- particleMCMC_death(x = xt, x0 = x0, n_particles = 1L,
                                  max_trials = 400L,
                                  target_sucess = 50L, dt = 1/10,
                                  ndpost = 50000L, nskip = 0L, filter = "franken",
                                  sd_prop = 0.4, printevery = 100L)
attr(res_franken, "AcceptRate")
coda::effectiveSize(coda::as.mcmc(res_franken))
coda::effectiveSize(coda::as.mcmc(res_franken[-(1:20000L)]))
plot(log(res_franken), type = "l")
plot(res_franken, type = "l");  abline(h = true_theta, col = "red")
mean(res_franken[-(1:1000L)]/true_theta)
plot(density(res_franken[-(1:1000L)]))

# Use exact likelhood
res_exact <- exactMCMC_death(y = xt, x0 = x0, niter = 50000L)
attr(res_exact, "AcceptRate")
coda::effectiveSize(coda::as.mcmc(res_exact))
coda::effectiveSize(coda::as.mcmc(res_exact[-(1:20000L)]))
plot(log(res_exact), type = "l")
plot(res_exact, type = "l");  abline(h = true_theta, col = "red")
plot(density(res_exact[-(1:1000L)]))
mean(res_exact[-(1:1000L)]/true_theta)


cbind(exact = mean(res_exact[-rmv]/true_theta),
      bf = mean(res_bf[-rmv]/true_theta),
      af = mean(res_alive[-rmv]/true_theta),
      ff = mean(res_franken[-rmv]/true_theta))
