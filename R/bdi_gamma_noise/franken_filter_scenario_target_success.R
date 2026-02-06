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

par(mfrow = c(1, 3))
true_alpha <- c(5.0, 10, 50)
y <- matrix(0.0, nrow = tmax, ncol = length(true_alpha))
for (k in seq_along(true_alpha)) {
  y[, k] <- stats::rgamma(n = tmax, shape = true_alpha[k], rate = true_alpha[k] / xt_retain)
  plot(0:29, y[, k], col = "red", ylim = c(0, max(max(y[, k]), max(xt))), pch = 19,
       ylab = "", xlab = "t", main = paste0("alpha = ", true_alpha[k]))
  lines(t_continuous, xt)
  lines(0:29, xt_retain, col = "blue", pch = 18)
}


# Run particle MCMC for different target ---------------------------------------
TARGET_SUCCESS <- seq.int(10, 100, by = 10)
MAX_TRIALS <- 1000L
NDPOST <- 50000L
NSKIP <- 50000L

list_sim <- vector("list", length = length(true_alpha))
for (k in seq_along(true_alpha)) {
  cat(k, "\n")
  y_cur <- y[, k]
  a_cur <- true_alpha[k]
  list_sim[[k]] <- mclapply(seq_along(TARGET_SUCCESS), function(i) {
    res_ff <- particleMCMC_bdi(y = y_cur, x0 = x0, dt = 1/10,
                               filter = "franken",
                               target_success = TARGET_SUCCESS[i],
                               max_trials = MAX_TRIALS,
                               shape = a_cur, ndpost = NDPOST, nskip = NSKIP,
                               theta_fix = true_theta[2L],
                               wh_fix = 2L,
                               sd_prior = c(1.0, 1.0),
                               S_prop = diag(1.0, 2),
                               printevery = 100L, theta_init = true_theta[-2L])
    return(res_ff)
  }, mc.cores = 10L)
}
saveRDS(object = list_sim, file = file.path("./R/bdi_gamma_noise/results/draws_scenarios_target_success.rds"))
list_sim <- readRDS(file = file.path("./R/bdi_gamma_noise/results/draws_scenarios_target_success.rds"))

m1 <- length(TARGET_SUCCESS)
m2 <- length(true_alpha)
efficiency <- matrix(nrow = m1 * m2, ncol = 4)
for (k in seq_along(true_alpha)) {
  for (j in seq_along(TARGET_SUCCESS)) {
    idx <- m1 * (k - 1) + j
    cat(idx, "\n")
    draws <- list_sim[[k]][[j]]
    ess <- min(coda::effectiveSize(coda::as.mcmc(draws)))
    efficiency[idx, 1L] <- TARGET_SUCCESS[j]
    efficiency[idx, 2L] <- true_alpha[k]
    efficiency[idx, 3L] <- ess
    efficiency[idx, 4L] <- ess / attr(draws, "time")[3]
  }
}

par(mfrow = c(1, 3))
plot(efficiency[1:m1, 1], efficiency[1:m1, 4], type = "o")
plot(efficiency[(m1+1):(2*m1), 1], efficiency[(m1+1):(2*m1), 4], type = "o")
plot(efficiency[(2*m1+1):(3*m1), 1], efficiency[(2*m1+1):(3*m1), 4], type = "o")

par(mfrow = c(1, 3))
plot(efficiency[1:m1, 1], efficiency[1:m1, 3], type = "o")
plot(efficiency[(m1+1):(2*m1), 1], efficiency[(m1+1):(2*m1), 3], type = "o")
plot(efficiency[(2*m1+1):(3*m1), 1], efficiency[(2*m1+1):(3*m1), 3], type = "o")




graphics.off()
TARGET_SUCCESS[10]
plot(ts(list_sim[[1]][[11]]))

