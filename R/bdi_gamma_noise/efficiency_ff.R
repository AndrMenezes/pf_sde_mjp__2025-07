rm(list = ls()); gc()
library(parallel)
library(ggplot2)
library(cowplot)
source("./R/utils.R")

# Paths
path_local <- "./R/bdi_gamma_noise/"
time_id <- format(Sys.time(), "%Y-%b-%d-%X") # "2026-Feb-13-09:33:44"#
path_res_old <- file.path(path_local, "results", "2026-Apr-17-11:04:29")
path_res <- file.path(path_local, "results", time_id)
if (!dir.exists(path_res)) dir.create(path = path_res)

# Source model and particleMCMC
Rcpp::sourceCpp(file = file.path(path_local, "bdi.cpp"))
source(file.path(path_local, "particleMCMC.R"))

# Plot theme
theme_set(theme_cowplot(font_size = 16) + background_grid())

# Seed for reproducibility
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

TRUE_ALPHA <- c(10.0, 20.0, 50.0)
LEN_ALPHA <- length(TRUE_ALPHA)
png(filename = file.path(path_res, "data.png"), width = 1200, height = 800)
par(mfrow = c(1, 3))
y <- matrix(0.0, nrow = tmax, ncol = length(TRUE_ALPHA))
for (k in seq_len(LEN_ALPHA)) {
  y[, k] <- stats::rgamma(n = tmax, shape = TRUE_ALPHA[k],
                          rate = TRUE_ALPHA[k] / xt_retain)
  plot(0:29, y[, k], col = "red", ylim = c(0, max(max(y[, k]), max(xt))), pch = 19,
       ylab = "", xlab = "t", main = paste0("alpha = ", TRUE_ALPHA[k]))
  lines(t_continuous, xt)
  lines(0:29, xt_retain, col = "blue", pch = 18)
}
graphics.off()



# Get an estimate of covariance matrix for the proposal -------------------
# list_pilot <- mclapply(seq_len(LEN_ALPHA), function(k) {
#   y_cur <- y[, k]
#   a_cur <- TRUE_ALPHA[k]
#   particleMCMC_bdi(y = y_cur, x0 = x0, dt = 1/10,
#                    filter = "franken",
#                    target_success = 50L,
#                    max_trials = 1000L,
#                    shape = a_cur, ndpost = 500000L, nskip = 20000L,
#                    theta_fix = true_theta[2L],
#                    wh_fix = 2L,
#                    sd_prior = c(1.0, 1.0),
#                    S_prop = diag(1.0, 2),
#                    printevery = 100L, theta_init = true_theta[-2L])
# }, mc.cores = LEN_ALPHA)
# saveRDS(object = list_pilot, file = file.path(path_res, "draws_pilot.rds"))
# list_pilot <- readRDS(file = file.path(path_res, "draws_pilot.rds"))
# S_prop <- vector(mode = "list", length = LEN_ALPHA)
# scale_opt <- 2.5^2/2
# for (k in seq_along(TRUE_ALPHA)) {
#   S_prop[[k]] <- scale_opt * cov(log(list_pilot[[k]]))
# }

# Import results from old the experiment to estimate common covariance
tmp <- readRDS(file = file.path(path_res_old, "draws.rds"))
length(tmp[[1]])
S_prop <- vector(mode = "list", length = LEN_ALPHA)
for (k in seq_along(TRUE_ALPHA)) {
  S_prop[[k]] <- 2.5^2/2*cov(log(do.call(rbind, tmp[[k]])))
}
rm(tmp);gc()

# Run particle MCMC for different target ---------------------------------------
min_s <- 10L
max_s <- 150
step <- 10L
TARGET_SUCCESS <- c(2L, 5L, 8L, seq.int(min_s, max_s, by = step))
LEN_S <- length(TARGET_SUCCESS)
MAX_TRIALS <- 1000L
NDPOST <- 80000L
NSKIP <- 0L

# remove the first iterations
to_rmv <- seq_len(30000L)
list_sim <- vector("list", length = length(TRUE_ALPHA))
for (k in seq_along(TRUE_ALPHA)) {
  y_cur <- y[, k]
  a_cur <- TRUE_ALPHA[k]
  S_prop_cur <- S_prop[[k]]
  list_sim[[k]] <- lapply(seq_len(LEN_S), function(i) {
    cat(k, "of", LEN_ALPHA, i, "of", LEN_S, "\n")
    res_ff <- particleMCMC_bdi(y = y_cur, x0 = x0, dt = 1/10,
                               filter = "franken",
                               target_success = TARGET_SUCCESS[i],
                               max_trials = MAX_TRIALS,
                               shape = a_cur, ndpost = NDPOST, nskip = NSKIP,
                               theta_fix = true_theta[2L],
                               wh_fix = 2L,
                               sd_prior = c(1.0, 1.0),
                               S_prop = S_prop_cur,
                               printevery = 100L, theta_init = true_theta[-2L])
    # res_ff[-to_rmv, ]
    return(res_ff)
  }) #, mc.cores = length(TARGET_SUCCESS)
}
saveRDS(object = list_sim, file = file.path(path_res, "draws.rds"))

list_sim <- readRDS(file = file.path(path_res, "draws.rds"))

# x= list_sim[[1]][[10]]
# plot(x[, 1L], x[, 2L])

# Append the results ------------------------------------------------------

m1 <- length(TARGET_SUCCESS)
m2 <- length(TRUE_ALPHA)
data_results <- matrix(nrow = m1 * m2, ncol = 9L)
for (k in seq_len(m2)) {
  # cat(k, "\n")
  for (j in seq_len(m1)) {
    idx <- m1 * (k - 1) + j
    cat(idx, "\n")
    draws <- list_sim[[k]][[j]]
    data_results[idx, 1L] <- TARGET_SUCCESS[j]
    data_results[idx, 2L] <- TRUE_ALPHA[k]
    if (is.matrix(draws)) {
      ess <- min(coda::effectiveSize(coda::as.mcmc(draws)))
      data_results[idx, 3L] <- ess
      data_results[idx, 4L] <- ess / attr(draws, "time")[3]
      data_results[idx, 5L] <- ess / TARGET_SUCCESS[j]
      data_results[idx, 6L] <- attr(draws, "avg_trials")
      data_results[idx, 7L] <- attr(draws, "acceptance_rate")
      data_results[idx, 8L] <- mean(attr(draws, "avg_pct_m_max_reached"))
      data_results[idx, 9L] <- attr(draws, "time")[3]
    }
  }
}

colnames(data_results) <- c("target_success", "alpha", "mESS", "mESS_s",
                            "mESS_target_s", "avg_trials", "accept_rate",
                            "pct_m_max_reached", "time")
data_results <- as.data.frame(data_results)
data_results$alpha <- paste0("alpha == ", data_results$alpha)
data_results$alpha <- forcats::fct_relevel(data_results$alpha, "alpha == 10",
                                           "alpha == 20")
#
p_ess <- ggplot(data = data_results, aes(x = target_success, y = mESS)) +
  facet_wrap(~alpha, labeller = label_parsed) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = tmax) +
  labs(x = "Target success (s)", y = "mESS") +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
save_plot(filename = file.path(path_res, "mESS.png"), plot = p_ess,
          base_height = 8, bg = "white")

p_ef <- ggplot(data = data_results, aes(x = target_success, y = mESS_s)) +
  facet_wrap(~alpha, labeller = label_parsed) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = tmax) +
  labs(x = "Target success (s)", y = "mESS/s") +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
save_plot(filename = file.path(path_res, "mESS_sec.png"), plot = p_ef,
          base_height = 8, bg = "white")

p_ef2 <- ggplot(data = data_results, aes(x = target_success, y = mESS_target_s)) +
  facet_wrap(~alpha, labeller = label_parsed) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = tmax) +
  labs(x = "Target success (s)", y = "mESS/target_success") +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
save_plot(filename = file.path(path_res, "mESS_targetS.png"), plot = p_ef2,
          base_height = 8, bg = "white")

p_m_max <- ggplot(data = data_results,
                  aes(x = target_success, y = pct_m_max_reached)) +
  facet_wrap(~alpha, labeller = label_parsed) +
  geom_point() +
  geom_line() +
  # geom_vline(xintercept = tmax) +
  labs(x = "Target success (s)", y = "mean(pct_m_max_reached)") +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
save_plot(filename = file.path(path_res, "pct_m_max_reached.png"),
          plot = p_m_max, base_height = 8, bg = "white")

p_ar <- ggplot(data = data_results, aes(x = target_success, y = accept_rate)) +
  facet_wrap(~alpha, labeller = label_parsed) +
  geom_point() +
  geom_line() +
  # geom_vline(xintercept = tmax) +
  labs(x = "Target success (s)", y = "Acceptance rate") +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
save_plot(filename = file.path(path_res, "acceptance_rate_targetS.png"),
          plot = p_ar, base_height = 8, bg = "white")

p_time <- ggplot(data = data_results, aes(x = target_success, y =time)) +
  facet_wrap(~alpha, labeller = label_parsed) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = tmax) +
  labs(x = "Target success (s)", y = "Time (sec)") +
  scale_x_continuous(breaks = scales::pretty_breaks(6))

# Trace plots for different \alpha's --------------------------------------

# k=1
dat_true <- lapply(seq_len(m1), function(i) {
  data.frame(parms = c("mu", "gamma"), value = true_theta[-2L],
             target_success = TARGET_SUCCESS[i])
})
dat_true <- do.call(rbind, dat_true)
for (k in seq_along(TRUE_ALPHA)) {
  cat(k, "\n")
  # if (!is.matrix(list_sim[[k]])) next
  dat <- as.data.frame(do.call(rbind, list_sim[[k]]))
  colnames(dat) <- c("mu", "gamma")
  dat$iter <- rep(seq_len(NDPOST), times = m1)
  dat$target_success <- rep(TARGET_SUCCESS, each = NDPOST)
  dat <- tidyr::pivot_longer(data = dat, cols = -c(iter, target_success),
                             names_to = "parms")
  p <- ggplot(data = dat, aes(x = iter, y = value)) +
    facet_wrap(parms~target_success, scales = "free_y", label = label_parsed) +
    geom_line() +
    geom_hline(data = dat_true, aes(yintercept = value), col = "red") +
    labs(x = "Iteration", y = "Draw")
  ff <- paste0("trace_", TRUE_ALPHA[k], ".png")
  save_plot(filename = file.path(path_res, ff), plot = p,
            base_height = 10, bg = "white")
}
