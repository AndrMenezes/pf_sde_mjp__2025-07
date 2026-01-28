rm(list = ls())
source("./R/utils.R")
source("./R/wquantile.R")
Rcpp::sourceCpp(file = "./R/bdi/bdi.cpp")

# Example: Birth-Death-Immigration (BDI) model
# X --> 2X (birth: \mu)
# X --> 0 (death: \lambda)
# 0 --> X (immigration: \lambda)
# theta = (mu[birth], lambda[death], gamma[immigration])

# Settings ---------------------------------------------------------------------
tmax <- 30L
dt <- 1 / 10 # Use 10 observations for each time period
x0 <- 500L
true_theta <- c(0.5, 0.8, 40)
hazard_bdi <- function(x, theta) c(x * theta[1], x * theta[2], theta[3])
S_bdi <- matrix(c(1, -1, 1), nrow = 1)
# S_bdi %*% c(4, 2 , 3)
# (1, -1, 1) [4, 2, 3] = 4 -2 + 3


# Simulating -------------------------------------------------------------------

sim_xt <- sim_poisson_leap(tmax = tmax, x0 = x0, theta = true_theta,
                           hazard = hazard_bdi, S = S_bdi, dt = dt)

# Adding a observational layer
xt <- sim_xt$x[, 1L]
t_continuous <- sim_xt$tt
n <- length(xt)
retain <- c(1, seq(11, length(t_continuous), by = 10))
t_retain <- t_continuous[retain]
xt_retain <- xt[retain]
yt <- stats::rpois(n = tmax, lambda = xt_retain)
plot(t_continuous, xt, type = "l")
lines(t_retain, xt_retain, type = "o", col = "blue")
lines(t_retain, yt, type = "o", col = "red")
graphics.off()



# Bootstrap filter -------------------------------------------------------------
Rcpp::sourceCpp(file = "./R/bdi/bdi.cpp")
n_particles <- 100L
mod <- new(BDI, yt, x0, n_particles, dt)
mod$RunBootstrapFilter(true_theta)
mod$lml

x_path <- mod$path

matplot(x_path, type = "l", lty = 1)
plot(xt, ylim = c(0, 500), cex = 3, pch = 20)
matlines(x_path, type = "l", lty = 1)

# plot(yt, ylim = c(0, 500))
# matlines(x_path, type = "l", lty = 1)
graphics.off()

particles <- matrix(mod$particles, nrow = tmax, ncol = n_particles, byrow = TRUE)
weights <- matrix(mod$weights, nrow = tmax, ncol = n_particles, byrow = TRUE)
dim(particles)
head(particles)
tail(particles)
head(weights)
rowSums(weights)
# Get the quantiles of the filtering distribution for each time
qs <- matrix(nrow = tmax, ncol = 3)
for (i in 1:tmax) {
  qs[i,] <- wquantile(particles[i, ], weights[i, ], prob = c(0.05, 0.5, 0.95))
}

yrange <- range(c(min(qs), max(qs), max(yt)))
plot(yt, type = "p", ylim = yrange)
lines(qs[, 1], col = "blue", lty = 2)
lines(qs[, 2], col = "blue")
lines(qs[, 3], col = "blue", lty = 2)

yrange <- c(min(c(0, qs)), max(qs))
plot(t_continuous, xt, type = "l", ylim = yrange)
lines(qs[, 1], col = "blue", lty = 2)
lines(qs[, 2], col = "blue")
lines(qs[, 3], col = "blue", lty = 2)

yrange <- c(min(c(0, qs)), max(qs))
plot(xt_retain, type = "o", ylim = yrange)
lines(qs[, 1], col = "blue", lty = 2)
lines(qs[, 2], col = "blue")
lines(qs[, 3], col = "blue", lty = 2)


plot(yt, pch = 19, cex = 2)
matlines(particles)

samp <- sample(1:nrow(weights), 4)
par(mfrow = c(2, 2))
for (i in seq_len(4)) hist(weights[samp[i], ])



# Plotting with ggplot2 --------------------------------------------------------

data_true_x <- data.frame(time = seq_len(tmax), x = xt_retain, y = yt)
data_filtering_x <- data.frame(time = seq_len(tmax), med = qs[, 2], lo = qs[, 1],
                               up = qs[, 3])
data_path_x <- data.frame(time = seq(1, tmax, length.out = nrow(x_path)), x_path)
data_path_x <- tidyr::pivot_longer(data_path_x, cols = -time)
head(data_path_x)

library(ggplot2)
theme_set(theme_bw(base_size = 16))
p1 <- ggplot(data_true_x, aes(x = time, y = x)) +
  geom_line(data = data_path_x, aes(x = time, y = value, group = name), alpha = 0.5,
            col = "dodgerblue", alpha = 0.5) +
  geom_point(size = 6, shape = 4)
# + geom_line()

p2 <- ggplot(data_true_x, aes(x = time)) +
  geom_line(data = data_filtering_x, aes(x = time, y = med), col = "dodgerblue") +
  # geom_line(data = data_filtering_x, aes(x = time, y = lo), col = "dodgerblue", linetype = 2) +
  # geom_line(data = data_filtering_x, aes(x = time, y = up), col = "dodgerblue", linetype = 2) +
  geom_ribbon(data = data_filtering_x, mapping = aes(x = time, ymin = lo, ymax = up),
              fill = "dodgerblue", alpha = 0.5) +
  geom_point(aes(y = x), size = 6, shape = 4)
# + geom_line(aes(y = x))

p_grid <- cowplot::plot_grid(p1 + ggtitle("Trajectories"), p2 + ggtitle("Filtering"))
p_grid
cowplot::save_plot(filename = "./figures/bdi_sim.png", plot = p_grid, base_height = 8,
                   bg = "white")
