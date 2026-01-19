rm(list = ls())
source("./R/bdi/utils.R")

# Example: Birth-Death-Immigration (BDI) model
# X --> 2X (birth: \mu)
# X --> 0 (death: \lambda)
# 0 --> X (immigration: \lambda)
# theta = (mu[birth], lambda[death], gamma[immigration])

# Settings ---------------------------------------------------------------------
tmax <- 30L
x0 <- 500L
true_theta <- c(0.5, 0.8, 40)
hazard_bdi <- function(x, theta) c(x * theta[1], x * theta[2], theta[3])
S_bdi <- matrix(c(1, -1, 1), nrow = 1)
# S_bdi %*% c(4, 2 , 3)
# (1, -1, 1) [4, 2, 3] = 4 -2 + 3

# Simulate using different methods ---------------------------------------------

out <- sim_gillespie(tmax = 30L, x0 = 500L, theta = true_theta,
                     hazard = hazard_bdi, S = S_bdi)

out_f <- sim_frm(tmax = 30, x0 = 500, theta = true_theta, hazard = hazard_bdi,
                 S = S_bdi)

# 10 obs per unit of time
out_p_01 <- sim_poisson_leap(tmax = 30, x0 = 500, theta = true_theta,
                             hazard = hazard_bdi, S = S_bdi, dt = 1/10)

# 20 obs per unit of time
out_p_005 <- sim_poisson_leap(tmax = 30, x0 = 500, theta = true_theta,
                              hazard = hazard_bdi, S = S_bdi, dt = 1/20)

# 30 obs per unit of time
out_p_003 <- sim_poisson_leap(tmax = 30, x0 = 500, theta = true_theta,
                              hazard = hazard_bdi, S = S_bdi, dt = 1/30)

# 100 obs per unit of time
out_p_001 <- sim_poisson_leap(tmax = 30, x0 = 500, theta = true_theta,
                              hazard = hazard_bdi, S = S_bdi, dt = 1/100)


# ODE solution for given parameter setting -------------------------------------
ode_solution_bdi <- function(tt, x0, theta) {
  mu <- theta[1]
  lambda <- theta[2]
  gama <- theta[3]
  gama / (lambda - mu) + (x0 - gama / (lambda - mu)) * exp((mu - lambda) * tt)
}
t_discrete <- seq.int(0, 30)
sol_bdi <- ode_solution_bdi(t_discrete, 500, true_theta)




# Plotting ---------------------------------------------------------------------

par(mfrow = c(2, 2))
plot(out$tt, out$x, type = "l", main = "Gillespie method")
lines(t_discrete, sol_bdi, col = "red")
points(t_discrete, sol_bdi, col = "red", pch = 19)

plot(out_f$tt, out_f$x, type = "l", main = "First reaction method")
lines(t_discrete, sol_bdi, col = "red")
points(t_discrete, sol_bdi, col = "red", pch = 19)

plot(out_p_01$tt, out_p_01$x, type = "l", main = "Poisson leap dt=0.1")
lines(t_discrete, sol_bdi, col = "red")
points(t_discrete, sol_bdi, col = "red", pch = 19)

plot(out_p_001$tt, out_p_001$x, type = "l", main = "Poisson leap dt=0.01")
lines(t_discrete, sol_bdi, col = "red")
points(t_discrete, sol_bdi, col = "red", pch = 19)

graphics.off()

par(mfrow = c(1, 2))
plot(out_p_005$tt, out_p_005$x, type = "l", main = "Poisson leap dt=1/20")
lines(t_discrete, sol_bdi, col = "red")
points(t_discrete, sol_bdi, col = "red", pch = 19)

plot(out_p_003$tt, out_p_003$x, type = "l", main = "Poisson leap dt=1/30")
lines(t_discrete, sol_bdi, col = "red")
points(t_discrete, sol_bdi, col = "red", pch = 19)



# Adding a observational layer -------------------------------------------------
xt <- out_p_01$x[, 1L]
t_continuous <- out_p_01$tt
n <- length(xt)

retain <- c(0, seq(11, length(t_continuous), by = 10))
t_ <- t_continuous[retain]
xt_ <- xt[retain]
yt <- stats::rpois(n = n, lambda = xt)
par(mfrow = c(1, 2))
# plot(t_continuous, yt, type = "l", main = "Observation")
# lines(t_discrete, sol_bdi, col = "red")
# points(t_discrete, sol_bdi, col = "red", pch = 19)

plot(t_, xt_, type = "o", main = "Observation")
lines(t_discrete, sol_bdi, col = "red")
points(t_discrete, sol_bdi, col = "red", pch = 19)

plot(t_continuous, xt, type = "l", main = "State")
lines(t_discrete, sol_bdi, col = "red")
points(t_discrete, sol_bdi, col = "red", pch = 19)
