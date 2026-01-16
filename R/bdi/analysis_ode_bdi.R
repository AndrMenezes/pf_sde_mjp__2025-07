#
ode_solution_bdi <- function(tt, x0, theta) {
  mu <- theta[1]
  lambda <- theta[2]
  gama <- theta[3]
  gama / (lambda - mu) + (x0 - gama / (lambda - mu)) * exp((mu - lambda) * tt)
}

# time points
tt <- seq.int(0, 30)

# initial value
x0 <- 500

# exponential decay for lambda > mu
theta_2 <- c(0.5, 0.8, 30) # lambda > mu
theta_2[3] / (theta_2[2] - theta_2[1])
sol_bdi_2 <- ode_solution_bdi(tt, x0, theta_2)
plot(tt, sol_bdi_2, ylim = c(0, x0))
abline(h = 100, lty = 2)


# theta_1 <- c(1, .5, 1) # lambda < mu
# theta_1[3] / (theta_1[2] - theta_1[1])
# sol_bdi_1 <- ode_solution_bdi(tt, x0, theta_1)
# plot(tt, ode_bdi_1, ylim = c(0, max(sol_bdi_1)))
