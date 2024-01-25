r <- terra::rast(xmin = -100, xmax = 100, 
                 ymin = -100, ymax = 100, res = 1)
r[] <- 0
r[80:100, ] <- 1
names(r) <- "x"
plot(r)


# ... Redistribution kernel ---------------------------------------------------

start <- make_start(c(0, 0), ta_ = pi/2)
m <- make_issf_model(
  coefs = c(x_end = 0), 
  sl = make_gamma_distr(shape = 2, scale = 2), 
  ta = make_vonmises_distr(kappa = 1))

rdk.1a <- redistribution_kernel(
  m, start = start, map = r, 
  stochastic = FALSE, as.rast = TRUE, max.dist = 5)
plot(rdk.1a$redistribution.kernel)

rdk.1a <- redistribution_kernel(
  m, start = start, map = r, 
  stochastic = TRUE, max.dist = 5, 
  n.control = 1e4)

# Repeat this for 50 animals
n <- 50
system.time(p1 <- replicate(n, simulate_path(rdk.1a, n = 15, 
                                             start = start), 
                            simplify = FALSE))

# Plot the trajectories
tibble(
  rep = 1:n, 
  path = p1
) |> unnest(cols = path) |> 
  ggplot(aes(x_, y_, group = rep)) + geom_path(alpha = 0.1) +
  coord_equal()

# Smooth output at different points in time
trks <- lapply(c(5, 10, 15), function(i) {
  tibble(
    rep = 1:n, 
    path = map(p1, ~ dplyr::slice(.x, i))
  ) |> unnest(cols = path) |> 
    make_track(x_, y_) |> hr_kde(trast = r)
})
plts <- terra::rast(lapply(trks, hr_ud))
names(plts) <- paste("n =", c("5", "10", "15"))
terra::plot(plts)


# Now with a random starting direction


n <- 50
system.time(p1 <- replicate(n, simulate_path(rdk.1a, n = 15, 
                                             start = make_start(c(0, 0), ta = runif(1, -pi, pi))), 
                            simplify = FALSE))

# Plot the trajectories
tibble(
  rep = 1:n, 
  path = p1
) |> unnest(cols = path) |> 
  ggplot(aes(x_, y_, group = rep)) + geom_path(alpha = 0.1) +
  coord_equal()

# Smooth output at different points in time
trks <- lapply(c(5, 10, 15), function(i) {
  tibble(
    rep = 1:n, 
    path = map(p1, ~ dplyr::slice(.x, i))
  ) |> unnest(cols = path) |> 
    make_track(x_, y_) |> hr_kde(trast = r)
})
plts <- terra::rast(lapply(trks, hr_ud))
names(plts) <- paste("n =", c("5", "10", "15"))
terra::plot(plts)
