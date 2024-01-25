#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#--------------- Module 08 -- Simulate ----------------X
#--------------- Last updated 2024-01-24 --------------X
#------------------ Code Walkthrough ------------------X
#######################################################X

library(lubridate)
library(tidyverse)
library(amt)
library(terra)
library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Simulate from scratch
# ... Landscape ------------------------------------------------------------

r <- terra::rast(xmin = -100, xmax = 100, 
                 ymin = -100, ymax = 100, res = 1)
r[] <- 0
r[80:100, ] <- 1
names(r) <- "x"
plot(r)


# ... Redistribution kernel ---------------------------------------------------

start <- make_start(c(0, 0), ta_ = pi/2)
m <- make_issf_model(
  coefs = c(x_end = 2), 
  sl = make_gamma_distr(shape = 2, scale = 2), 
  ta = make_vonmises_distr(kappa = 1))

rdk.1a <- redistribution_kernel(
  m, start = start, map = r, 
  landscape = "discrete", as.rast = TRUE, max.dist = 5)
plot(rdk.1a$redistribution.kernel)

rdk.1a <- redistribution_kernel(
  m, start = start, map = r, 
  max.dist = 5, 
  n.control = 1e4)



# Now simulate one path (i.e. one animal)
p1 <- simulate_path(rdk.1a, n.steps = 20, start = start)
plot(r)
lines(p1$x_, p1$y_)

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
  coord_equal() +
  xlim(-100, 100) +
  ylim(-100, 100)

# Smooth output at different points in time
trks <- lapply(c(2, 5, 10, 15), function(i) {
  tibble(
    rep = 1:n, 
    path = map(p1, ~ dplyr::slice(.x, i))
  ) |> unnest(cols = path) |> 
    make_track(x_, y_) |> hr_kde(trast = r)
})
plts <- terra::rast(lapply(trks, hr_ud))
names(plts) <- paste("n =", c("2", "5", "10", "15"))
terra::plot(plts)

# Contrast this with unconstrained movement

m <- make_issf_model(
  coefs = c(x_end = 0), # No more habitat selection
  sl = make_gamma_distr(shape = 2, scale = 2), 
  ta = make_vonmises_distr(kappa = 1))

rdk.1a <- redistribution_kernel(
  m, start = start, map = r, 
  landscape = "discrete", as.rast = TRUE, max.dist = 5)
plot(rdk.1a$redistribution.kernel)

rdk.1a <- redistribution_kernel(
  m, start = start, map = r, 
  max.dist = 5, 
  n.control = 1e4)

# Now simulate one path (i.e. one animal)
p1 <- simulate_path(rdk.1a, n.steps = 20, start = start)
plot(r)
lines(p1$x_, p1$y_)

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
  coord_equal() +
  xlim(-100, 100) +
  ylim(-100, 100)

# Smooth output at different points in time
trks <- lapply(c(2, 5, 10, 15), function(i) {
  tibble(
    rep = 1:n, 
    path = map(p1, ~ dplyr::slice(.x, i))
  ) |> unnest(cols = path) |> 
    make_track(x_, y_) |> hr_kde(trast = r)
})
plts1 <- terra::rast(lapply(trks, hr_ud))
names(plts) <- paste("n =", c("2", "5", "10", "15"))

terra::plot(plts1)
terra::plot(plts) # constrained movement as comparison from before‚


# African Buffalo example ----

set.seed(123333)

# Preparing the data ----

# We will use an example data set of a African buffalo here. This data set has
# been used as a case study in many other papers.
cilla <- read_rds("data/cilla.rds")
env <- rast(read_rds("data/env_covar.rds"))

# The first step is to prepare the data to fit an iSSF model. `cilla` is already
# and `amt` track.
ssf_cilla <- cilla %>% steps_by_burst() %>% random_steps(n_control = 100) %>% 
  extract_covariates(env, where = "both") %>% 
  mutate(hour = hour(t1_) + minute(t1_) / 60) %>% 
  filter(complete.cases(.))

# Modelling ----
# ... First model ----
# Next, we fit a model with all movement relevant covariates for full movement
# kernel and only one habitat covariate (the distance to water).
m_0 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + 
                    water_dist_end + strata(step_id_))
summary(m_0)

# ... Simulation ---- 
# We have to define a start step (i.e., we need the start coordinates, a
# direction and sampling rate, which will be used for all consecutive steps).
start <- make_start(ssf_cilla[1, ])

# We can then create a redistribution kernel (i.e., a probability surface of
# where the animal will move next).
k1 <- redistribution_kernel(m_0, map = env, start = start, as.rast = TRUE, 
                            landscape = "discrete")
plot(k1$redistribution.kernel)

# There are a number of options that can be set here: 
# See table 1 here: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14263
k1 <- redistribution_kernel(m_0, map = env, start = start,
                            landscape = "continuous", tolerance.outside = 0.2, as.rast = FALSE, 
                            n.control = 1e3)
# Or we could draw 1000 100 steps 
k1 <- redistribution_kernel(m_0, map = env, start = start,
                            tolerance.outside = 0.2, as.rast = FALSE, 
                            n.control = 1e4, n.sample = 1000)

# Lets visualize this
plot(k1$redistribution.kernel[, 1:2], asp = 1)
points(start$x_, start$y_, col = "red", pch = 20)

# Finally, we can repeat this step many times to simulate a path (i.e., a series
# of consecutive steps). We will simulate 500 steps. It should take approx 10
# seconds to run.
tictoc::tic()
s1 <- simulate_path(k1, n.steps = 500)
tictoc::toc()

# The resulting object is a `track_xyt`, like an object that we obtained from
# the `make_track()`.
s1


# Next we can visualize the results and compare it with the real trajectory.
raster::plot(env[["elev"]])
with(cilla[1:501, ], lines(x_, y_))
lines(s1$x_, s1$y_, col = "red")

# We see that the trajectory shows some kind of home-range behavior.

# Second model ----- 
# We will make the model now slightly more complex and add home ranging behavior
# (credits for this way to consider home-ranging behavior are to Chris
# Flemming). We include x2_, y2_ and their squared product in the model and this
# will pull the animal always back.
m_1 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + 
                    water_dist_end + x2_ + y2_ + I(x2_^2 + y2_^2) + strata(step_id_))
summary(m_1)

k2 <- redistribution_kernel(m_1, map = env, start = make_start(ssf_cilla[1, ]),
                            tolerance.outside = 0.01, as.rast = FALSE, n.sample = 100, 
                            n.control = 1e3)
s2 <- simulate_path(k2, n.steps = 200)

raster::plot(env[["elev"]])
lines(cilla$x_, cilla$y_)
lines(s2$x_, s2$y_, col = "red")

# This looks much better, but there is still a rive that the animal crosses in
# the simulated data set but does not do so in reality.
raster::plot(env$water_dist < 100)
lines(cilla$x_, cilla$y_)
lines(s2$x_, s2$y_, col = "red")

# We can include the river as a barrier

# ... Third Model ----

# First lets create a new covariate, that distinguished the two sides of the
# river.
water <- env[["water_dist"]] > 100
water <- crop(water, water - 5000)
plot(water)

# We can use the connected components algorigthm that is implemented in the
# raster package with the `clump()` function.
water <- env[["water_dist"]] > 100
water <- crop(water, ext(water) - 5000)
plot(water)
ww <- patches(water, zeroAsNA = TRUE)
plot(ww)
ww[ww == 3] <- 2

# And give the new layer a meaningful name and add it to the raster stack.
names(ww) <- "water_crossed"
env1 <- c(ww, crop(env, ww))

# We have to extract the covariates again, to also have the side of the river,
# where a step started and ended.
ssf_cilla <- cilla %>% steps_by_burst() %>% random_steps() %>% 
  # Note that we use `where = "both"` here, in order to get the river side at
  # the satart and at the end of each step.
  extract_covariates(env1, where = "both") 

# We can coulbe check this again
ssf_cilla %>% select(starts_with("water"))

# In the model we now create a new variable on the fly, that checks if the start
# and end position are on the same side of the river.
m_3 <- fit_clogit(ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + 
                    water_dist_end + x2_ + y2_ + I(x2_^2 + y2_^2) +
                    I(water_crossed_end != water_crossed_start) + strata(step_id_))

# And see the model summary
summary(m_3)

# Note, we get an Inf large upper confidence interval for the river crossing.
# This is because we did not observe a single crossing in the observed data.
ssf_cilla %>% mutate(cross = water_crossed_end != water_crossed_start) %>% 
  filter(case_, cross)
ssf_cilla %>% mutate(cross = water_crossed_end != water_crossed_start) %>% 
  filter(!case_, cross)

# Again we to create a redistribution kernel
k3 <- redistribution_kernel(
  m_3, map = env1, start = make_start(ssf_cilla[1, ]),
  stochastic = TRUE, tolerance.outside = 0.1, as.rast = FALSE, 
  n.sample = 100, n.control = 1e3)

# And we can then simulate a trajectory from this kernel.
s3 <- simulate_path(k3, n.steps = 500)

# And finally plot the result
plot(env1[["water_crossed"]])
lines(cilla$x_, cilla$y_)
lines(s3$x_, s3$y_, col = "red")

# More complex model ---- 
# The function `redistribution_kernel()` has a further argument that is called
# `fun`. By default, fun will extract all covariates provided by `map` at the
# beginning and at the end. But we could pass another function. Lets revisit the
# previous example and do the calculation -- whether or not the buffalo crossed
# the river -- before. And save this result in a new variable called `crossed`.
ssf_cilla <- ssf_cilla %>% 
  mutate(crossed = water_crossed_end != water_crossed_start)

# We can fit the model with the new variable
m_3a <- fit_clogit(
  ssf_cilla, case_ ~ cos(ta_) + sl_ + log(sl_) + 
    water_dist_end + x2_ + y2_ + I(x2_^2 + y2_^2) +
    crossed + strata(step_id_))

# And compare coefficients (which are the same)
coef(m_3)
coef(m_3a)

coef(m_3) == coef(m_3a)

# Calculating the redistribution kernel again, will result in an error:
k3 <- redistribution_kernel(
  m_3a, map = env1, start = make_start(ssf_cilla[1, ]),
  tolerance.outside = 0.1, as.rast = FALSE, 
  n.sample = 100, n.control = 1e3)

# Because the new variable crossed is not found. We have to adjust the function
# argument and calculate it, as we did before.
k3 <- redistribution_kernel(
  m_3a, 
  fun = function(x, map) extract_covariates(x, map, where = "both") %>% 
    mutate(crossed = water_crossed_end != water_crossed_start),
  map = env1, start = make_start(ssf_cilla[1, ]),
  stochastic = TRUE, tolerance.outside = 0.1, as.rast = FALSE, 
  n.sample = 100, n.control = 1e3)

