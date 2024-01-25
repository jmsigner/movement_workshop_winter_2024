#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Module 08 -- iSSF pt 2----------------X
#----------------Last updated 2023-03-16---------------X
#-------------------Exercise Solution------------------X
#######################################################X

# Using UT cougar data

# Load packages ----
library(tidyverse)
library(amt)
library(lubridate)
library(terra)

# Load data ----
# Location data
dat <- read_csv("data/coyote_cougar.csv") |> 
  # Subset to just cougar F64
  filter(id == "F64")

# Set the timezone
tz(dat$t_) <- "US/Mountain"

# Habitat data
hab <- rast("data/coyote_cougar_habitat.tif")
names(hab) <- c("elevation", "trees", "biomass", "dist_to_road")

# Format as track_xyt
trk <- dat |> 
  make_track(x_, y_, t_, crs = 32612)

# Do we have regular steps?
summarize_sampling_rate(trk)

# The median sampling rate is 4h, but the track isn't perfectly regular.
# Let's separate into bursts of 4h steps, with a tolerance of 15 minutes.
# Only keep bursts that have at least 3 relocations. Then we can turn those
# locations into steps.

stp <- track_resample(trk, rate = hours(4), tolerance = minutes(15)) |> 
  filter_min_n_burst(min_n = 3) |> 
  steps_by_burst()

# What do our step durations look like now?
hist(as.numeric(stp$dt_)) #good

# Fit iSSF ----
issf_dat <- stp |> 
  random_steps(n_control = 50) |> 
  extract_covariates(hab) |>
  # step_id_ needs to account for burst_
  mutate(step_id_ = paste(burst_, step_id_, sep = "_"))

# Simple model
m <- issf_dat |> 
  fit_issf(case_ ~ 
             # Habitat
             elevation + #I(elevation^2) +
             # Movement
             sl_ + log(sl_) + cos(ta_) +
             # Don't forget the strata
             strata(step_id_),
           model = TRUE)

summary(m) 

# Visualize selection ----
# Selection for elevation
x1 <- data.frame(elevation = seq(1600, 3400, length.out = 100),
                 sl_ = 100,
                 ta_ = 0)
x2 <- data.frame(elevation = 2000, 
                 sl_ = 100,
                 ta_ = 0)
lr <- log_rss(m, x1, x2)
# Default plot:
plot(lr)

# Simulate ----
# Starting place and time
start <- make_start(trk[1, ], dt = hours(4))

# Redistribution kernel
k <- redistribution_kernel(m, map = hab, start = start, 
                            fun = function(xy, map) {
                              extract_covariates(xy, map, where = "end")
                            })

# Simulate
s <- simulate_path(k, n.steps = 500)

# Plot
plot(hab$elevation)
lines(s, col = "red")
