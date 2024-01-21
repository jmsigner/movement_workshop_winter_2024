#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#--------------Module 02 -- Data Cleaning--------------X
#----------------Last updated 2023-03-13---------------X
#----------------Bonus Script: aniMotum----------------X
#######################################################X

# We talked about the new R package 'aniMotum'
# Paper describing the package here: https://doi.org/10.1111/2041-210X.14060

# You need compiler tools on Windows (RTools) or Mac (Xcode) to
# install. If you need more detailed installation instructions, 
# check here: https://ianjonsen.github.io/aniMotum/index.html#installation

# If you need to install from GitHub:
# remotes::install_github("ianjonsen/aniMotum")

# Load packages ----
library(tidyverse)
library(mvtnorm)
library(amt)
library(aniMotum)

# Generate example data ----
# This is the exact same simulation of erroneous data as in '02b_walkthrough.R',
# but without all the extra comments.

# The curly braces ({}) run the whole code chunk at once:
{
  start <- c("x" = 446589, "y" = 4625899)
  n_loc <- 500
  
  # Initialize matrix of locations
  locs <- matrix(NA, nrow = n_loc, ncol = 2, dimnames = list(NULL, c("x", "y")))
  locs[1, ] <- start
  
  Sigma <- diag(2) * 1000^2
  
  # Set random seed
  set.seed(20221114)
  
  # Loop to generate new locations
  for (t in 2:n_loc) {
    locs[t, ] <- rmvnorm(n = 1, mean = locs[(t-1), ], sigma = Sigma)
  }
  
  # Convert to data.frame
  real_dat <- as.data.frame(locs)
  
  # Assign actual dates and times (start on 20 June, 2021)
  real_dat$t <- as.POSIXct("2021-06-20 12:00:00") +
    1:n_loc * 60 * 60
  
  # Randomly generate some DOPs
  # These are good locations, so keep DOP small
  real_dat$dop <- runif(n = n_loc,
                        min = 1.2,
                        max = 3.5) %>% 
    round(1)
  
  # Add a bit of noise to our otherwise perfect timestamps
  real_dat$t <- real_dat$t + runif(n = nrow(real_dat),
                                   min = 55, max = 80)
  
  # Add a column for the truth
  real_dat$truth <- "real"
  
  # Convert to track_xyt
  real_trk <- real_dat %>% 
    make_track(x, y, t, truth = truth, crs = 32612)
  
  # ... 1. Locations before deployment ----
  office <- c(432702, 4621579)
  
  # Times
  office_times <- seq(as.POSIXct("2021-06-15 12:00"),
                      as.POSIXct("2021-06-17 10:00"),
                      by = "1 hour")
  # Random locations (with some GPS error)
  office_locs <- rmvnorm(n = length(office_times),
                         mean = office,
                         sigma = diag(2) * 5^2)
  # Convert to data.frame and add time
  office_dat <- as.data.frame(office_locs)
  names(office_dat) <- c("x", "y")
  office_dat$t <- office_times
  
  # Add DOP
  # These are probably indoors and may not get good reception
  office_dat$dop <- runif(n = nrow(office_dat),
                          min = 6,
                          max = 8) %>% 
    round(1)
  
  # Add a column for the truth
  office_dat$truth <- "pre_deployment"
  
  # ... 2. Capture effect ----
  deployment <- start - c(250, 250)
  
  # Times
  deployment_times <- seq(as.POSIXct("2021-06-17 11:00"),
                          as.POSIXct("2021-06-20 11:00"),
                          by = "1 hour")
  
  # Random locations (with some small movements)
  deployment_locs <- rmvnorm(n = length(deployment_times),
                             mean = deployment,
                             sigma = diag(2) * 50^2)
  # Convert to data.frame and add time
  deployment_dat <- as.data.frame(deployment_locs)
  names(deployment_dat) <- c("x", "y")
  deployment_dat$t <- deployment_times
  
  # Add DOP
  # We'll assume these are decent locations
  deployment_dat$dop <- runif(n = nrow(deployment_dat),
                              min = 2,
                              max = 5) %>% 
    round(1)
  
  # Add column for the truth
  deployment_dat$truth <- "capture_effect"
  
  # ... 4. Low quality fixes ----
  # Pick 10 rows
  lq_rows <- sample.int(n = nrow(real_dat), size = 10, replace = FALSE)
  
  # Pull them out
  lq_dat <- real_dat[lq_rows, ]
  
  # Remove them from the real data
  real_dat2 <- real_dat[-lq_rows, ]
  
  # Make the DOP large
  lq_dat$dop <- runif(n = nrow(lq_dat),
                      min = 10,
                      max = 12) %>% 
    round(1)
  
  # Add location error
  lq_dat$x <- rnorm(n = nrow(lq_dat),
                    mean = lq_dat$x, 
                    # Note that DOP is a multiplier of the SD
                    sd = (lq_dat$dop * 25))
  
  lq_dat$y <- rnorm(n = nrow(lq_dat),
                    mean = lq_dat$y, 
                    # Note that DOP is a multiplier of the SD
                    sd = (lq_dat$dop * 25))
  
  # Replace the truth
  lq_dat$truth <- "low_quality"
  
  
  # ... 5. Low quality duplicates ----
  # Pick 20 rows
  dup_rows <- sample.int(n = nrow(real_dat2), size = 20, replace = FALSE)
  
  # Pull them out
  dup_dat <- real_dat2[dup_rows, ]
  
  # Make the DOP larger than the real data (but not unacceptable)
  dup_dat$dop <- runif(n = nrow(dup_dat),
                       min = 6,
                       max = 8) %>% 
    round(1)
  
  # Add location error
  dup_dat$x <- rnorm(n = nrow(dup_dat),
                     mean = dup_dat$x, 
                     # Note that DOP is a multiplier of the SD
                     # Also note that we're using SD, not variance here
                     sd = (dup_dat$dop * 25))
  
  dup_dat$y <- rnorm(n = nrow(dup_dat),
                     mean = dup_dat$y, 
                     # Note that DOP is a multiplier of the SD
                     # Also note that we're using SD, not variance here
                     sd = (dup_dat$dop * 25))
  
  # Jitter the timestamp by a few seconds
  dup_dat$t <- dup_dat$t + 
    round(runif(n = nrow(dup_dat),
                min = -30,
                max = 15))
  
  # Replace the truth
  dup_dat$truth <- "duplicate"
  
  # ... 6. Unreasonably fast steps ----
  # Let's make our SDR cutoff 3000.
  delta <- 3000
  
  fast_step_rows <- which(sdr(real_trk) > 1700 & 
                            sdr(real_trk) < 2500)
  # Are any of these already low-quality replacements
  sum(fast_step_rows %in% lq_rows)
  # Or low-quality duplicates?
  sum(fast_step_rows %in% dup_rows)
  ## Yep, get rid of it
  fast_step_rows <- fast_step_rows[which(!fast_step_rows %in% dup_rows)]
  
  # Get those rows
  fast_step_dat <- real_dat[fast_step_rows, ]
  
  # Loop through those rows and add large distance to the previous location
  for (i in 1:length(fast_step_rows)) {
    # Get the current row number
    r <- fast_step_rows[i]
    # Get the current location
    c <- fast_step_dat[i, ]
    # Get the previous location
    p <- real_dat[r-1, ]
    # Figure out if real step is +/-
    dir_x <- (c$x - p$x)/abs(c$x - p$x)
    dir_y <- (c$y - p$y)/abs(c$y - p$y)
    
    # Template for new location
    n <- c
    
    # Add 5 km to each location by adding +/- 3k to x and +/- 4k to y
    n$x <- c$x + dir_x * 3000
    n$y <- c$y + dir_y * 4000
    # New DOP
    n$dop <- runif(n = 1, min = 6, max = 8) %>% 
      round(1)
    # New truth
    n$truth <- "fast_step"
    
    # Replace new row into data
    fast_step_dat[i, ] <- n[1,]
  }
  
  # ... 7. Unreasonably fast roundtrips ----
  # Pick out 5 locations.
  
  rt_rows <- sample(x = 2:nrow(real_dat),
                    size = 5, replace = FALSE)
  
  # Check that none are already manipulated
  # Are any of these already low-quality replacements
  sum(rt_rows %in% lq_rows)
  # Or low-quality duplicates?
  sum(rt_rows %in% dup_rows)
  # Or fast steps?
  sum(rt_rows %in% fast_step_rows)
  
  # Get the rows
  rt_dat <- real_dat[rt_rows, ]
  
  # Add 3k to the x- and y-coords
  rt_dat$x <- rt_dat$x + 3000
  rt_dat$y <- rt_dat$y + 3000
  
  # Lower the DOP
  rt_dat$dop <- runif(n = nrow(rt_dat),
                      min = 6, max = 8) %>% 
    round(1)
  
  # Update the truth
  rt_dat$truth <- "fast_rt"
  
  # ... 8. Defunct clusters ----
  drop <- as.matrix(real_dat[nrow(real_dat), c("x", "y")])
  
  # Times
  drop_times <- seq(real_dat$t[nrow(real_dat)],
                    real_dat$t[nrow(real_dat)] + days(7),
                    by = "1 hour")
  
  # Random locations (with some location error)
  drop_locs <- rmvnorm(n = length(drop_times),
                       mean = drop,
                       sigma = diag(2) * 10^2)
  # Convert to data.frame and add time
  drop_dat <- as.data.frame(drop_locs)
  names(drop_dat) <- c("x", "y")
  drop_dat$t <- drop_times
  
  # Add DOP
  # We'll assume these are decent locations
  drop_dat$dop <- runif(n = nrow(drop_dat),
                        min = 2,
                        max = 5) %>% 
    round(1)
  
  # Add column for the truth
  drop_dat$truth <- "defunct"
  
  # ... combine ----
  real_dat3 <- real_dat[-c(fast_step_rows, lq_rows, rt_rows), ]
  
  # Combine
  dat <- rbind(real_dat3,
               office_dat,
               deployment_dat,
               lq_dat, 
               dup_dat,
               fast_step_dat,
               rt_dat,
               drop_dat) %>% 
    # Sort by timestamp
    arrange(t)
  
}

# Preprocess data ----
# aniMotum doesn't deal with GPS DOP metrics or other sources of error. 
# We can do some pre-filtering using the usual 'amt' workflow to get started.
# (I.e., up through step 5, then step 8; SSM takes care of 6 and 7)
pre <- dat %>%
  make_track(x, y, t, crs = 32612, all_cols = TRUE) %>% 
  # Step 1
  tracked_from_to(from = as.POSIXct("2021-06-17 11:00"),
                  to = as.POSIXct("2021-07-18 08:00")) %>% 
  # Step 2
  remove_capture_effect(start = days(3)) %>% 
  # Step 4
  filter(dop < 10) %>% 
  # Step 5
  flag_duplicates(gamma = minutes(5)) %>% 
  filter(!duplicate_) %>% 
  # Step 8
  flag_defunct_clusters(zeta = 50,
                        eta = 24,
                        theta = hours(24)) %>% 
  filter(!defunct_cluster_)

# One data format that 'aniMotum' will accept is an sf tibble.
# amt can convert a track to this format:
pre_sf <- as_sf(pre)

# The package requires some columns we don't have
# Individual ID
pre_sf$id <- "animal_1"
# Location class (G for GPS)
pre_sf$lc <- "G"

# Format data for 'aniMotum'
pre_ani <- format_data(pre_sf, date = "t_")

# Fit SSM ----
# Now fit the state-space model
# Using an ordinary random walk because that's how we simulated the data!
m <- fit_ssm(pre_ani, model = "rw", time.step = 1,
             # Max speed (sqrt of our delta)
             vmax = sqrt(3000))

# Predict ----
# Plot fitted model
plot(m, what = "fitted", type = 1)
# Note the black x's are points that were pre-filtered

# Map track (this is a ggplot, so we can add our own layers)
map(m, what = "predicted")

# Add our real points on top of the SSM points
map(m, what = "predicted") +
  geom_point(data = real_dat, aes(x = x, y = y), inherit.aes = FALSE,
             size = 0.1)

# Not bad!