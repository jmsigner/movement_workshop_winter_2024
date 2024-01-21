#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#--------------Module 02 -- Data Cleaning--------------X
#----------------Last updated 2023-03-12---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(mvtnorm)
library(amt)

# Generate real data ----
# We need to generate some location data for analysis here. We will use a 
# pure random walk, where each new location is drawn from a bivariate normal
# distribution with the previous location as its mean.

# We'll generate data for 1 individuals in the Bear River mountains 
# just east of Logan, UT. We'll assume our spatial coordinates are in 
# UTMs (zone 12).

start <- c("x" = 446589, "y" = 4625899)

# We'll simulate until we have 500 locations.
n_loc <- 500

# Initialize matrix of locations
locs <- matrix(NA, nrow = n_loc, ncol = 2, dimnames = list(NULL, c("x", "y")))
locs[1, ] <- start
head(locs)

# We need to define a variance-covariance matrix for the bivariate normal
# distribution. We'll just use a diagonal matrix with no covariance.
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

# See what we got
head(real_dat)
tail(real_dat)

# Format as 'track_xyt' for `amt`
real_trk <- real_dat %>% 
  make_track(x, y, t, truth = truth, crs = 32612)

# Plot locations
plot(real_trk, main = "True Track")
lines(real_trk)
# Start location
points(real_trk[1,], pch = 16, cex = 2, col = "blue")
# End location
points(real_trk[n_loc,], pch = 17, cex = 2, col = "orange")
# Legend
legend("topright", legend = c("Start", "End"),
       pch = c(16, 17), col = c("blue", "orange"),
       pt.cex = 2)

# Add some erroneous data ----

# ... 1. Locations before deployment ----

# My office is in the Natural Resources building on the campus of Utah State 
# University. Here are the coordinates for the building:
office <- c(432702, 4621579)

# Imagine that I activated the collar on 15 June, 2021, although the collar
# wasn't deployed until 17 June at 11:00. It would have taken these locations.

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

# See what we got
head(office_dat)

# ... 2. Capture effect ----
# We deployed our collar on 17 June at 11:00, but we know from previous
# research that our study animals show signs of capture stress for about 3
# days.

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

# See what we got
head(deployment_dat)

# ... 3. Unreachable habitat ----
# This one is tricky to simulate and easy to imagine, so we'll skip it for now.
# We'll look at this again later when we cover habitat selection.

# ... 4. Low quality fixes ----
# Let's pick 10 of our real locations to delete and substitute with low-quality
# versions.

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

# See what we got
print(lq_dat)


# ... 5. Low quality duplicates ----
# This is very similar to #4, but this time the locations are duplicated
# rather than replaced.

# Let's pick 20 of our real locations to duplicate.

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

# See what we got
head(dup_dat)

# ... 6. Unreasonably fast steps ----
# Now we need to define what we mean by "unreasonably fast".

# We will want to define this based on SDR. We can do this two ways:
#   (1.) Use the biology of our organism.
#   (2.) Find a break in the data

# Let's look at quick examples:
#   (1.) Use the biology of our organism.
#         E.g., a cheetah can run 100 km/h, but only for 60 seconds
(cheetah_sdr <- calculate_sdr(speed = 100, 
                              time = minutes(1), 
                              speed_unit = "km/h"))
# What does this imply for a 1-hour step?
# (divide by 1000 to get km)
get_displacement(cheetah_sdr, hours(1))/1000

#   (2.) Find a break in the data
# We can use the function 'sdr()' to calculate the SDR for a track_xyt object.
sdr(real_trk)
# Let's look at a histogram
hist(sdr(real_trk), breaks = 30)

# We can see a break after an SDR of 2500 m^2/s 
# (even though there is a real data point with ~ 4000 m^2/s).

# Let's make our SDR cutoff 3000.
delta <- 3000
# What's a reasonable 1-hour displacement?
get_displacement(delta, hours(1)) # >3 km

# We want to create a fast step without necessarily making a fast roundtrip
# (that's our next step). So we want to bias the location in the direction
# it was already headed -- even better if it was already a long step.

# Let's pick the 10 steps that already have an SDR between 1800 and 2500
# BUT aren't already low quality replacements or duplicates
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

# See what we got
print(fast_step_dat)


# ... 7. Unreasonably fast roundtrips ----
# These will be easier to construct, since we do want them to be far from
# both positions, p(t-1) and p(t+1). Let's pick out 5 locations.

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
# Imagine that after our last real location, the animal dropped its collar,
# and it remained there for 5 days.

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

# See what we got
head(drop_dat)

# Combine data ----
# Now we put it all together.

# Remove lq_rows, fast_step_rows, rt_rows
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

# Summarize
# How many of each location type are we left with?
table(dat$truth)

# Clean data ----

# We're finally ready to start cleaning! From here on is what the process
# looks like with real data.

# First, convert our real data to a track_xyt.
trk <- dat %>% 
  make_track(x, y, t, dop = dop, truth = truth, crs = 32612)

# Now we plot it
plot(trk, main = "Simulated Dataset")
lines(trk)

# Because we know what the "real" data look like, we can see some of the 
# errors here already. But note that they really aren't that obvious most
# of the time!

# We can also take a look on a map.
trk %>% 
  inspect()

# Now we start the cleaning process.

# ... 1. pre-deployment ----
# Recall we deployed our tag on 17 June at 11:00 and picked it up on 
# 18 July at 08:00.
trk2 <- trk %>% 
  tracked_from_to(from = as.POSIXct("2021-06-17 11:00"),
                  to = as.POSIXct("2021-07-18 08:00"))

# ... 2. capture effect ----
# Recall that we said that previous research showed the capture effect can last
# up to 3 days in our study species.
trk3 <- trk2 %>% 
  remove_capture_effect(start = days(3))

# How did we do?
head(trk3)

# Well, we actually missed the very last "capture_effect" point, but that's okay!

# ... 3. Unreachable habitat ----
# Skipping this for now

# ... 4. Low-quality fixes ----
# We want to remove any fix with very high DOP. A good cutoff for this is 10,
# but if you have a very low incidence of lower qualities, you can go lower.
hist(trk3$dop)

# Notice that most of our fixes are 2 - 4. We see a slight bump from 6 - 8,
# and then another at 10 - 12. For this illustration, we will just remove > 10,
# but on a real dataset, I would think about removing > 6.

# Recall, there is no amt function here. All we need is dplyr::filter().
trk4 <- trk3 %>% 
  filter(dop < 10)

# ... 5. Low-quality duplicates ----
# Let's check our fix rates for where there might be duplicates
summarize_sampling_rate(trk4)

# We're pretty consistent with our 1 hour fixes, but the min is very small.
sort(diff(trk4$t_))

# If we set gamma to 1 minute, we'd catch all of our duplicates. For 
# hourly data, much larger values are still safe.

# Note that this operation takes a long time on large datasets.
trk5 <- trk4 %>% 
  flag_duplicates(gamma = minutes(5))

# Recall that this function just flags duplicates -- you may decide you didn't
# like your gamma parameter and want to try again. Let's look.
trk5 %>% 
  filter(duplicate_)

# You can see we did a good job of finding our "true" duplicates here.
# Let's remove them.
trk5 <- trk5 %>% 
  filter(!duplicate_)

# ... 6. Unreasonably fast steps ----
# This is where we need to decide what "unreasonably fast" means.
# Let's look at the data to help decide.
hist(sdr(trk5), breaks = 30)

# We see a big drop off by the time we get to SDR = 5000. Let's use that.
trk6 <- trk5 %>% 
  flag_fast_steps(delta = 5000)

# How did we do?
trk6 %>% 
  filter(fast_step_)

# We're going to remove a couple of "real" locations and a "fast_rt" location,
# but we mostly got it right. Here, we can afford to lose a few "real" locations
# for the sake of cleaning out bad ones.

# Get rid of them
trk6 <- trk6 %>% 
  filter(!fast_step_)

# ... 7. Fast roundtrips ----
# We can use the delta (=5000) from before, but now we also have to decide on
# epsilon. A good way to do this is to try a few different values and see
# how many steps you flag.

# epsilon = 1
trk7a <- trk6 %>% 
  flag_roundtrips(delta = 5000, epsilon = 1)
# How many?
trk7a %>% 
  filter(fast_roundtrip_)

# epsilon = 3
trk7b <- trk6 %>% 
  flag_roundtrips(delta = 5000, epsilon = 3)
# How many?
trk7b %>% 
  filter(fast_roundtrip_)

# epsilon = 5
trk7c <- trk6 %>% 
  flag_roundtrips(delta = 5000, epsilon = 5)
# How many?
trk7c %>% 
  filter(fast_roundtrip_)

# epsilon = 10
trk7d <- trk6 %>% 
  flag_roundtrips(delta = 5000, epsilon = 10)
# How many?
trk7d %>% 
  filter(fast_roundtrip_)

# epsilon = 3 was the best for us (because we know the truth)
# If we didn't, we would need to do a lot more examining of our tracks to know
# that we would be removing a lot of real points with larger epsilons

trk7 <- trk7b %>% 
  filter(!fast_roundtrip_)

# ... 8. Defunct clusters ----
# Since we know we picked up a dropped collar, we know we want this step.

# By looking at the plot, we can see that the cluster at the end lasted for a
# while. Say we have some idea that it was there for at least two days (maybe
# from a mortality signal or a scheduled drop-off). What do the last 48
# points look like?
lasts <- tail(trk7, 48)
plot(lasts)
lines(lasts)

# This cluster is about 50 m across, so maybe we should consider zeta to be
# 25. Again, we can check what really gets flagged and adjust if necessary.

# Try zeta = 25
trk8a <- trk7 %>% 
  flag_defunct_clusters(zeta = 25,
                        eta = 24,
                        theta = hours(24))

trk8a %>% 
  filter(defunct_cluster_) %>% 
  plot()

# Try zeta = 50
trk8b <- trk7 %>% 
  flag_defunct_clusters(zeta = 50,
                        eta = 24,
                        theta = hours(24))

trk8b %>% 
  filter(defunct_cluster_) %>% 
  plot()

# Try zeta = 100
trk8c <- trk7 %>% 
  flag_defunct_clusters(zeta = 100,
                        eta = 24,
                        theta = hours(24))

trk8c %>% 
  filter(defunct_cluster_) %>% 
  plot()

# Try zeta = 200
trk8d <- trk7 %>% 
  flag_defunct_clusters(zeta = 200,
                        eta = 24,
                        theta = hours(24))

trk8d %>% 
  filter(defunct_cluster_) %>% 
  plot()

# Try zeta = 4000
trk8e <- trk7 %>% 
  flag_defunct_clusters(zeta = 4000,
                        eta = 24,
                        theta = hours(24))

trk8e %>% 
  filter(defunct_cluster_) %>% 
  plot()

# We can see that zeta = 25 is too small.
# When zeta goes from 50 -- 200, we don't see much change.
# But if we make zeta = 4000, we start getting some movement data

# 50 was fine, we'll keep it.
final <- trk8b %>% 
  filter(!defunct_cluster_)

# ... how did we do? ----
table(final$truth)

# We should have been left with 476 real points (nrow(real_dat3))
# and no other classes.

# We lost just a few real locations and kept just a couple of erroneous
# locations. 

# This is a very likely outcome with real data. The point is to do our best
# to try to make our cleaning process reproducible and our answers robust
# to that cleaning.

# Piped workflow ----
# Lastly, let's see this all at once.

final2 <- trk %>% 
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
  # Step 6
  flag_fast_steps(delta = 5000) %>% 
  filter(!fast_step_) %>% 
  # Step 7
  flag_roundtrips(delta = 5000, epsilon = 3) %>% 
  filter(!fast_roundtrip_) %>% 
  # Step 8
  flag_defunct_clusters(zeta = 50,
                        eta = 24,
                        theta = hours(24)) %>% 
  filter(!defunct_cluster_)

# Identical?
identical(final, final2)  

# Now we're ready to move on to data analysis!

