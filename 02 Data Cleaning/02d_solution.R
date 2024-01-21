#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#--------------Module 02 -- Data Cleaning--------------X
#----------------Last updated 2022-11-13---------------X
#-------------------Exercise Solution------------------X
#######################################################X

# Using amt_fisher data

# Load packages ----
library(tidyverse)
library(amt)

# Part 1 ----
# Load in your data or the `amt_fisher` dataset. If you're working with your 
# own data, format it as a `track_xyt` object using `amt::make_track()`.
dat <- amt_fisher %>% 
  filter(name == "Ricky T")

# Part 2 ----
# Clean these data! 

# ... Quick explore ----
# Plot locations
plot(dat, main = "Ricky T -- Raw")
lines(dat)
# Start location
points(dat[1,], pch = 16, cex = 2, col = "blue")
# End location
points(dat[nrow(dat),], pch = 17, cex = 2, col = "orange")
# Legend
legend("topright", legend = c("Start", "End"),
       pch = c(16, 17), col = c("blue", "orange"),
       pt.cex = 2)

# Could be a capture effect! Large cluster amidst the start points.
# I'll throw out 3 days to be safe.

# We have no DOP column. Would be useful to look at a histogram here.
# We'll pretend all locations have the same DOP. Whenever there is a tie,
# the functions will simply keep the first location (in row order).
dat$dop <- 1

# What about fix rate?
summarize_sampling_rate(dat)
# The fix rate is very high! Would be worth flagging duplicates within 1 minute
# of each other.

# How about SDRs?
hist(sdr(dat), breaks = 30)

# Really hard to see a break here, so let's get rid of the SDRs less than 250
# and look again.
big_sdr <- dat %>% 
  mutate(sdr = sdr(.)) %>% 
  filter(sdr > 250)
hist(big_sdr$sdr, breaks = 30)

# Much easier to see a break a dip just before 500 there. This could be a
# separate behavior, but for now, we'll plan to filter it.
delta <- 500

# So how far is reasonable for Ricky T to go in an hour?
get_displacement(delta, hours(1)) # Over 1.3 km.

# ... clean ----
clean <- dat %>% 
  # Step 2
  remove_capture_effect(start = days(3)) %>% 
  # Step 5
  flag_duplicates(gamma = minutes(1), DOP = "dop") %>% 
  filter(!duplicate_) %>% 
  # Step 6
  flag_fast_steps(delta = delta) %>% 
  filter(!fast_step_) %>% 
  # Step 7
  flag_roundtrips(delta = delta, epsilon = 3) %>% 
  filter(!fast_roundtrip_)

# How many points did we clean?
nrow(dat) - nrow(clean)

# And what does the track look like now?
# Plot locations
plot(clean, main = "Ricky T -- Cleaned")
lines(clean)
# Start location
points(clean[1,], pch = 16, cex = 2, col = "blue")
# End location
points(clean[nrow(clean),], pch = 17, cex = 2, col = "orange")
# Legend
legend("topright", legend = c("Start", "End"),
       pch = c(16, 17), col = c("blue", "orange"),
       pt.cex = 2)
