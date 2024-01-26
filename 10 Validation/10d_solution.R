#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Module 10 -- Validation---------------X
#----------------Last updated 2021-01-27---------------X
#-------------------Exercise Solution------------------X
#######################################################X

# Refitting the cougar iSSF from module 8.

# Load packages ----
library(tidyverse)
library(amt)
library(lubridate)
library(terra)
library(survMisc)

# 1. Load data ----
# Location data
dat <- read_csv("data/coyote_cougar.csv") %>% 
  # Subset to just cougar F53
  filter(id == "F53")

# Set the timezone
tz(dat$t_) <- "US/Mountain"

# Habitat data
hab <- rast("data/coyote_cougar_habitat.tif")
names(hab) <- c("elevation", "trees", "biomass", "dist_to_road")

# Prep data
set.seed(1)
mod_dat <- dat %>% 
  # Make track
  make_track(x_, y_, t_, crs = 32612) %>% 
  # Find regular bursts
  track_resample(rate = hours(4), tolerance = minutes(15)) %>% 
  filter_min_n_burst(min_n = 3) %>% 
  # Convert points to steps
  steps_by_burst() %>% 
  # Generate available steps
  random_steps(n_control = 20) %>% 
  # Attach habitat variables
  extract_covariates(hab) %>% 
  # Add additional movement covariates
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_))

head(mod_dat)

# 2. Separate into training and testing ----
set.seed(2)
step_df <- data.frame(stratum = sort(unique(mod_dat$step_id_))) %>% 
  mutate(train = as.logical(rbinom(n = nrow(.), size = 1, prob = 0.8)))

train <- mod_dat %>% 
  filter(step_id_ %in% step_df$stratum[step_df$train])

test <- mod_dat %>% 
  filter(step_id_ %in% step_df$stratum[!step_df$train]) 

# Drop the strata with NA for turn angle
NA_strata <- test %>% 
  filter(is.na(ta_)) %>% 
  pull(step_id_)

test <- test %>% 
  filter(!step_id_ %in% NA_strata)

nrow(train)/nrow(mod_dat)
nrow(test)/nrow(mod_dat)

# 3. Fit models ----
m1 <- train %>% 
  fit_issf(case_ ~ 
             # Habitat selection main effects
             elevation + I(elevation^2) + 
             trees + I(trees^2) + dist_to_road +
             # Movement main effects
             sl_ + log_sl_ + cos_ta_ +
             # Don't forget the strata
             strata(step_id_),
           # And include model = TRUE so we can use 'log_rss()' later
           model = TRUE)

summary(m1)

# Drop trees
m2 <- train %>% 
  fit_issf(case_ ~ 
             # Habitat selection main effects
             # (removed elevation^2 and trees + trees^2)
             elevation + dist_to_road +
             # Movement main effects
             sl_ + log_sl_ + cos_ta_ +
             # Don't forget the strata
             strata(step_id_),
           # And include model = TRUE so we can use 'log_rss()' later
           model = TRUE)

summary(m2)

# Evaluate model ----

# ... concordance ----
# Concordance is already in the summary
# Analogous to AUC
summary(m1)$concordance
summary(m2)$concordance

# We can also get it with a specific function
survival::concordance(m1$model)$concordance
survival::concordance(m2$model)$concordance

# Overall, fairly low

# ... MER ----
# Measure of explained randomness
survMisc::rsq(m1$model)
survMisc::rsq(m2$model)

# Quite low

# ... UHC plots ----

uhc1 <- prep_uhc(m1, 
                 test_dat = dplyr::select(test, -burst_), 
                 n_samp = 500) # n = 500 for speed, but should use more
uhc2 <- prep_uhc(m2, 
                 test_dat = dplyr::select(test, -burst_),
                 n_samp = 500) # n = 500 for speed, but should use more

# Plot
uhc1 %>% 
  as.data.frame() %>% 
  conf_envelope(levels = c(0.9, 0.95, 1)) %>% 
  plot()

# Notice that even without trees in the model, the UHC plot for trees
# looks fine. This is a good indication that we don't need it!
uhc2 %>% 
  as.data.frame() %>% 
  conf_envelope(levels = c(0.9, 0.95, 1)) %>% 
  plot()
