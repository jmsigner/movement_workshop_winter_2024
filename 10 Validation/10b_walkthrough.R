#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------------Module 10 -- Validation---------------X
#----------------Last updated 2023-03-16---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(terra)
library(tidyverse)
library(amt)
library(sf)
library(rcompanion)
library(pROC)
library(survMisc)

# Simulated data ----
# We saved the dataset we simulated in module 5.
# We can load it here.
gps <- read.csv("05 HSF/module05_HSF_data.csv")

# Load habitat rasters
hab <- rast("05 HSF/geo/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")

# Cover is a factor, let's code it that way
levels(hab[[4]]) <- data.frame(id = 1:3,
                               cover = c("grassland", "forest", "wetland"))

# Split into training and testing ----
# We want to do some out-of-sample validation, so let's withhold ~20% of our
# data for testing our model, using the remaining 80% to fit the model.
set.seed(20220128)
gps$train <- rbinom(n = nrow(gps), size = 1, prob = 0.8)

train <- gps |> 
  filter(train == 1) |>
  make_track(x, y, crs = 32612)

test <- gps |> 
  filter(train == 0) |>
  make_track(x, y, crs = 32612)

nrow(train)/nrow(gps)
nrow(test)/nrow(gps)

# Let's format both our training and testing data for an HSF.
# We want to generate available points across our entire raster.
# We'll use a polygon for that.

r_poly <- st_bbox(hab) |> 
  st_as_sfc() |> 
  st_sf()

set.seed(123456)

# Format training
train <- r_poly |>
  random_points(n = nrow(train) * 100, presence = train) |>
  extract_covariates(hab) |> 
  # Assign large weights to available points
  mutate(weight = ifelse(case_, 1, 1e5))

# Format testing
test <- r_poly |>
  random_points(n = nrow(test) * 100, presence = test) |>
  extract_covariates(hab) |> 
  # Assign large weights to available points
  mutate(weight = ifelse(case_, 1, 1e5))

# Fit model ----
# Let's fit two models. One that is correctly specified, and one that is
# missing the quadratic term for temperature, predators, and cover.

# Correct
m1 <- glm(case_ ~ forage + temp + I(temp^2) + predator + cover,
          family = binomial(), weights = weight, data = train)
# Incorrect
m2 <- glm(case_ ~ forage + temp,
          family = binomial(), weights = weight, data = train)


# We're ready to evaluate our model

# Model evaluation ----
# ... Pseudo R-squared ----

# Before we jump into measures of external validity, let's take a 
# quick look at a measure of internal validity, i.e., how well does
# our model fit to the data we actually used to fit it.

?rcompanion::nagelkerke

(m1r2 <- nagelkerke(m1)) # Quite low
(m2r2 <- nagelkerke(m2)) # Even lower

cbind(m1r2$Pseudo.R.squared.for.model.vs.null,
      m2r2$Pseudo.R.squared.for.model.vs.null)

# ... AUC ----

# One method that's commonly used to assess models for binary responses 
# is the Area Under the Curve metric, or AUC. It measures the area under the 
# receiver-operator curve (ROC), which is a good measure of the trade-off 
# between sensitivity and specificity in a binary predictor.

# Note that you can calculate AUC as a measure of internal validity (without
# a testing set) or as a measure of external validity (on a testing set).

# We'll use the package `pROC` here, but there are *many* R packages that
# can calculate this statistic for you.

# We will use the function 'roc()' to produce the ROC. We can then plot it
# or calculate the area under it.

# We need to provide 'roc()' two arguments:
#   1. 'response' the binary data used to fit the model
#   2. 'predictor' the predicted probability under the model
#       In our case this is the linear prediction from our model,
#       backtransformed using the inverse-logit.

m1_roc <- roc(response = as.numeric(test$case_),
              predictor = predict(m1, 
                                  newdata = test,
                                  type = "response"))

# Now we can plot the curve
plot(m1_roc)

# We can also calculate AUC
auc(m1_roc)

# Not bad!

m2_roc <- roc(response = as.numeric(test$case_),
              predictor = predict(m2, 
                                  newdata = test,
                                  type = "response"))
auc(m2_roc) # lower

# ... Boyce method ----
# The general approach with binning is to aggregate the testing data into
# bins of:
#   1. equal area (geographic space)
#   2. equal count (geographic space)
#   3. equal interval (environmental space)

# Here, we'll use equal-interval bins.

# We also need to know which cell of the landscape each point falls in
# so that we can adjust our expected number of points by the area covered.

# Let's also break our habitat variables into 5 bins. Note that cover is
# already in discrete classes, so we don't need to bin.

# The base R function 'cut()' can bin our data for us, but the tidyverse
# function 'cut_interval()' will make it easier to control the details
?cut_interval

# Let's bin and summarize:
test_bins <- test |>
  # Keep only the observed points
  filter(case_) %>%
  # Assign cell number
  mutate(cell = cellFromXY(hab, as.matrix(.[, c("x_", "y_")]))) |> 
  # Bin variables
  mutate(forage_bin = cut_interval(forage, n = 5),
         temp_bin = cut_interval(temp, n = 5),
         predator_bin = cut_interval(predator, n = 5)) |> 
  # Now group by bins
  dplyr::group_by(forage_bin, temp_bin, predator_bin, cover) |> 
  # Summarize observations
  summarize(cells = n_distinct(cell),
            obs = sum(case_),
            obs_per_cell = obs/cells) |> 
  # Ungroup
  ungroup()

# What did we get?
test_bins

# You can see we have the number of observed points by habitat 
# bin, as well as the density of points. Now how can we compare this to 
# our model?

# The basic idea is that we want to predict the value of each habitat and
# elevation value using our HSF, and then see how strong the correlation
# is between the HSF and the density of observed points.

# Let's convert our text label for each bin into the mean value
# for that bin. Since we don't have nice, round numbers from 'cut_interval()', 
# we need some string manipulation.

# Here's a function that will do it for us.
get_mean <- function(x) {
  # Get rid of parentheses
  x <- gsub(pattern = "(", replacement = "", x, fixed = TRUE)
  x <- gsub(pattern = ")", replacement = "", x, fixed = TRUE)
  # Get rid of square brackets
  x <- gsub(pattern = "[", replacement = "", x, fixed = TRUE)
  x <- gsub(pattern = "]", replacement = "", x, fixed = TRUE)
  # Split by comma
  y <- strsplit(x, ",")
  # Average
  z <- unlist(lapply(y, function(yy) {
    mean(as.numeric(yy))
  }))
  # Return
  return(z)
}

# Example
levels(test_bins$forage_bin)
get_mean(levels(test_bins$forage_bin))

test_bins <- test_bins |> 
  mutate(forage = get_mean(forage_bin),
         temp = get_mean(temp_bin),
         predator = get_mean(predator_bin))

# Now that we have our habitat variables, we can use 'predict()' to
# calculate w(x). Remember, we can get the linear prediction (g(x)), 
# subtract the intercept, and exponentiate to get w(x).

test_bins <- test_bins |> 
  # Linear predictor
  mutate(g = predict(m1, newdata = test_bins, type = "link")) |> 
  # Subtract intercept and exp()
  mutate(w = exp(g - coef(m1)[1]))

# Done! Now we can evaluate our model by:
#   1. plotting
#   2. calculating the correlation

# Plot
ggplot(test_bins, aes(x = w, y = obs_per_cell)) +
  geom_point() +
  geom_smooth() +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw()

# Correlation
cor(test_bins$w, test_bins$obs_per_cell, method = "spearman")

# The Boyce method took us significantly more coding than the pR2 or
# AUC approaches, but hopefully you can see much clearer ties to our
# inhomogeneous Poisson point process model here.

# ... UHC plots ----

# We have implemented UHC plots in 'amt'. There is also an accompanying
# vignette demonstrating how they work with both HSFs and iSSFs.

# You can do the bootstrap resampling using 'prep_uhc()'.
?prep_uhc

# And there is a default 'plot()' method for making the actual plots.
?plot.uhc_data

# Let's see how it works with our data.
# Model 1 (correct)
uhc1 <- prep_uhc(m1, 
                 test_dat = select(test, -weight), 
                 n_samp = 500)

# Model 2 (incorrect)
uhc2 <- prep_uhc(m2,
                 test_dat = select(test, -weight), 
                 n_samp = 500)

# Plot all variables
plot(uhc1)

# Plot one at a time
plot(uhc1["predator"])

# The default plots can be slow because they draw a line (with ~512 vertices)
# for every bootstrap sample (plus the two from the observed data).

# One way to speed that up is to convert the lines to confidence
# envelopes using amt::conf_envelope().

# We must first convert the object from 'prep_uhc()' to a data.frame
# if we want to make our own plots.

df1 <- as.data.frame(uhc1)
df2 <- as.data.frame(uhc2)

# We can use these data.frames directly if we want to make our own plots:
head(df1)

df2 |> 
  filter(var %in% c("forage", "predator")) |>
  mutate(dist_sort = factor(dist, levels = c("S", "U", "A"))) |>
  ggplot(aes(x = x, y = y, color = dist_sort, linetype = dist_sort)) +
  facet_wrap(~ var, scales = "free") +
  geom_line() +
  scale_color_manual(name = "Distribution",
                     breaks = c("S", "U", "A"),
                     labels = c("Sampled", "Used", "Avail"),
                     values = c("gray70", "black", "red")) +
  scale_linetype_manual(name = "Distribution",
                        breaks = c("S", "U", "A"),
                        labels = c("Sampled", "Used", "Avail"),
                        values = c("solid", "solid", "dashed")
  ) +
  xlab("Forage (g/m2)") +
  ylab("Density") +
  theme_bw()

# Or we can convert to envelopes
# The default envelopes are the 95% and 100% envelopes
# Let's add the 90% envelope
env1 <- conf_envelope(df1, levels = c(0.9, 0.95, 1))
env2 <- conf_envelope(df2, levels = c(0.9, 0.95, 1))

# We can either plot them all
plot(env1)

# Or plot them one at a time:
env2 |>
  filter(var == "predator") |>
  plot()

env2 |>
  filter(var == "cover") |>
  plot()

env2 |>
  filter(var == "temp") |>
  plot()
