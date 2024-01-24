#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#---------- Module 04 -- Multiple Instances -----------X
#----------------Last updated 2024-01-18---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(amt)
library(lubridate) # to deal with date and time

# Question ----

# We will use the `amt_fisher` data set. And calculate for each fisher the daily
# home-range size using two different estimators (MCP and KDE).

# Data preparation ----

# ... Loading the data ----
# Let's first load the data (the data is shipped with the amt package). 

data("amt_fisher")
dat <- amt_fisher
head(dat)

# We have four individuals
unique(dat$name)

# ... Obtain additional variables ---- In order to calculate the daily
# home-range size, we need an other covariate (the day when a relocation was
# observed). We can use the function `floor_date(unit = "day")` from the
# `lubridate` package. This will floor (i.e. the times stamp is rounded down to
# the day).

# Here is an example
floor_date(ymd_hms("2022-01-20 05:12:10"), unit = "day")
floor_date(ymd_hms("2022-01-20 20:12:10"), unit = "day")

round_date(ymd_hms("2022-01-20 05:12:10"), unit = "day")
round_date(ymd_hms("2022-01-20 20:12:10"), unit = "day")

ceiling_date(ymd_hms("2022-01-20 05:12:10"), unit = "day")
ceiling_date(ymd_hms("2022-01-20 20:12:10"), unit = "day")

# Note, the argument `unit` could also be set to "hour", "week" and other
# commonly used units. 
floor_date(ymd_hms("2022-01-20 05:12:10"), unit = "hour")
floor_date(ymd_hms("2022-01-20 20:12:10"), unit = "hour")

# Also multiples are possible. 
floor_date(ymd_hms("2022-01-20 05:12:10"), unit = "3 hours")
floor_date(ymd_hms("2022-01-20 20:12:10"), unit = "3 hours")

# But see `?floor_date` for more information.

# Here we add a new column called `day` to our data set.
dat <- dat |> mutate(day = floor_date(t_, unit = "day"))
dat

# To get a brief overview of the data, it makes sense to count how many
# observation we have per animal
dat |> count(name)

# and per animal and day.
dat |> count(name, day) |> print(n = Inf)

# We have 114 instances. I use the word instance here to refer to one unit of
# interest. Here a unit of interest are all relocations of on animal at a given
# day.

# Note, for a proper analysis we would have to get the sampling right corrected
# here, but we will skip this for this illustration.

# Approaches ----

# ... Naive approach ----

# The "naive" approach is characterized by calculating the a home-range size for
# each animal and day by subsetting. The first instance is "Leroy" on the 11th
# of February 2009.
leroy.1 <- filter(dat, day == ymd("2009-02-11"), name == "Leroy")

# Now we can calculate the home-range size
hr.1 <- leroy.1 |> hr_mcp() |> hr_area()

# We would we will now have to save the results in a new data frame with a
# column for the `name` and the `day`.
res <- distinct(dat, name, day)
res$area <- NA

# Instead of manually typing the `name` and `day`, we can retrieve it from the
# `res` data.frame.
filter(dat, day == res$day[1], name == res$name[1])

# With this data we can calculate again a home range, extract its size, and save
# it in the results data.frame.
res$area[1] <- filter(dat, day == res$day[1], name == res$name[1]) |> 
  hr_mcp() |> hr_area() |> pull(area)

# We now have to repeat this for every instance
res$area[1] <- filter(dat, day == res$day[1], name == res$name[1]) |> 
  hr_mcp() |> hr_area() |> pull(area)
res$area[2] <- filter(dat, day == res$day[2], name == res$name[2]) |> 
  hr_mcp() |> hr_area() |> pull(area)

res
# And so on ...

# ... Loop approach ----
# If we investigate the the code from above again, we could substitute the row
# number with a variable, lets say `i`.
i <- 1
res$area[i] <- filter(dat, day == res$day[i], name == res$name[i]) |> 
  hr_mcp() |> hr_area() |> pull(area)

# For the second animal, we simply have to change the value of `i` to 2
i <- 2
res$area[i] <- filter(dat, day == res$day[i], name == res$name[i]) |> 
  hr_mcp() |> hr_area() |> pull(area)

# And we could continue like this. This is very tedious and does not save us a lot
# of time. Fortunately, it easy to automatize this process using `for`-loops.

# For loops work like this: for each iteration (each time the loop is called) a
# variable (for example `i`) changes it's value.
for (i in 1:3) {
  print(i)
}

# Instead of printing the value of `i`, we can actually use a loop to calculate
# the home-range size for different animals.
for (i in 1:10) {
  res$area[i] <- filter(dat, day == res$day[i], name == res$name[i]) |> 
    hr_mcp() |> hr_area() |> pull(area)
}
res |> print(n = 12)

# In order to iterate over all individuals, we can replace `10` with the number
# of instances (= `nrow(res)`).
for (i in 1:nrow(res)) {
  res$area[i] <- filter(dat, day == res$day[i], name == res$name[i]) |> 
    hr_mcp() |> hr_area() |> pull(area)
}
res |> print(n = Inf)

# ... Nesting approach ---- 
# A third approach is to use list columns. This slightly more challenging, but a
# lot more elegant.

# The idea is as follows: We take all relocations and group them by one or more
# variables (this is equivalent to the `filter` part in the for loop).

# With the tidyverse, we can use the function `nest()` from the `purrr` package. 
dat1 <- nest(dat, my.fisher.data = -name)
dat1 <- nest(dat, data = -name)
dat1

# `data` is a new list column (i.e., a list) with the relocation for each
# animal.
dat1$data

# Or all relocations for Leroy
dat1
dat1$data[[1]]

# There are two ways to tell `nest()` the grouping variables. Either by
# explicitly naming all variables that should be used to group (using the `-`).
# Or by naming all variables that should be grouped.

# We have to redo the nesting, because we have two grouping variables (the
# individual and the day).

# Using the first approach (naming all grouping variables) looks like this: 
dat |> nest(data = -c(name, day))

# With the second approach (naming all variables that should be grouped and the
# remaining variables are then used as grouping variables) looks like this:
dat |> nest(data = x_:id)

# The result is the same. I usually prefer the first approach, because it is
# more explicit. But that's a matter of taste.
dat1 <- dat |> nest(data = -c(name, day))

# Note, I always used `data = ....` to name the new list column `data`. But of a
# different name could also be used.
dat |> nest(daily.gps.points = -c(name, day))

# We now want to apply the MCP home-range estimator to each instance. We could
# use a loop over this column, but it is much easier to use an apply-like
# approach here. The base R `lapply()` function could be used, or the `map()`
# function from the `tidyverse`.
map(dat1$data[1:5], hr_mcp) 

# Note, I used `dat1$data[1:5]` to only work with the first 5 instances. 

# This returns yet another list. Each entry of this list a MCP home-range
# estimate. Staying within the tidyverse, we can create a new list column with
# `mutate()`.
dat1 <- dat1 |> mutate(hr.mcp = map(data, hr_mcp))
dat1

# In the column `hr.mcp` the MCP home-range estimates for each instance are
# saved. We can now extract attributes from this estimates. We are interested in
# the area. Thus we call the function `hr_area()` on each instance. We could use
# the `map()` function again, but this would return a list, but we would like to
# obtain a single number. Here the function `map_dbl()` is the choice.
# `map_dbl()` work exactly as `map()` with the key difference that `map()`
# returns a list and `map_dbl()` a vector with double values (i.e.,
# decimal numbers).
map(dat1$hr.mcp[1:5], ~ hr_area(.x)$area) # This is a list

# We don't want this.
map_dbl(dat1$hr.mcp[1:5], ~ hr_area(.x)$area) # This is a vector

# And we calculate the area for all animals
dat1 <- dat1 |> mutate(area.mcp = map_dbl(hr.mcp, ~ hr_area(.x)$area))
dat1

dat1 |> select(-data, -hr.mcp)


# Do the number of relocations matter? ----
# We could now ask, if the number of relocations influence the home-range size?
# Using the current setup this is easy to answer. We just need one new column
# `n` (i.e., the number of data points per instance). We can use again a
# `map`-approach. But since the number of relocation will always be an integer
# (i.e., a number without decimals), we should use the function `map_int()`.
dat1 <- dat1 |> mutate(n = map_int(data, nrow))
dat1

# Now we can calculate a correlation coefficient
with(dat1, cor(n, area.mcp))

# or just plot it.
ggplot(dat1, aes(n, area.mcp)) + geom_point() + geom_smooth()

# There seems to be an other pattern, it is likely that the individual is more
# important.
ggplot(dat1, aes(n, area.mcp, col = name)) + geom_point() + geom_smooth()

# Let's try to calculate some numbers to support this pattern. 

# We start with a correlation coefficient again, but this time per individual.
dat1 |> group_by(name) |> 
  summarize(cor = cor(n, area.mcp))

# Or we can use a simple linear model here:
lm(area.mcp ~ n + name, data = dat1) |> summary()

# Subsampling would clearly be a good idea. Let's rerun everything with
# subsampling. We can just add one step to into our tidy workflow.
dat2 <- dat |> nest(data = -c(name, day)) |> 
  mutate(
    # Subsample to 1 relocation per hour
    data = map(data, ~ track_resample(.x, rate = hours(1), tolerance = minutes(5))), 
    # Calculate MCP
    hr.mcp = map(data, hr_mcp),
    # Calculate the area
    area.mcp = map_dbl(hr.mcp, ~ hr_area(.x)$area),
    # n
    n = map_int(data, nrow)
  )
dat2


# Different sampling rates for different animals
dat2a <- dat |> nest(data = -name) 
dat2a$sampling_rate <- 1:4

dat2a <- mutate(dat2a, data1 = map2(data, sampling_rate, ~ 
                                      track_resample(.x, rate = hours(.y), tolerance = minutes(5))))

# And we can redo the same summary as before
ggplot(dat2, aes(n, area.mcp, col = name)) + geom_point() + geom_smooth()

# And the linear model
lm(area.mcp ~ n + name, data = dat2) |> summary()

# Maybe we should only consider days with at least 10 relocations? This is easy
# to achieve, since we already have a column called `n`, so we just need to
# filter by this column. 

ggplot(dat2 |> filter(n >= 10), aes(n, area.mcp, col = name)) + 
  geom_point() + geom_smooth()

# And the linear model
lm(area.mcp ~ n + name, data = filter(dat2, n >= 10)) |> summary()

# So we can conclude that sample size does not matter, since the coefficient for
# `n` is not significantly different from zero.

# Does estimator choice matter?   ---- 
# Lastly, we can check how the estimated
# daily space use differs between two different estimators.
dat2 <- dat |> nest(data = -c(name, day)) |> 
  mutate(
    # Subsample to 1 relocation per hour
    data = map(data, ~ track_resample(.x, rate = hours(1), tolerance = minutes(5))), 
    n = map_int(data, nrow)) |> 
  filter(n >= 10) |> 
  mutate(
    # Calculate MCP
    hr.mcp = map(data, hr_mcp),
    # Calculate KDE
    hr.kde = map(data, hr_kde)
  )

dat2

# Since all estimators implemented in `amt` have an area method implemented, we
# can simplify the workflow by changing from the wide format to the long format.
dat3 <- dat2 |> 
  pivot_longer(hr.mcp:hr.kde, names_to = "estimator", values_to = "hr") |> 
  # Now it is easy to extract the home-range size
  mutate(hrs = map_dbl(hr, ~ hr_area(.x)$area))

dat3

# Again we can inspect the result visually
dat3 |> ggplot(aes(name, hrs, fill = estimator)) + geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  theme_light()

# Overlap between MCP and KDE ----

dat2 <- dat2[1:10, ] |> mutate(overlap = map2(hr.mcp, hr.kde, ~ hr_overlap(.x, .y)))
dat2
dat2$overlap[[2]]
dat2 |> unnest(overlap)

‚

