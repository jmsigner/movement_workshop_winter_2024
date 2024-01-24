# When fitting KDE home ranges:
# What's the difference between using a single bandwidth for the population
# or estimating a separate one for each individual?

library(tidyverse)
library(amt)

# Use the fisher data and calculate the reference bandwidth for each individual
dat <- amt_fisher %>% 
  nest(data = -name) %>% 
  mutate(bw = map_dbl(data, ~ hr_kde_ref(.x)[1]))

# What's the population-level mean bandwidth?
mean_bw <- mean(dat$bw)

# Fit KDEs with each bandwidth
kdes <- dat %>% 
  mutate(kde_ind = map2(data, bw, ~ hr_kde(.x, h = .y)),
         kde_pop = map(data, ~ hr_kde(.x, h = mean_bw))) %>% 
  select(name, bw, kde_ind, kde_pop) %>% 
  # Pivot longer
  pivot_longer(cols = kde_ind:kde_pop,
               names_to = "type", 
               values_to = "kde") %>% 
  # Get polygons
  mutate(poly = map(kde, hr_isopleths)) %>% 
  # Unnest
  unnest(poly)

# Plot to compare
ggplot(kdes, aes(geometry = geometry)) +
  facet_grid(name ~ type) +
  geom_sf()


