#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#-----------Module 04 -- Multiple Instances------------X
#----------------Last updated 2023-03-07---------------X
#-------------------Exercise Solution------------------X
#######################################################X

library(tidyverse)
library(amt)
library(lubridate)

## Elephants

ele <- read_csv(here::here("data/elephants.csv"))
tr <- make_track(ele, long, lat, timestamp, crs = 4326, all_cols = TRUE)
tr1 <- tr %>% mutate(year = year(t_), week = week(t_)) %>% 
  nest(data = -c(id, year, week)) %>% 
  mutate(n = map_int(data, nrow)) %>% 
  mutate(mean.temp = map_dbl(data, ~ .x %>% pull(temperature) %>% mean()))

tr2 <- filter(tr1, year == 2010)
tr3 <- tr2 %>% mutate(mcp = map(data, hr_mcp), 
               mcp.area = map_dbl(mcp, ~ hr_area(.x)$area))


ggplot(tr3, aes(week, mean.temp)) + geom_point()
ggplot(tr3, aes(week, mcp.area, col = factor(year))) + geom_point()
ggplot(tr3, aes(week, mcp.area, col = factor(year))) + geom_point() + 
  scale_y_continuous(trans = "log10")
ggplot(tr3, aes(mean.temp, mcp.area)) + geom_point() +
  scale_y_continuous(trans = "log10") +
  geom_smooth(method = "lm")

m1 <- lm(log(mcp.area) ~ mean.temp, data = tr3) 
summary(m1) # Weak evidence for an effect. 

m2 <- lm(log(mcp.area) ~ poly(mean.temp, 2), data = tr3) 
summary(m2) # Weak evidence for an effect. 

anova(m1, m2)