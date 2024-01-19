required.packages <- c(
  "amt", 
  "tidyverse", 
  "lubridate", 
  "moveHMM", 
  "CircStats",
  "sf", 
  "here", 
  "raster", 
  "terra",
  "move", 
  "crawl", 
  "here",
  "glmmTMB", 
  "ctmm", 
  "conflicted",
  "mvtnorm",
  "scales",
  "rcompanion",
  "pROC",
  "wrswoR", 
  "ResourceSelection", 
  "terra", 
  "suncalc", 
  "remotes"
)

# Suggested packages

# Used to create slides
slides <- c(
  "knitr",
  "xaringan",
  "RefManageR"
)

# Used to create figures
figs <- c(
  "ragg"
)

all.pkgs <- c(required.packages, slides, figs)

install.packages(all.pkgs[!all.pkgs %in% row.names(installed.packages())])

# Make sure you have the latest version of all packages.
update.packages(ask = FALSE)

# Add amt development version
remotes::install_github("jmsigner/amt", upgrade = "never")

