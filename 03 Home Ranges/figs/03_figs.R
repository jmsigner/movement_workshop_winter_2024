# Figures for lecture 03 -- Home Ranges
# Making home range figures with same data as Fig. 2 in Signer and Fieberg 2021.

# Load packages----
library(tidyverse)
library(amt)
library(raster)
library(ragg)
library(patchwork)
library(here)
library(adehabitatHR)
library(sp)
library(sf)
library(tlocoh)

# Range vs Occurence ----

# Source function for BCRW
source("fun/bcrw.R")

# Centroids
centroids <- list("A01" = c("x" = 446589, "y" = 4625899),
                  "A02" = c("x" = 440284, "y" = 4618197),
                  "A03" = c("x" = 448796, "y" = 4613795),
                  "A04" = c("x" = 442367, "y" = 4627042)
)

# Set random seed
set.seed(20220126)

# Generate random walk for range distribution
rng <- lapply(centroids, function(cent) {
  # The basic BCRW
  x <- bcrw(start_loc = cent + rnorm(n = 2, mean = 0, sd = 500),
            centroid = cent,
            n_steps = 1000,
            sl_distr = c('shape' = 1, 'scale' = 300),
            rho = 0.50, # stronger correlation (was 0.25)
            beta = 0.15) # stronger bias (was 0.1)
  
  # Assign actual dates and times (start on June 20, 2021)
  x$date <- as.POSIXct("2021-06-20 00:00:00") +
    x$t * 60 * 60
  
  # Return
  return(x)
}) %>% 
  bind_rows(.id = "ID") %>% 
  # Keep just individual 4
  filter(ID == "A04") %>% 
  # Make track
  make_track(x, y, date, crs = 32612)

# Data for occurrence distribution
occ <- rng[1:300,]

# Fit KDE to occurrence
occ_kde <- occ %>% 
  hr_kde()

# Fit aKDE to occurrence
occ_akde <- occ %>% 
  hr_akde(model = fit_ctmm(occ, model = "ou"))

# Get polygons
kde_poly <- hr_isopleths(occ_kde) %>% 
  mutate(what = "Est", type = "KDE") %>% 
  dplyr::select(what, type)
akde_poly <- hr_isopleths(occ_akde) %>% 
  mutate(what = c("95% CI", "Est", "95% CI"), type = "aKDE") %>% 
  dplyr::select(what, type)
polys <- rbind(kde_poly, akde_poly)

# Make plots
(occ_plot <- ggplot() +
    geom_sf(data = kde_poly, color = "gray70", linewidth = 1, fill = NA) +
    # Just to set the limits (100% transparent)
    geom_point(data = rng, aes(x = x_, y = y_), color = "#00000000") +
    geom_point(data = occ, aes(x = x_, y = y_), color = "blue", size = 1.25) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw())

(occ_plot_long <- occ_plot +
    geom_point(data = rng, aes(x = x_, y = y_), color = "orange", size = 1) +
    geom_point(data = occ, aes(x = x_, y = y_), color = "blue", size = 1.25))

(rng_plot <- ggplot() +
    geom_sf(data = polys, aes(color = type, linetype = what), 
            linewidth = 1, fill = NA) +
    geom_point(data = rng, aes(x = x_, y = y_), color = "orange", size = 1) +
    geom_point(data = occ, aes(x = x_, y = y_), color = "blue", size = 1.25) +
    xlab(NULL) +
    ylab(NULL) +
    scale_linetype_manual(name = "CI",
                          breaks = c("Est", "95% CI"), 
                          values = c("solid", "dashed")) +
    scale_color_manual(name = "Distribution",
                       breaks = c("aKDE", "KDE"),
                       labels = c("Range", "Occurrence"),
                       values = c("black", "gray70")) +
    theme_bw())

(dist_fig <- occ_plot + rng_plot)

ggsave(here::here("03 Home Ranges", "figs/occ_fig.png"), plot = occ_plot, 
       width = 800, height = 500, units = "px", dpi = 110)
ggsave(here::here("03 Home Ranges", "figs/occ_fig_long.png"), plot = occ_plot_long, 
       width = 800, height = 500, units = "px", dpi = 110)
ggsave(here::here("03 Home Ranges", "figs/rng_fig.png"), plot = rng_plot, 
       width = 1000, height = 500, units = "px", dpi = 110)
ggsave(here::here("03 Home Ranges", "figs/distr_fig.png"), plot = dist_fig, 
       width = 1000, height = 350, units = "px", dpi = 110)

# FISHER FIGURES ####
# Load data ----
dat <- amt_fisher %>% 
  filter(id == "F1")

# Fit home ranges ----
mcp <- hr_mcp(dat, levels = c(1, 0.95, 0.5))
locoh <- hr_locoh(dat, levels = seq(0.2, 1, by = 0.2),
                  type = "a", n = max(dist(dat[, c("x_", "y_")])))
kde <- hr_kde(dat, levels = c(0.95, 0.5))
akde <- hr_akde(dat, levels = c(0.95, 0.5), model = fit_ctmm(dat, model = "ou"))

# Figures ----

# ... locations ----
loc_fig <- dat %>% 
  ggplot(aes(x = x_, y = y_)) +
  geom_point(alpha = 0.3) +
  coord_sf(crs = 5070) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw()

ggsave(here::here("03 Home Ranges", "figs/locs.png"), plot = loc_fig, device = agg_png,
       width = 800, height = 750, units = "px", dpi = 150)

# ... MCP ----
mcp_fig <- ggplot() +
  geom_sf(data = mcp$mcp, aes(color = factor(level, 
                                             levels = c("1", "0.95", "0.5"),
                                             labels = c("100%", "95%", "50%"))), 
          fill = NA, linewidth = 2) +
  geom_point(data = mcp$data, aes(x = x_, y = y_),
             color = "gray30", alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(name = "MCP Level",
                     palette = "Set1") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("03 Home Ranges", "figs/mcp.png"), plot = mcp_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)

# ... LoCoH ----
locoh_fig <- ggplot() +
  geom_sf(data = locoh$locoh, aes(fill = fct_rev(factor(level)),
                                  color = fct_rev(factor(level))), 
          size = 1) +
  geom_point(data = locoh$data, aes(x = x_, y = y_),
             color = "gray", size = 0.5, alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_brewer(name = "Level", palette = "Reds") +
  scale_color_brewer(name = "Level", palette = "Reds") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("03 Home Ranges", "figs/locoh.png"), plot = locoh_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)

# ... KDE ----

# ... ... KDE UD ----
ud_fig <- kde$ud %>% 
  as.data.frame(xy = TRUE) %>% 
  ggplot(aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Probability\nDensity") +
  coord_sf(crs = 5070,
           xlim = range(dat$x_),
           ylim = range(dat$y_)) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("03 Home Ranges", "figs/ud.png"), plot = ud_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)

# ... ... KDE isopleths ----
kde_fig <- hr_isopleths(kde) %>% 
  ggplot() +
  geom_sf(aes(color = factor(level, 
                             levels = c("0.95", "0.5"),
                             labels = c("95%", "50%"))), 
          fill = NA, linewidth = 2) +
  geom_point(data = mcp$data, aes(x = x_, y = y_),
             color = "gray30", alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(name = "KDE Isopleth",
                     palette = "Set2") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("03 Home Ranges", "figs/kde.png"), plot = kde_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)

# ... aKDE ----

# ... ... plot ----
akde_fig <- hr_isopleths(akde) %>% 
  filter(what == "estimate") %>% 
  mutate(level = factor(level, 
                        levels = c("0.95", "0.5"),
                        labels = c("95%", "50%"))) %>% 
  ggplot() +
  geom_sf(aes(color = level), 
          fill = NA, linewidth = 2) +
  geom_point(data = akde$data, aes(x = x_, y = y_),
             color = "gray30", alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_brewer(name = "aKDE Isopleth",
                     palette = "Set2") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("03 Home Ranges", "figs/akde.png"), plot = akde_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)

# Other home ranges ----
# ... BBMM ----
# Inspired by figure from 'adehabitatHR' vignette
cents <- data.frame(x = c(-10, 0, 10), y = c(-10, 10, 0))

set.seed(20230312)
pts <- lapply(1:100, function(i) {
  res <- cents
  res$x <- cents$x + rnorm(3, mean = 0, sd = 1)
  res$y <- cents$y + rnorm(3, mean = 0, sd = 1)
  
  # Order
  if ((i %% 2) == 0) {
    return(res)
  } else {
    res <- res[3:1,]
    return(res)
  }
}) %>% 
  bind_rows()

pts <- pts %>% 
  mutate(x = x * 100 +  432702,
         y = y * 100 + 4621579, 
         t = as.POSIXct("2021-01-01 00:00:00") + (1:300) * 60 * 60)

# plot(pts$x, pts$y)
# lines(pts$x, pts$y)

# Format as track
trk <- make_track(pts, x, y, t, crs = 32612)
# Format as ltraj
lt <- as_ltraj(trk)

# Raster for estimation
r <- raster(extent(c(430500, 435000, 
                     4619500, 4623500)), res = 25, crs = 32612)
# SpatialPixels
spr <- as(r, "SpatialPixels")

# Ordinary kernel density
kud <- hr_kde(trk, trast = r)$ud

# Brownian bridge UD
bb <- raster(kernelbb(lt, sig1 = 7.5, sig2 = 10, grid = spr))

# Plot
kud_df <- as.data.frame(kud, xy = TRUE)
bb_df <- as.data.frame(bb, xy = TRUE)

kud_panel <- ggplot(kud_df, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Ordinary KDE") +
  coord_sf(crs = 32612,
           xlim = range(kud_df$x),
           ylim = range(kud_df$y)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank())

bb_panel <- ggplot(bb_df, aes(x = x, y = y, fill = ud)) +
  geom_raster() +
  scale_fill_viridis_c() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brownian Bridge KDE") +
  coord_sf(crs = 32612,
           xlim = range(kud_df$x),
           ylim = range(kud_df$y)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank())

loc_panel <- trk %>% 
  steps() %>% 
  ggplot(aes(x = x1_, xend = x2_, y = y1_, yend = y2_)) +
  geom_segment(color = "gray70") +
  geom_point() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Track") +
  coord_sf(crs = 32612) +
  theme_bw() +
  theme(legend.position = "none")

bb_fig <- (kud_panel + bb_panel)

ggsave(here::here("03 Home Ranges", "figs/bb_locs.png"), plot = loc_panel, 
       device = agg_png, width = 800, height = 400, units = "px", dpi = 100)
ggsave(here::here("03 Home Ranges", "figs/bb_ex.png"), plot = bb_fig, 
       device = agg_png, width = 800, height = 400, units = "px", dpi = 100)

# ... ... fisher BBMM ----
dat2 <- as_ltraj(dat)
# Determine parameters
liker(dat2, rangesig1 = c(0.5, 5), sig2 = 25)
# Fit HR
bb_fisher <- kernelbb(dat2, sig1 = 3, sig2 = 25, grid = 200, extent = 0.05)
plot(raster(bb_fisher))
# Get polygons
bb_95 <- st_as_sf(getverticeshr(bb_fisher, percent = 95)) %>% 
  mutate(level = "0.95")
bb_50 <- st_as_sf(getverticeshr(bb_fisher, percent = 50)) %>% 
  mutate(level = "0.5")
bb_poly <- rbind(bb_95, bb_50)
st_crs(bb_poly) <- st_crs(5070)

# Figure
bb_ud <- ggplot() + 
  geom_sf(data = bb_poly, aes(color = factor(level, 
                             levels = c("0.95", "0.5"),
                             labels = c("95%", "50%"))), 
          fill = NA, linewidth = 2) + 
  coord_sf(default_crs = 5070,
           crs = 5070,
           xlim = range(dat$x_),
           ylim = range(dat$y_)) +
  scale_fill_viridis_c(name = "Probability Density") +
  scale_color_brewer(name = "BBMM Isopleth",
                     palette = "Set2") +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", linewidth = 0.8))

ggsave(here::here("03 Home Ranges", "figs/bb.png"), plot = bb_ud, 
       device = agg_png, width = 800, height = 800, units = "px", dpi = 150)

# Figure
bb_ud2 <- as.data.frame(raster(bb_fisher), xy = TRUE) %>% 
  ggplot() + 
  geom_raster(aes(x = x, y = y, fill = ud), show.legend = FALSE)  +
  geom_sf(data = bb_poly, aes(color = factor(level, 
                                             levels = c("0.95", "0.5"),
                                             labels = c("95%", "50%"))), 
          fill = NA, linewidth = 2) + 
  coord_sf(default_crs = 5070,
           crs = 5070,
           xlim = range(dat$x_),
           ylim = range(dat$y_)) +
  scale_fill_viridis_c(name = "Probability Density") +
  scale_color_brewer(name = "BBMM Isopleth",
                     palette = "Set2") +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", linewidth = 0.8))

ggsave(here::here("03 Home Ranges", "figs/bb2.png"), plot = bb_ud2, 
       device = agg_png, width = 800, height = 800, units = "px", dpi = 150)

# ... T-LoCoH ----
# Convert data to 'move' object
dat_mov <- as_move(dat)

# Format for tlocoh
dat_lxy <- move.lxy(dat_mov)

# Decide on parameter s
# (Read the y-axis)
# Say we're interested in locations close in space on the order of a few hours
lxy.plot.sfinder(dat_lxy, delta.t = 60 * c(30, 60, 120, 180, 240, 480))
# Looks like 0.02 is good if we're interested in 3, 4, or 8h relocations
abline(h = 0.02, col = "red", lty = 2)

# Add nearest neighbors
nn <- lxy.nn.add(dat_lxy, k = 30, s = 0.02)
# Create hull set
hs <- lxy.lhs(nn, s = 0.02, k = 30, iso.add = TRUE)

# Plot
plot(hs, iso = TRUE)

# Get the polygons
polys <- hs$F1.pts1349.k30.s0.02.kmin0$isos$`iso.srt-area.iso-q.h1349.i5`$polys %>% 
  st_as_sf()

# Make the plot

tlocoh_fig <- ggplot() +
  geom_sf(data = polys, aes(fill = fct_rev(factor(iso.level)),
                                  color = fct_rev(factor(iso.level))), 
          size = 1) +
  geom_point(data = dat, aes(x = x_, y = y_),
             color = "gray", size = 0.5, alpha = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_brewer(name = "Level", palette = "Reds") +
  scale_color_brewer(name = "Level", palette = "Reds") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.box.background = element_rect(colour = "black", size = 0.8))

ggsave(here::here("03 Home Ranges", "figs/tlocoh.png"), plot = tlocoh_fig, device = agg_png,
       width = 800, height = 800, units = "px", dpi = 150)
