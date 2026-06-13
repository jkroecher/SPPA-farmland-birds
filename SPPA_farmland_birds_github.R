##########################
# Author:      Jenny Kroecher
# Description: Data pre-processing – filtering, spatial analysis, and
#              covariate modelling of bird observations (breeding season 2021)
# Date:        2026-06-13
##########################


# Clear workspace
rm(list = ls())

# Set working directory (adjust path if needed)
setwd("D:/BirdMon")





# --- Load packages -----------------------------------------------------------
# Spatial data (vector & raster)
library(sp)          # SpatialPoints, SpatialPolygons
library(sf)          # Simple Features (modern vector data handling)
library(terra)       # Raster processing (successor to raster)
library(raster)      # Raster (still needed for as.im.RasterLayer)

# Point pattern analysis
library(spatstat)    # Kernel density, envelopes, PPM models

# Data management & import
library(openxlsx)    # Read/write Excel files
library(lubridate)   # Date/time handling
library(dplyr)       # Data transformation (arrange, etc.)

# Visualisation
library(ggplot2)
library(RColorBrewer)
library(viridis)

# ============================================================================
# USER SETTINGS – adjust these file paths before running the script
# ============================================================================

# Bird observation data 
file_bird_obs     <- "path/to/your/bird_observations.xlsx"

# Flower strip vertices 
file_flstrip      <- "path/to/your/flstrip_vertices.xlsx"

# Vegetation cover and weed cover 
file_veg_cover    <- "path/to/your/vegetation_cover_weed_cover.xlsx"

# Vegetation height 
file_veg_height   <- "path/to/your/vegetation_height.xlsx"

# Weed control activity 
file_weed_control <- "path/to/your/weed_control.xlsx"

# Pesticide activity 
file_psm          <- "path/to/your/psm_control.xlsx"

# Land use intensity and crop types 
file_land_use     <- "path/to/your/land_use_intensity_crop_types.xlsx"

# spatial geometry of observation windows
file_obs_window   <- "path/to/your/observation_window.shp"

# spatial geometry of patch polygons
file_patches      <- "path/to/your/patches.shp"


# ============================================================================
# 1. LOAD BIRD OBSERVATION DATA
# ============================================================================

# Load raw bird observation data
# Path: replace with a relative path once the project folder is finalised
#       e.g. "data/1_Bird_observations.xlsx"
df <- read.xlsx(file_bird_obs)

# Excel stores dates as integers (days since 1899-12-30); convert to Date
df$DATE <- as.Date(df$DATE_TIME, origin = "1899-12-30")


# ============================================================================
# 2. FILTER: KEEP ONLY RELEVANT BEHAVIOURAL OBSERVATIONS
# ============================================================================

# Retain only records that can be spatially assigned to plots:
# SINGING, PRESENT, FEEDING
df_a <- df[df$OBSERVATION_DETAIL %in% c("SINGING", "PRESENT", "FEEDING"), ]


# ============================================================================
# 3. FILTER: KEEP ONLY OBSERVATIONS INSIDE THE MONITORING WINDOWS
# ============================================================================

# Load observation windows (100 m buffers around transect routes) as polygons
obs_win <- st_read(file_obs_window,
                   "obs_win")
obs_win <- as(obs_win, "Spatial")   # Convert to sp format for downstream use

# Ensure coordinates are numeric
df_a$COORD_LON <- as.numeric(df_a$COORD_LON)
df_a$COORD_LAT <- as.numeric(df_a$COORD_LAT)

# Create a SpatialPointsDataFrame in WGS84
xy       <- df_a[, c("COORD_LON", "COORD_LAT")]
df_a_sp <- SpatialPointsDataFrame(
  coords      = xy,
  data        = df_a,
  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
)

# Project observation windows to WGS84 to match the point CRS
obs_win_t <- spTransform(
  obs_win,
  CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
)

# Visual check of the observation window geometry
plot(obs_win_t)

# Remove points that fall outside the observation windows
df_sp_sub <- df_a_sp[obs_win_t, ]

# Reproject points back to the original CRS of beob_win
df_sp_sub_t <- spTransform(df_sp_sub, CRS(obs_win@proj4string@projargs))

# Convert result to a plain data frame
df_obs <- data.frame(df_sp_sub_t)


# ============================================================================
# 4. CREATE OBSERVATION WINDOW OBJECT FOR SPATSTAT (owin)
# ============================================================================

# Extract polygons without attributes
obs_spol <- as(obs_win, "SpatialPolygons")

# Convert to spatstat observation window (owin)
# The individual polygons correspond to the four field areas;
# points are later assigned to their respective polygon by spatial location.
mon_win <- as(obs_spol, "owin")


# ============================================================================
# 5. CREATE POINT PATTERN OBJECT (all years, all field areas)
# ============================================================================

# Columns 17 and 18 hold the projected X and Y coordinates
pp_obs <- ppp(df_obs[, 17], df_obs[, 18], window = mon_win)

# Verify class
class(pp_obs)

# Set coordinate unit
unitname(pp_obs) <- c("meter", "meter")

# Plot full point pattern
par(mfrow = c(1, 1))
plot(pp_obs)


# ============================================================================
# 6. KERNEL DENSITY PLOTS (per field area, three bandwidth values)
# ============================================================================
# Field areas are processed in the order 3, 2, 4, 1 to match the intended
# plot layout. For each area three sigma values are tested:
#   - 100 m  : fixed bandwidth
#   - bw.ppl(): data-driven bandwidth (likelihood cross-validation)
#   - mean(nndist()): mean nearest-neighbour distance

par(mfrow = c(4, 3))

for (i in c(3, 2, 4, 1)) {
  
  # Subset polygon and observation window for field area i
  obs_spol_sub_fl <- obs_spol[i]
  obs_win_fl      <- as(obs_spol_sub_fl, "owin")
  
  # Restrict point pattern to field area i
  win_in     <- inside.owin(pp_obs, w = obs_win_fl)
  pp_obs_sub <- pp_obs[win_in]
  Window(pp_obs_sub) <- obs_win_fl
  
  # Remove duplicate locations
  pp_obs_sub <- unique(pp_obs_sub)
  unitname(pp_obs_sub) <- c("meter", "meter")
  
  # Bandwidth values and corresponding colour scale limits
  ranges <- data.frame(
    j   = c(100, bw.ppl(pp_obs), mean(nndist(pp_obs))),
    min = c(0, 0, 0),
    max = c(35, 60, 120)
  )
  
  for (j in ranges$j) {
    
    # Adjust plot margins per field area (different shapes/sizes)
    if (i == 3) par(mar = c(1, 1, 1, 3))
    if (i == 2) par(mar = c(5.0, 5.0, 5.0, 5.5))
    if (i == 4) par(mar = c(1, 1, 1, 3))
    if (i == 1) par(mar = c(2.0, 2.0, 2.0, 3.0))
    
    # Multiply density by 10 000 to convert from obs/m² to obs/ha
    plot(
      density(pp_obs_sub, sigma = j, kernel = "gaussian",
              edge = TRUE, diggle = (i == 3)) * 10000,
      col      = colourmap(viridis(256),
                           range = c(ranges$min[ranges$j == j],
                                     ranges$max[ranges$j == j])),
      main     = "",
      ribbon   = TRUE,
      ribscale = 1,
      box      = FALSE,
      riblab   = list("Density [observation/ha]",
                      line = if (i == 2) -0.5 else if (i == 1) -1 else 2)
    )
  }
}


# ============================================================================
# 7. SPATIAL POINT PATTERN ANALYSIS – SECOND-ORDER ANALYSIS: CSR NULL MODEL
#    (Complete Spatial Randomness)
# ============================================================================
# Tests whether the observed point pattern deviates from CSR using:
#   J(r): nearest-neighbour / empty-space function (Jest)
#   L(r): Ripley's L function (Lest)
# Results are visualised as colour-coded tiles per field area:
#   Green      = Aggregation  (points cluster together)
#   Grey       = Random       (no significant deviation)
#   Yellow-green = Segregation (points repel each other)
#
# Sign convention:
#   J function: obs < envelope → aggregation; obs > envelope → segregation
#   L function: obs > envelope → aggregation; obs < envelope → segregation

yaxis <- c("J(r)", "L(r)")

par(mfrow = c(4, 3))

for (i in c(3, 2, 4, 1)) {
  
  # Subset polygon, window and points for field area i
  obs_spol_sub_fl <- obs_spol[i]
  obs_win_fl      <- as(obs_spol_sub_fl, "owin")
  win_in           <- inside.owin(pp_obs, w = obs_win_fl)
  pp_obs_sub       <- pp_obs[win_in]
  Window(pp_obs_sub) <- obs_win_fl
  pp_obs_sub       <- unique(pp_obs_sub)
  unitname(pp_obs_sub) <- c("meter", "meter")
  
  # J function: simulation envelope under CSR (199 simulations, rank 5)
  env_j <- envelope(pp_obs_sub, fun = Jest,
                    funargs     = list(correction = "rs"),
                    alternative = "two.sided",
                    nrank = 5, nsim = 199, verbose = FALSE)
  
  # L function: global simulation envelope under CSR (199 simulations, rank 10)
  env_l <- envelope(pp_obs_sub, fun = Lest,
                    funargs     = list(correction = "border"),
                    alternative = "two.sided",
                    nrank = 10, nsim = 199, verbose = FALSE, global = TRUE)
  
  # Classify each distance r: is the observed value outside the envelope?
  env_j$dif_lo <- env_j$obs - env_j$lo
  env_j$dif_hi <- env_j$obs - env_j$hi
  env_j$result <- ifelse(env_j$dif_hi < 0 & env_j$dif_lo < 0, "aggregation",
                         ifelse(env_j$dif_hi > 0 & env_j$dif_lo > 0, "segregation",
                                "regular"))
  env_j <- env_j[!is.na(env_j$result), ]
  
  env_l$dif_lo <- env_l$obs - env_l$lo
  env_l$dif_hi <- env_l$obs - env_l$hi
  env_l$result <- ifelse(env_l$dif_hi > 0 & env_l$dif_lo > 0, "aggregation",
                         ifelse(env_l$dif_hi < 0 & env_l$dif_lo < 0, "segregation",
                                "regular"))
  
  # Tile plot: J(r) at y = 0.5, L(r) at y = 1.0
  ggp <- ggplot() +
    theme_light() +
    theme(
      axis.title.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90"),
      legend.title       = element_blank(),
      panel.border       = element_blank(),
      axis.text          = element_text(size = 14),
      axis.title.x       = element_text(size = 14),
      plot.margin        = unit(c(0, 0.4, 0, 0), "cm"),
      legend.position    = "none",
      axis.line          = element_line(colour = "grey75")
    ) +
    geom_tile(aes(x      = env_j$r,
                  y      = rep(0.5, length(env_j$r)),
                  width  = rep(0.5, length(env_j$r)),
                  height = rep(0.3, length(env_j$r)),
                  fill   = env_j$result)) +
    geom_tile(aes(x      = env_l$r,
                  y      = rep(1.0, length(env_l$r)),
                  width  = rep(0.5, length(env_l$r)),
                  height = rep(0.3, length(env_l$r)),
                  fill   = env_l$result)) +
    scale_fill_manual(
      values = c("aggregation" = "#3AAE6C", "regular" = "grey95", "segregation" = "#9DAE21"),
      labels = c("aggregation" = "Aggregation", "regular" = "Regular", "segregation" = "Segregation")
    ) +
    scale_y_continuous(breaks = c(0.5, 1), labels = yaxis) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
    xlab("Distance r")
  
  print(ggp)
}


# ============================================================================
# 8. SPATIAL POINT PATTERN ANALYSIS – SECOND-ORDER ANALYSIS: 
#    INHOMOGENEOUS POISSON NULL MODEL (HP)
# ============================================================================
# Same structure as Section 7, but using inhomogeneous versions of the
# summary functions (Jinhom / Linhom) to account for spatially varying
# intensity. Intensity is estimated via kernel density (bandwidth = bw.ppl).

yaxis <- c("J(r)", "L(r)")

par(mfrow = c(4, 3))

for (i in c(3, 2, 4, 1)) {
  
  obs_spol_sub_fl <- obs_spol[i]
  obs_win_fl      <- as(obs_spol_sub_fl, "owin")
  win_in           <- inside.owin(pp_obs, w = obs_win_fl)
  pp_obs_sub       <- pp_obs[win_in]
  Window(pp_obs_sub) <- obs_win_fl
  pp_obs_sub       <- unique(pp_obs_sub)
  unitname(pp_obs_sub) <- c("meter", "meter")
  
  # Estimate intensity at observation points (required by Jinhom)
  kern_dens <- density(pp_obs_sub, sigma = bw.ppl(pp_obs),
                       kernel = "gaussian", at = "points")
  
  # Inhomogeneous J function (Jinhom)
  r_sub_j <- seq(0, 50, by = 0.1)
  env_j   <- envelope(pp_obs_sub, fun = Jinhom, lambda = kern_dens,
                      funargs     = list(correction = "rs"),
                      r           = r_sub_j,
                      alternative = "two.sided",
                      nrank = 5, nsim = 199, verbose = FALSE, fix.n = TRUE)
  
  # Inhomogeneous L function (Linhom)
  r_sub_l <- seq(0, 100, by = 0.5)
  env_l   <- envelope(pp_obs_sub, fun = Linhom,
                      funargs     = list(correction = "border"),
                      r           = r_sub_l,
                      alternative = "two.sided",
                      nrank = 5, nsim = 199, verbose = FALSE,
                      global = TRUE, fix.n = TRUE)
  
  # Classification (same sign convention as Section 7)
  env_j$dif_lo <- env_j$obs - env_j$lo
  env_j$dif_hi <- env_j$obs - env_j$hi
  env_j$result <- ifelse(env_j$dif_hi < 0 & env_j$dif_lo < 0, "aggregation",
                         ifelse(env_j$dif_hi > 0 & env_j$dif_lo > 0, "segregation",
                                "regular"))
  
  env_l$dif_lo <- env_l$obs - env_l$lo
  env_l$dif_hi <- env_l$obs - env_l$hi
  env_l$result <- ifelse(env_l$dif_hi > 0 & env_l$dif_lo > 0, "aggregation",
                         ifelse(env_l$dif_hi < 0 & env_l$dif_lo < 0, "segregation",
                                "regular"))
  
  ggp <- ggplot() +
    theme_light() +
    theme(
      axis.title.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90"),
      legend.title       = element_blank(),
      panel.border       = element_blank(),
      axis.text          = element_text(size = 14),
      axis.title.x       = element_text(size = 14),
      plot.margin        = unit(c(0, 0.4, 0, 0), "cm"),
      legend.position    = "none",
      axis.line          = element_line(colour = "grey75")
    ) +
    geom_tile(aes(x      = env_j$r,
                  y      = rep(0.5, length(env_j$r)),
                  width  = rep(0.5, length(env_j$r)),
                  height = rep(0.3, length(env_j$r)),
                  fill   = env_j$result)) +
    geom_tile(aes(x      = env_l$r,
                  y      = rep(1.0, length(env_l$r)),
                  width  = rep(0.5, length(env_l$r)),
                  height = rep(0.3, length(env_l$r)),
                  fill   = env_l$result)) +
    scale_fill_manual(
      values = c("aggregation" = "#3AAE6C", "regular" = "grey95", "segregation" = "#9DAE21"),
      labels = c("aggregation" = "Aggregation", "regular" = "Regular", "segregation" = "Segregation")
    ) +
    scale_y_continuous(breaks = c(0.5, 1), labels = yaxis) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
    xlab("Distance r")
  
  print(ggp)
}


# ============================================================================
# 9. COVARIATE MODELS (breeding season 2021, PatchCROP main field)
# ============================================================================

# --- 9.1 Restrict observations to breeding season 2021 ---------------------

df_obs21 <- df_obs[df_obs$DATE >= as.Date("2021-04-01") &
                     df_obs$DATE <= as.Date("2021-07-31"), ]

pp_obs21 <- ppp(df_obs21[, 17], df_obs21[, 18], window = mon_win)
class(pp_obs21)

par(mfrow = c(1, 1))
plot(pp_obs21)


# --- 9.2 Restrict to PatchCROP main field (polygon 3) ----------------------
obs_win_fl      <- as(obs_spol[c(3)], "owin")
win_in           <- inside.owin(pp_obs21, w = obs_win_fl)
pp_obs_sub       <- pp_obs21[win_in]
Window(pp_obs_sub) <- obs_win_fl

# Species composition as percentage of total observations
table(pp_obs_sub$marks) / length(pp_obs_sub$marks) * 100

plot(pp_obs_sub)


# --- 9.3 Prepare patch geometry ---------------------------------------------

# Load polygon boundaries for the 30 crop patches and flower strips
patch <- st_read(file_patches,
                 "patches")

# Reconstruct flower strip polygons from vertex table (Excel)
# Path: replace with a relative path once the project folder is finalised
fs <- read.xlsx(file_flstrip)
fs <- fs %>% arrange(ID, VERTEX_INDEX)

poly_list <- lapply(split(fs, fs$ID), function(x) {
  coords <- as.matrix(x[, c("X_W84", "Y_W84")])
  # Close the polygon if the first and last vertex differ
  if (!all(coords[1, ] == coords[nrow(coords), ])) {
    coords <- rbind(coords, coords[1, ])
  }
  st_polygon(list(coords))
})

fs_poly <- st_sf(
  ID       = as.numeric(names(poly_list)),
  geometry = st_sfc(poly_list, crs = 4326)
)
fs_poly <- st_transform(fs_poly, 25833)
plot(st_geometry(fs_poly))

# Combine patches and flower strips; add residual area (PatchID = 400)
beob_sf <- st_as_sf(obs_spol[c(3)])
patch   <- st_transform(patch, st_crs(beob_sf))

patch2          <- patch[, "PatchID"]
fs_poly$PatchID <- fs_poly$ID
fs_poly2        <- fs_poly[, "PatchID"]
patch_fs        <- rbind(patch2, fs_poly2)

# Area inside the monitoring window but outside all defined patches
rest   <- st_difference(st_union(beob_sf), st_union(patch_fs))
rest   <- st_sf(PatchID = 400, geometry = rest)
result <- rbind(patch_fs, rest)

plot(result)


# --- 9.4 Null model (Poisson, no trend, no interaction) --------------------

m_pois <- ppm(pp_obs_sub ~ 1, interaction = NULL)
print(m_pois)
stats::BIC(m_pois)


# ============================================================================
# 10. PREPARE COVARIATES
# ============================================================================
# For each covariate, the patch-level mean (or count) is calculated,
# joined to the patch geometry, and rasterised at 5 m resolution (EPSG:25833).
# The raster is then converted to a spatstat image (im) for use in PPM models.


# --- 10.1 Total vegetation cover -------------------------------------------
# Path: replace with a relative path once the project folder is finalised
veg <- read.xlsx(file_veg_cover)
veg$DATE <- as.Date(veg$DATE, origin = "1899-12-30")
veg_21   <- veg[veg$DATE >= as.Date("2021-05-01") & veg$DATE <= as.Date("2021-08-31"), ]

# Mean seasonal deviation from the daily field-wide mean per patch
veg_21_a <- aggregate(TOTAL_VEGETATION_COVER ~ PATCH_ID + DATE, data = veg_21, mean)
for (i in seq_len(nrow(veg_21_a))) {
  veg_21_a$dif_mean[i] <- veg_21_a$TOTAL_VEGETATION_COVER[i] -
    mean(veg_21_a$TOTAL_VEGETATION_COVER[veg_21_a$DATE == veg_21_a$DATE[i]])
}
veg_21_dif <- aggregate(dif_mean ~ PATCH_ID, data = veg_21_a, mean)
result_veg <- merge(result, veg_21_dif, by.x = "PatchID", by.y = "PATCH_ID", all.x = TRUE)

r        <- rast(ext(result_veg), resolution = 5, crs = "epsg:25833")
veg_cover    <- raster(rasterize(vect(result_veg), r, field = "dif_mean"))
veg_cover_im <- as.im.RasterLayer(veg_cover)

plot(veg_cover_im, col = colorRampPalette(brewer.pal(9, "Greens"))(50))
plot(pp_obs_sub, add = TRUE, pch = 16, cex = 0.6, col = "grey20")


# --- 10.2 Weed cover -------------------------------------------------------
veg_21_w <- aggregate(WEED_COVER ~ PATCH_ID + DATE, data = veg_21, mean)
for (i in seq_len(nrow(veg_21_w))) {
  veg_21_w$dif_mean[i] <- veg_21_w$WEED_COVER[i] -
    mean(veg_21_w$WEED_COVER[veg_21_w$DATE == veg_21_w$DATE[i]])
}
veg_21_w_dif <- aggregate(dif_mean ~ PATCH_ID, data = veg_21_w, mean)
result_weed  <- merge(result, veg_21_w_dif, by.x = "PatchID", by.y = "PATCH_ID", all.x = TRUE)

r       <- rast(ext(result_weed), resolution = 5, crs = "epsg:25833")
weed_cover   <- raster(rasterize(vect(result_weed), r, field = "dif_mean"))
weed_cover_im <- as.im.RasterLayer(weed_cover)

plot(weed_cover_im, col = colorRampPalette(brewer.pal(9, "Purples"))(50))
plot(pp_obs_sub, add = TRUE, pch = 16, cex = 0.6, col = "grey20")


# --- 10.3 Vegetation height ------------------------------------------------
# Path: replace with a relative path once the project folder is finalised
vegh <- read.xlsx(file_veg_height)
vegh$DATE <- as.Date(vegh$DATE, origin = "1899-12-30")
vegh_21   <- vegh[vegh$DATE >= as.Date("2021-05-01") & vegh$DATE <= as.Date("2021-08-31"), ]

vegh_21_a <- aggregate(HEIGHT_AVG ~ PATCH_ID + DATE, data = vegh_21, mean)
for (i in seq_len(nrow(vegh_21_a))) {
  vegh_21_a$dif_mean[i] <- vegh_21_a$HEIGHT_AVG[i] -
    mean(vegh_21_a$HEIGHT_AVG[vegh_21_a$DATE == vegh_21_a$DATE[i]])
}
vegh_21_dif <- aggregate(dif_mean ~ PATCH_ID, data = vegh_21_a, mean)
result_vegh <- merge(result, vegh_21_dif, by.x = "PatchID", by.y = "PATCH_ID", all.x = TRUE)

r       <- rast(ext(result_vegh), resolution = 5, crs = "epsg:25833")
veg_h    <- raster(rasterize(vect(result_vegh), r, field = "dif_mean"))
veg_h_im <- as.im.RasterLayer(veg_h)

plot(veg_h_im, col = colorRampPalette(brewer.pal(9, "BuGn"))(50))
plot(pp_obs_sub, add = TRUE, pch = 16, cex = 0.6, col = "grey20")


# --- 10.4 Weed control activity --------------------------------------------
# Path: replace with a relative path once the project folder is finalised
wc <- read.xlsx(file_weed_control)
wc$DATE <- as.Date(wc$DATE, origin = "1899-12-30")
wc_21   <- wc[wc$DATE >= as.Date("2021-04-01") & wc$DATE <= as.Date("2021-07-31"), ]

# Count number of weed control events per patch
weed_count_21 <- data.frame(table(wc_21$PATCH_ID))
result_wc     <- merge(result, weed_count_21, by.x = "PatchID", by.y = "Var1", all.x = TRUE)
result_wc$Freq[is.na(result_wc$Freq)] <- 0   # Patches with no events → 0

r     <- rast(ext(result_wc), resolution = 5, crs = "epsg:25833")
weed_control    <- raster(rasterize(vect(result_wc), r, field = "Freq"))
weed_control_im <- as.im.RasterLayer(weed_control)

plot(weed_control_im, col = colorRampPalette(brewer.pal(9, "Greens"))(50))
plot(pp_obs_sub, add = TRUE, pch = 16, cex = 0.6, col = "grey20")


# --- 10.5 Pesticide activity -----------------------------------------------
# Path: replace with a relative path once the project folder is finalised
psm <- read.xlsx(file_psm)
psm$DATE <- as.Date(psm$DATE, origin = "1899-12-30")
psm_21   <- psm[psm$DATE >= as.Date("2021-04-01") & psm$DATE <= as.Date("2021-07-31"), ]

# Count number of pesticide application events per patch
psm_21    <- data.frame(table(psm_21$PATCH_ID))
result_psm <- merge(result, psm_21, by.x = "PatchID", by.y = "Var1", all.x = TRUE)
result_psm$Freq[is.na(result_psm$Freq)] <- 0

r      <- rast(ext(result_psm), resolution = 5, crs = "epsg:25833")
psm_control   <- raster(rasterize(vect(result_psm), r, field = "Freq"))
psm_control_im <- as.im.RasterLayer(psm_control)

plot(psm_control_im, col = colorRampPalette(brewer.pal(9, "Greens"))(50))
plot(pp_obs_sub, add = TRUE, pch = 16, cex = 0.6, col = "grey20")


# --- 10.6 Land use intensity (LUI) -----------------------------------------
# Path: replace with a relative path once the project folder is finalised
lu <- read.xlsx(file_land_use)

result_lui <- merge(result, lu[, c(1, 2)], by.x = "PatchID", by.y = "PATCH_ID", all.x = TRUE)

r      <- rast(ext(result_lui), resolution = 5, crs = "epsg:25833")
lui    <- raster(rasterize(vect(result_lui), r, field = "LUI"))
lui_im <- as.im.RasterLayer(lui)

plot(lui_im, col = colorRampPalette(brewer.pal(9, "Greens"))(50))
plot(pp_obs_sub, add = TRUE, pch = 16, cex = 0.6, col = "grey20")


# --- 10.7 Crop type --------------------------------------------------------
# Encode crop type as integer factor for raster representation
lu$CROP_TYPE_BREEDING_SEASON_2021 <- as.integer(
  factor(lu$CROP_TYPE_BREEDING_SEASON_2021)
)
result_ct <- merge(result, lu[, c(1, 4)], by.x = "PatchID", by.y = "PATCH_ID", all.x = TRUE)

r     <- rast(ext(result_ct), resolution = 5, crs = "epsg:25833")
ct_r  <- rasterize(vect(result_ct), r, field = "CROP_TYPE_BREEDING_SEASON_2021")
crop_type    <- raster(ct_r)
crop_type_im <- as.im.RasterLayer(crop_type)

plot(crop_type_im, col = colorRampPalette(brewer.pal(10, "Paired")))
plot(pp_obs_sub, add = TRUE, pch = 16, cex = 0.6, col = "grey20")


# --- 10.8 Spatial crop diversity -------------------------------------------
# For each raster cell: count the number of distinct crop types within
# a 100 m radius (local diversity index)

# Get cell indices with valid values (i.e. inside the field)
ncell  <- cells(ct_r)
coords <- st_sfc(
  st_multipoint(cbind(xyFromCell(crop_type, ncell)), dim = "XY"),
  crs = 25833
)
coords_points <- st_cast(coords, "POINT")

# 100 m buffers around each cell centre
b <- st_cast(st_as_sf(st_buffer(coords_points, 100)), "POLYGON")

# Count distinct crop types within each buffer
f <- lapply(seq_len(nrow(b)), function(i) {
  freq(mask(crop(ct_r, ext(b[i, ])), vect(b[i, ])))
})
b$crop <- unlist(lapply(f, function(x) length(x$value)))

# Write diversity values back into a raster
crop_div_r <- ct_r
values(crop_div_r)[ncell] <- b$crop

crop_div    <- raster(crop_div_r)
crop_div_im <- as.im.RasterLayer(crop_div)

plot(crop_div_im, col = colorRampPalette(brewer.pal(9, "Greens"))(50))
plot(pp_obs_sub, add = TRUE, pch = 16, cex = 0.6, col = "grey20")


# --- 10.9 Distance to flower strips ----------------------------------------
fs_r   <- rasterize(vect(fs_poly), r)
dist_r <- distance(fs_r)
dist   <- raster(dist_r)
dist_m <- mask(crop(dist, obs_spol[c(3)]), obs_spol[c(3)])

dist_im <- as.im.RasterLayer(dist_m)
plot(dist_im, col = colorRampPalette(brewer.pal(9, "Greens"))(50))
plot(pp_obs_sub, add = TRUE, pch = 16, cex = 0.6, col = "grey20")


# ============================================================================
# 11. COLLINEARITY CHECK
# ============================================================================

# Estimate observation density as a continuous surface (reference layer)
dens_obs     <- density(pp_obs_sub, sigma = 1, kernel = "gaussian", eps = 5)
dens_obs_res <- resample(raster(dens_obs), dist, method = "bilinear")
dens_obs_res <- projectRaster(dens_obs_res, dist, method = "bilinear")

# Stack all covariates as a terra SpatRaster
rast_all <- c(
  rast(dist), rast(veg_h), rast(weed_cover), rast(veg_cover), rast(psm_control),
  rast(weed_control), rast(lui), rast(crop_div), rast(crop_type),
  rast(dens_obs_res)
)
names(rast_all) <- c("dist", "veg_h", "weed_cover", "veg_cover", "psm_control",
                     "weed_control", "lui", "crop_div", "crop_type", "dens_obs")

# Pearson correlation matrix across all covariate layers
coll <- layerCor(rast_all, "pearson", na.rm = TRUE)
print(coll$pearson)

# Draw a random sample of 1 000 pixels (without replacement) for VIF estimation
smp <- spatSample(rast_all, 1000, replace = FALSE, na.rm = TRUE, as.df = TRUE)

# Linear model of observation density as a function of all covariates
mod <- lm(dens_obs ~ ., data = smp)

# Variance Inflation Factor (VIF) – values > 10 indicate multicollinearity
usdm::vif(rast_all)


# ============================================================================
# 12. PPM MODEL FITTING (Point Process Model)
# ============================================================================

# --- 12.1 Full model (all covariates) ---------------------------------------
m_all <- ppm(
  pp_obs_sub ~
    crop_div + crop_type + dist + weed_cover + 
    veg_cover + veg_h + weed_control + lui,
  # + psm   # Pesticides: currently excluded
  covariates = list(
    crop_div   = crop_div_im,
    crop_type  = crop_type_im,
    dist  = dist_im,
    weed_cover  = weed_cover_im,
    veg_cover = veg_cover_im,
    veg_h  = veg_h_im,
    weed_control    = weed_control_im,
    lui   = lui_im
    # psm = psm_im
  )
)

# --- 12.2 Backward stepwise selection via BIC --------------------------------
# k = log(n) corresponds to the BIC penalty
step(m_all, k = log(pp_obs_sub$n), direction = "backward")


# --- 12.3 Final model (after selection) -------------------------------------
# Selected covariates: spatial crop diversity (cd),
#                      vegetation height (vegh),
#                      land use intensity (lui)
m_21 <- ppm(
  pp_obs_sub ~ crop_div + veg_h + lui,
  covariates = list(
    crop_div   = crop_div_im,
    veg_h = veg_h_im,
    lui  = lui_im
  )
)

# BIC comparison: null model vs. final model
stats::BIC(m_pois, m_21)
