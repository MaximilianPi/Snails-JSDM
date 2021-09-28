################################################################################
#### Rhone River
################################################################################
# Description: Collect and download some data on the Rhone river with the
# primary aim of generating a plot of the study area

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)             # To handle spatial data
library(rgdal)              # To handle spatial data
library(rgeos)              # To buffer
library(tmap)               # For beautiful spatial plots
library(elevatr)            # Download elevation data
library(cleangeo)           # To make geometries valid
library(rnaturalearth)      # To download shapefile of the earth

################################################################################
#### Greater Extent
################################################################################
# Specify area of interest
aoi <- as(extent(c(4, 8, 43, 46.5)), "SpatialPolygons")
crs(aoi) <- CRS("+init=epsg:4326")

# Get shapefile of country borders
map <- ne_countries(scale = "medium", returnclass = "sp")

# Make crop and make valid
map <- crop(map, extent(-50, 50, 0, 100))
map <- gBuffer(map, byid = T, width = 0)

# Download shapefile of mayor european rivers
dir.create("DownloadedFiles", showWarnings = F)
file <- "https://www.eea.europa.eu/data-and-maps/data/wise-large-rivers-and-large-lakes/zipped-shapefile-with-wise-large-rivers-vector-line/zipped-shapefile-with-wise-large-rivers-vector-line/at_download/file"
file_name <- file.path("DownloadedFiles", basename(file))
download.file(file, destfile = file_name)
unzip(zipfile = file_name, exdir = "DownloadedFiles")

# Identify downloaded shapefile
rhone <- dir(path = "DownloadedFiles", pattern = "Large.*.shp$", full.names = T)

# Load it
rhone <- readOGR(rhone)
rhone <- subset(rhone, NAME == "Rhone")
rhone <- spTransform(rhone, crs(map))

# Visualize
p1 <- tm_shape(map, bbox = c(-10, 35, 25, 50)) +
    tm_polygons(
        col = "gray95"
      , lwd = 1.5
    ) +
    # tm_text("ADM0_A3") +
  tm_shape(aoi) +
    tm_polygons(
        col        = "blue"
      , alpha      = 0.1
      , border.col = "blue"
      , lwd        = 1.5
    ) +
  tm_shape(rhone) +
    tm_lines(col = "red", lwd = 2) #+
  # tm_graticules(
  #     n.y                 = 10
  #   , n.x                 = 10
  #   , labels.inside.frame = F
  #   , lines               = T
  #   , ticks               = T
  #   , lwd                 = 0.1
  #   , alpha               = 0.3
  # ) +
  # tm_compass(
  #     position    = c("left", "bottom")
  #   , color.light = "black"
  # ) +
  # tm_scale_bar(width = 0.2)

################################################################################
#### Detailed Extent (With Elevation + Hillshade)
################################################################################
# Download sampling sites
sit <- read.csv("https://raw.githubusercontent.com/pennekampster/Bayesian_thinking_UZH/main/Group_projects/SDM_Rhone_snails/data/coordinates.csv")
coordinates(sit) <- c("Lon_dd", "Lat_dd")
crs(sit) <- crs(rhone)

# Get the extent of these sites
ext_sites <- extent(sit) + 0.3
ext_sites <- as(ext_sites, "SpatialPolygons")
crs(ext_sites) <- crs(sit)

# Download high resolution country borders
che <- getData("GADM", level = 0, country = "CHE")
fra <- getData("GADM", level = 0, country = "FRA")
ita <- getData("GADM", level = 0, country = "ITA")
map2 <- rbind(che, fra, ita)
map2_crop <- crop(map2, aoi)

# Let's also download some elevation data (will take a moment)
ele <- get_elev_raster(aoi, z = 8)
ele <- crop(ele, aoi)

# Create hillshade
aspect  <- terrain(ele, opt = "aspect")
slope   <- terrain(ele, opt = "slope")
hill    <- hillShade(aspect = aspect, slope = slope)

# Download river data from diva gis and a shapefile of large rivers from the
# WISE dataset
files <- c(
    "https://biogeo.ucdavis.edu/data/diva/wat/FRA_wat.zip"
  , "https://biogeo.ucdavis.edu/data/diva/wat/CHE_wat.zip"
  , "https://biogeo.ucdavis.edu/data/diva/wat/ITA_wat.zip"
)

# Download files
for (i in seq_along(files)){
  file_name <- file.path("DownloadedFiles", basename(files[i]))
  download.file(files[i], destfile = file_name)
  unzip(zipfile = file_name, exdir = "DownloadedFiles")
}

# Identify downloaded shapefiles
river_files <- dir(path = "DownloadedFiles", pattern = "lines.*.shp", full.names = T)
lake_files <- dir(path = "DownloadedFiles", pattern = "areas.*.shp", full.names = T)

# Load them
rivers <- lapply(river_files, readOGR)
rivers <- do.call(rbind, rivers)
lakes <- lapply(lake_files, readOGR)
lakes <- do.call(rbind, lakes)

# Clean the "Rhone" line object
geneva <- lakes[grepl(lakes$NAME, pattern = "GENEVA"), ]
geneva <- gBuffer(geneva, width = 1.3 / 111)
rhone <- erase(rhone, geneva)

# Mask out raster values at see
hill <- mask(hill, map2)
ele <- mask(ele, map2)

# Visualite
p2 <- tm_shape(hill) +
    tm_raster(
        palette     = gray(60:100/100)
      , style       = "cont"
      , legend.show = F
    ) +
  tm_shape(ele) +
    tm_raster(
        palette     = terrain.colors(n = 100)
      , style       = "cont"
      , legend.show = F
      , alpha       = 0.3
  ) +
  tm_shape(map2) +
    tm_borders(col = "black") +
  # tm_shape(sea) +
  #   tm_polygons(col = "cornflowerblue", border.alpha = 0) +
  tm_shape(rivers) +
    tm_lines(col = "cornflowerblue", lwd = 0.3) +
  tm_shape(lakes) +
    tm_polygons(col = "cornflowerblue", border.alpha = 0) +
  tm_shape(rhone) +
    tm_lines(col = "red", lwd = 2) +
  tm_shape(sit) +
    tm_dots(col = "blue") +
  tm_shape(ext_sites) +
    tm_polygons(col = "blue", alpha = 0.2, border.col = "blue") +
  tm_graticules(
      n.y                 = 10
    , n.x                 = 10
    , labels.inside.frame = F
    , lines               = T
    , ticks               = T
    , lwd                 = 0.1
    , alpha               = 0.5
  ) +
  tm_compass(
      position    = c("left", "bottom")
    , color.light = "black"
  ) +
  tm_scale_bar(width = 0.2) +
  tm_layout(bg.color = "cornflowerblue")

# Show the two plots
p1
p2

# Store the maps to file
tmap_save(p1, filename = "Rhone_01.png", scale = 1.3)
tmap_save(p2, filename = "Rhone_02.png", scale = 1.3)
