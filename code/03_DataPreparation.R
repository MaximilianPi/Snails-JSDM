################################################################################
#### Rhone River
################################################################################
# Description: Preparing the Rhone dataset for fitting with TMB

# Load required packages
library(tidyverse)
library(parzer)
library(rgdal)
library(sp)

# Clear R's brain
rm(list = ls())

################################################################################
#### Prepare Data
################################################################################
# Load keys
url <- "https://raw.githubusercontent.com/pennekampster/Bayesian_thinking_UZH/main/Group_projects/SDM_Rhone_snails/data/"
fauna_key <- read_csv2(file.path(url, "gastero_fauna_key.csv"), col_names = F)
envir_key <- read_csv2(file.path(url, "gastero_environment_key.csv"), col_names = F)

# Make proper column names
names(fauna_key) <- c("family", "genus", "species", "shortname")

# Make proper species names
fauna_key$name <- paste(fauna_key$genus, fauna_key$species)

# Load datasets
envir <- read.csv2(file.path(url, "gastero_environment.csv"), row.names = 1, stringsAsFactors = F)
fauna <- read.csv2(file.path(url, "gastero_fauna.csv"), row.names = 1, stringsAsFactors = F)
sampl <- read_csv2(file.path(url, "gastero_samples.csv"))

# We also need the restoration information
resto <- read_csv2(file.path(url, "restoration_info.csv"))

# Parse the latitude and longitudes
resto$Latitude <- parse_lat(resto$Latitude)
resto$Longitude <- parse_lon(resto$Longitude)

# Reproject to utm
coordinates(resto) <- ~ Longitude + Latitude
proj4string(resto) <- CRS("+init=epsg:4326")
resto <- spTransform(resto, CRS("+proj=utm +zone=31 ellps=WGS84"))
resto <- as.data.frame(resto)
resto <- rename(resto, Northing = Latitude, Easting = Longitude)

# Transpose
envir <- as.data.frame(t(envir), stringsAsFactors = F)
fauna <- as.data.frame(t(fauna), stringsAsFactors = F)

# Replace NAs
fauna[is.na(fauna)] <- 0
fauna <- as.data.frame(fauna)

# Bind data together and convert covariates to numeric
dat <- cbind(fauna, envir, sampl) %>%
  mutate(across(depth:veg_sub_1, as.numeric))

# Correct year
dat$year <- dat$year + 2000
dat$site <- gsub("AV", "DO", dat$site)
dat$site <- gsub("AM", "UP", dat$site)

# Remove some sites
dat <- subset(dat, !grepl("N", dat$site))

# Add information on restoration make column names lower_case
info <- left_join(dat, resto, by = c("channel" = "Channel", "site" = "Site")) %>%
  dplyr::select(
      restoration_year = rest_year
    , restoration_type = Type
    , northing         = Northing
    , easting          = Easting
  )
dat <- cbind(dat, info)

# Indicate if something is pre or post restoration
dat$restoration <- ifelse(dat$year <= dat$restoration_year, "pre", "post")

# Everything that has never been restored is here assumed to be "pre-restoration"
dat$restoration[is.na(dat$restoration)] <- "pre"

# Recode the seasons to more intuitive names
dat <- mutate(dat, season = case_when(
    season == "P" ~ "Spring"
  , season == "E" ~ "Summer"
  , season == "A" ~ "Autumn"
))

# Make some columns factorial
dat$channel     <- as.factor(dat$channel)
dat$site        <- as.factor(dat$site)
dat$season      <- as.factor(dat$season)
dat$restoration <- as.factor(dat$restoration)

# Remove undesired columns
dat$number <- NULL
dat$restoration_year <- NULL
dat$restoration_type <- NULL

# Make unique channel names for upstream, downstream etc.
dat$channel <- paste0(dat$channel, "_", dat$site)
dat$site <- NULL

# Scale land cover covariates to values between 0 and 1
dat <- dat %>%
  mutate(across(subs_lime:veg_sub_1, function(x) {
    (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
  }))

# Aggregate data across "q" samples
dat <- dat %>%
  group_by(channel, year, season, restoration, northing, easting) %>%
  nest() %>%
  mutate(data = map(data, function(x) {
    x %>%
      colMeans() %>%
      t() %>%
      as.data.frame() %>%
      mutate(across(Acr_la:Viv_sp, function(y) {
        ifelse(y > 0, 1, 0)
      })) %>%
      dplyr::select(-q)
  })) %>%
  unnest(data) %>%
  ungroup()

# Remove autumn samples
dat <- subset(dat, season != "Autumn")

# Arrange data
dat <- arrange(dat, channel, year) %>% select(Acr_la:Viv_sp, everything())

# Find species that are present before and after the restoration
keep <- dat %>%
  dplyr::select(restoration, Acr_la:Viv_sp) %>%
  pivot_longer(2:ncol(.), names_to = "species", values_to = "count") %>%
  group_by(restoration, species) %>%
  summarize(present = ifelse(max(count) > 0, T, F), .groups = "drop") %>%
  pivot_wider(names_from = restoration, values_from = present) %>%
  subset(post & pre) %>%
  pull(species)

# Subset to those species
dat <- dat %>% dplyr::select(all_of(keep), channel:veg_sub_1)

# Store data to file
save(dat, file = "SnailData.Rdata")
