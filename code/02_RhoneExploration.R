################################################################################
#### Rhone River
################################################################################
# Description: Exploring the Rhone River data

# Load required packages
library(tidyverse)
library(parzer)

# Clear R's brain
rm(list = ls())

################################################################################
#### Prepare Data
################################################################################
# Load keys
url <- "https://raw.githubusercontent.com/pennekampster/Bayesian_thinking_UZH/main/Group_projects/SDM_Rhone_snails/data/"
fauna_key <- read_csv2(file.path(url, "gastero_fauna_key.csv"))
envir_key <- read_csv2(file.path(url, "gastero_environment_key.csv"))

# Make proper column names
names(fauna_key) <- c("Family", "Genus", "Species", "Short")

# Make proper species names
fauna_key$Name <- paste(fauna_key$Genus, fauna_key$Species)

# Load datasets
envir <- read.csv2(file.path(url, "gastero_environment.csv"), row.names = 1, stringsAsFactors = F)
fauna <- read.csv2(file.path(url, "gastero_fauna.csv"), row.names = 1, stringsAsFactors = F)
sampl <- read_csv2(file.path(url, "gastero_samples.csv"))

# We also need the restoration information
resto <- read_csv2(file.path(url, "restoration_info.csv"))

# Parse the latitude and longitudes
resto$Latitude <- parse_lat(resto$Latitude)
resto$Longitude <- parse_lon(resto$Longitude)

# Transpose
envir <- as.data.frame(t(envir))
fauna <- as.data.frame(t(fauna))

# Replace NAs with 0s (THIS IS A BIG ASSUMPTION)
fauna[is.na(fauna)] <- 0
fauna <- as.data.frame(fauna)

# Bind data together
dat <- cbind(fauna, envir, sampl)
print(names(dat))

# Correct year
dat$year <- dat$year + 2000
dat$site <- gsub("AV", "DO", dat$site)
dat$site <- gsub("AM", "UP", dat$site)

# Remove double downlstream sites (why though?)
table(dat$site)
dat <- subset(dat, !grepl("N", dat$site))

# Add information on restoration
info <- left_join(dat, resto, by = c("channel" = "Channel", "site" = "Site")) %>%
  select(RestorationYear = rest_year, RestorationType = Type, Latitude, Longitude)
dat <- cbind(dat, info)

# Indicate if something is pre or post restoration
dat$restoration <- ifelse(dat$year <= dat$RestorationYear, "pre", "post")
dat$restoration[is.na(dat$restoration)] <- "not_restored"

# Let's recode the seasons a bit nicer
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

################################################################################
#### Explore Data
################################################################################
# Let's visualize the different locations
locs <- as.matrix(unique(dat[, c("Longitude", "Latitude")]))
plot(locs, pch = 20)

# Are there any NA's
sum(is.na(envir))
sum(is.na(fauna))
sum(is.na(sampl))

# Let's see in which columns they are
print(apply(envir, 2, function(x){!all(!is.na(x))}))
print(apply(fauna, 2, function(x){!all(!is.na(x))}))
print(apply(sampl, 2, function(x){!all(!is.na(x))}))

# Check dimensions
dim(dat)   # Whole dataset
dim(envir) # 19 Covariates
dim(fauna) # 26 Species
dim(sampl) # 6 Variables describing the sampling site

# What kind of data is there?
summary(dat$number)        # Unique ID for each sampling occasion
length(unique(dat$number)) # Not all samples remained in the dataset
count(dat, channel)  # Different channels (19 in total)
count(dat, site)     # Different sites (3 in total) -> Up Down Center?
count(dat, year)     # Year of sampling
count(dat, season)   # Season of sampling -> Spring (P), Summer (E), Autumn (A)
count(dat, restoration)       # Most samples collected after the restoration
count(dat, RestorationType)   # Different restoration methods

# Are all species present before and after the restoration?
dat %>%
  select(restoration, Acr_la:Viv_sp) %>%
  pivot_longer(2:ncol(.), names_to = "Species", values_to = "Count") %>%
  group_by(restoration, Species) %>%
  summarize(Present = ifelse(max(Count) > 0, T, F), .groups = "drop") %>%
  ggplot(aes(x = Species, y = restoration, fill = Present)) +
    geom_tile(col = "black") +
    scale_fill_manual(values = c("lightgray", "cornflowerblue")) +
    theme_minimal() +
    coord_equal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.8)) +
    ylab("Restoration")

# How does the sampling intensity vary over time?
dat %>%
  group_by(year, .drop = F) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = year, y = n)) +
    geom_point() +
    geom_line() +
    ylab("Number of Samples")
dat %>%
  group_by(year, season, .drop = F) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = year, y = n, col = season)) +
    geom_point() +
    geom_line() +
    ylab("Number of Samples by Season")
dat %>%
  group_by(year, channel, .drop = F) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = year, y = n)) +
    geom_point() +
    geom_line() +
    facet_wrap(~ channel) +
    ylab("Number of Samples by Site")

# Generate groups to separate channels from the upstream and downstream
dat$grp <- paste0(dat$channel, "_", dat$site)

# How many samples are there per group?
table(dat$grp, dat$year)

# Visualize timing of restoration
ggplot(dat, aes(x = year, y = grp, col = factor(restoration))) +
  geom_jitter(size = 0.5, width = 0.1) +
  theme_minimal() +
  scale_color_manual(values = c("gray30", "orange", "cornflowerblue"), name = "Restoration") +
  ylab("Group") +
  xlab("Year")

# Visualize with tiles
dat %>%
  dplyr::select(channel, year) %>%
  distinct() %>%
  ggplot(aes(x = year, y = channel)) +
    geom_tile() +
    coord_equal() +
    theme_minimal()

# Let's check the number of observations at different seasons
table(dat$channel, dat$season)      # Season "A" only rarely sampled
table(dat$channel, dat$year)        # Up to 16 samples per year
table(dat$year, dat$season)         # Season "A" only sampled in year 2

# Which species are present the most?
dat %>%
  dplyr::select(Acr_la:Viv_sp) %>%
  colSums() %>%
  as.tibble(rownames = "Species") %>%
  arrange(desc(value)) %>%
  left_join(fauna_key, by = c("Species" = "Short"))
