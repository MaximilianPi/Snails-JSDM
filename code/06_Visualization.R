################################################################################
#### Visualizations
################################################################################
# Description: Visualizing results from TMB models

# Load required packages
library(tidyverse)
library(raster)

# Clear R's brain
rm(list = ls())

# Set working directory
setwd("/home/david/Schreibtisch/Snails-JSDM")

# Reload model results
load("SnailData.Rdata")
source("code/model.R")
models <- readRDS("model_results.RDS")
models$before <- models$before$before_CAR
models$after <- models$after$after_CAR


url <- "https://raw.githubusercontent.com/pennekampster/Bayesian_thinking_UZH/main/Group_projects/SDM_Rhone_snails/data/"
fauna_key <- read_csv2(file.path(url, "gastero_fauna_key.csv"), col_names = F)
names(fauna_key) <- c("family", "genus", "species", "shortname")
fauna_key$name <- paste(fauna_key$genus, fauna_key$species)
fauna_key <- fauna_key[, c("shortname", "name")]

################################################################################
#### Species Matrix
################################################################################
# Using fields
fields::image.plot(vcov(models$before), main = "Before Restoration")
fields::image.plot(vcov(models$after), main = "After Restoration")

# For ggplot, reshape data
par(mfrow = c(1, 2))
species <- colnames(models$before$data$Y)
species <- factor(species, levels = species)
covariances <- expand_grid(species1 = species, species2 = species)
covariances$before <- as.vector(vcov(models$before))
covariances$after <- as.vector(vcov(models$after))
covariances <- pivot_longer(covariances, 3:4, names_to = "when", values_to = "covar")
covariances$when <- factor(covariances$when, levels = c("before", "after"))

# Compute total covariance before and after
sum(vcov(models$before))
sum(vcov(models$after))

# # Compute difference in matrices
# fields::image.plot(vcov(models$before) - vcov(models$after))

# Plot
p1 <- ggplot(covariances, aes(x = species1, y = species2, fill = covar)) +
  geom_tile() +
  facet_wrap("when") +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = levels(covariances$species1)) +
  scale_fill_gradientn(colors = rev(rainbow(20, end = 0.7))) +
  xlab("") +
  ylab("")
ggsave("CovariancePlot.png", plot = p1, width = 8, height = 5)

################################################################################
#### Weigths
################################################################################
# Reshape
betas <- as.data.frame(rbind(models$before$W, models$after$W))
names(betas) <- colnames(models$before$data$Y)
betas$covariate <- rep(colnames(models$before$data$X), 2)
betas$model <- rep(c("before", "after"), each = nrow(models$before$W_SD))
betas <- pivot_longer(betas, Anc_fl:Val_pi
  , names_to  = "species"
  , values_to = "beta"
)

sds <- as.data.frame(rbind(models$before$W_SD, models$after$W_SD))
names(sds) <- colnames(models$before$data$Y)
sds$covariate <- rep(colnames(models$before$data$X), 2)
sds$model <- rep(c("before", "after"), each = nrow(models$before$W_SD))
sds <- pivot_longer(sds, Anc_fl:Val_pi
  , names_to  = "species"
  , values_to = "sd"
)

# Put together
dat <- cbind(betas, sds[, c("sd")])
dat <- mutate(dat, across(covariate:species, as.factor))

# We don't care about the intercept
dat <- subset(dat, covariate != "(Intercept)")

# Remove "scale"
dat$covariate <- gsub(dat$covariate, pattern = "scale\\(", replacement = "")
dat$covariate <- gsub(dat$covariate, pattern = "\\)", replacement = "")

# Find species with lowest correlation of niche before and after restoration
correlations <- dat %>%
  dplyr::select(-sd) %>%
  pivot_wider(names_from = model, values_from = beta) %>%
  group_by(species) %>%
  nest() %>%
  mutate(correlation = map(data, function(x) {
    sub <- subset(x, covariate != "(Intercept)")
    correlation <- cor(sub$before, sub$after)
    return(correlation)
  })) %>%
  dplyr::select(-data) %>%
  unnest(correlation) %>%
  arrange(correlation) %>%
  head(n = 4)

# Subset to those species and remove the intercept
dat_lowcorr <- subset(dat, species %in% correlations$species & covariate != "(Intercept)")

# For each remaining species, identify the strongest niches before the
# restoration
dat_lowcorr_strongest <- dat_lowcorr %>%
  group_by(species) %>%
  nest() %>%
  mutate(data = map(data, function(x) {
    sub <- subset(x, covariate != "(Intercept)")
    sub <- arrange(sub, desc(model), desc(abs(beta)))
    sub <- head(sub$covariate, n = 5)
    sub <- subset(x, covariate %in% sub)
    return(sub)
  })) %>%
  unnest(data)

# Visualize
dat_lowcorr_strongest <- left_join(dat_lowcorr_strongest, fauna_key, by = c("species" = "shortname"))
dat_lowcorr_strongest$model <- factor(dat_lowcorr_strongest$model, levels = c("before", "after"))
p2 <- dat_lowcorr_strongest %>%
  ggplot(aes(x = beta, y = covariate, col =  model)) +
    geom_vline(xintercept = 0, col = "gray50", size = 0.2) +
    geom_point(position = position_dodge(width = 0.7)) +
    geom_errorbarh(aes(xmin = beta - sd, xmax = beta + sd), height = 0, position = position_dodge(width = 0.7)) +
    facet_wrap("name", scales = "free", nrow = 1) +
    theme_minimal() +
    scale_color_manual(values = c("cornflowerblue", "orange"), name = "Restoration") +
    xlab(expression(beta*-"Coefficient")) +
    ylab("Covariate")
ggsave("BetaPlot.png", plot = p2, bg = "white", height = 5, width = 12)
