# Fitting JSDMs before and after restoration ------------------------------


# Load required packages
library(TMB)


# Clear R's brain
rm(list = ls())

source("code/model.R")

# Load and prepare data
load("SnailData.Rdata")
summary(dat) # 1 NA in shade
dat_no_na = dat[complete.cases(dat),]

# Separate into X (environment) and Y (Species)
Y = dat_no_na[,1:16]
X = dat_no_na[,-(1:16)]

# Fit model
D = as.matrix(dist(scale(cbind(X$northing[!duplicated(X$channel)], X$easting[!duplicated(X$channel)])), diag=TRUE, upper=TRUE))


model = JSDM_TMB(as.matrix(Y), X, 
                 formula = ~scale(depth)+scale(current), 
                 l = 2L, 
                 randomEffects = X$year,
                 spatialRandomEffects = list(X$channel, D))

AIC(model)
fields::image.plot( vcov(model1) )
fields::image.plot( vcov(model2) )

cor(as.vector(model1$W), as.vector(model2$W))
