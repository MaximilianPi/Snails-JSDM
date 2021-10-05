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
Y = dat_no_na[,1:21]
X = dat_no_na[,-(1:21)]



# Fit model


## Before restoration
dd = c(101:105)
#dd = 1:117
Y_before = Y[X$restoration=="pre" & X$season == "Summer",][dd,]
X_before = X[X$restoration=="pre"& X$season == "Summer",][dd,]
X_before$year = as.factor(X_before$year)
print(X[X$restoration=="pre"& X$season == "Summer",][dd,]$northing[1], digits = 15 )
print(X[X$restoration=="pre"& X$season == "Summer",][dd,]$easting[1], digits = 15)
# 5058701.60065336 701700.312977955 not working!

X_removed = X[abs(X$northing-5058701.60065336 )>0.01 & abs(X$easting-701700.312977955 )>0.01 , ]
Y_removed = Y[abs(X$northing-5058701.60065336 )>0.01 & abs(X$easting-701700.312977955 )>0.01 , ]


# Pre
Y_before = Y_removed[X_removed$restoration=="pre",]
X_before = X_removed[X_removed$restoration=="pre",]
X_before$year = as.factor(X_before$year)


coords = cbind(X_before$northing[!duplicated(X_before$channel[X_before$restoration=="pre"])], 
               X_before$easting[!duplicated(X_before$channel[X_before$restoration=="pre"])])
D_before = as.matrix(dist(scale(coords)), diag=TRUE, upper=TRUE)
X_before= X_before

before = JSDM_TMB(Y = as.matrix(Y_before), 
                 X = X_before, 
                 formula = ~ season + scale(depth) + scale(temperature)+
                   scale(current) + scale(conductivity) + scale(temperature) + scale(shade) +
                   scale(subs_lime) + scale(subs_sand) + scale(subs_grav) + scale(subs_rocks) + 
                   scale(org_rip) + scale(org_helo) + scale(org_wood) + scale(veg_sub) + scale(veg_float)+
                   scale(veg_emer) + scale(veg_sub_4)+ scale(veg_sub_3)+ scale(veg_sub_2)+ scale(veg_sub_1), 
                 l = 2L, 
                 randomEffects = X_before$year,
                 spatialRandomEffects = list(as.factor(X_before$channel), D_before))

AIC(before)
fields::image.plot( vcov(before) )


## After
Y_after = Y_removed[X_removed$restoration=="post",]
X_after = X_removed[X_removed$restoration=="post",]
X_after$year = as.factor(X_after$year)


coords = cbind(X_after$northing[!duplicated(X_after$channel[X_after$restoration=="post"])], 
               X_after$easting[!duplicated(X_after$channel[X_after$restoration=="post"])])
D_after = as.matrix(dist(scale(coords)), diag=TRUE, upper=TRUE)

after = JSDM_TMB(Y = as.matrix(Y_after), 
                  X = X_after, 
                  formula = ~ season + scale(depth) + scale(temperature)+
                    scale(current) + scale(conductivity) + scale(temperature) + scale(shade) +
                    scale(subs_lime) + scale(subs_sand) + scale(subs_grav) + scale(subs_rocks) + 
                    scale(org_rip) + scale(org_helo) + scale(org_wood) + scale(veg_sub) + scale(veg_float)+
                    scale(veg_emer) + scale(veg_sub_4)+ scale(veg_sub_3)+ scale(veg_sub_2)+ scale(veg_sub_1), 
                  l = 2L, 
                  randomEffects = X_after$year,
                  spatialRandomEffects = list(as.factor(X_after$channel), D_after))

AIC(after)
fields::image.plot( vcov(after) )

saveRDS(list(before = before, after = after), file = "model_results.RDS")

