# Script for Generating Spatial Eigen Values.
# generateSpatialEv function was adopted from the sjSDM package

# Clear R's brain
rm(list = ls())
load("SnailData.Rdata") #load the clean dataset
# source("code/generateSpatialEV.R")

generateSpatialEV = function(coords = NULL, threshold = 0.0) {
  ## create dist ##
  dist = as.matrix(stats::dist(coords))
  zero = diag(0.0, ncol(dist))
  
  ## create weights ##
  if (threshold > 0) dist[dist < threshold] = 0
  
  distW = 1/dist
  distW[is.infinite(distW)] = 1
  diag(distW) <- 0
  rowSW =  rowSums(distW)
  rowSW[rowSW == 0] = 1
  distW <- distW/rowSW
  
  ## scale ##
  rowM = zero + rowMeans(distW)
  colM = t(zero + colMeans(distW))
  distC = distW - rowM - colM + mean(distW)
  
  eigV = eigen(distC, symmetric = TRUE)
  values = eigV$values / max(abs(eigV$values))
  SV = eigV$vectors[, values>0]
  colnames(SV) = paste0("SE_", 1:ncol(SV))
  return(SV)
}

#For Pre-Restoration 
pre_resData <- dat[dat$restoration=="pre",] #subset the data to pre-restoration
pre_resCoord <- data.frame(Easting=pre_resData$easting, Northing=pre_resData$northing)  #extract the coordinates for the pre-restoration sites
head(pre_resCoord)
str(pre_resCoord)
pre_resCoord_SEV <- generateSpatialEV(coords = pre_resCoord, threshold = 0.2)

# For Post restoration
post_resDat <- dat[dat$restoration=="post",]
post_resCoord <- data.frame(Easting=post_resDat$easting, Northing=post_resDat $northing) 
post_resCoord_SEV <- generateSpatialEV(coords = post_resCoord, threshold = 0.2)

save(pre_resCoord_SEV, file="Pre_restoration_spatial_eigenValues.RData")
save(post_resCoord_SEV, file="Post_restoration_spatial_eigenValues.RData")

# class(pre_resCoord)

barplot(pre_resCoord_SEV)
barplot(post_resCoord_SEV)
                                       