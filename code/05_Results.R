library(TMB)
library(sjSDM)
source("code/model.R")
models = readRDS( file = "model_results.RDS")


create_variance_partitioning = function(object) {
  indices = grep("SE_",colnames(object$data$X))
  sp = object$W[indices, ]
  covSP = stats::cov(object$data$X[, indices])
  
  beta = rbind(object$W[-c( 1, indices), ], object$LF)
  covX = stats::cov(cbind(object$data$X[, -c( 1, indices)], object$LV))
  vp = fa$importance(beta = beta, betaSP = sp, sigma = diag(1.0, 21), covX = covX, covSP = covSP)
  vp$biotic = vp$env[,22:23]
  vp$env = vp$env[,-(22:23)]

  res = list(split = vp, total = list(env = rowSums(vp$env), spatial = rowSums(vp$spatial), biotic = rowSums(vp$biotic)))
  out = list()
  out$names = colnames(object$data$Y)
  out$res = res
  out$spatial = TRUE
  class(out) = "sjSDMimportance"
  return(out)
}
vp_before = create_variance_partitioning(models$before$before_Eigen)
vp_after = create_variance_partitioning(models$after$after_Eigen)
par(mfrow = c(1,2))
colnames(models$before$before_Eigen$data$Y)[vp_before$res$total$spatial > 0.25]
plot(vp_before, col.points = ifelse(vp_before$res$total$spatial > 0.25, "red", "blue"))
plot(vp_after,  col.points = ifelse(vp_before$res$total$spatial > 0.25, "red", "blue"))
colnames(models$after$after_Eigen$data$Y)[vp_after$res$total$spatial > 0.25]

plot(models$before$before_CAR$LV)
plot(models$after$after_CAR$LV)
