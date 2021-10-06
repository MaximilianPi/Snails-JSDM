library(TMB)
source("code/model.R")
models = readRDS( file = "model_results.RDS")
before = models$before
after = models$after

par(mfrow = c(1, 2))
fields::image.plot( vcov(before), main = "Before restoration")
fields::image.plot( vcov(after), main = "After restoration")


as.vector(before$W[1:5,])
eff_names = colnames(before$data$X)[1:5]
expand.grid(eff_names, paste0(colnames(before$data$Y), ":") )
effect = data.frame( Estimate=as.vector(before$W[1:5,]),
                     Std.Err=as.vector(before$W_SD[1:5,]),P=0,rownames=)

