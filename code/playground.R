library(TMB)
library(gllvm)

################ LVM-JSDM example ################

# Simulate data
l = 2L # number of latent variables
e = 2L # number of environmental variables
s = 20L # number of species
n = 200L # number of sites

X = matrix(runif(e*n), n, e)
W = matrix(rnorm(e*s), e, s)
LF = matrix(rnorm(l*s), l, s)
LV = mvtnorm::rmvnorm(n, c(0.0, 0.0))
YR = X%*% W + LV%*%(LF)
Y = apply(1/(1+exp(-YR)), 1:2, function(p) rbinom(1, 1, p))

# Fit model via TMB
## Compile and load model 
compile("TMB/LVM.cpp")
dyn.load(dynlib("TMB/LVM"))

## Define parameters
parameters = list(W=matrix(0.0, e, s) , 
                  LF = matrix(rnorm(l*s), l, s), 
                  LV = mvtnorm::rmvnorm(n, c(0.0, 0.0)))

## Make loss func
model = MakeADFun(list(x=X, y = Y), parameters, DLL = "LVM", random = "LV")

## Fit model
fit = nlminb(model$par, model$fn, model$gr)
W_TMB = matrix(fit$par, e, s)
LF_TMB = matrix(fit$par[names(fit$par) %in% "LF"], l, s)

# Compare with gllvm
gl = gllvm(Y, data.frame(X), family = binomial(), num.lv = 2L)
W_GLLVM = t(coef(gl)$Xcoef)
LF_GLLVM = t(coef(gl)$theta)

# W/Beta:
sqrt(mean((W-W_TMB)**2))
sqrt(mean((W-W_GLLVM)**2))

# LF:
plot_covariance = function(L) fields::image.plot( cov2cor(t(L) %*% L + diag(1.0, ncol(L))))
par(mfrow = c(1,3))
plot_covariance(LF)
plot_covariance(LF_TMB)
plot_covariance(LF_GLLVM)
