library(TMB)
library(gllvm)

################ LVM-JSDM example ################

# Simulate data
l = 2L # number of latent variables
e = 3L # number of environmental variables
s = 25L # number of species
n = 200L # number of sites

X = matrix(runif(e*n), n, e) # Env
W = matrix(rnorm(e*s), e, s) # Env response -> estimated
LF = matrix(rnorm(l*s), l, s) # Factor loadings -> estimated
LV = mvtnorm::rmvnorm(n, rep(0, l)) # Latent variables -> estimated
YR = X%*% W + LV %*% LF
Y = apply(1/(1+exp(-YR)), 1:2, function(p) rbinom(1, 1, p))

# Fit model via TMB
## Compile and load model 
compile("TMB/LVM.cpp")
dyn.load(dynlib("TMB/LVM"))

## Define parameters
parameters = list(W=matrix(0.0, e, s) , 
                  LF = as.vector(matrix(rnorm(l*s), l, s))[-(1:sum(1:(l-1)))], 
                  LV = mvtnorm::rmvnorm(n, rep(0, l)))

## Make loss func
model = MakeADFun(list(x=X, y = Y), parameters, DLL = "LVM", random = c("LV"))

## Fit model
fit = nlminb(model$par, model$fn, model$gr)
res = sdreport(model)

W_TMB = matrix(fit$par[1:(e*s)], e, s)
LF_TMB = matrix(0, s, l)
counter = (e*s)+1
for(j in 1:l) {
  for(i in 1:s) {
    if(j>i) { LF_TMB[i, j] = 0 }
    else if (j == i) { 
      LF_TMB[i, j] = pnorm(fit$par[counter])
      counter= counter+1
    } else {
      LF_TMB[i, j] = -1+2*pnorm(fit$par[counter])
      counter= counter+1
    }
  }
}

LF_TMB = t(LF_TMB)

# Compare with gllvm
gl = gllvm(Y, data.frame(X), family = binomial(), num.lv = l, method = "VA")
W_GLLVM = t(coef(gl)$Xcoef)
LF_GLLVM = t(coef(gl)$theta)

# W/Beta:
sqrt(mean((W-W_TMB)**2))
sqrt(mean((W-W_GLLVM)**2))

# LF:
plot_covariance = function(L) fields::image.plot( cov2cor(t(L) %*% L + diag(1.0, ncol(L))))
return_cov = function(L) cov2cor(t(L) %*% L + diag(1.0, ncol(L)))
par(mfrow = c(1,3))
plot_covariance(LF)
plot_covariance(LF_TMB)
plot_covariance(LF_GLLVM)

cor(as.vector(return_cov(LF_TMB)), as.vector(return_cov(LF_GLLVM)))
