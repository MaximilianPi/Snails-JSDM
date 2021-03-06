#' JSDM model function
#' 
#' @param Y sites by species matrix (must be a matrix)
#' @param X environment data.frame or matrix 
#' @param formula formula
#' @param l number of latent variables
#' @param randomEffects random intercept column
#' @param spatialRandomEffects list of two objects, 1. spatial covariates and 2. positive semidefinite matrix to inform the spatial CAR structure
#' 
#' JSDM model based on the LVM model (Warton, 2015). 
#' 
#' @import TMB
#' @import checkmate

JSDM_TMB = function(Y, X, formula = NULL, l = 2L, randomEffects, spatialRandomEffects) {
  
  library(TMB)
  library(checkmate)
  
  # check inputs
  assert(checkMatrix(Y))
  assert(checkList(spatialRandomEffects))
  assert(checkDataFrame(X), checkMatrix(X))
  qassert(l, "i1[0,)")
  
  # prepare input
  if(is.data.frame(X)) {
    if(!is.null(formula)){
      mf = match.call()
      m = match("formula", names(mf))
      if(class(mf[3]$formula) == "name") mf[3]$formula = eval(mf[3]$formula, envir = parent.env(environment()))
      formula = stats::as.formula(mf[m]$formula)
      X = stats::model.matrix(formula, X)
    } else {
      formula = stats::as.formula("~.")
      X = stats::model.matrix(formula, X)
    }
  } else {
    if(!is.null(formula)) {
      mf = match.call()
      m = match("formula", names(mf))
      if(class(mf[3]$formula) == "name") mf[3]$formula = eval(mf[3]$formula, envir = parent.env(environment()))
      formula = stats::as.formula(mf[m]$formula)
      X = data.frame(X)
      X = stats::model.matrix(formula, X)
    } else {
      formula = stats::as.formula("~.")
      X = stats::model.matrix(formula,data.frame(X))
    }
  }
  Z = model.matrix(~-1+as.factor(randomEffects)) # random effects w/o an intercept
  SP = model.matrix(~-1+as.factor(droplevels(spatialRandomEffects[[1]])))
  
  compile("TMB/LVM.cpp")
  dyn.load(dynlib("TMB/LVM"))
  
  # Define parameters
  s = ncol(Y)
  e = ncol(X)
  n = nrow(X)
  print(rep(0, ncol(SP)))
  parameters = list(W=matrix(0.0, e, s) , 
                    LF = as.vector(matrix(rnorm(l*s, sd = 0.1), l, s))[-(1:sum(1:(l-1)))], 
                    LV = mvtnorm::rmvnorm(n, rep(0, l), sigma = diag(0.1, l)),
                    dev=rep(0, ncol(Z)),
                    spatial=rep(0, ncol(SP)),
                    log_sd_dev = 0,
                    lambda = c(0,0))
  
  # Make loss func
  model = MakeADFun(list(x=X, 
                         y = Y, 
                         z = Z,
                         sp = SP,
                         D = spatialRandomEffects[[2]]), 
                    parameters, DLL = "LVM", random = c("LV", "dev", "spatial"))
  
  # Fit model
  fit = nlminb(model$par, model$fn, model$gr,control = list(eval.max = 400, iter.max = 300) )
  res = sdreport(model)
  
  # prepare output
  W = matrix(res$par.fixed[1:(e*s)], e, s)
  W_SD = matrix(res$sd[1:(e*s)], e, s)
  LF = matrix(0, s, l)
  counter = (e*s)+1
  for(j in 1:l) {
    for(i in 1:s) {
      if(j>i) { LF[i, j] = 0 }
      else if (j == i) { 
        LF[i, j] = pnorm(res$par.fixed[counter])
        counter= counter+1
      } else {
        LF[i, j] = -1+2*pnorm(res$par.fixed[counter])
        counter= counter+1
      }
    }
  }
  LF = t(LF)
  LV = matrix(res$par.random[1:(n*l)], n, l)
  RE = res$par.random[grep("dev", names(res$par.random))]
  
  # create predictions?
  pred = X%*%W + LV%*%LF
  ll = sum(dbinom(Y, 1, 1/(1+exp(-pred)), log = TRUE))
  
  # return results
  out = list()
  out$result = res
  out$model = model
  out$pred = pred
  out$logLik = ll
  out$W = W
  out$W_SD = W_SD
  out$LF = LF
  out$LV = LV
  out$RE = RE
  out$RE_SD = exp(res$par.fixed[["log_sd_dev"]])
  out$spatial = res$par.random[grep("spatial", names(res$par.random))]
  out$lambda = res$par.fixed[grep("lambda", names(res$par.fixed))]
  out$data = list(X = X, Y = Y, 
                  randomEffects = randomEffects, 
                  spatialRandomEffects = spatialRandomEffects)
  class(out) = "JSDM_TMB"
  return(out)
}

#' Calculate AIC
AIC.JSDM_TMB = function(object, ..., k = 2) {
  p = length(object$model$par)
  return(-2*object$logLik+2*p)
}

#' Calculate co-occurrence matrix
vcov.JSDM_TMB = function(object, ...) {
  return(cov2cor(t(object$LF) %*% object$LF + diag(1, ncol(object$LF))))
}
