#' JSDM model function
#' 
#' @param Y sites by species matrix (must be a matrix)
#' @param X environment data.frame or matrix 
#' @param formula formula
#' @param l number of latent variables
#' @param randomEffects random intercept column
#' 
#' JSDM model based on the LVM model (Warton, 2015). 
#' 
#' @import TMB
#' @import checkmate

JSDM_TMB = function(Y, X, formula = NULL, l = 2L, randomEffects) {
  
  library(TMB)
  library(checkmate)
  
  # check inputs
  assert(checkMatrix(Y))
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
  
  compile("TMB/LVM.cpp")
  dyn.load(dynlib("TMB/LVM"))
  
  # Define parameters
  s = ncol(Y)
  e = ncol(X)
  n = nrow(X)
  parameters = list(W=matrix(0.0, e, s) , 
                    LF = as.vector(matrix(rnorm(l*s, sd = 0.1), l, s))[-(1:sum(1:(l-1)))], 
                    LV = mvtnorm::rmvnorm(n, rep(0, l), sigma = diag(0.1, l)),
                    dev=rep(0, ncol(Z)),
                    log_sd_dev = 0)
  
  # Make loss func
  model = MakeADFun(list(x=X, y = Y, z = Z), parameters, DLL = "LVM", random = c("LV", "dev"))
  
  # Fit model
  fit = nlminb(model$par, model$fn, model$gr)
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
  RE = res$par.random[(n*l+1):length(res$par.random)]
  
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
  out$data = list(X = X, Y = Y)
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
