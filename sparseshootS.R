sparseshooting <- function(x, y, k = 3.420, maxIteration = 100, tol = 10^-2, 
                           betaEst = NULL, intercept = NULL, scaleVar = NULL, xhat = NULL, xtilde = NULL,
                           maxituniv = 1, maxitscale = 100, wvalue = 3, shoot_order = "default",
                           nlambda = 100, post = TRUE, standardize = TRUE, lambda_grid = NULL, predset = NULL){ 
  #### Function to compute sparse shooting S ####
  
  # x: nxp matrix of predictors
  # y: n-dimensional response vector
  # k: tuning parameter for rho-function
  # maxIteration: maximum number of iterations shooting loop
  # tol: numerical covergence shooting loop
  # betaEst: initial regression estimates
  # intercept: initial intercept 
  # scaleVar: initial scale
  # xhat : nxp matrix of expected predictors
  # xtilde: nxp matrix of pre-processed predictors
  # maxituniv: number of iterations IRLS (1: one-step reweighting)
  # maxitscale: maximum number of iterations fixed point algorithm for computing residual scale
  # wvalue : cut-off for outlier flagging 
  # shoot_order : if "default" : loop through variables in order 1 ...p if "robcor" loop through variables in order of correlation with response
  # nlambda: number of values in lambda grid
  # post : TRUE for post-lasso estimation, FALSe for lasso estimation
  # standardize : TRUE standardize predictors, FALSE do not standardize predictors
  # lambda_grid : sequence of regularization parameters
  
  #### Preliminaries ####
  n <- nrow(x) # sample size
  p <- ncol(x) # number of predictors 
  delta <- (1 - 3/k^2 + 5/k^4 - k^2/3) *pnorm(k) + (4/(3*k) - k/3 - 5/k^3)*dnorm(k) - 1/2 + 3/(2*k^2) - 5/(2*k^4) + k^2/3
  
  if(standardize){
    Xest <- Xestimfast(x, value = wvalue)
    Xclean <- apply(x, 2, Xinitftc, Xmatrix = x, Xest = Xest, value = wvalue) 
    sx <- apply(Xclean, 2, sd)
    mx <- apply(Xclean, 2, mean)
    x <-  scale(x, center = mx, scale = sx)
  }
  
  if(is.null(xhat)){
    xhat <- Xestimfast(Xmatrix = x, value = wvalue)
  }
  if(is.null(xtilde)){
    xtilde <- apply(X = x, 2, Xinitftc, Xmatrix = x, Xest = xhat)
  }
  xhat0 <- xhat
  xtilde0 <- xtilde
  start_flag <- xtilde!=x
  
  if(shoot_order == "default"){
    order_variables <- 1:p
  }
  
  robcor_fit <- NULL
  if(shoot_order == "robcor"){
    robcor_fit <- robcorr(X = x, Y = y)
    order_variables <- robcor_fit$predcor
  }
  
  startfit <- NULL
  if(is.null(betaEst)){ # Initialization
    startfit <- startvalue_MM(X = x, Y = y, value = wvalue, robcor_fit = robcor_fit, Xinit = xtilde, predset = predset)
    betaEst <- as.matrix(startfit$betaEst)
    intercept <- as.matrix(startfit$intercept)
    scaleVar <- startfit$scaleVar
  }
  
  betaEst <- as.matrix(betaEst[order_variables, 1])
  intercept <- as.matrix(intercept[order_variables, 1])
  scaleVar <- scaleVar[order_variables]
  xtilde <- xtilde[, order_variables]
  ytilde <- y - xtilde%*%betaEst
  x <- x[, order_variables]
  xhat <- xhat[, order_variables]
  
  # Initialization
  wt <- matrix(NA, ncol=p, nrow=n)
  y <- as.matrix(y)
  
  # Get lambda grid
  if(is.null(lambda_grid)){
    lambda_max <- get_lambda_max(xtilde,y,ytilde,x,xhat,betaEst,intercept,scaleVar,k,delta,maxIteration,tol,maxitscale,maxituniv, 
                                 post, value = wvalue)
    lambda_grid <- getLambdaGrid(lambda_max,n,p,nlambda)
  }
  
  
  fit_regression <- selectModel(lambda_grid = lambda_grid, ytilde = ytilde, xtilde = xtilde, xhat = xhat, x = x, wt = wt,
                                betaEst = betaEst, intercept = intercept, scaleVar = scaleVar, k = k, delta = delta,
                                maxIteration = maxIteration, y = y, maxitscale = maxitscale, maxituniv = maxituniv, tol = tol, 
                                post = post,  value = wvalue)
  
  # Put the coefficients and the weights back in the order of the columns X provided by the user
  betahat <- betahat_ln <- rep(NA, p)
  
  # Results for BIC with sigma and ln_sigma
  betahat[order_variables] <- fit_regression$coef_opt
  betahat_ln[order_variables] <- fit_regression$coef_opt_ln
  betahats <- lapply(fit_regression$get_coefs$betaEst, p = p, order_variables = order_variables, function(X, p, order_variables){get_beta <- rep(NA, p); get_beta[order_variables] <- X; return(get_beta)})
  
  alphahat <- fit_regression$alpha_opt
  alphahat_ln <- fit_regression$alpha_opt_ln
  alphahats <- fit_regression$get_coefs$alpha
  weights <- weights_ln <- matrix(NA, n, p)
  weights[, order_variables] <- fit_regression$weights_opt
  weights_ln[, order_variables] <- fit_regression$weights_opt_ln
  
  if(standardize){
    betahat <-  betahat/sx 
    betahat_ln <- betahat_ln/sx
    alphahat <- alphahat - sum(betahat*mx)
    alphahat_ln <- alphahat_ln - sum(betahat_ln*mx)
    
    betahats <- lapply(betahats, sx = sx, function(X, sx){X/sx})
    alphahats <- lapply(fit_regression$lambdas, lambdas = fit_regression$lambdas, alphas = alphahats, betas = betahats, 
                        mx = mx, function(X, lambdas, alphas, betas, mx){
                          index <- which(lambdas==X)
                          constant <- alphas[[index]] - sum(betas[[index]]*mx)
                          return(constant)
                        })
  }
  
  out <- list("coef" = c(alphahat, betahat), "weights" = weights, "iter" = fit_regression$iter_opt, 
              "fits" = fit_regression, "betahats" = betahats, "alphahats" = alphahats , 
              "lambda_opt" = fit_regression$lambda_opt, "lambdas" = fit_regression$lambdas, "start_flag" = start_flag,
              "coef_ln" = c(alphahat_ln, betahat_ln), "weights_ln" = weights_ln, "startfit" = startfit, "xhat" = xhat0, "xtilde" = xtilde0)
}

shooting <- function(x, y, k = 3.420, maxIteration = 100, tol = 10^-2, 
                     betaEst = NULL, intercept = NULL, scaleVar = NULL, xhat = NULL, xtilde = NULL,
                     maxituniv = 1, maxitscale = 100, wvalue = 3, shoot_order = "default", standardize = TRUE, predset = NULL){
  #### Function to compute shooting S ####
  
  # x: nxp matrix of predictors
  # y: n-dimensional response vector
  # k: tuning parameter for rho-function
  # maxIteration: maximum number of iterations shooting loop
  # tol: numerical covergence shooting loop
  # betaEst: initial regression estimates
  # intercept: initial intercept 
  # scaleVar: initial scale
  # xhat : nxp matrix of expected predictors
  # xtilde: nxp matrix of pre-processed predictors
  # maxituniv: number of iterations IRLS (1: one-step reweighting)
  # maxitscale: maximum number of iterations fixed point algorithm for computing residual scale
  # wvalue : cut-off for outlier flagging 
  # shoot_order : if "default" : loop through variables in order 1 ...p if "robcor" loop through variables in order of correlation with response
  # standardize : TRUE standardize predictors, FALSE do not standardize predictors
  
  
  #### Preliminaries ####
  n <- nrow(x) # sample size
  p <- ncol(x) # number of predictors 
  delta <- (1 - 3/k^2 + 5/k^4 - k^2/3) *pnorm(k) + (4/(3*k) - k/3 - 5/k^3)*dnorm(k) - 1/2 + 3/(2*k^2) - 5/(2*k^4) + k^2/3
  
  if(standardize){
    Xest <- Xestimfast(x, value = wvalue)
    Xclean <- apply(x, 2, Xinitftc, Xmatrix = x, Xest = Xest, value = wvalue) 
    sx <- apply(Xclean, 2, sd)
    mx <- apply(Xclean, 2, mean)
    x <-  scale(x, center = mx, scale = sx)
  }
  
  if(is.null(xhat)){
    xhat <- Xestimfast(Xmatrix = x, value = wvalue)
  }
  if(is.null(xtilde)){
    xtilde <- apply(X = x, 2, Xinitftc, Xmatrix = x, Xest = xhat)
  }
  xhat0 <- xhat
  xtilde0 <- xtilde
  start_flag <- xtilde!=x 
  
  if(shoot_order == "default"){
    order_variables <- 1:p
  }
  
  robcor_fit <- NULL
  if(shoot_order == "robcor"){
    robcor_fit <- robcorr(X = x, Y = y)
    order_variables <- robcor_fit$predcor
  }
  
  startfit <- NULL
  if(is.null(betaEst)){ # Initialization
    startfit <- startvalue_MM(X = x, Y = y, value = wvalue, robcor_fit = robcor_fit, Xinit = xtilde, predset = predset)
    betaEst <- as.matrix(startfit$betaEst)
    intercept <- as.matrix(startfit$intercept)
    scaleVar <- startfit$scaleVar
  }
  
  betaEst <- as.matrix(betaEst[order_variables, 1])
  intercept <- as.matrix(intercept[order_variables, 1])
  scaleVar <- scaleVar[order_variables]
  xtilde <- xtilde[, order_variables]
  ytilde <- y - xtilde%*%betaEst
  x <- x[, order_variables]
  xhat <- xhat[, order_variables]
  
  # Initialization
  wt <- matrix(NA, ncol=p, nrow=n)
  y <- as.matrix(y)
  
  #### 1. Shooting loop ####
  shootfit <- shootloop(ytilde = ytilde, xtilde = xtilde, x = x, wt = wt, Xexp = xhat,  
                        betaEst = betaEst, intercept = intercept, scaleVar = scaleVar, scaleVarOLD = scaleVar,
                        k = k, delta = delta, maxIteration = maxIteration, tol=mad(y)*10^-6,
                        tolout=tol, tolscale=10^-5, maxitscale = maxitscale, y = y, maxituniv = 1 , wcut = wvalue) # implemented in C++
  
  # Put the coefficients and the weights back in the order of the columns X provided by the user
  betahat <- rep(NA, p)
  betahat[order_variables] <- shootfit$betaEst
  alphahat <- shootfit$alpha
  weights <- matrix(NA, n, p)
  weights[, order_variables] <- shootfit$weights
  
  if(standardize){
    betahat <-  betahat/sx 
    alphahat <- alphahat - sum(betahat*mx)
  }
  
  out <- list("coef" = c(alphahat, betahat), "weights" = weights,  "iter" = shootfit$iter, "startfit" = startfit,
              "xhat" = xhat0, "xtilde" = xtilde0, "start_flag" = start_flag) 
}



#### AUXILIARY FUNCTIONS ####

robcorr <- function(X, Y){
  # Kendall's correlations
  robcor <- apply(X, 2, cor.fk, y=Y)
  predcor <- order(abs(robcor), decreasing = T)
  
  out <- list("robcor" = robcor, "predcor" = predcor)
}

Xestimfast <- function(Xmatrix, value = 3){ 
  # Function to obtain xhat
  
  cork <- abs(cor.fk(Xmatrix))
  diag(cork) <- 0
  cormax <- apply(cork, 2, which.max)
  
  MMuniv <- function(U, indexmax, Xdata){
    i.variable <- which(apply(U==Xdata, 2, all)==T)
    i.pred <- indexmax[i.variable]
    Xstand <- robStandardize(U)
    mx <- attr(Xstand, "center")
    sx <- attr(Xstand, "scale")
    
    # Univariate Robust Regression
    fit <- suppressWarnings(lmrob(U ~ Xdata[, i.pred], setting = "KS2011"))
    rs <- fit$residuals/fit$scale
    Xstar <- fit$fitted.values # Fitted values
    xstars <- (Xstar - mx)/sx
    
    flag <- abs(xstars) > value
    Xhat <- Xstar*(1-flag) + median(U)*flag 
  }
  
  Xhat <- apply(Xmatrix, 2, MMuniv, indexmax = cormax, Xdata = Xmatrix)
  
  return(Xhat)
}

Xinitftc <- function(U, Xmatrix, Xest, value = 3){
  # Function to obtain pre-processed x-values
  i.variable <- which(apply(U==Xmatrix, 2, all)==T)
  
  # Check original observation with univariate outlier detection tool
  xx <- robStandardize(U)
  mx <- attr(xx, "center")
  sx <- attr(xx, "scale")
  critvalue <- value*sx
  flagX <- 1*(abs(U-mx)>critvalue) # flag as outlier
  Xinit <- (1-flagX)*U + flagX*(Xest[,i.variable])
  
}

get_lambda_max <- function(xtilde, y, ytilde, x, xhat, betaEst, intercept, scaleVar, k, delta, maxIteration, tol, maxitscale,
                           maxituniv, post, value){
  
  # Get estimate of smallest sparsity parameter that sets everything equal to zero.
  
  # xtilde: nxp matrix
  # y: nx1 response 
  # ytilde: nx1 new response 
  # x: nxp matrix predictor matrix
  # xhat: nxp matrix expected predictors
  # betaEst: px1 matrix
  # intercept:px1 matrix
  # scaleVar:p-dimensional vector
  # k: tuning parameter for rho-function
  # delta: constant of consistency of M-scale (E[Z]=delta)
  # maxIteration: maximum number of iterations shooting loop
  # tol: numerical tolerance for convergence IRWLS
  # maxitscale: maximum number of iterations fixed point iterations
  # maxituniv : maximum number of iterations in reweighted least squares
  # post: 1 if post-lasso approach should be used, 0 else
  # value: cut-off value for outlier flagging  
  
  # Robust correlations
  cork <- apply(xtilde, 2, cor.fk, y = y)
  max <- which.max(abs(cork))
  lambda_init <- 30*abs(cork[max])
  n <- nrow(x)
  p <- ncol(x)
  
  lambda_max <- lambda_init
  lambda_last_max <- lambda_init
  last_frac <- 0.5
  count <- 1
  for(i in 1:10){
    wt <- matrix(NA, ncol=p, nrow=n)
    
    ft <- shootloop_sparse(ytilde=ytilde, xtilde=xtilde, x=x, wt=wt, Xexp=xhat,  
                           betaEst=betaEst,intercept=intercept, scaleVar=scaleVar, scaleVarOLD=scaleVar,
                           k=k, delta=delta, maxIteration=maxIteration,tol= mad(y)*10^-6,
                           tolout=tol, tolscale=10^-5, maxitscale=maxitscale, y=y, maxituniv=maxituniv, lambda = lambda_max, 
                           post = post, wcut = value) 
    if(sum(ft$betaEst == 0) == p){
      lambda_last_max = lambda_max
      count = 1
      last_frac = 0.5
      lambda_max = (0.5)*lambda_max
    } else{
      if(i ==1){
        lambda_max = lambda_last_max = 10*lambda_max
      } else{
        count = count + 1
        last_frac = last_frac + 1/2^count
        lambda_max = (last_frac)*lambda_last_max
        if(last_frac >=1){
          break
        }
      }
    }
  }
  return(lambda_last_max)
}


getLambdaGrid <- function(lambda_max, n,p, size){
  # Ger sparsity grid
  
  lmin = lambda_max/10^5
  
  grid <- c(exp(seq(log(lambda_max),log(lmin),length = size)))
  if(p <n){
    grid = c(grid,0)
  }
  return(grid) 
}

selectModel <- function(lambda_grid,ytilde,xtilde,xhat,x,wt,betaEst,intercept,scaleVar,k,delta,maxIteration,y,maxitscale,
                        maxituniv, tol, post, value){
  # Sparse shooting S estimates along the sparsity grid
  
  # lambda_grid : vector of sparsity parameters
  # ytilde: nx1 new response 
  # xtilde: nxp matrix
  # xhat: nxp matrix expected predictors
  # x: nxp matrix predictor matrix
  # betaEst: px1 matrix
  # intercept:px1 matrix
  # scaleVar:p-dimensional vector
  # k: tuning parameter for rho-function
  # delta: constant of consistency of M-scale (E[Z]=delta)
  # maxIteration: maximum number of iterations shooting loop
  # y: nx1 response 
  # maxitscale: maximum number of iterations fixed point iterations
  # maxituniv : maximum number of iterations in reweighted least squares
  # tol: numerical tolerance for convergence IRWLS
  # post: 1 if post-lasso approach should be used, 0 else
  # value: cut-off value for outlier flagging  
  
  n <- length(ytilde)
  p <- ncol(xtilde)
  
  get_coefs <- shootloop_sparse_lambdas(ytilde = ytilde, xtilde = xtilde, x = x, wt = wt, Xexp = xhat,  
                                        betaEst = betaEst, intercept = intercept, scaleVar = scaleVar, scaleVarOLD = scaleVar,
                                        k = k, delta = delta, maxIteration = maxIteration, tol = mad(y)*10^-6,
                                        tolout = tol, tolscale = 10^-5, maxitscale = maxitscale, y = y, maxituniv = maxituniv, 
                                        lambda = lambda_grid, post = post,  wcut = value)
  
  dfs <- unlist(lapply(get_coefs$betaEst, function(U){length(which(U!=0))}))
  sigmahat <- unlist(get_coefs$sigmahat)
  bic_sigmas <- sigmahat^2 + dfs*(log(n)/n)
  bic_ln_sigmas <- log(sigmahat^2) + dfs*(log(n)/n)
  minimum_sigmas <- which.min(bic_sigmas)
  minimum_ln_sigmas <- which.min(bic_ln_sigmas)
  
  out <- list("get_coefs" = get_coefs, "lambdas" = lambda_grid,
              "coef_opt" = get_coefs$betaEst[[minimum_sigmas]], "lambda_opt" = lambda_grid[minimum_sigmas], "BIC_opt" = bic_sigmas[minimum_sigmas],
              "weights_opt" = get_coefs$weights[[minimum_sigmas]], "alpha_opt" = get_coefs$alpha[[minimum_sigmas]], "iter_opt" = get_coefs$iter[[minimum_sigmas]],
              
              "coef_opt_ln" = get_coefs$betaEst[[minimum_ln_sigmas]], "lambda_opt_ln" = lambda_grid[minimum_ln_sigmas], "BIC_opt_ln" = bic_sigmas[minimum_ln_sigmas],
              "weights_opt_ln" = get_coefs$weights[[minimum_ln_sigmas]], "alpha_opt_ln" = get_coefs$alpha[[minimum_ln_sigmas]], "iter_opt_ln" = get_coefs$iter[[minimum_ln_sigmas]]
  )
}

startvalue_MM <- function(X, Y, value, robcor_fit = NULL, Xinit = NULL, predset = NULL){
  # Function to obtain starting value
  # X : nxp matrix of regressors
  # Y: n-dimensional response vector
  # value : cut-off value for outlier flagging 
  # robcor_fit : robust correlations
  # Xinit : Xclean 
  
  # 0. Initialization # 
  p <- ncol(X)
  n <- nrow(X)
  betaEst <- rep(0,p) 
  
  # 1. Select initial predictor set. Try out: Most robustly correlated one # 
  if(is.null(predset)){
    kpred <- min(round(n/2),p) # Set of predictors
    if(is.null(robcor_fit)){
      robcor_fit <- robcorr(X = X, Y = Y)
      
    }
    robcor <- robcor_fit$robcor
    predcor <- robcor_fit$predcor
    predset <- predcor[1:kpred]
  }

  
  # MM-FIT
  fitinit <- suppressWarnings(lmrob(Y~., data = cbind(Y, data.frame(Xinit[, predset])), setting = "KS2011"))
  betaEst[predset] <- fitinit$coef[-1]
  intercept <- rep.int(fitinit$coef[1], p)
  scaleVar <- rep.int(fitinit$scale, p)
  
  out <- list("betaEst" = as.matrix(betaEst), "intercept" = as.matrix(intercept), "scaleVar" = scaleVar, "Xinit" = Xinit)
  
}