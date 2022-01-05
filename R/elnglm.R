#' Data generator under the framework of generalized linear model
#'
#' @param n integer. The number of observations.
#' @param d integer. The number of predictors.
#' @param family string. The value of starting inertia weight in PSO updating procedure. The default is 1.2.
#' @param trueb0 double.
#' @param trueb vector. 
#' @param s double. 
#' @param seed integer. The random seed.
#' @return An List.
#' \itemize{
#' \item{x}{ the matrix of the predictors.}
#' \item{y}{ the vector of the response variable. 
#' For \code{family = 'multinomial'}, the output is a matrix of size $n\times$#(categories).}
#' }
#' @examples
#' # Intercept
#' trueb0 <- 1
#' # Regression Coefficients (the first 3 are active)
#' trueact <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
#' trueb <- runif(10, -1, 1)*10
#' trueb[which(trueact == 0)] <- 0 
#' 
#' # Generate data of continuous response
#' df <- glmDataGen(n = 500, d = 10, family = "gaussian", trueb0, trueb, s = 0.5, seed = 1)
#' 
#' # Generate data of binary response
#' dfb <- glmDataGen(n = 500, d = 10, family = "binomial", trueb0, trueb, seed = 1)
#' 
#' # Generate data of multi-categorical response (not run)
#' #
#' # Intercepts of the three log-linear models 
#' # trueb0 <- c(1, 1, 1)
#' # Regression Coefficients of the three log-linear models 
#' # trueact <- cbind(
#' #   c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0),
#' #   c(0, 0, 0, 1, 1, 1, 1, 0, 0, 0),
#' #   c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0)
#' # )
#' # trueb <- matrix(runif(10*3, -1, 1)*10, 10, 3)
#' # for (m in 1:3) { trueb[which(trueact[,m] == 0),m] <- 0 }
#' # 
#' # Generate data. The response y is a matrix of 3 columns.
#' # dfm <- glmDataGen(n = 500, d = 10, family = "multinomial", trueb0, trueb, seed = 1)
#' 
#' @name glmDataGen
#' @rdname glmDataGen
#' @export
glmDataGen <- function(n, d, family = c("gaussian", "binomial", "multinomial"), trueb0, trueb, s = 0.5, seed = NULL) {

  set.seed(seed)
  x <- matrix(rnorm(n*d), n, d)
  if (family == "gaussian") {
    y <- trueb0 + x %*% trueb + rnorm(n, 0, s)
  } else if (family == "binomial") {
    eta <- trueb0 + x %*% trueb; yp <- exp(eta)/(1 + exp(eta))
    y <- rbinom(n, 1, yp)
  } else if (family == "multinomial") {
    eta <- trueb0 + x %*% trueb
    yp <- exp(eta)/matrix(rowSums(exp(eta)), n, ncol(trueb))
    y <- matrix(0, n, ncol(trueb))
    for (i in 1:n) {
      y[i,] <- rmultinom(1, 1, yp[i,])
    }
  }
  return(list(x = x, y = y))
}

#' Run cross-validation for the penalty parameter of the generalized linear model.
#' 
#' @param y the vector of the response variable.
#' @param x the matrix of the predictors.
#' @param family string. One of the response families, "gaussian", "binomial" or "multinomial".
#' @param lambdaLength integer. The number of tuning penalty parameters. The default is 100.
#' @param minLambdaRatio double. The ratio of the minimal value to the maximal value 
#' of the penalty parameter. The default is \code{1e-3}.
#' @param lambdaVec vector. The optional input of the tuning penalty parameters.
#' The default is \code{NULL} that the function automatically 
#' computes the maximal value of the penalty parameter and generates a sequence of 
#' penalty parameter values of length \code{lambdaLength}.
#' @param alpha double. The elastic net parameter between 0 and 1. The default value is 0.5.
#' @param standardize boolean. If \code{TRUE}, the function first standardizes the predictor matrix.
#' @param maxit integer. The number of maximal iterations of the coordinate descent algorithm.
#' The default is 100.
#' @param tol double. The value of the convergence tolerance of the coordinate descent algorithm.
#' The default is \code{1e-4}.
#' @param nfolds integer. The number of folds. The default value is 3.
#' @param ver string. The version of the coordinate descent engine, "r": R codes or
#'  "arma": C++ codes with armadillo library. 
#' @return An List.
#' \itemize{  
#' \item{nfolds}{ the matrix of the predictors.}
#' \item{foldid}{ the vector of the response variable. }
#' \item{lambdaBestId}{ the vector of the response variable. }
#' \item{cvscore}{ the vector of the response variable. }
#' }
#' @examples
#' # Generate data of continuous response
#' trueb0 <- 1
#' trueact <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
#' trueb <- runif(10, -1, 1)*10
#' trueb[which(trueact == 0)] <- 0 
#' df <- glmDataGen(n = 500, d = 10, family = "gaussian", trueb0, trueb, s = 0.5, seed = 1)
#' 
#' # Run cross-validation
#' mdlcv <- glmPenaltyCV(y = df$y, x = df$x, family = "gaussian", lambdaLength = 200,
#'                       minLambdaRatio = 1e-3, maxit = 1e5, tol = 1e-7, alpha = 0.5, nfolds = 10, ver = "arma")
#' # Best Lambda value
#' mdlcv$lambda[mdlcv$lambdaBestId]
#' plot(log(mdlcv$lambda), mdlcv$cvscore, type = "l", xlab = "log(lambda)", ylab = "rmse")
#' # Estimated intercept of the best model
#' print(mdlcv$b0[mdlcv$lambdaBestId])
#' # Estimated coefficients of the best model
#' print(mdlcv$b[,mdlcv$lambdaBestId])
#' 
#' @name glmPenaltyCV
#' @rdname glmPenaltyCV
#' @export
glmPenaltyCV <- function(
  y, x, family = c("gaussian", "binomial", "multinomial"), lambdaLength = 100, minLambdaRatio = 1e-3, lambdaVec = NULL,
  alpha = 0.5, standardize = TRUE, maxit = 100, tol = 1e-4, nfolds = 3, ver = c("r", "arma"))
{
  if (is.vector(x)) x <- as.matrix(x, length(x), 1)
  n <- nrow(x)
  if (is.vector(y)) y <- as.matrix(y, n, 1)
  if (!(family %in% c("gaussian", "binomial", "multinomial"))) {
    stop("Illegal input of family.")
  }
  if (!(ver %in% c("r", "arma"))) {
    message("Illegal input of ver. Use ver = 'arma' instead.")
    ver <- "arma"
  }
  if (lambdaLength <= 0) {
    stop("lambdaLength must be positive integer.")
  }
  if ((minLambdaRatio <= 0) | (minLambdaRatio >= 1)) {
    stop("minLambdaRatio should be in (0, 1).")
  }
  if ((alpha <= 0) | (alpha >= 1)) {
    stop("alpha should be in (0, 1).")
  }
  #
  #
  if (standardize) {
    xm <- colMeans(x)
    xs <- sapply(1:ncol(x), function(j) sd(x[,j]))
    xs[which(xs == 0)] <- 1
    xn <- (x - matrix(xm, nrow(x), ncol(x), byrow = TRUE))/xs
  } else {
    xn <- x
  }
  #
  if (is.null(lambdaVec)) {
    lambdaMax <- calcMaxLambda(family, x, y, alpha)
    lambdaVec <- exp(seq(log(lambdaMax), log(lambdaMax*minLambdaRatio), length = lambdaLength))
  }

  mdl <- glmPenaltyFit(y, x, family = family, lambdaVec = lambdaVec,
                       alpha = alpha, standardize = FALSE, maxit = maxit, tol = tol, ver = ver)

  if (nfolds > 1) {
    foldid <- (rep(1:nfolds, ceiling(n/nfolds))[sample(1:(ceiling(n/nfolds)*nfolds))])[1:n]
  } else {
    foldid <- rep(1, n)
  }

  if (nfolds > 1) {
    if (family %in% c("gaussian", "binomial")) {
      yp <- matrix(0, n, lambdaLength)
    } else if (family == "multinomial") {
      yp <- array(0, c(n, ncol(y), lambdaLength))
    }

    for (ifold in 1:nfolds) {
      trainid <- which(foldid != ifold)
      testid <- which(foldid == ifold)
      x0 <- xn[trainid,]
      y0 <- y[trainid,]
      x1 <- xn[testid,]
      cvmdl <- glmPenaltyFit(y0, x0, family = family, lambdaVec = lambdaVec,
                             alpha = alpha, standardize = FALSE, maxit = maxit, tol = tol, ver = ver)
      
      cvpred <- glmPenaltyPred(cvmdl, x1, type = "response")
      if (family %in% c("gaussian", "binomial")) {
        yp[testid,] <- cvpred
      } else if (family == "multinomial") {
        yp[testid,,] <- cvpred
      }
    }

    if (family == "gaussian") {
      # Compute CV-RMSE
      cvscore <- sqrt(colMeans((matrix(y, nrow(x), length(lambdaVec)) - yp)^2))
    } else if (family == "binomial") {
      # Compute -log(acc)
      cvscore <- sapply(1:length(lambdaVec), function(l) {
        -log(sum(diag(confusion(y, yp[,l])))/n)
      })
    } else if (family == "multinomial") {
      # Compute -log(acc)
      ycate <- max.col(y) - 1
      cvscore <- sapply(1:length(lambdaVec), function(l) {
        ypl <- max.col(yp[,,l]) - 1
        -log(sum(diag(confusion(ycate, ypl)))/n)
      })
    }
    lambdaBestId <- which.min(cvscore)
  } else {
    cvscore <- 0; lambdaBestId <- 1
  }
  mdl$nfolds <- nfolds
  mdl$foldid <- foldid
  mdl$lambdaBestId <- lambdaBestId
  mdl$cvscore <- cvscore
  return(mdl)
}


#' Fit the generalized linear model with a vector of penalty parameters.
#' 
#' @param y the vector of the response variable.
#' @param x the matrix of the predictors.
#' @param family string. One of the response families, "gaussian", "binomial" or "multinomial".
#' @param lambdaLength integer. The number of tuning penalty parameters. The default is 100.
#' @param minLambdaRatio double. The ratio of the minimal value to the maximal value 
#' of the penalty parameter. The default is \code{1e-3}.
#' @param lambdaVec vector. The optional input of the tuning penalty parameters.
#' The default is \code{NULL} that the function automatically 
#' computes the maximal value of the penalty parameter and generates a sequence of 
#' penalty parameter values of length \code{lambdaLength}.
#' @param alpha double. The elastic net parameter between 0 and 1. The default value is 0.5.
#' @param standardize boolean. If \code{TRUE}, the function first standardizes the predictor matrix.
#' @param maxit integer. The number of maximal iterations of the coordinate descent algorithm.
#' The default is 100.
#' @param tol double. The value of the convergence tolerance of the coordinate descent algorithm.
#' The default is \code{1e-4}.
#' @param ver string. The version of the coordinate descent engine, "r": R codes or
#' "arma": C++ codes with armadillo library. 
#' @examples
#' # Generate data of continuous response
#' trueb0 <- 1
#' trueact <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
#' trueb <- runif(10, -1, 1)*10
#' trueb[which(trueact == 0)] <- 0 
#' df <- glmDataGen(n = 500, d = 10, family = "gaussian", trueb0, trueb, s = 0.5, seed = 1)
#' 
#' # Run cross-validation
#' mdl <- glmPenaltyFit(y = df$y, x = df$x, family = "gaussian", lambdaLength = 100,
#'                      minLambdaRatio = 1e-3, maxit = 1e5, tol = 1e-7, alpha = 0.5, ver = "arma")
#' # Estimated intercept of the best model
#' print(mdl$b0)
#' # Estimated coefficients of the best model
#' print(mdl$b)
#' 
#' @importFrom Rcpp cppFunction sourceCpp
#' @useDynLib elnglm
#' @export
glmPenaltyFit <- function(
  y, x, family = c("gaussian", "binomial", "multinomial"), lambdaLength = 100L, minLambdaRatio = 1e-3, lambdaVec = NULL,
  alpha = 0.5, standardize = TRUE, maxit = 100L, tol = 1e-4, ver = c("r", "arma"))
{
  
  if (is.vector(x)) x <- as.matrix(x, length(x), 1)
  n <- nrow(x)
  if (is.vector(y)) y <- as.matrix(y, n, 1)
  if (!(family %in% c("gaussian", "binomial", "multinomial"))) {
    stop("Illegal input of family.")
  }
  if (!(ver %in% c("r", "arma"))) {
    message("Illegal input of ver. Use ver = 'arma' instead.")
    ver <- "arma"
  }
  if (lambdaLength <= 0) {
    stop("lambdaLength must be positive integer.")
  }
  if ((minLambdaRatio <= 0) | (minLambdaRatio >= 1)) {
    stop("minLambdaRatio should be in (0, 1).")
  }
  if ((alpha <= 0) | (alpha >= 1)) {
    stop("alpha should be in (0, 1).")
  }
  #
  if (standardize) {
    xm <- colMeans(x)
    xs <- sapply(1:ncol(x), function(j) sd(x[,j]))
    xs[which(xs == 0)] <- 1
    xn <- (x - matrix(xm, nrow(x), ncol(x), byrow = TRUE))/matrix(xs, nrow(x), ncol(x), byrow = TRUE)
  } else {
    xn <- x; xm <- NULL; xs <- NULL
  }

  if (is.null(lambdaVec)) {
    lambdaMax <- calcMaxLambda(family, x, y, alpha)
    lambdaVec <- exp(seq(log(lambdaMax), log(lambdaMax*minLambdaRatio), length = lambdaLength))
  }

  if (family == "gaussian") {
    if (ver == "r") {
      mdl <- gaussian_elastic(y, xn, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
    } else if (ver == "arma") {
      mdl <- gaussian_elastic_arma(y, xn, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
    }

  } else if (family == "binomial") {
    if (ver == "r") {
      mdl <- binomial_elastic(y, xn, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
    } else if (ver == "arma") {
      mdl <- binomial_elastic_arma(y, xn, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
    }
  } else if (family == "multinomial") {
    if (ver == "r") {
      mdl <- multinomial_elastic(y, xn, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
    } else if (ver == "arma") {
      mdl <- multinomial_elastic_arma(y, xn, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)
    }
  }

  mdl$family <- family
  mdl$standardize <- standardize
  mdl$xm <- xm
  mdl$xs <- xs
  class(mdl) <- "glmPenalty"
  return(mdl)
}



#' Predict the generalized linear model with a vector of penalty parameters.
#' 
#' @param object the vector of the response variable.
#' @param xnew the matrix of the predictors.
#' @param type string. One of the response families, "response", "probability" or "link".
#' @examples
#' # Generate data of continuous response
#' trueb0 <- 1
#' trueact <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
#' trueb <- runif(10, -1, 1)*10
#' trueb[which(trueact == 0)] <- 0
#' df <- glmDataGen(n = 500, d = 10, family = "gaussian", trueb0, trueb, s = 0.5, seed = 1)
#' 
#' # Run cross-validation
#' mdlcv <- glmPenaltyCV(y = df$y, x = df$x, family = "gaussian", lambdaLength = 200,
#'                       minLambdaRatio = 1e-3, maxit = 1e5, tol = 1e-7, alpha = 0.5, nfolds = 10, ver = "arma")
#' # Predict for new data
#' xnew <- matrix(rnorm(10), 1, 10)
#' yp <- glmPenaltyPred(mdlcv, xnew)
#' yp[,mdlcv$lambdaBestId]
#' 
#' @export
glmPenaltyPred <- function(object, xnew, type = c("response", "probability", "link"))
{
  #
  if (is.vector(xnew)) {
    xnew <- matrix(xnew, 1, length(xnew))
  }
  if (ncol(xnew) != length(object$xm)) {
    stop("Dimension of the data does not match.")
  }
  if (class(mdl) != "glmPenalty") {
    stop("Wrong class of the input object.")
  }
  #
  if (object$standardize) {
    xn <- (xnew - matrix(object$xm, nrow(xnew), ncol(xnew), byrow = TRUE))/matrix(object$xs, nrow(xnew), ncol(xnew), byrow = TRUE)
  } else {
    xn <- xnew
  }

  if (object$family == "gaussian") {
    pred <- gaussian_predict(object, xnew)
  } else if (object$family == "binomial") {
    pred <- binomial_predict(object, xnew, type)
  } else if (object$family == "multinomial") {
    pred <- multinomial_predict(object, xnew, type)
  }
  return(pred)
}




