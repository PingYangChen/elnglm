
#' @export
glmDataGen <- function(n, d, family, trueb0, trueb, s = 0.5, seed = NULL) {

  set.seed(seed)
  x <- matrix(rnorm(n*d), n, d)
  if (family == "gaussian") {
    y <- trueb0[1] + x %*% trueb[,1] + rnorm(n, 0, s)
  } else if (family == "binomial") {
    eta <- trueb0[1] + x %*% trueb[,1]; yp <- exp(eta)/(1 + exp(eta))
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

#' @export
glmPenaltyCV <- function(
  y, x, family, lambdaLength = 100, minLambdaRatio = 1e-3, lambdaVec = NULL,
  alpha = 0.5, standardize = TRUE, maxit = 100, tol = 1e-4, nfold = 3)
{
  n <- nrow(x)
  if (is.vector(y)) y <- as.matrix(y, n, 1)
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
                       alpha = alpha, standardize = FALSE, maxit = maxit, tol = tol)

  if (nfold > 1) {
    foldid <- (rep(1:nfold, ceiling(n/nfold))[sample(1:(ceiling(n/nfold)*nfold))])[1:n]
  } else {
    foldid <- rep(1, n)
  }

  if (nfold > 1) {
    if (family %in% c("gaussian", "binomial")) {
      yp <- matrix(0, n, lambdaLength)
    } else if (family == "multinomial") {
      yp <- array(0, c(n, ncol(y), lambdaLength))
    }

    for (ifold in 1:nfold) {
      trainid <- which(foldid != ifold)
      testid <- which(foldid == ifold)
      x0 <- xn[trainid,]
      y0 <- y[trainid,]
      x1 <- xn[testid,]
      cvmdl <- glmPenaltyFit(y0, x0, family = family, lambdaVec = lambdaVec,
                             alpha = alpha, standardize = FALSE, maxit = maxit, tol = tol)
      cvpred <- glmPenaltyPred(cvmdl, x1)
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
  mdl$nfold <- nfold
  mdl$foldid <- foldid
  mdl$lambdaBestId <- lambdaBestId
  mdl$cvscore <- cvscore
  return(mdl)
}

#' @export
glmPenaltyFit <- function(
  y, x, family, lambdaLength = 100, minLambdaRatio = 1e-3, lambdaVec = NULL,
  alpha = 0.5, standardize = TRUE, maxit = 100, tol = 1e-4)
{
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

    mdl <- gaussian_elastic(y, xn, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)

  } else if (family == "binomial") {

    mdl <- binomial_elastic(y, xn, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)

  } else if (family == "multinomial") {

    mdl <- multinomial_elastic(y, xn, lambdaVec, alpha = 0.5, maxit = 100, tol = 1e-4)

  }

  mdl$family <- family
  mdl$standardize <- standardize
  mdl$xm <- xm
  mdl$xs <- xs
  class(mdl) <- "glmPenalty"
  return(mdl)
}

#' @export
glmPenaltyPred <- function(object, xnew)
{
  if (is.vector(xnew)) {
    xnew <- matrix(xnew, 1, length(xnew))
  }
  if (object$standardize) {
    xn <- (xnew - matrix(object$xm, nrow(xnew), ncol(xnew), byrow = TRUE))/matrix(object$xs, nrow(x), ncol(x), byrow = TRUE)
  } else {
    xn <- xnew
  }

  if (object$family == "gaussian") {
    pred <- gaussian_predict(object, xnew)
  } else if (object$family == "binomial") {
    pred <- binomial_predict(object, xnew)
  } else if (object$family == "multinomial") {
    pred <- multinomial_predict(object, xnew)
  }
  return(pred)
}


# glmPenaltyPredCpp <- function(y, x, engine = "arma") {
#
# }
#

