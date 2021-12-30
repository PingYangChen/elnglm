


path <- "D:\\rProject\\elnglm\\src"

library(Rcpp)
library(RcppArmadillo)

sourceCpp(file.path(path, "elnglm.cpp"))



alpha = 0.5
maxit = 100
lambdaLength = 100
minLambdaRatio = 1e-3
lambdaVec = NULL
tol = 1e-4
nfold = 3
standardize = TRUE


trueb0 <- c(0,0,0)
trueact <- cbind(
  c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0),
  c(0, 0, 0, 1, 1, 1, 1, 0, 0, 0),
  c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0)
)
trueb <- matrix(runif(10*3, -1, 1)*10, 10, 3)
for (m in 1:3) { trueb[which(trueact[,m] == 0),m] <- 0 }
family <- "multinomial" # "gaussian" # "binomial" #      #

df <- glmDataGen(500, 10, family, trueb0, trueb, s = 0.5, NULL)
#y = df$y; x = df$x
out <- glmPenaltyCV(y = df$y, x = df$x, family, lambdaLength = 100,
                    alpha = 0.5, minLambdaRatio = 1e-3, nfolds = 5,
                    ver = "arma")

out$lambdaBestId
out$cvscore[out$lambdaBestId]
