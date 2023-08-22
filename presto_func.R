dir_main <- getwd()
dir_ordnet <- paste(dir_main, "/ordinalNet modified", sep="")
setwd(dir_ordnet)

source("cdIn.R")
source("cdOut.R")
source("links.R")
source("mirlsNet.R")
source("misc.R")
source("ordinalNet-methods.R")
source("ordinalNet.R")
source("ordinalNetCV.R")
source("ordinalNetTune.R")

setwd(dir_main)

# Function that selects a tuning parameter lambda and fits PRESTO
#
# X: A design matrix, as in the output of model.matrix(). (X should not contain
# an intercept column.)
# y: A response vector; should be an ordered factor, and should have the same
# length as the number of rows of X.
# printIter: Logical. If TRUE, the optimization routine progress is printed to
# the terminal. Default is TRUE.
# alpha_param: The elastic net mixing parameter, with 0 <= alpha_param <= 1.
# alpha_param=1 corresponds to the lasso penalty, and alpha_param=0 corresponds
# to the ridge penalty. Default is 1.
#
# Output: A list with two elements.
# alpha_hat: A numeric vector of length K - 1 containing the estimated
# intercepts (bias terms) for the K - 1 estimated decision boundaries.
# Beta_hat: A numeric p x (K - 1) matrix. Column k contains the coefficients for
# the kth decision boundary.
presto <- function(X, y, printIter=TRUE, alpha_param=1){

	n <- length(y)
	K <- length(unique(y))

	stopifnot("ordered" %in% class(y))
	stopifnot("factor" %in% class(y))
	stopifnot(K < n)
	if(K < 3){
		stop("Response variable must have at least three ordered classes.")
	}
	stopifnot(nrow(X) == n)

	# Use data to get lambda
	if(printInter){
		print("Getting optimal lambda...")
	}
	
	lambda_star <- getLambdaFusedPolr(X=X, y=y, printProgress=printIter)

	if(printInter){
		print("Done! Fitting model on selected lambda...")
	}

	ret <- getParamEstsMod(X=X, y=y, lambdaVals=lambda_star, K=K,
		printIter=printIter, alpha_param=alpha_param)

	if(printIter){
		print("Done!")
	}
	
	return(ret)
}

# Function that gets final estimated class probabilties using fitted PRESTO
# model
#
# res: The output from presto.
# X_test: An unlabeled set of covariates with the same structure as the X used
# to fit PRESTO.
# 
# Output: A n_test x K matrix of probabilities, where n_test is the number of
# row of X_test and K is the number of responses. X_test[i, j] contains the
# estimated probability that observation i will lie in class j.
predictPrestoProbs <- function(res, X_test){
	# Confirm that paramter estimates needed for calculating probability
	# estimates are available in method output, as expected

	stopifnot("alpha_hat" %in% names(res))
	stopifnot("Beta_hat" %in% names(res))

	alpha_hat <- res$alpha_hat
	Beta_hat <- res$Beta_hat

	stopifnot(all(!is.na(alpha_hat)))
	stopifnot(all(!is.na(Beta_hat)))

	K <- length(alpha_hat) + 1
	p <- nrow(Beta_hat)

	stopifnot(K >= 3)
	stopifnot(all(dim(Beta_hat) == c(p, K - 1)))

	n <- nrow(X_test)

	stopifnot(n >= 1)
	stopifnot(ncol(X_test) == p)

	p_hat <- matrix(rep(as.numeric(NA), n*K), nrow=n, ncol=K)

	if(p > 1){
		p_hat[, 1] <- plogis(alpha_hat[1] + X_test %*% Beta_hat[, 1])[, 1]
	} else{
		p_hat[, 1] <- plogis(alpha_hat[1] + Beta_hat[, 1]*X_test)
	}

	for(k in 2:(K - 1)){
		if(p > 1){
			p_hat[, k] <- plogis(alpha_hat[k] + X_test %*% Beta_hat[, k])[, 1] -
				plogis(alpha_hat[k - 1] + X_test %*% Beta_hat[, k - 1])[, 1]
		} else{
			p_hat[, k] <- plogis(alpha_hat[k] + Beta_hat[, k]*X_test) -
				plogis(alpha_hat[k - 1] + Beta_hat[, k - 1]*X_test)
		}
	}
	
	if(p > 1){
		p_hat[, K] <- 1 - plogis(alpha_hat[K - 1] +
			X_test %*% Beta_hat[, K - 1])[, 1]
	} else{
		p_hat[, K] <- 1 - plogis(alpha_hat[K - 1] +
			Beta_hat[, K - 1]*X_test)
	}

	# Verify predicted probabilities either way
	stopifnot(all(!is.na(p_hat)))
	stopifnot(is.numeric(p_hat) | is.integer(p_hat))
	stopifnot(nrow(p_hat) == n)
	stopifnot(ncol(p_hat) == K)
	stopifnot(all(p_hat >= 0))
	stopifnot(all(p_hat <= 1))

	return(p_hat)
}



#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#


# Other "under-the-hood" functions:

getLambdaFusedPolr <- function(X, y, printProgress=FALSE, nLambda=20,
	lambdaMinRatio=.01, alpha_param=1){
	# Check input

	stopifnot(nrow(X) == length(y))

	fit_tune <- ordinalNetTuneMod(X, y, family="cumulative", link="logit",
		printProgress=printProgress, nLambda=nLambda,
		lambdaMinRatio=lambdaMinRatio, alpha=alpha_param, warn=FALSE)

	briers <- rowMeans(fit_tune$brier)

	lambda_star <- max(fit_tune$fit$lambdaVals[briers == min(briers)])

	# Check output

	stopifnot(is.numeric(lambda_star) | is.integer(lambda_star))
	stopifnot(!is.na(lambda_star))
	stopifnot(length(lambda_star) == 1)
	stopifnot(lambda_star >= 0)

	return(lambda_star)
}

getParamEstsMod <- function(X, y, lambdaVals, K, printIter=FALSE, alpha_param=1){

	fit <- ordinalNetMod(X, y, family="cumulative", link="logit",
		printIter=printIter, lambdaVals=lambdaVals, alpha=alpha_param,
		warn=FALSE)

	p <- ncol(X)

	# Get estimates

	alpha_hat <- fit$coefs[1, 1:(K - 1)]
	names(alpha_hat) <- names(fit$coefs[1, 1:(K - 1)])

	psi <- fit$coefs[1, K:length(fit$coefs)]

	stopifnot(length(psi) == p*(K - 1))

	Beta_hat <- matrix(as.numeric(NA), nrow=p, ncol=K - 1)

	Beta_hat[, 1] <- psi[1:p]

	for(k in 2:(K - 1)){
	     Beta_hat[, k] <- Beta_hat[, k - 1] + psi[((k-1)*p + 1):(k*p)]
	}

	stopifnot(nrow(Beta_hat) == ncol(X))
	rownames(Beta_hat) <- colnames(X)
	colnames(Beta_hat) <- names(fit$coefs[1, 1:(K - 1)])

	stopifnot(all(!is.na(alpha_hat)))
	stopifnot(all(!is.na(Beta_hat)))


	return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat))
}

