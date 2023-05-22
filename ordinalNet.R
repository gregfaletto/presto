#' Ordinal regression models with elastic net penalty
#'
#' Fits ordinal regression models with elastic net penalty by coordinate descent.
#' Supported model families include cumulative probability, stopping ratio, continuation ratio,
#' and adjacent category. These families are a subset of vector glm's which belong to a model
#' class we call the elementwise link multinomial-ordinal (ELMO) class. Each family
#' in this class links a vector of covariates to a vector of class probabilities.
#' Each of these families has a parallel form, which is appropriate for ordinal response
#' data, as well as a nonparallel form that is appropriate for an unordered categorical
#' response, or as a more flexible model for ordinal data. The parallel model
#' has a single set of coefficients, whereas the nonparallel model has a set of coefficients
#' for each response category except the baseline category. It is also possible
#' to fit a model with both parallel and nonparallel terms, which we call the semi-parallel
#' model. The semi-parallel model has the flexibility of the nonparallel model,
#' but the elastic net penalty shrinks it toward the parallel model.
#'
#' @param x Covariate matrix. It is recommended that categorical covariates are
#' converted to a set of indicator variables with a variable for each category
#' (i.e. no baseline category); otherwise the choice of baseline category will
#' affect the model fit.
#' @param y Response variable. Can be a factor, ordered factor, or a matrix
#' where each row is a multinomial vector of counts. A weighted fit can be obtained
#' using the matrix option, since the row sums are essentially observation weights.
#' Non-integer matrix entries are allowed.
#' @param alpha The elastic net mixing parameter, with \code{0 <= alpha <= 1}.
#' \code{alpha=1} corresponds to the lasso penalty, and \code{alpha=0} corresponds
#' to the ridge penalty.
#' @param standardize If \code{standardize=TRUE}, the predictor variables are
#' scaled to have unit variance. Coefficient estimates are returned on the
#' original scale.
#' @param penaltyFactors Optional nonnegative vector of penalty factors with
#' length equal to the number of columns in \code{x}. If this argument is used,
#' then the penalty for each variable is scaled by its corresponding factor.
#' If \code{NULL}, the penalty factor is one for each coefficient.
#' @param positiveID Logical vector indicating whether each coefficient should
#' be constrained to be non-negative. If \code{NULL}, the default value is \code{FALSE}
#' for all coefficients.
#' @param family Specifies the type of model family. Options are "cumulative"
#' for cumulative probability, "sratio" for stopping ratio, "cratio" for continuation ratio,
#' and "acat" for adjacent category.
#' @param reverse Logical. If \code{TRUE}, then the "backward" form of the model
#' is fit, i.e. the model is defined with response categories in reverse order.
#' For example, the reverse cumulative model with \eqn{K+1} response categories
#' applies the link function to the cumulative probabilities \eqn{P(Y \ge 2),
#' \ldots, P(Y \ge K+1)}, rather then \eqn{P(Y \le 1), \ldots, P(Y \le K)}.
#' @param link Specifies the link function. The options supported are logit,
#' probit, complementary log-log, and cauchit. Only used if \code{customLink=NULL}.
#' @param customLink Optional list containing a vectorized link function \code{g},
#' a vectorized inverse link \code{h}, and the Jacobian function of the inverse link
#' \code{getQ}. The Jacobian should be defined as \eqn{\partial h(\eta) / \partial \eta^T}
#' (as opposed to the transpose of this matrix).
#' @param parallelTerms Logical. If \code{TRUE}, then parallel coefficient terms
#' will be included in the model. \code{parallelTerms} and \code{nonparallelTerms}
#' cannot both be \code{FALSE}.
#' @param nonparallelTerms Logical. if \code{TRUE}, then nonparallel coefficient terms
#' will be included in the model. \code{parallelTerms} and \code{nonparallelTerms}
#' cannot both be \code{FALSE}.
#' @param parallelPenaltyFactor Nonnegative numeric value equal to one by
#' default. The penalty on all parallel terms is scaled by this factor (as well
#' as variable-specific \code{penaltyFactors}). Only used if
#' \code{parallelTerms=TRUE}.
#' @param lambdaVals An optional user-specified lambda sequence (vector). If \code{NULL},
#' a sequence will be generated based on \code{nLambda} and \code{lambdaMinRatio}.
#' In this case, the maximum lambda is the smallest value that sets all penalized
#' coefficients to zero, and the minimum lambda is the maximum value multiplied
#' by the factor \code{lambdaMinRatio}.
#' @param nLambda Positive integer. The number of lambda values in the solution path.
#' Only used if \code{lambdaVals=NULL}.
#' @param lambdaMinRatio A factor greater than zero and less than one. Only used
#' if \code{lambdaVals=NULL}.
#' @param includeLambda0 Logical. If \code{TRUE}, then zero is added to the end
#' of the sequence of \code{lambdaVals}. This is not done by default because
#' it can significantly increase computational time. An unpenalized saturated model
#' may have infinite coefficient solutions, in which case the fitting algorithm
#' will still terminate when the relative change in log-likelihood becomes small.
#' Only used if \code{lambdaVals=NULL}.
#' @param alphaMin \code{max(alpha, alphaMin)} is used to calculate the starting
#' lambda value when \code{lambdaVals=NULL}. In this case, the default lambda
#' sequence begins with the smallest lambda value such that all penalized
#' coefficients are set to zero (i.e. the value where the first penalized
#' coefficient enters the solution path). The purpose of this argument is to
#' help select a starting value for the lambda sequence when \code{alpha = 0},
#' because otherwise it would be infinite. Note that \code{alphaMin} is only
#' used to determine the default lamba sequence and that the model is always fit
#' using \code{alpha} to calculate the penalty.
#' @param pMin Value greater than zero, but much less than one. During the optimization
#' routine, the Fisher information is calculated using fitted probabilities. For
#' this calculation, fitted probabilities are capped below by this value to prevent
#' numerical instability.
#' @param stopThresh In the relative log-likelihood change between successive
#' lambda values falls below this threshold, then the last model fit is used for all
#' remaining lambda.
#' @param threshOut Convergence threshold for the coordinate descent outer loop.
#' The optimization routine terminates when the relative change in the
#' penalized log-likelihood between successive iterations falls below this threshold.
#' It is recommended to set \code{theshOut} equal to \code{threshIn}.
#' @param threshIn Convergence threshold for the coordinate descent inner loop. Each
#' iteration consists of a single loop through each coefficient. The inner
#' loop terminates when the relative change in the penalized approximate
#' log-likelihood between successive iterations falls below this threshold.
#' It is recommended to set \code{theshOut} equal to \code{threshIn}.
#' @param maxiterOut Maximum number of outer loop iterations.
#' @param maxiterIn Maximum number of inner loop iterations.
#' @param printIter Logical. If \code{TRUE}, the optimization routine progress is
#' printed to the terminal.
#' @param printBeta Logical. If \code{TRUE}, coefficient estimates are printed
#' after each coordinate descent outer loop iteration.
#' @param warn Logical. If \code{TRUE}, the following warning message is displayed
#' when fitting a cumulative probability model with \code{nonparallelTerms=TRUE}
#' (i.e. nonparallel or semi-parallel model).
#' "Warning message: For out-of-sample data, the cumulative probability model
#' with nonparallelTerms=TRUE may predict cumulative probabilities that are not
#' monotone increasing."
#' The warning is displayed by default, but the user may wish to disable it.
#' @param keepTrainingData Logical. If \code{TRUE}, then \code{x} and \code{y}
#' are saved with the returned "ordinalNet" object. This allows
#' \code{predict.ordinalNet} to return fitted values for the training data
#' without passing a \code{newx} argument.

#' @details
#' The \code{ordinalNet} function fits regression models for a categorical response
#' variable with \eqn{K+1} levels. Conditional on the covariate vector \eqn{x_i}
#' (the \eqn{i^{th}} row of \code{x}), each observation has a vector of \eqn{K+1}
#' class probabilities \eqn{(p_{i1}, \ldots, p_{i(K+1)})}. These probabilities
#' sum to one, and can therefore be parametrized by \eqn{p_i = (p_{i1}, \ldots, p_{iK})}.
#' The probabilities are mapped to a set of \eqn{K} quantities
#' \eqn{\delta_i = (\delta_{i1}, \ldots, \delta_{iK}) \in (0, 1)^K}, which depends on the choice
#' of model \code{family}. The elementwise \code{link} function maps
#' \eqn{\delta_i} to a set of \eqn{K} linear predictors. Together, the \code{family}
#' and \code{link} specifiy a link function between \eqn{p_i} and \eqn{\eta_i}.
#'
#' \strong{\emph{Model families:}}
#'
#' Let \eqn{Y} denote the random response variable for a single observation,
#' conditional on the covariates values of the observation. The random variable
#' \eqn{Y} is discrete with support \{\eqn{1, \ldots, K+1}\}. The following model
#' families are defined according to these mappings between the class
#' probabilities and the values \eqn{\delta_1, \ldots, \delta_K}:
#' \describe{
#'   \item{Cumulative probability}{\eqn{\delta_j = P(Y \le j)}}
#'   \item{Reverse cumulative probability}{\eqn{\delta_j = P(Y \ge j + 1)}}
#'   \item{Stopping ratio}{\eqn{\delta_j = P(Y = j | Y \ge j)}}
#'   \item{Reverse stopping ratio}{\eqn{\delta_j = P(Y=j + 1 | Y \le j + 1)}}
#'   \item{Continuation ratio}{\eqn{\delta_j = P(Y > j | Y \ge j)}}
#'   \item{Reverse continuation ratio}{\eqn{\delta_j = P(Y < j | Y \le j)}}
#'   \item{Adjacent category}{\eqn{\delta_j = P(Y = j + 1 | j \le Y \le j+1)}}
#'   \item{Reverse adjacent category}{\eqn{\delta_j = P(Y = j | j \le Y \le j+1)}}
#' }
#'
#' \strong{\emph{Parallel, nonparallel, and semi-parallel model forms:}}
#'
#' Models within each of these families can take one of three forms, which have
#' different definitions for the linear predictor \eqn{\eta_i}. Suppose each
#' \eqn{x_i} has length \eqn{P}. Let \eqn{b} be a length \eqn{P} vector of
#' regression coefficients. Let \eqn{B} be a \eqn{P \times K} matrix of regression
#' coefficient. Let \eqn{b_0} be a vector of \eqn{K} intercept terms.
#' The three model forms are the following:
#' \describe{
#'   \item{Parallel}{\eqn{\eta_i = b_0 + b^T x_i} (\code{parallelTerms=TRUE}, \code{nonparallelTerms=FALSE})}
#'   \item{Nonparallel}{\eqn{\eta_i = b_0 + B^T x_i} (\code{parallelTerms=FALSE}, \code{nonparallelTerms=TRUE})}
#'   \item{Semi-parallel}{\eqn{\eta_i = b_0 + b^T x_i + B^T x_i} (\code{parallelTerms=TRUE}, \code{nonparallelTerms=TRUE})}
#' }
#' The parallel form has the defining property of ordinal models, which is that
#' a single linear combination \eqn{b^T x_i} shifts the cumulative class probabilities
#' \eqn{P(Y \le j)} in favor of either higher or lower categories. The linear predictors
#' are parallel because they only differ by their intercepts (\eqn{b_0}). The nonparallel form
#' is a more flexible model, and it does not shift the cumulative probabilities together.
#' The semi-parallel model is equivalent to the nonparallel model, but the
#' elastic net penalty shrinks the semi-parallel coefficients toward a common
#' value (i.e. the parallel model), as well as shrinking all coefficients toward zero.
#' The nonparallel model, on the other hand, simply shrinks all coefficients toward zero.
#' When the response categories are ordinal, any of the three model forms could
#' be applied. When the response categories are unordered, only the nonparallel
#' model is appropriate.
#'
#' \strong{\emph{Elastic net penalty:}}
#'
#' The elastic net penalty is defined for each model form as follows. \eqn{\lambda}
#' and \eqn{\alpha} are the usual elastic net tuning parameters, where \eqn{\lambda}
#' determines the degree to which coefficients are shrunk toward zero, and \eqn{\alpha}
#' specifies the amound of weight given to the L1 norm and squared L2 norm penalties.
#' Each covariate is allowed a unique penalty factor \eqn{c_j}, which is specified with the
#' \code{penaltyFactors} argument. By default \eqn{c_j = 1} for all \eqn{j}.
#' The semi-parallel model has a tuning parameter \eqn{\rho} which determines the degree to
#' which the parallel coefficients are penalized. Small values of \eqn{\rho} will
#' result in a fit closer to the parallel model, and large values of \eqn{\rho}
#' will result in a fit closer to the nonparallel model.
#' \describe{
#'   \item{Parallel}{\eqn{\lambda \sum_{j=1}^P c_j \{ \alpha |b_j| +
#'                        \frac{1}{2} (1-\alpha) b_j^2 \}}}
#'   \item{Nonparallel}{\eqn{\lambda \sum_{j=1}^P c_j \{ \sum_{k=1}^K \alpha |B_{jk}| +
#'                           \frac{1}{2} (1-\alpha) B_{jk}^2 \}}}
#'   \item{Semi-parallel}{\eqn{\lambda [ \rho \sum_{j=1}^P c_j \{ \alpha |b_j| +
#'                             \frac{1}{2} (1-\alpha) b_j^2 \} +
#'                             \sum_{j=1}^P c_j \{ \sum_{k=1}^K \alpha |B_{jk}| +
#'                             \frac{1}{2} (1-\alpha) B_{jk}^2 \}]}}
#' }
#'
#' \code{ordinalNet} minimizes the following objective function. Let \eqn{N} be
#' the number of observations, which is defined as the sum of the \code{y} elements
#' when \code{y} is a matrix.
#' \deqn{objective = -1/N*loglik + penalty}
#'
#' @return An object with S3 class "ordinalNet".  Model fit information can be accessed
#' through the \code{coef}, \code{predict}, and \code{summary} methods.
#' \describe{
#'   \item{coefs}{Matrix of coefficient estimates, with each row corresponding to a lambda value.
#'   (If covariates were scaled with \code{standardize=TRUE}, the coefficients are
#'   returned on the original scale).}
#'   \item{lambdaVals}{Sequence of lambda values. If user passed a sequence to the
#'   \code{lambdaVals}, then it is this sequence. If \code{lambdaVals} argument
#'   was \code{NULL}, then it is the sequence generated.}
#'   \item{loglik}{Log-likelihood of each model fit.}
#'   \item{nNonzero}{Number of nonzero coefficients of each model fit, including intercepts.}
#'   \item{aic}{AIC, defined as \code{-2*loglik + 2*nNonzero}.}
#'   \item{bic}{BIC, defined as \code{-2*loglik + log(N)*nNonzero}.}
#'   \item{devPct}{Percentage deviance explained, defined as \eqn{1 - loglik/loglik_0},
#'   where \eqn{loglik_0} is the log-likelihood of the null model.}
#'   \item{iterOut}{Number of coordinate descent outer loop iterations until
#'   convergence for each lambda value.}
#'   \item{iterIn}{Number of coordinate descent inner loop iterations on last outer loop
#'   for each lambda value.}
#'   \item{dif}{Relative improvement in objective function on last outer loop
#'   for each lambda value. Can be used to diagnose convergence issues. If \code{iterOut}
#'   reached \code{maxiterOut} and \code{dif} is large, then \code{maxiterOut} should
#'   be increased. If \code{dif} is negative, this means the objective did not improve
#'   between successive iterations. This usually only occurs when the model is
#'   saturated and/or close to convergence, so a small negative value is not of concern.
#'   (When this happens, the algorithm is terminated for the current lambda value,
#'   and the coefficient estimates from the previous outer loop iteration are returned.)}
#'   \item{nLev}{Number of response categories.}
#'   \item{nVar}{Number of covariates in \code{x}.}
#'   \item{xNames}{Covariate names.}
#'   \item{args}{List of arguments passed to the \code{ordinalNet} function.}
#' }
#'
#' @examples
#' # Simulate x as independent standard normal
#' # Simulate y|x from a parallel cumulative logit (proportional odds) model
#' set.seed(1)
#' n <- 50
#' intercepts <- c(-1, 1)
#' beta <- c(1, 1, 0, 0, 0)
#' ncat <- length(intercepts) + 1  # number of response categories
#' p <- length(beta)  # number of covariates
#' x <- matrix(rnorm(n*p), ncol=p)  # n x p covariate matrix
#' eta <- c(x %*% beta) + matrix(intercepts, nrow=n, ncol=ncat-1, byrow=TRUE)
#' invlogit <- function(x) 1 / (1+exp(-x))
#' cumprob <- t(apply(eta, 1, invlogit))
#' prob <- cbind(cumprob, 1) - cbind(0, cumprob)
#' yint <- apply(prob, 1, function(p) sample(1:ncat, size=1, prob=p))
#' y <- as.factor(yint)
#'
#' # Fit parallel cumulative logit model
#' fit1 <- ordinalNet(x, y, family="cumulative", link="logit",
#'                    parallelTerms=TRUE, nonparallelTerms=FALSE)
#' summary(fit1)
#' coef(fit1)
#' coef(fit1, matrix=TRUE)
#' predict(fit1, type="response")
#' predict(fit1, type="class")
#'
#' # Fit nonparallel cumulative logit model
#' fit2 <- ordinalNet(x, y, family="cumulative", link="logit",
#'                    parallelTerms=FALSE, nonparallelTerms=TRUE)
#' fit2
#' coef(fit2)
#' coef(fit2, matrix=TRUE)
#' predict(fit2, type="response")
#' predict(fit2, type="class")
#'
#' # Fit semi-parallel cumulative logit model (with both parallel and nonparallel terms)
#' fit3 <- ordinalNet(x, y, family="cumulative", link="logit",
#'                    parallelTerms=TRUE, nonparallelTerms=TRUE)
#' fit3
#' coef(fit3)
#' coef(fit3, matrix=TRUE)
#' predict(fit3, type="response")
#' predict(fit3, type="class")
#'
#' @export
ordinalNetMod <- function(x, y, alpha=1, standardize=TRUE, penaltyFactors=NULL, positiveID=NULL,
                       family=c("cumulative", "sratio", "cratio", "acat"), reverse=FALSE,
                       link=c("logit", "probit", "cloglog", "cauchit"), customLink=NULL,
                       parallelTerms=FALSE, nonparallelTerms=TRUE, parallelPenaltyFactor=1,
                       lambdaVals=NULL, nLambda=20, lambdaMinRatio=0.01, includeLambda0=FALSE, alphaMin=0.01,
                       pMin=1e-8, stopThresh=1e-8, threshOut=1e-8, threshIn=1e-8, maxiterOut=100, maxiterIn=100,
                       printIter=FALSE, printBeta=FALSE, warn=TRUE, keepTrainingData=TRUE)
{
    family <- match.arg(family)
    link <- match.arg(link)
    args <- as.list(environment())  # list of arguments to return
    if (!keepTrainingData) args$x <- args$y <- NULL

    # Initial argument checks
    if (!is.matrix(x))
        stop("x should be a matrix.")
    if (any(is.na(x)))
        stop("x must not contain missing values.")
    if (any(abs(x) == Inf))
        stop("x must not contain infinite values.")
    if (!is.factor(y) && !is.matrix(y))
        stop("y should be a factor or matrix.")
    if (any(is.na(y)))
        stop("y must not contain missing values.")

    # Variable definitions
    yMat <- if (is.matrix(y)) y else yFactorToMatrix(y)
    wts <- attr(yMat, "wts") <- rowSums(yMat)
    wtsum <- attr(yMat, "wtsum") <- sum(wts)
    nVar <- ncol(x)
    nLev <- ncol(yMat)
    if (reverse) yMat <- yMat[, nLev:1]

    # Other argument checks
    if (nrow(x) != nrow(yMat))
        stop("x and y dimensions do not match.")
    if (alpha<0 || alpha>1)
        stop("alpha should be a number such that 0 <= alpha <= 1.")
    if (!is.null(penaltyFactors) && length(penaltyFactors)!=nVar)
        stop(paste0("penaltyFactors should be a numeric vector of length equal to the number ",
                    "of variables in x. Set penaltyFactor=NULL to penalize each variable equally."))
    if (!is.null(penaltyFactors) && any(penaltyFactors < 0))
        stop("penaltyFactors values should be nonnegative.")
    if (!is.null(positiveID) && length(positiveID)!=nVar)
        stop(paste0("positiveID should be a logical vector of length equal to the number ",
                    "of variables in x. Set positiveID=NULL for no positive restrictions."))
    if (!is.null(positiveID) && any(!is.logical(positiveID)))
        stop("positiveID values should be logical.")
    if (!is.null(lambdaVals) && any(lambdaVals<0))
        stop("lambdaVals values should be nonnegative.")
    if (is.null(lambdaVals) && nLambda < 1)
        stop("nLambda should be >= 1.")
    if (is.null(lambdaVals) && lambdaMinRatio<=0)
        stop("lambdaMinRatio should be strictly greater than zero.")
    if (alpha<alphaMin && alphaMin<=0)
        stop("alphaMin should be strictly greater than zero.")
    if (length(parallelPenaltyFactor) > 1)
        stop("parallelPenaltyFactor should be a single value.")
    if (parallelTerms && parallelPenaltyFactor<0)
        stop("parallelPenaltyFactor should be >= 0.")
    if (!parallelTerms && parallelPenaltyFactor!=1)
        warning("parallelPenaltyFactor is not used when parallelTerms=FALSE.")
    if (!parallelTerms && !nonparallelTerms)
        stop("parallelTerms and nonparallelTerms cannot both be FALSE.")
    if (warn && family=="cumulative" && nonparallelTerms) {
        warning(paste0("For out-of-sample data, the cumulative probability model with ",
                       "nonparallelTerms=TRUE may predict cumulative probabilities that are not ",
                       "monotone increasing."))
    }
    if (!is.null(customLink)) {
        link <- customLink
        message("customLink should be a list containing:\n
                \ \ $linkfun := vectorized link function with domain (0, 1)\n
                \ \ $linkinv := vectorized inverse link with domain (-Inf, Inf)\n
                \ \ $mu.eta  := vectorized Jacobian of linkinv\n
                The customLink argument is not checked, so user should be cautious
                using it.")
    }
    if ((parallelPenaltyFactor != 1) && !parallelTerms)
        warning("parallelPenaltyFactor is not used when parallelTerms = FALSE.")

    #
    #
    #
    #
    #
    #
    # MODIFICATIONS FOR OUR IMPLEMENTATION OF PRESTO
    #
    #
    #
    #
    if(!nonparallelTerms){
        stop("Only non-parallel model currently supported (nonparallelTerms must be TRUE)")
    }
    if(parallelTerms){
        stop("Only non-parallel model currently supported (parallelTerms must be FALSE)")
    }
    if(nLev < 3){
    	stop("Minimum of three levels in response supported")
    }

    # Create linkfun
    linkfun <- makeLinkfun(family, link)

    # Center x and create xList; also scale x if standardize=TRUE
    xMeans <- colMeans(x)
    if (standardize) {
        xSD <- sqrt(rowSums(wts*(t(x)-xMeans)^2) / wtsum)
        xSD[xSD==0] <- 1
        xStd <- t((t(x) - xMeans) / xSD)
    } else {
        xStd <- t(t(x) - xMeans)
    }
    xList <- makexList(xStd, nLev, parallelTerms, nonparallelTerms)

    # Augment penaltyFactors to include all model coefficients
    if (is.null(penaltyFactors)) penaltyFactors <- rep(1, nVar)
    penaltyFactorsIntercept <- rep(0, nLev-1)
    penaltyFactorsParallel <- if (parallelTerms) penaltyFactors * parallelPenaltyFactor else NULL
    #
    #
    #
    #
    #
    #
    # MODIFICATIONS FOR OUR IMPLEMENTATION OF PRESTO
    #
    #
    #
    #
    # BELOW MODIFIED AGAIN FOR NEW VERSION OF PRESTO
    # (COMMENTED OUT PREVIOUS VERSION)
    # penaltyFactorsNonparallel <- if(nonparallelTerms) c(rep(0, nVar), rep(penaltyFactors, nLev-2)) else NULL
    penaltyFactorsNonparallel <- if(nonparallelTerms) c(rep(penaltyFactors, nLev-1)) else NULL
    stopifnot(length(penaltyFactorsNonparallel) == nVar*(nLev - 1))
    # penaltyFactorsNonparallel <- if(nonparallelTerms) rep(penaltyFactors, nLev-1) else NULL
    penaltyFactors <- c(penaltyFactorsIntercept, penaltyFactorsParallel, penaltyFactorsNonparallel)

    # Augment positiveID to include all model coefficients
    if (is.null(positiveID)) positiveID <- rep(FALSE, nVar)
    positiveID <- c(rep(FALSE, nLev-1), rep(positiveID, parallelTerms + nonparallelTerms*(nLev-1)))

    # Initialize coefficient values to intercept-only model
    yFreq <- colSums(yMat) / wtsum
    interceptStart <- linkfun$g(yFreq[-nLev])
    interceptStart <- pmin(100, pmax(-100, interceptStart))
    noninterceptStart <- rep(0, nVar*(parallelTerms + nonparallelTerms*(nLev-1)))
    betaStart <- c(interceptStart, noninterceptStart)

    # Fit solution path
    mirlsNetFit <- mirlsNet(xList, yMat, alpha, penaltyFactors, positiveID, linkfun, betaStart,
                            lambdaVals, nLambda, lambdaMinRatio, includeLambda0, alphaMin,
                            pMin, stopThresh, threshOut, threshIn, maxiterOut, maxiterIn,
                            printIter, printBeta)
    betaHat <- mirlsNetFit$betaHat
    lambdaVals <- mirlsNetFit$lambdaVals
    loglik <- mirlsNetFit$loglik
    iterOut <- mirlsNetFit$iterOut
    iterIn <- mirlsNetFit$iterIn
    dif <- mirlsNetFit$dif

    # Change coefficient estimates back to original scale if standardize=TRUE
    intercepts0 <- betaHat[ , 1:(nLev-1), drop=FALSE]
    nonintercepts0 <- betaHat[ , -(1:(nLev-1)), drop=FALSE]
    unscaleFact <- if (standardize) xMeans / xSD else xMeans
    intAdjust <- matrix(0, nrow=nrow(betaHat), ncol=nLev-1)
    if (parallelTerms) intAdjust <- intAdjust +
        (nonintercepts0[ , 1:nVar, drop=FALSE] %*% unscaleFact)[ , rep(1, nLev-1), drop=FALSE]
    if (nonparallelTerms) intAdjust <- intAdjust + sapply(1:(nLev-1), function(i) {
        nonintercepts0[ , (nVar*(i-1+parallelTerms)+1):(nVar*(i+parallelTerms)), drop=FALSE] %*% unscaleFact
    })
    intercepts <- intercepts0 - intAdjust
    nonintercepts <- if (standardize) t(t(nonintercepts0) / xSD) else nonintercepts0
    coefs <- cbind(intercepts, nonintercepts)

    # Create coefficient column names
    catOrder <- if (reverse) nLev:2 else 1:(nLev-1)
    interceptNames <- paste0("(Intercept):", catOrder)
    xNames <- if (is.null(colnames(x))) paste0("X", 1:nVar) else colnames(x)
    parallelNames <- nonparallelNames <- NULL
    if (parallelTerms) parallelNames <- xNames
    if (nonparallelTerms) nonparallelNames <- paste0(rep(xNames, nLev-1), ":", rep(catOrder, each=nVar))
    colnames(coefs) <- c(interceptNames, parallelNames, nonparallelNames)

    # Calculate approximate AIC, BIC, and %deviance
    nNonzero <- apply(coefs, 1, function(b) sum(b!=0))
    aic <- -2 * loglik + 2 * nNonzero
    bic <- -2 * loglik + log(wtsum) * nNonzero
    loglikNull <- getLoglikNull(yMat)
    devPct <- 1 - loglik / loglikNull

    fit <- list(coefs=coefs, lambdaVals=lambdaVals, loglik=loglik,
                nNonzero=nNonzero, aic=aic, bic=bic, devPct=devPct,
                iterOut=iterOut, iterIn=iterIn, dif=dif,
                nLev=nLev, nVar=nVar, xNames=xNames, args=args)
    class(fit) <- "ordinalNet"
    fit
}
