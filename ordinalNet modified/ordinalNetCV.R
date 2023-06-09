#' Uses K-fold cross validation to obtain out-of-sample log-likelihood and
#' misclassification rates. Lambda is tuned within each cross validation fold.
#'
#' The data is divided into K folds. \code{ordinalNet} is fit \eqn{K} times, each time
#' leaving out one fold as a test set. For each of the \eqn{K} model fits, lambda
#' can be tuned by AIC or BIC, or cross validation. If cross validation is used,
#' the user can choose whether to user the best average out-of-sample log-likelihood,
#' misclassification rate, Brier score, or percentage of deviance explained.
#' The user can also choose the number of cross validation folds to use for tuning.
#' Once the model is tuned, the out of sample log-likelihood,
#' misclassification rate, Brier score, and percentage of deviance explained
#' are calculated on the held out test set.
#'
#' @param x Covariate matrix.
#' @param y Response variable. Can be a factor, ordered factor, or a matrix
#' where each row is a multinomial vector of counts. A weighted fit can be obtained
#' using the matrix option, since the row sums are essentially observation weights.
#' Non-integer matrix entries are allowed.
#' @param lambdaVals An optional user-specified lambda sequence (vector). If \code{NULL},
#' a sequence will be generated using the model fit to the full training data.
#' This default sequence is based on \code{nLambda} and \code{lambdaMinRatio},
#' which can be passed as additional arguments (otherwise \code{ordinalNet} default
#' values are used). The maximum lambda is the smallest value that sets all penalized
#' coefficients to zero, and the minimum lambda is the maximum value multiplied
#' by the factor \code{lambdaMinRatio}.
#' @param folds An optional list, where each element is a vector of row indices
#' corresponding to a different cross validation fold. Indices correspond to rows
#' of the \code{x} matrix. Each index number should be used in exactly one fold.
#' If \code{NULL}, the data will be randomly divided into equally-sized partitions.
#' It is recommended to call \code{set.seed} before calling \code{ordinalNetCV}
#' for reproducibility.
#' @param nFolds Numer of cross validation folds. Only used if \code{folds=NULL}.
#' @param nFoldsCV Number of cross validation folds used to tune lambda for each
#' training set (i.e. within each training fold). Only used of \code{tuneMethod} is
#' "cvLoglik", "cvMisclass", "cvBrier", or "cvDevPct.
#' @param tuneMethod Method used to tune lambda for each training set (ie. within
#' each training fold). The "cvLoglik", "cvMisclass", "cvBrier", and "cvDevPct"
#' methods use cross validation with \code{nFoldsCV} folds and select the
#' lambda value with the best average out-of-sample performance. The "aic" and "bic"
#' methods are less computationally intensive because they do not require the
#' model to be fit multiple times.
#' Note that for the methods that require cross validation, the fold splits are
#' determined randomly and cannot be specified by the user. The \code{set.seed()}
#' function should be called prior to \code{ordinalNetCV} for reproducibility.
#' @param printProgress Logical. If \code{TRUE} the fitting progress is printed
#' to the terminal.
#' @param warn Logical. If \code{TRUE}, the following warning message is displayed
#' when fitting a cumulative probability model with \code{nonparallelTerms=TRUE}
#' (i.e. nonparallel or semi-parallel model).
#' "Warning message: For out-of-sample data, the cumulative probability model
#' with nonparallelTerms=TRUE may predict cumulative probabilities that are not
#' monotone increasing."
#' The warning is displayed by default, but the user may wish to disable it.
#' @param ... Other arguments (besides \code{x}, \code{y}, \code{lambdaVals}, and \code{warn})
#' passed to \code{ordinalNet}.
#'
#' @details
#' \itemize{
#'   \item The fold partition splits can be passed by the user via the \code{folds}
#'   argument. By default, the data are randomly divided into equally-sized partitions.
#'   Note that if lambda is tuned by cross validation, the fold splits are
#'   determined randomly and cannot be specified by the user. The \code{set.seed}
#'   function should be called prior to \code{ordinalNetCV} for reproducibility.
#'   \item A sequence of lambda values can be passed by the user via the
#'   \code{lambdaVals} argument. By default, the sequence is generated by first
#'   fitting the model to the full data set (this sequence is determined by the
#'   \code{nLambda} and \code{lambdaMinRatio} arguments of \code{ordinalNet}).
#'   \item The \code{standardize} argument of \code{ordinalNet} can be modified through
#'   the additional arguments (...). If \code{standardize=TRUE}, then the data are scaled
#'   within each cross validation fold. If \code{standardize=TRUE} and lambda is tuned by
#'   cross validation, then the data are also scaled within each tuning sub-fold.
#'   This is done because scaling is part of the statistical procedure and should
#'   be repeated each time the procedure is applied.
#' }
#'
#' @return
#' An S3 object of class "ordinalNetCV", which contains the following:
#' \describe{
#'   \item{loglik}{Vector of out-of-sample log-likelihood values.
#'   Each value corresponds to a different fold.}
#'   \item{misclass}{Vector of out-of-sample misclassificaton rates.
#'   Each value corresponds to a different fold.}
#'   \item{brier}{Vector of out-of-sample Brier scores.
#'   Each value corresponds to a different fold.}
#'   \item{devPct}{Vector of out-of-sample percentages of deviance explained.
#'   Each value corresponds to a different fold.}
#'   \item{bestLambdaIndex}{The index of the value within the lambda sequence
#'   selected for each fold by the tuning method.}
#'   \item{lambdaVals}{The sequence of lambda values used for all cross validation folds.}
#'   \item{folds}{A list containing the index numbers of each fold.}
#'   \item{fit}{An object of class "ordinalNet", resulting from fitting
#'   \code{ordinalNet} to the entire dataset.}
#' }
#'
#' @examples
#' \dontrun{
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
#' # Evaluate out-of-sample performance of the  cumulative logit model
#' # when lambda is tuned by cross validation (best average out-of-sample log-likelihood)
#' cv <- ordinalNetCV(x, y, tuneMethod="cvLoglik")
#' summary(cv)
#' }
#'
#' @export
ordinalNetCVMod <- function(x, y, lambdaVals=NULL, folds=NULL, nFolds=5, nFoldsCV=5,
                         tuneMethod=c("cvLoglik", "cvMisclass", "cvBrier", "cvDevPct",
                                      "aic", "bic"),
                         printProgress=TRUE, warn=TRUE, ...)
{
    tuneMethod <- match.arg(tuneMethod)
    # cvID := indicator to use cross validation within folds
    cvID <- tuneMethod %in% c("cvLoglik", "cvMisclass", "cvBrier", "cvDevPct")
    # tuneMethod := argument passed to ordinalNetTune
    if (tuneMethod == "cvLoglik") cvCriterion <- "loglik"
    if (tuneMethod == "cvMisclass") cvCriterion <- "misclass"
    if (tuneMethod == "cvBrier") cvCriterion <- "brier"
    if (tuneMethod == "cvDevPct") cvCriterion <- "devPct"

    # Argument checks
    if (is.matrix(y) && any(rowSums(y)!=1))
        warning(paste0("Data is split by row for cross validation, but note that ",
                       "y matrix rows have different weights. Be sure this is what you want."))
    if (!is.null(folds) && length(folds)<2)
        stop(paste0("\'folds\' should be a list of at least two vectors. ",
                    "Each vector should contain indices of a cross validation fold. ",
                    "Each index from 1:nrow(x) should be used exactly once."))
    if (!is.null(folds) && !setequal(unlist(folds), 1:nrow(x)))
        stop("\'folds\' should include each index from 1:nrow(x) exactly once.")

    yMat <- if (is.matrix(y)) y else yFactorToMatrix(y)  # for computing log-likelihood
    if (printProgress) cat("Fitting ordinalNet on full training data\n")
    fit <- ordinalNetMod(x, y, lambdaVals=lambdaVals, warn=warn, ...)
    if (is.null(lambdaVals)) lambdaVals <- fit$lambdaVals

    if (is.null(folds))
    {
        n <- nrow(x)
        randIndex <- sample(n)
        folds <- split(randIndex, rep(1:nFolds, length.out=n))
    } else
    {
        nFolds <- length(folds)
    }

    nLambda <- length(lambdaVals)
    loglik <- misclass <- brier <- devPct <- bestLambdaIndex <- rep(NA, nFolds)
    names(loglik) <- names(misclass) <- names(brier) <- names(devPct) <-
        names(bestLambdaIndex) <- paste0("fold", 1:nFolds)
    for (i in 1:nFolds)
    {
        testFold <- folds[[i]]
        xTrain <- x[-testFold, , drop=FALSE]
        xTest <- x[testFold, , drop=FALSE]
        yTrain <- if (is.matrix(y)) y[-testFold, , drop=FALSE] else y[-testFold]
        yMatTest <- yMat[testFold, , drop=FALSE]
        if (printProgress) cat("Fitting ordinalNet on fold", i, "of", nFolds, '\n')

        if (cvID)
        {
            fitTrainCV <- ordinalNetTuneMod(xTrain, yTrain, lambdaVals=lambdaVals, folds=NULL,
                               nFolds=5, printProgress=FALSE, warn=FALSE, ...)
            fitTrain <- fitTrainCV$fit
            if (cvCriterion %in% c("loglik", "devPct"))
                wm <- which.max
            if (cvCriterion %in% c("misclass", "brier"))
                wm <- which.min
            bestLambdaIndex[[i]] <- wm(rowMeans(fitTrainCV[[cvCriterion]]))
        } else  # tuneMethod is either "aic" or "bic"
        {
            fitTrain <- ordinalNetMod(xTrain, yTrain, lambdaVals=lambdaVals, warn=FALSE, ...)
            bestLambdaIndex[[i]] <- which.min(fitTrain[[tuneMethod]])
        }

        pHatFull <- predict.ordinalNetMod(fitTrain, newx=xTest, type="response", whichLambda=bestLambdaIndex[[i]])
        pHat <- pHatFull[, -ncol(pHatFull), drop=FALSE]
        loglik[i] <- getLoglik(pHat, yMatTest)
        misclass[i] <- getMisclass(pHat, yMatTest)
        brier[i] <- getBrier(pHat, yMatTest)
        loglikNull <- getLoglikNull(yMatTest)
        devPct[i] <- 1 - loglik[i] / loglikNull
    }

    if (printProgress) cat("Done\n")

    out <- list(loglik=loglik, misclass=misclass, brier=brier, devPct=devPct,
                bestLambdaIndex=bestLambdaIndex, lambdaVals=lambdaVals, folds=folds, fit=fit)
    class(out) <- "ordinalNetCV"
    out
}

#' Summary method for an "ordinalNetCV" object.
#'
#' Provides a data frame which summarizes the cross validation results, which
#' can be used as an estimate of the out-of-sample performance of a model tuned
#' by a particular method.
#'
#' @param object An "ordinalNetCV" S3 object
#' @param ... Not used. Additional summary arguments.
#'
#' @return A data frame containing a record for each cross validation fold.
#' Each record contains the following: lambda value, log-likelihood,
#' misclassification rate, Brier score, and percentage of deviance explained.
#'
#' @seealso
#' \code{\link{ordinalNetCV}}
#'
#' @examples
#' # See ordinalNetCV() documentation for examples.
#'
#' @export
summary.ordinalNetCVMod <- function(object, ...)
{
    lambda <- object$lambdaVals[object$bestLambdaIndex]
    loglik <- object$loglik
    misclass <- object$misclass
    brier <- object$brier
    devPct <- object$devPct
    data.frame(lambda=lambda, loglik=loglik, misclass=misclass, brier=brier, devPct=devPct)
}

#' Print method for an "ordinalNetCV" object.
#'
#' Prints the data frame returned by the \code{summary.ordinalNetCVMod()} method.
#'
#' @param x An "ordinalNetCV" S3 object
#' @param ... Not used. Additional print arguments.
#'
#' @seealso
#' \code{\link{ordinalNetCV}}
#'
#' @examples
#' # See ordinalNetCV() documentation for examples.
#'
#' @export
print.ordinalNetCVMod <- function(x, ...)
{
    cat("\nCross validation summary:\n\n")
    print(summary.ordinalNetCVMod(x))
    cat("\n")
    invisible(x)
}
