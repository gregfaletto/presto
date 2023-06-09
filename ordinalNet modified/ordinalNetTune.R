#' Uses K-fold cross validation to obtain out-of-sample log-likelihood and
#' misclassification rates for a sequence of lambda values.
#'
#' The data is divided into K folds. \code{ordinalNet} is fit \eqn{K} times (\code{K=nFolds}),
#' each time leaving out one fold as a test set. The same sequence of lambda values is used
#' each time. The out-of-sample log-likelihood, misclassification rate, Brier score,
#' and percentage of deviance explained are obtained for each lambda value from
#' the held out test set. It is up to the user to determine how to tune the model
#' using this information.
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
#' If \code{NULL}, the data will be randomly divided into equal-sized partitions.
#' It is recommended to use \code{set.seed} before calling this function to make
#' results reproducible.
#' @param nFolds Numer of cross validation folds. Only used if \code{folds=NULL}.
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
#'   The \code{set.seed} function should be called prior to \code{ordinalNetCV} for reproducibility.
#'   \item A sequence of lambda values can be passed by the user via the
#'   \code{lambdaVals} argument. By default, the sequence is generated by first
#'   fitting the model to the full data set (this sequence is determined by the
#'   \code{nLambda} and \code{lambdaMinRatio} arguments of \code{ordinalNet}).
#'   \item The \code{standardize} argument of \code{ordinalNet} can be modified through
#'   the additional arguments (...). If \code{standardize=TRUE}, then the data are scaled
#'   within each cross validation fold. This is done because scaling is part of
#'   the statistical procedure and should be repeated each time the procedure is applied.
#' }
#'
#' @return
#' An S3 object of class "ordinalNetTune", which contains the following:
#' \describe{
#'   \item{loglik}{Matrix of out-of-sample log-likelihood values.
#'   Each row corresponds to a lambda value, and each column corresponds to a fold.}
#'   \item{misclass}{Matrix of out-of-sample misclassificaton rates.
#'   Each row corresponds to a lambda value, and each column corresponds to a fold.}
#'   \item{brier}{Matrix of out-of-sample Brier scores. Each row corresponds
#'   to a lambda value, and each column corresponds to a fold.}
#'   \item{devPct}{Matrix of out-of-sample percentages of deviance explained.
#'   Each row corresponds to a lambda value, and each column corresponds to a fold.}
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
#' # Fit parallel cumulative logit model; select lambda by cross validation
#' tunefit <- ordinalNetTune(x, y)
#' summary(tunefit)
#' plot(tunefit)
#' bestLambdaIndex <- which.max(rowMeans(tunefit$loglik))
#' coef(tunefit$fit, whichLambda=bestLambdaIndex, matrix=TRUE)
#' predict(tunefit$fit, whichLambda=bestLambdaIndex)
#' }
#'
#' @export
ordinalNetTuneMod <- function(x, y, lambdaVals=NULL, folds=NULL, nFolds=5, printProgress=TRUE, warn=TRUE, ...)
{
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
    loglik <- misclass <- brier <- devPct <- matrix(nrow=nLambda, ncol=nFolds)
    colnames(loglik) <- colnames(misclass) <- colnames(brier) <- colnames(devPct) <-
        paste0("fold", 1:nFolds)
    rownames(loglik) <- rownames(misclass) <- rownames(brier) <- rownames(devPct) <-
        paste0("lambda", 1:nLambda)
    for (i in 1:nFolds)
    {
        testFold <- folds[[i]]
        xTrain <- x[-testFold, , drop=FALSE]
        xTest <- x[testFold, , drop=FALSE]
        yTrain <- if (is.matrix(y)) y[-testFold, , drop=FALSE] else y[-testFold]
        yMatTest <- yMat[testFold, , drop=FALSE]
        if (printProgress) cat("Fitting ordinalNet on fold", i, "of", nFolds, '\n')
        fitTrain <- ordinalNetMod(xTrain, yTrain, lambdaVals=lambdaVals, warn=FALSE, ...)
        nLambdaTrain <- length(fitTrain$lambdaVals)

        for (j in 1:nLambda)
        {
            if (j > nLambdaTrain) {
                loglik[j, i] <- NA
                misclass[j, i] <- NA
            } else {
                pHatFull <- predict.ordinalNetMod(fitTrain, newx=xTest, type="response", whichLambda=j)
                pHat <- pHatFull[, -ncol(pHatFull), drop=FALSE]
                loglik[j, i] <- getLoglik(pHat, yMatTest)
                misclass[j, i] <- getMisclass(pHat, yMatTest)
                brier[j, i] <- getBrier(pHat, yMatTest)
                loglikNull <- getLoglikNull(yMatTest)
                devPct[j, i] <- 1 - loglik[j, i] / loglikNull
            }
        }
    }

    if (printProgress) cat("Done\n")

    out <- list(loglik=loglik, misclass=misclass, brier=brier, devPct=devPct,
                lambdaVals=lambdaVals, folds=folds, fit=fit)
    class(out) <- "ordinalNetTune"
    out
}






#' Summary method for an "ordinalNetTune" object.
#'
#' Provides a data frame which summarizes the cross validation results and may
#' be useful for selecting an appropriate value for the tuning parameter lambda.
#'
#' @param object An "ordinalNetTune" S3 object.
#' @param ... Not used. Additional summary arguments.
#'
#' @return A data frame containing a record for each lambda value in the solution
#' path. Each record contains the following: lambda value, average log-likelihood,
#' average misclassification rate, average Brier score, and average percentage
#' of deviance explained. Averages are taken across all cross validation folds.
#'
#' @seealso
#' \code{\link{ordinalNetTune}}
#'
#' @examples
#' # See ordinalNetTune() documentation for examples.
#'
#' @export
summary.ordinalNetTuneMod <- function(object, ...)
{
    lambda <- object$lambdaVals
    loglik_avg <- unname(rowMeans(object$loglik))
    misclass_avg <- unname(rowMeans(object$misclass))
    brier_avg <- unname(rowMeans(object$brier))
    devPct_avg <- unname(rowMeans(object$devPct))
    data.frame(lambda=lambda, loglik_avg=loglik_avg, misclass_avg=misclass_avg,
               brier_avg=brier_avg, devPct_avg=devPct_avg)
}

#' Print method for an "ordinalNetTune" object.
#'
#' Prints the data frame returned by the \code{summary.ordinalNetTuneMod()} method.
#'
#' @param x An "ordinalNetTune" S3 object.
#' @param ... Not used. Additional print arguments.
#'
#' @seealso
#' \code{\link{ordinalNetTune}}
#'
#' @examples
#' # See ordinalNetTune() documentation for examples.
#'
#' @export
print.ordinalNetTuneMod <- function(x, ...)
{
    cat("\nCross validation summary:\n\n")
    print(summary.ordinalNetTuneMod(x))
    cat("\n")
    invisible(x)
}

#' Plot method for "ordinalNetTune" object.
#'
#' Plots the average out-of-sample log-likelihood, misclassification rate,
#' Brier score, or percentage of deviance explained for each lambda value in the
#' solution path. The averae is taken over all cross validation folds.
#'
#' @param x An "ordinalNetTune" S3 object.
#' @param type Which performance measure to plot. Either "loglik", "misclass",
#' "brier", or "devPct".
#' @param ... Additional plot arguments.
#'
#' @seealso
#' \code{\link{ordinalNetTune}}
#'
#' @examples
#' # See ordinalNetTune() documentation for examples.
#'
#'@export
plot.ordinalNetTuneMod <- function(x, type=c("loglik", "misclass", "brier", "devPct"), ...)
{
    type <- match.arg(type)
    y <- rowMeans(x[[type]])
    loglambda <- log(x$lambdaVals)
    if (type == "misclass")
        ylab <- "avg misclass rate"
    if (type == "loglik")
        ylab <- "avg loglik"
    if (type == "brier")
        ylab <- "avg brier score"
    if (type == "devPct")
        ylab <- "avg pct deviance explained"
    graphics::plot(y ~ loglambda, ylab=ylab, xlab="log(lambda)", ...)
}
