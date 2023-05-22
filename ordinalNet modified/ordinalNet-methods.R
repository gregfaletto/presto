# Function to get column names for coefficient matrices and predictions.
getDeltaNames <- function(family, reverse, nLev)
{
    index <- if (reverse) nLev:2 else 1:(nLev-1)
    deltaNames <- sapply(index, function(i)
    {
        if (family=="cumulative") {
            if (reverse) {
                paste0("P[Y>=", i, "]")  # P(Y>=i)
            } else {
                paste0("P[Y<=", i, "]")  # P(Y<=i)
            }
        } else if (family=="sratio") {
            if (reverse) {
                paste0("P[Y=", i, "|Y<=", i, "]")  # P(Y=i|Y<=i)
            } else {
                paste0("P[Y=", i, "|Y>=", i, "]")  # P(Y=i|Y>=i)
            }
        } else if (family=="cratio") {
            if (reverse) {
                paste0("P[Y<", i, "|Y<=", i, "]")  # P(Y<i|Y<=i)
            } else {
                paste0("P[Y>", i, "|Y>=", i, "]")  # P(Y>i|Y>=i)
            }
        } else if (family=="acat") {
            if (reverse) {
                paste0("P[Y=", i, "|", i, "<=Y<=", i+1, "]")  # P(Y=i|i<=Y<=i+1)
            } else {
                paste0("P[Y=", i+1, "|", i, "<=Y<=", i+1, "]")  # P(Y=i+1|i<=Y<=i+1)
            }
        }
    })

    deltaNames
}

#' Summary method for an "ordinalNet" object.
#'
#' Provides a data frame which summarizes the model fit at each lambda value in
#' the solution path.model fit summary as a data frame.
#'
#' @param object An "ordinalNet" S3 object
#' @param ... Not used. Additional summary arguments.
#'
#' @return A data frame containing a record for each lambda value in the solution
#' path. Each record contains the following fields: lambda value, degrees of freedom
#' (number of nonzero parameters), log-likelihood, AIC, BIC, and percent deviance explained.
#'
#' @seealso
#' \code{\link{ordinalNet}}
#'
#' @examples
#' # See ordinalNet() documentation for examples.
#'
#' @export
summary.ordinalNetMod <- function(object, ...)
{
    with(object, data.frame(lambdaVals, nNonzero, loglik, devPct, aic, bic))
}

#' Print method for an "ordinalNet" object.
#'
#' Prints the data frame returned by the \code{summary.ordinalNetMod()} method.
#'
#' @param x An "ordinalNet" S3 object
#' @param ... Not used. Additional plot arguments.
#'
#' @seealso
#' \code{\link{ordinalNet}}
#'
#' @examples
#' # See ordinalNet() documentation for examples.
#'
#' @export
print.ordinalNetMod <- function(x, ...)
{
    cat("\nSummary of fit:\n\n")
    print(summary.ordinalNetMod(x))
    cat("\n")
    invisible(x)
}

#' Method to extract fitted coefficients from an "ordinalNet" object.
#'
#' @param object An "ordinalNet" S3 object.
#' @param matrix Logical. If \code{TRUE}, coefficient estimates are returned in
#' matrix form. Otherwise a vector is returned.
#' @param whichLambda Optional index number of the desired \code{lambda} within
#' the sequence of \code{lambdaVals}. By default, the solution with the best AIC
#' is returned.
#' @param criteria Selects the best \code{lambda} value by AIC or BIC. Only used
#' if \code{whichLambda=NULL}.
#' @param ... Not used. Additional coef arguments.
#'
#' @return The object returned depends on \code{matrix}.
#'
#' @seealso
#' \code{\link{ordinalNet}}
#'
#' @examples
#' # See ordinalNet() documentation for examples.
#'
#' @export
coef.ordinalNetMod <- function(object, matrix=FALSE, whichLambda=NULL, criteria=c("aic", "bic"), ...)
{
    if (!is.null(whichLambda) && length(whichLambda)!=1)
        stop("whichLambda should be a single value, or NULL.")
    criteria <- match.arg(criteria)
    if (is.null(whichLambda)) whichLambda <- which.min(object[[criteria]])
    betaHat <- object$coefs[whichLambda, ]
    if (!matrix) return(betaHat)
    # Extract variables from ordinalNet object
    nLev <- object$nLev
    nVar <- object$nVar
    parallelTerms <- object$args$parallelTerms
    nonparallelTerms <- object$args$nonparallelTerms
    family <- object$args$family
    reverse <- object$args$reverse
    xNames <- object$xNames
    link <- if (!is.null(object$args$customLink)) "g" else object$args$link
    # Create coefficient matrix
    intercepts <- betaHat[1:(nLev-1)]
    nonintercepts <- matrix(0, nrow=nVar, ncol=nLev-1)
    if (parallelTerms) nonintercepts <- nonintercepts + betaHat[nLev:(nLev-1+nVar)]
    if (nonparallelTerms) nonintercepts <- nonintercepts + betaHat[-(1:(nLev-1+nVar*parallelTerms))]
    betaMat <- rbind(intercepts, nonintercepts)
    rownames(betaMat) <- c("(Intercept)", xNames)
    deltaNames <- getDeltaNames(family, reverse, nLev)
    colnames(betaMat) <- paste0(link, "(", deltaNames, ")")
    betaMat
}


#' Predict method for an "ordinalNet" object
#'
#' Obtains predicted probabilities, predicted class, or linear predictors.
#'
#' @param object An "ordinalNet" S3 object.
#' @param newx Optional covariate matrix. If NULL, fitted values will be obtained
#' for the training data, as long as the model was fit with the argument
#' \code{keepTrainingData=TRUE}.
#' @param whichLambda Optional index number of the desired lambda value within
#' the solution path sequence.
#' @param criteria Selects the best lambda value by AIC or BIC. Only used
#' if \code{whichLambda=NULL}.
#' @param type The type of prediction required.  Type "response" returns a
#' matrix of fitted probabilities. Type "class" returns a vector containing the
#' class number with the highest fitted probability. Type "link" returns a
#' matrix of linear predictors.
#' @param ... Not used. Additional predict arguments.
#'
#' @return The object returned depends on \code{type}.
#'
#' @seealso
#' \code{\link{ordinalNet}}
#'
#' @examples
#' # See ordinalNet() documentation for examples.
#'
#' @export
predict.ordinalNetMod <- function(object, newx=NULL, whichLambda=NULL, criteria=c("aic", "bic"),
                                  type=c("response", "class", "link"), ...)
{
    criteria <- match.arg(criteria)
    type <- match.arg(type)
    if (is.null(newx) && !object$args$keepTrainingData)
        stop(paste0("Model was fit with keepTrainingData=FALSE, so training data was not saved. ",
             "A newx argument is required."))
    # Extract variables from ordinalNet object
    nLev <- object$nLev
    xNames <-  object$xNames
    parallelTerms <- object$args$parallelTerms
    nonparallelTerms <- object$args$nonparallelTerms
    family <- object$args$family
    link <- object$args$link
    reverse <- object$args$reverse
    linkfun <- if (is.null(object$args$customLink)) makeLinkfun(family, link) else object$args$customLink
    x <- if (is.null(newx)) object$args$x else newx
    betaMat <- coef.ordinalNetMod(object, matrix=TRUE, whichLambda=whichLambda, criteria=criteria)
    # Compute prediction values
    etaMat <- cbind(1, x) %*% betaMat
    deltaNames <- getDeltaNames(family, reverse, nLev)
    colnames(etaMat) <- colnames(betaMat)
    if (type == "link") return(etaMat)
    probMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
    probMat <- cbind(probMat, 1-rowSums(probMat))
    if (reverse) probMat <- probMat[, nLev:1]
    colnames(probMat) <- paste0("P[Y=", 1:nLev, "]")
    if (type == "response") return(probMat)
    class <- c(apply(probMat, 1, which.max))
    class
}
