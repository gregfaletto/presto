## General optimization algorithm for multinomial regression models via coordinate descent
mirlsNet <- function(xList, yMat, alpha, penaltyFactors, positiveID, linkfun, betaStart,
                     lambdaVals, nLambda, lambdaMinRatio, includeLambda0, alphaMin,
                     pMin, stopThresh, threshOut, threshIn, maxiterOut, maxiterIn,
                     printIter, printBeta)
{
    wtsum <- if (is.null(attr(yMat, "wtsum"))) sum(yMat) else attr(yMat, "wtsum")
    nObs <- nrow(yMat)
    xMat <- do.call(rbind, xList)
    if (!is.null(lambdaVals)) lambdaVals <- sort(lambdaVals, decreasing=TRUE)
    lambdaNum <- if (is.null(lambdaVals)) nLambda + includeLambda0 else length(lambdaVals)
    fits <- vector("list", length=lambdaNum)

    # If lambdaVals=NULL, need to determine the minimum lambda value that sets all penalized coefficients to zero
    if (is.null(lambdaVals))
    {
        lambdaMod <- ifelse(penaltyFactors==0, 0, Inf)  # to find solution with only unpenalized terms
        fits[[1]] <- cdOut(betaHat=betaStart, lambdaIndex=1, lambdaNum, lambdaMod,
                           xList, xMat, yMat, max(alpha, alphaMin), positiveID, linkfun,
                           pMin, threshOut, threshIn, maxiterOut, maxiterIn,
                           printIter, printBeta)
        betaStart <- fits[[1]]$betaHat

        # Calculate starting lambda value
        etaMat <- matrix(xMat %*% betaStart, nrow=nObs, byrow=TRUE)
        pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
        si <- getScoreInfo(xList, yMat, pMat, pMin, linkfun)
            # betaHat is zero for all penalized terms, so the soft threshold argument is just the score function
        penID <- penaltyFactors != 0
        lambdaMaxVals <- si$score[penID] / (wtsum * max(alpha, alphaMin) * penaltyFactors[penID])
        lambdaMaxVals[positiveID[penID]] <- pmax(0, lambdaMaxVals[penID & positiveID])
        lambdaMaxVals <- abs(lambdaMaxVals)
        lambdaMax <- max(lambdaMaxVals)
        lambdaMin <- lambdaMax * lambdaMinRatio
        lambdaVals <- exp(seq(log(lambdaMax), log(lambdaMin), length.out=nLambda))
        if (includeLambda0) lambdaVals <- c(lambdaVals, 0)
    }

    # If alpha < alphaMin, the model needs to be re-fit the first lambda value
    # using alpha instead of alphaMin
    if (alpha < alphaMin)
        fits[1] <- list(NULL)

    # fits[[1]] is NULL if alpha < alphaMin or if lambdaVals is specified by user.
    llik <- if (is.null(fits[[1]])) -Inf else fits[[1]]$loglik
    for (i in (1+!is.null(fits[[1]])):lambdaNum)
    {
        # If relative change in loglik is < stopThresh, then use the current fit
        # for all remaining lambda. Do not stop if loglik stays the same, because
        # this can happen if the first several lambda values produce null models,
        # e.g. in cross validation.
        if ((i > 2) && llikOld != llik && (abs((llikOld - llik) / llikOld) < stopThresh))
        {
            fits[[i]] <- fits[[i-1]]
        } else
        {
            lambdaMod <- lambdaVals[i] * penaltyFactors
            lambdaMod <- ifelse(penaltyFactors==0, 0, lambdaVals[i] * penaltyFactors)
            fits[[i]] <- cdOut(betaHat=betaStart, lambdaIndex=i, lambdaNum, lambdaMod,
                               xList, xMat, yMat, alpha, positiveID, linkfun,
                               pMin, threshOut, threshIn, maxiterOut, maxiterIn,
                               printIter, printBeta)
            betaStart <- fits[[i]]$betaHat
            llikOld <- llik
            llik <- fits[[i]]$loglik
        }
    }

    iterOut <- sapply(fits, function(x) x$iterOut)
    iterIn <- sapply(fits, function(x) x$iterIn)
    dif <- sapply(fits, function(x) x$dif)

    if (any(iterOut==maxiterOut))
        warning(paste0("Reached outer loop iteration limit before convergence ",
                       "for at least one lambda value. Consider increasing maxiterOut."))
    # if (any(iterIn==maxiterIn & dif>0))
    #     warning(paste0("Outer loop converged upon inner loop reaching iteration ",
    #                    "limit for at least one lambda value. Consider increasing maxiterIn."))

    betaHat <- t(sapply(fits, function(f) f$betaHat))
    loglik <- sapply(fits, function(f) f$loglik)
    list(lambdaVals=lambdaVals, betaHat=betaHat, loglik=loglik, iterOut=iterOut, iterIn=iterIn, dif=dif)
}
