# Coordinate descent outer loop function
cdOut <- function(betaHat, lambdaIndex, lambdaNum, lambdaMod,
                  xList, xMat, yMat, alpha, positiveID, linkfun,
                  pMin, threshOut, threshIn, maxiterOut, maxiterIn,
                  printIter, printBeta)
{
    nObs <- nrow(yMat)
    wts <- if (is.null(attr(yMat, "wts"))) rowSums(yMat) else attr(yMat, "wts")
    wtsum <- if (is.null(attr(yMat, "wtsum"))) sum(wts) else attr(yMat, "wtsum")
    nLev <- ncol(yMat)
    nVar <- ncol(xMat)

    # Could carry over epsMat, pMat, and loglik from previous cdOut iteration
    betaNonzeroIndex <- which(betaHat != 0)
    etaMat <- matrix(xMat[, betaNonzeroIndex, drop=FALSE] %*% betaHat[betaNonzeroIndex], nrow=nObs, byrow=TRUE)
    pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
    loglik <- getLoglik(pMat, yMat)
    penalty <- getPenalty(betaHat, lambdaMod, alpha)
    obj <- -loglik/wtsum + penalty

    conv <- FALSE
    iterOut <- 0
    while (!conv && iterOut<maxiterOut)
    {
        iterOut <- iterOut + 1

        # Update score and info
        si <- getScoreInfo(xList, yMat, pMat, pMin, linkfun)
        score <- si$score
        info <- si$info

        # Run coordinate descent inner loop
        betaHatOld <- betaHat
        cdInResult <- cdIn(wtsum, betaHat, score, info, alpha, lambdaMod, positiveID, threshIn, maxiterIn)
        betaHat <- cdInResult$betaHat
        iterIn <- cdInResult$iterIn
        betaNonzeroIndex <- which(betaHat != 0)
        etaMat <- matrix(xMat[, betaNonzeroIndex, drop=FALSE] %*% betaHat[betaNonzeroIndex], nrow=nObs, byrow=TRUE)
        pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))

        # Update log-likelihood and objective
        loglikOld <- loglik
        objOld <- obj
        loglik <- getLoglik(pMat, yMat)
        penalty <- getPenalty(betaHat, lambdaMod, alpha)
        obj <- -loglik / wtsum + penalty

        # Take half steps if obj does not improve. Loglik is set to -Inf
        # if any fitted probabilities are negative, which can happen for
        # the nonparallel or semiparallel cumulative probability model.
        nhalf <- 0
        while (obj > objOld && nhalf < 10) {
            nhalf <- nhalf + 1
            betaHat <- (betaHat + betaHatOld) / 2
            betaNonzeroIndex <- which(betaHat != 0)
            etaMat <- matrix(xMat[, betaNonzeroIndex, drop=FALSE] %*% betaHat[betaNonzeroIndex], nrow=nObs, byrow=TRUE)
            pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
            loglik <- getLoglik(pMat, yMat)
            penalty <- getPenalty(betaHat, lambdaMod, alpha)
            obj <- -loglik / wtsum + penalty
        }
        dif <- (objOld - obj) / (abs(objOld) + 1e-100)
        conv <- dif < threshOut

        # Convergence is declared if objective worsens. In this case, set betaHat
        # to previous value. (Typically means model is saturated.)
        objImproved <- obj <= objOld
        if (!objImproved)
        {
            betaHat <- betaHatOld
            loglik <- loglikOld
        }

        # Print iteration info if printIter=TRUE
        if (printIter)
        {
            if (iterOut == 1) cat("\nLambda", lambdaIndex, " of ", lambdaNum, '\n')
            cat("outer iteration ", iterOut, ":  ", iterIn,
                " inner iterations, relative change in objective: ",
                signif(dif, 2), '\n', sep='')
        }

        # Print betaHat if printBeta=TRUE
        if (printBeta)
        {
            if (!printIter)
            {
                if (iterOut==1) cat("\nLambda", lambdaIndex, " of ", lambdaNum, '\n', sep='')
                cat("outer iteration ", iterOut, " ", '\n', sep='')
            }
            cat("betaHat: ", signif(betaHat, 2), '\n\n')
        }

    }  # end while (!conv && iterOut<maxiterOut)

    # Opting not to return penalty or obj because they depend on covariate scaling.
    list(betaHat=betaHat, loglik=loglik, iterOut=iterOut, iterIn=iterIn, dif=dif)
}
