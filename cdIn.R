## Coordinate descent inner loop function
# Note: should not need to check for improvement because each coordinate step necessarily improves the approximate objective
cdIn <- function(wtsum, betaHat, score, info, alpha, lambdaMod, positiveID, threshIn, maxiterIn)
{
    # Update Active Set
    activeSet <- which((betaHat!=0 | lambdaMod==0) & diag(info)!=0)
        # check lambdaMod==0 because if data are balanced, an intercept could be initialized exactly to zero
    betaHat[-activeSet] <- 0
    betaHatActive <- betaHat[activeSet]
    scoreActive <- score[activeSet]
    infoActive <- info[activeSet, activeSet, drop=FALSE]
    infoInactive <- info[-activeSet, activeSet, drop=FALSE]
    lambdaModActive <- lambdaMod[activeSet]
    lambdaModInactive <- lambdaMod[-activeSet]
    positiveIDActive <- positiveID[activeSet]
    positiveIDInactive <- positiveID[-activeSet]

    # softThreshTerms vector does not change during the inner loop, even if active set changes
    softThreshTerms <- c(info[, activeSet, drop=FALSE] %*% betaHatActive) + score
        # softThreshTerms = I(beta^(r)) %*% beta^(r) + U(beta^(r)) = I(beta^(r)) %*% beta^(r+1)
    softThreshTermsActive <- softThreshTerms[activeSet]
    softThreshTermsInactive <- softThreshTerms[-activeSet]

    # Initialize quadratic approximation to the log-likelihood and objective
    loglikApprox <- getLoglikApprox(betaHatActive, scoreActive, infoActive)
    penalty <- getPenalty(betaHat, lambdaMod, alpha)
    objApprox <- -loglikApprox/wtsum + penalty

    iterIn <- 0
    kktAll <- FALSE
    while (!kktAll && iterIn<maxiterIn)
    {
        conv <- FALSE
        while (!conv && iterIn<maxiterIn)
        {
            iterIn <- iterIn + 1
            for (i in seq_along(activeSet))
            {
                numTerm <- softThreshTermsActive[i] - sum(infoActive[i, -i, drop=FALSE] * betaHatActive[-i])
                denTerm <- infoActive[i, i]
                penTerm <- wtsum * lambdaModActive[i]
                betaHatActive[i] <- softThresh(numTerm, penTerm*alpha) / (denTerm + penTerm*(1-alpha))
                if (positiveIDActive[i]) betaHatActive[i] <- max(0, betaHatActive[i])
            }

            betaHatOld <- betaHat
            betaHat[activeSet] <- betaHatActive
            loglikApprox <- getLoglikApprox(betaHatActive, scoreActive, infoActive)
            penalty <- getPenalty(betaHat, lambdaMod, alpha)
            objApproxOld <- objApprox
            objApprox <- -loglikApprox / wtsum + penalty
            dif <- abs((objApproxOld - objApprox) / (abs(objApproxOld) + 1e-100))
            conv <- dif < threshIn

        }  # end while (!conv && iterIn<maxiterIn)

        kkt <- rep(TRUE, length(betaHat))
        kktInactiveTerms <- softThreshTermsInactive - c(infoInactive %*% betaHatActive)
        kktInactiveTerms[!positiveIDInactive] <- abs(kktInactiveTerms[!positiveIDInactive])
        kktInactive <- kktInactiveTerms <= wtsum * lambdaModInactive * alpha
        kkt[-activeSet] <- kktInactive
        kktAll <- all(kkt)
        if (!kktAll)
        {
            iterIn <- iterIn - 1  # repeat the iteration if kkt conditions are not satisfied
            activeSet <- union(activeSet, which(!kkt))
            betaHatActive <- betaHat[activeSet]
            scoreActive <- score[activeSet]
            infoActive <- info[activeSet, activeSet, drop=FALSE]
            infoInactive <- info[-activeSet, activeSet, drop=FALSE]
            softThreshTermsActive <- softThreshTerms[activeSet]
            softThreshTermsInactive <- softThreshTerms[-activeSet]
            lambdaModActive <- lambdaMod[activeSet]
            lambdaModInactive <- lambdaMod[-activeSet]
            positiveIDActive <- positiveID[activeSet]
            positiveIDInactive <- positiveID[-activeSet]
        }

    }  # end while (!kktAll && iterIn<maxiterIn)

    list(betaHat=betaHat, iterIn=iterIn)
}  # end cdIn
