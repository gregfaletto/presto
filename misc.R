# returns wtsum * lambda * (alpha*|betaHat|_1 + (1-alpha)/2*|betaHat|_2^2)
getPenalty <- function(betaHat, lambdaMod, alpha)
{
    lambdaMod[betaHat == 0] <- 0  # 0*Inf=0
    l1 <- sum(abs(betaHat) * lambdaMod)
    l22 <- sum(betaHat^2 * lambdaMod)
    pen1 <- alpha * l1
    pen2 <- .5 * (1-alpha) * l22
    pen <- pen1 + pen2
    pen
}

getLoglik <- function(pMat, yMat)
{
    pkplusone <- 1 - rowSums(pMat)
    pMatFull <- cbind(pMat, pkplusone)
    if (any(pMatFull < 0)) return(-Inf)
    llMat <- yMat * log(pMatFull)
    llMat[yMat==0] <- 0  # -Inf*0 = 0
    llik <- sum(llMat)
    llik
}

getMisclass <- function(pMat, yMat)
{
    pkplusone <- 1 - rowSums(pMat)
    pMatFull <- cbind(pMat, pkplusone)
    predClass <- apply(pMatFull, 1, which.max)
    nMisclass <- sapply(1:nrow(yMat), function(i) sum(yMat[i, -predClass[i]]))
    misclass <- sum(nMisclass) / sum(yMat)
    misclass
}

getBrier <- function(pMat, yMat)
{
    pkplusone <- 1 - rowSums(pMat)
    pMatFull <- cbind(pMat, pkplusone)
    n <- rowSums(yMat)
    brier <- sum(yMat * (1 - pMatFull)^2 + (n - yMat) * pMatFull^2) / sum(n)
    brier
}

# Returns approximate log-likelihood (as a function of beta, up to a constant)
getLoglikApprox <- function(betaHatActive, scoreActive, infoActive)
{
    -sum(betaHatActive * (infoActive %*% betaHatActive)) + sum(scoreActive * betaHatActive)
}

getLoglikNull <- function(yMat)
{
    pHatNull <- colSums(yMat) / sum(yMat)
    llmat0 <- yMat * rep(log(pHatNull), each=nrow(yMat))
    llmat0[yMat == 0] <- 0  # 0*-Inf=0
    loglik0 <- sum(llmat0)
    loglik0
}

invLogit <- function(x) 1/(1+exp(-x))

softThresh <- function(z, g) sign(z)*max(0, abs(z)-g)

yFactorToMatrix <- function(y)
{
    nObs <- length(y)
    nLev <- length(levels(y))
    yMat <- matrix(0, nrow=nObs, ncol=nLev, dimnames=list(NULL, levels(y)))
    yInt <- as.integer(y)
    yMat[cbind(1:nObs, yInt)] <- 1
    yMat
}

# Subroutine for makexMat function
makeNonparallelBlock <- function(x, nLev)
{
    nVar <- length(x)
    zeros <- rep(0, nVar)
    block <- do.call(rbind, lapply(1:(nLev-1), function(i)
    {
        # MODIFICATION HERE FOR OUR MODEL
        # c(rep(zeros, i-1), x, rep(zeros, nLev-1-i))
        c(rep(x, i-1), x, rep(zeros, nLev-1-i))
    }))
    block
}

# at least one of parallelTerms and nonparellelTerms should be TRUE
makexList <- function(x, nLev, parallelTerms, nonparallelTerms)
{
    if(!nonparallelTerms){
        stop("Only non-parallel model currently supported (nonparallelTerms must be TRUE)")
    }
    if(parallelTerms){
        stop("Only non-parallel model currently supported (parallelTerms must be FALSE)")
    }
    x1 <- diag(nLev-1)  # intercept columns
    xList <- lapply(1:nrow(x), function(i)
    {
        xi <- x[i, ]
        x2 <- if (!parallelTerms) NULL else rbind(xi)[rep(1, nLev-1), , drop=FALSE]
        x3 <- if (!nonparallelTerms) NULL else makeNonparallelBlock(xi, nLev)
        xListi <- cbind(x1, x2, x3)
        rownames(xListi) <- NULL
        xListi
    })
    xList
}

getScoreInfo <- function(xList, yMat, pMat, pMin, linkfun)
{
    nObs <- length(xList)
    nCoef <- ncol(xList[[1]])
    nLev <- ncol(yMat)

    # Fitted probabilities less than pMin are set to pMin; everything is rescaled to sum to 1
    pMatFull <- cbind(pMat, 1-rowSums(pMat))
    pMatFull[pMatFull < pMin] <- pMin
    pMatFull <- pMatFull / rowSums(pMatFull)
    pkplusone <- pMatFull[, nLev]  # vector of fitted probabilities for class K+1
    pMat <- pMatFull[, -nLev, drop=FALSE]

    d <- yMat[, nLev] / pkplusone
    d[yMat[, nLev] == 0] <- 0  # 0/0 = 0
    uMat <- yMat[, -nLev, drop=FALSE] / pMat
    uMat[yMat[, -nLev] == 0] <- 0  # 0/0 = 0
    uminusdMat <- uMat - d

    wts <- if (is.null(attr(yMat, "wts"))) rowSums(yMat) else attr(yMat, "wts")
    wpMat <- wts / pMat
    wpkplusone <- wts / pkplusone

    score <- rep(0, nCoef)
    info <- matrix(0, nrow=nCoef, ncol=nCoef)
    for (i in 1:nObs)
    {
        # Compute score term
        x <- xList[[i]]
        uminusd <- uminusdMat[i, ]
        eta <- linkfun$g(pMat[i, ])  # compute eta from p, capped at pMin
        q <- linkfun$getQ(eta)
        score <- score + crossprod(x, crossprod(q, uminusd))

        # Compute info term
        sigInv <- diag(wpMat[i, ], nrow=nLev-1) + wpkplusone[i]
        w <- crossprod(q, sigInv) %*% q
        info <- info + crossprod(x, w %*% x)
    }

    score <- c(score)
    list(score=score, info=info)
}
