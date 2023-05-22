# Each of these functions creates a list of functions consisting of
# g (link), h (inverse link), and getQ (Jacobian of inverse link, dh/ddeta^T).
# delta is a transformation of p to which the element-wise link function is applied (e.g. logit)
# tt(p) = delta, ttinv(delta) = p, and ttinvprime(delta) = dttinv/ddelta^T

# Wrapper for link function construction functions
makeLinkfun <- function(family, link)
{
    if (family == "cumulative") {
        linkfun <- makeLinkCumulative(link)
    } else if (family == "sratio") {
        linkfun <- makeLinkRatio(link, stopping=TRUE)
    } else if (family == "cratio") {
        linkfun <- makeLinkRatio(link, stopping=FALSE)
    } else if (family == "acat") {
        linkfun <- makeLinkACAT(link)
    }

    linkfun
}

# Cumulative probability family
makeLinkCumulative <- function(link)
{
    lf <- stats::make.link(link)

    tt <- function(p) cumsum(p)

    ttinv <- function(delta) delta - c(0, delta[-length(delta)])

    ttinvprime <- function(delta)
    {
        k <- length(delta)
        ttipDiag <- diag(1, nrow=k)
        ttipOffDiag <- cbind(-ttipDiag[, -1], 0)
        ttip <- ttipDiag + ttipOffDiag
        ttip
    }

    g <- function(p) lf$linkfun(tt(p))
    h <- function(eta) ttinv(lf$linkinv(eta))
    getQ <- function(eta) rep(lf$mu.eta(eta), each=length(eta)) * ttinvprime(lf$linkinv(eta))
    list(g=g, h=h, getQ=getQ)
}

# (Forward) stopping and continuation ratio families
# stopping is logical argument, if FALSE, then make continuation ratio link
makeLinkRatio <- function(link, stopping)
{
    lf <- stats::make.link(link)

    tt <- function(p)
    {
        k <- length(p)
        cp <- c(1, 1-cumsum(p[-k]))  # 1-cumulative probabilites
        delta <- p / cp
        delta
    }

    ttinv <- function(delta)
    {
        k <- length(delta)
        p <- rep(NA, k)
        p[1] <- delta[1]
        cp <- 1 - p[1]  # 1-cumulative probability dummy variable
        if (k >= 2) {
            for (i in 2:k) {
                p[i] <- delta[i] * cp
                cp <- cp - p[i]
            }
        }
        p
    }

    ttinvprime <- function(delta)
    {
        k <- length(delta)
        p <- ttinv(delta)
        cp <- c(1, 1-cumsum(p))  # 1-cumulative probabilities
        ttip <- matrix(nrow=k, ncol=k)
        ttip[1, ] <- c(1, rep(0, k-1))
        cs <- ttip[1, ]  # cumulative row sum dummy variable
        if (k >= 2) {
            for (i in 2:k) {
                ttip[i, ] <- -delta[i] * cs
                ttip[i, i] <- cp[i]
                cs <- cs + ttip[i, ]
            }
        }
        ttip
    }

    if (stopping)
    {
        g <- function(p) lf$linkfun(tt(p))
        h <- function(eta) ttinv(lf$linkinv(eta))
        getQ <- function(eta) rep(lf$mu.eta(eta), each=length(eta)) * ttinvprime(lf$linkinv(eta))
    } else
    {
        g <- function(p) lf$linkfun(1 - tt(p))
        h <- function(eta) ttinv(1 - lf$linkinv(eta))
        getQ <- function(eta) rep(lf$mu.eta(eta), each=length(eta)) * -ttinvprime(1-lf$linkinv(eta))
    }
    list(g=g, h=h, getQ=getQ)
}

## Adjacent category
makeLinkACAT <- function(link)
{
    lf <- stats::make.link(link)

    tt <- function(p)
    {
        pp <- c(p[-1], 1-sum(p))
        delta <- pp / (p + pp)
        delta
    }

    ttinv <- function(delta)
    {
        k <- length(delta)
        ar <- delta / (1-delta)  # adjacent ratios p2/p1, p3/p2, ... , pk+1/pk
        r1 <- cumprod(ar)  # p relative to p1:  p2/p1, p3/p1, ..., pk+1/p1
        p1 <- 1 / (1 + sum(r1))
        p <- c(p1, p1*r1[-k])
        p
    }

    ttinvprime <- function(delta)
    {
        k <- length(delta)
        p <- ttinv(delta)
        p1 <- p[1]
        cp <- cumsum(p)  # cumulative probabilities
        ar <- delta / (1-delta)  # adjacent ratios p2/p1, p3/p2, ... , pk+1/pk
        dttdar <- matrix(nrow=k, ncol=k)  # jacobian of ttinv with respect to ar
        dttdar[1, ] <- -p1 * (1-cp) / ar
        if (k >= 2) {
            for (i in 2:k) {
                dttdar[i, ] <- ar[i-1] * dttdar[i-1, ]
                dttdar[i, i-1] <- dttdar[i, i-1] + p[i-1]
            }
        }
        darddelta <- 1 / (1-delta)^2  # jacobian of ar with respect to delta (diagonal matrix)
        ttip <- rep(darddelta, each=k) * dttdar  # multiply each row of dttdar by darddelta
        ttip
    }

    g <- function(p) lf$linkfun(tt(p))
    h <- function(eta) ttinv(lf$linkinv(eta))
    getQ <- function(eta) rep(lf$mu.eta(eta), each=length(eta)) * ttinvprime(lf$linkinv(eta))
    list(g=g, h=h, getQ=getQ)
}
