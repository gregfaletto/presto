# setwd("/Users/gregfaletto/Dropbox/Jacob and Greg/AISTATS 2023")
# Figure 8 (from Remark D.1) in the appendix of the version of the paper
# submitted to ICML

rm(list=ls())

library(ggplot2)

nsims <- 25

n <- 10^6
p <- 10

alpha <- c(0, 20)
beta <- rep(1, p)
K <- 3
stopifnot(length(alpha) == K - 1)

prob <- function(x, k){
  stopifnot(k %in% 1:K)
  stopifnot(length(x) == p)
  if(k == 1){
    return(plogis(alpha[k] + beta %*% x))
  } else if(k == K){
    return(1 - plogis(alpha[K - 1] + beta %*% x))
  }
  return(plogis(alpha[k] + beta %*% x) - plogis(alpha[k - 1] + beta %*% x))
}

t0 <- Sys.time()

mins <- rep(as.numeric(NA), nsims)

set.seed(105723)

for(j in 1:nsims){
  I_beta_beta <- array(as.numeric(NA), dim=c(p, p, n))

  I_beta_alpha_1 <- matrix(as.numeric(NA), nrow=p, ncol=n)

  I_alpha_1_alpha_1 <- rep(as.numeric(NA), n)

  rare_probs_j <- rep(as.numeric(NA), n)

  # Generate X
  X <- matrix(runif(n*p, -1, 1), nrow=n, ncol=p)

  for(i in 1:n){
    I_beta_beta[, , i] <- outer(X[i, ], X[i, ])*as.numeric((1 - prob(X[i, ], 1))*
                                                             (1 - prob(X[i, ], 2))*(1 - prob(X[i, ], 3)))

    I_beta_alpha_1[, i] <- X[i, ]*as.numeric((1 - prob(X[i, ], 3))*prob(X[i, ], 1)*
                                               (1 - prob(X[i, ], 1)))

    I_alpha_1_alpha_1[i] <- (1 - prob(X[i, ], 3))*prob(X[i, ], 1)*
      (1 - prob(X[i, ], 1))

    rare_probs_j <- prob(X[i, ], 3)
  }

  stopifnot(all(!is.na(I_beta_beta)))
  stopifnot(all(!is.na(I_beta_alpha_1)))
  stopifnot(all(!is.na(I_alpha_1_alpha_1)))
  stopifnot(all(!is.na(rare_probs_j)))
  stopifnot(all(I_alpha_1_alpha_1 > 0))
  stopifnot(all(rare_probs_j >= 0))
  stopifnot(all(rare_probs_j <= 1))

  I_beta_beta_mean <- apply(I_beta_beta, c(1, 2), mean)
  stopifnot(all(dim(I_beta_beta_mean) == c(p, p)))

  I_beta_alpha_1_mean <- rowMeans(I_beta_alpha_1)
  stopifnot(length(I_beta_alpha_1_mean) == p)

  I_alpha_1_alpha_1_mean <- mean(I_alpha_1_alpha_1)

  mat <- I_beta_beta_mean - 2*outer(I_beta_alpha_1_mean,
                                    I_beta_alpha_1_mean)/I_alpha_1_alpha_1_mean
  min <- min(eigen(mat)$values)

  mins[j] <- min

}


print("Total time:")
print(Sys.time() - t0)

print("mean(mins):")
print(mean(mins))
print("sum(mins > 0)/nsims:")
print(sum(mins > 0)/nsims)

df <- data.frame(lambda_min=mins)

# Figure 10
plot <- ggplot(df, aes(x=lambda_min)) + geom_boxplot() + coord_flip() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Minimum Eigenvalue")

print(plot)

print("Standard error:")

se <- sd(mins)/sqrt(nsims)

print(se)

print("95% confidence interval:")

print(round(c(mean(mins) - qnorm(.975)*se, mean(mins) + qnorm(.975)*se), 5))

print("pi_rare:")

pi_rare <- prob(rep(-1, p), K)

print(pi_rare)

print("Condition met?")

M <- sqrt(p)

print(pi_rare <= mean(mins)/(3*M^2*(2+M)))




