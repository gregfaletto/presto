rm(list=ls())

library(simulator)
library(MASS)
library(ggplot2)
library(parallel)

# Ran on 7 cores in figures in paper
n_cores <- detectCores() - 1

# for each setting, takes about one hour per simulation
nsims <- 1

n <- 10^6
p <- 10

beta <- rep(1, p)
K <- 3
rare_prob <- 0
bound <- 3

getRareIntcpt <- function(rare_prob, beta){
    # qlogis is inverse of plogis
    stopifnot(length(beta) == p)
    param <- qlogis(1 - rare_prob)
    # on distribution bounded to [-bound, bound]^p, largest beta_max can be is
    # beta %*% rep(bound, p)
    return(as.numeric(beta %*% rep(bound, p)) + param)
}

prob <- function(x, k, alpha, beta){
  stopifnot(k %in% 1:K)
  stopifnot(length(x) == p)
  stopifnot(length(alpha) == K - 1)
  if(k == 1){
    return(plogis(alpha[k] + beta %*% x))
  } else if(k == K){
    return(1 - plogis(alpha[K - 1] + beta %*% x))
  }
  return(plogis(alpha[k] + beta %*% x) - plogis(alpha[k - 1] + beta %*% x))
}

t0 <- Sys.time()

set.seed(105723)


eigensim_func <- function(n, p, rare_prob, beta, K, corr, nsim){

    # Get intercept
    stopifnot(K == 3)
    stopifnot(length(beta) == p)
    alpha <- c(0, 0)

    
    if(rare_prob == !0){
        alpha[2] <- getRareIntcpt(rare_prob, beta)
    }

    ret_list <- list()

    I_beta_beta <- array(as.numeric(NA), dim=c(p, p, n))

    I_beta_alpha_1 <- matrix(as.numeric(NA), nrow=p, ncol=n)

    I_alpha_1_alpha_1 <- rep(as.numeric(NA), n)

    rare_probs_j <- rep(as.numeric(NA), n)

    Sigma <- matrix(corr, p, p)
    diag(Sigma) <- 1

    for(j in 1:nsim){
        
        # Generate X
        X <- matrix(runif(n*p, -1, 1), nrow=n, ncol=p)
        X <- MASS::mvrnorm(n=n, mu=rep(0, p), Sigma=Sigma)
        # Truncate so X is bounded
        if(any(X > bound)){
            X[X > bound] <- bound
        }
        if(any(X < -1*bound)){
            X[X < -1*bound] <- -1*bound
        }
        
        for(i in 1:n){

            prob_1 <- prob(X[i, ], 1, alpha, beta)
            if(rare_prob == 0){
                prob_2 <- 1 - prob_1
                prob_3 <- 0
            } else{
                prob_2 <- prob(X[i, ], 2, alpha, beta)
                prob_3 <- prob(X[i, ], 3, alpha, beta)
            }

            stopifnot(abs(prob_1 + prob_2 + prob_3 - 1) < 10^(-8))

            I_beta_beta[, , i] <- outer(X[i, ], X[i, ])*as.numeric((1 - prob_1)*
                (1 - prob_2)*(1 - prob_3))

            I_beta_alpha_1[, i] <- X[i, ]*as.numeric((1 - prob_3)*prob_1*
                (1 - prob(X[i, ], 1, alpha, beta)))

            I_alpha_1_alpha_1[i] <- (1 - prob_3)*prob_1*(1 - prob_1)

            rare_probs_j[i] <- prob_3
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

        min_eigen <- min(eigen(mat)$values)

        stopifnot(length(min_eigen) == 1)
        stopifnot(is.numeric(min_eigen))

        # Minimum eigenvalue
        ret_list[[j]] <- min_eigen

        stopifnot(length(ret_list) == j)

    }

    stopifnot(length(ret_list) == nsim)

    return(ret_list)
}

eigenmodel <- function(n, p, rare_prob, beta, K, corr) {

    stopifnot(length(beta) == p)
    
    mod <- new_model(name = "eigenmodel"
        , label = "Eigenvalue Simulation Study"
        , params = list(n=n, p=p, rare_prob=rare_prob,
            beta=beta, K=K, corr=corr)
        , simulate = eigensim_func
    )
    return(mod)
}

eigenmethod <- new_method("eigenmethod", "Eigenvalue method",
    method = function(model, draw) {
    # Just return minimum eigenvalue
    stopifnot(length(draw) == 1)
    stopifnot(is.numeric(draw))
    return(draw)
    }
)

eigeneval <- new_metric("eigeneval", "Minimum Eigenvalue",
    metric = function(model, out) {
        stopifnot(length(out$out) == 1)
        stopifnot(is.numeric(out$out))
        return(out$out)
    }
)

eigeneval2 <- new_metric("eigeneval2", "Required Condition",
    metric = function(model, out) {
        stopifnot(length(out$out) == 1)
        stopifnot(is.numeric(out$out))

        min_eigen <- out$out

        pi_rare <- model$rare_prob

        M <- bound*sqrt(p)

        return(as.integer(pi_rare <= min_eigen/(3*M^2*(2+M))))
    }
)

eigensim2 <- new_simulation("eigensim2", "Theorem 2.3 Assumption Simulation") |>
        generate_model(eigenmodel, n=n, p=p,
            rare_prob=as.list(c(0, 10^(seq(-7, -5, by=1)))),
            beta=rep(1, p), K=K
            , corr=as.list(c(0, 0.25, 0.5, 0.75)),
            vary_along=c("rare_prob", "corr")
            ) |> 
        simulate_from_model(nsim = nsims,
            , index = 1:n_cores
        ) |>
        run_method(list(eigenmethod)
            # , parallel = list(socket_names=n_cores)
            ) |> 
        evaluate(list(eigeneval, eigeneval2))

save_simulation(eigensim2)

# Figure 13
fig_13 <- subset_simulation(eigensim2, corr==0) |>
    plot_eval_by("eigeneval", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0")

print(fig_13)

# Figure 14
fig_14 <- subset_simulation(eigensim2, corr==0.25) |>
    plot_eval_by("eigeneval", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.25")

print(fig_14)

# Figure 15
fig_15 <- subset_simulation(eigensim2, corr==0.5) |>
    plot_eval_by("eigeneval", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.5")

print(fig_15)

# Figure 16
fig_16 <- subset_simulation(eigensim2, corr==0.75) |>
    plot_eval_by("eigeneval", varying="rare_prob",
        include_zero=TRUE) +
    scale_x_log10() +
    geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.75")

print(fig_16)

# Plots for condition being satisfied

# Figure 17
fig_17 <- subset_simulation(eigensim2, corr==0) |>
    plot_eval_by("eigeneval2", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0")

print(fig_17)

# Figure 18
fig_18 <- subset_simulation(eigensim2, corr==0.25) |>
    plot_eval_by("eigeneval2", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.25")

print(fig_18)

# Figure 19
fig_19 <- subset_simulation(eigensim2, corr==0.5) |>
    plot_eval_by("eigeneval2", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.5")

print(fig_19)

# Figure 20
fig_20 <- subset_simulation(eigensim2, corr==0.75) |>
    plot_eval_by("eigeneval2", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.75 ")

print(fig_20)

print("Total time:")
print(Sys.time() - t0)




