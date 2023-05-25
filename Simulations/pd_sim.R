# setwd("/Users/gregfaletto/Dropbox/Jacob and Greg/AISTATS 2023")
# Figure 8 (from Remark D.1) in the appendix of the version of the paper
# submitted to ICML

rm(list=ls())

library(simulator)
library(MASS)
library(ggplot2)
library(parallel)

n_cores <- detectCores() - 1

# for all 6 settings, takes about one hour per simulation
nsims <- 1

n <- 10^6
p <- 10

# alpha <- c(0, 20)
beta <- rep(1, p)
K <- 3
# corr <- 0
# rare_prob <- 4*10^(-5)
rare_prob <- 0
bound <- 3

getRareIntcpt <- function(rare_prob, beta){
    # qlogis is inverse of plogis
    stopifnot(length(beta) == p)
    param <- qlogis(1 - rare_prob)
    # on distribution bounded to [-bound, bound]^p, largest beta_max can be is
    # beta %*% rep(5, p)
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

conf_int_agg_func <- function(ev){
     stopifnot(is.list(ev))

     values <- unlist(ev)
     stopifnot(is.numeric(values))
     n <- sum(!is.na(values))

     se <- sd(values, na.rm=TRUE)/sqrt(n)

     # print("max(values):")
     # print(max(values))

     # print("min(values):")
     # print(min(values))
     # print("se:")
     # print(se)

     margin <- se*qt(.975, df=n - 1)

     # print("margin:")
     # print(margin)

     # return(margin)
     return(se)
}

conf_int_agg <- new_aggregator("95% Margin of Error", conf_int_agg_func)

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

eigeneval3 <- new_metric("eigeneval3", "Positive Definite",
    metric = function(model, out) {
        stopifnot(length(out$out) == 1)
        stopifnot(is.numeric(out$out))

        return(as.integer(out$out > 0))
    }
)


# eigensim2
# eigensim

eigensim2 <- new_simulation("eigensim2", "Theorem 2.3 Assumption Simulation") |>
        generate_model(eigenmodel, n=n, p=p,
            rare_prob=as.list(c(0, 10^(seq(-7, -5, by=1)))),
            beta=rep(1, p), K=K
            , corr=as.list(c(0, 0.25, 0.5, 0.75)),
            vary_along=c("rare_prob", "corr")
            # , corr= 0.75,
            # vary_along=c("rare_prob")
            ) |> 
        simulate_from_model(nsim = nsims,
            , index = 1:n_cores
        ) |>
        run_method(list(eigenmethod)
            # , parallel = list(socket_names=n_cores)
            ) |> 
        evaluate(list(eigeneval, eigeneval2, eigeneval3))

save_simulation(eigensim2)

# eigensim <- load_simulation("eigensim")

# plot_eval(eigensim2, "eigeneval")

eigenplot4 <- plot_eval_by(eigensim2, "eigeneval", varying="rare_prob",
    include_zero=TRUE
    # , spread_aggregator=conf_int_agg
    ) +
    scale_x_log10() +
    # scale_y_log10() + 
    geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.75")

print(eigenplot4)

eigenplot1 <- subset_simulation(eigensim, corr==0) |>
    plot_eval_by("eigeneval", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    # scale_y_log10() + 
    geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0")

print(eigenplot1)

eigenplot2 <- subset_simulation(eigensim, corr==0.25) |>
    plot_eval_by("eigeneval", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    # scale_y_log10() + 
    geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.25")

print(eigenplot2)

eigenplot3 <- subset_simulation(eigensim, corr==0.5) |>
    plot_eval_by("eigeneval", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    # scale_y_log10() + 
    geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.5")

print(eigenplot3)





# Plots for condition being satisfied
cond_plot_1 <- subset_simulation(eigensim, corr==0) |>
    plot_eval_by("eigeneval2", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    # scale_y_log10() + 
    # geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0")

print(cond_plot_1)

cond_plot_2 <- subset_simulation(eigensim, corr==0.25) |>
    plot_eval_by("eigeneval2", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    # scale_y_log10() + 
    # geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.25")

print(cond_plot_2)

cond_plot_3 <- subset_simulation(eigensim, corr==0.5) |>
    plot_eval_by("eigeneval2", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    # scale_y_log10() + 
    # geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.5")

print(cond_plot_3)

cond_plot_4 <- subset_simulation(eigensim, corr==0.75) |>
    plot_eval_by("eigeneval2", varying="rare_prob",
    include_zero=TRUE) +
    scale_x_log10() +
    # scale_y_log10() + 
    # geom_hline(yintercept=0, color="red", linetype="dashed") +
    xlab("Rare Class Probability") + theme(legend.position="none") +
    ggtitle("Correlation = 0.75 ")

print(cond_plot_4)

print("Total time:")
print(Sys.time() - t0)

# print("mean(mins):")
# print(mean(mins))
# print("sum(mins > 0)/nsims:")
# print(sum(mins > 0)/nsims)

# df <- data.frame(lambda_min=mins)

# plot <- ggplot(df, aes(x=lambda_min)) + geom_boxplot() + coord_flip() +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   ylab("Minimum Eigenvalue")

# print(plot)

# print("Standard error:")

# se <- sd(mins)/sqrt(nsims)

# print(se)

# print("95% confidence interval:")

# print(round(c(mean(mins) - qnorm(.975)*se, mean(mins) + qnorm(.975)*se), 5))

# print("pi_rare:")

# pi_rare <- prob(rep(-1, p), K)

# print(pi_rare)






