rare_prob_mse_gen <- new_metric("rare_prob_mse_gen",
     "Rare Probability MSE", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          ret <- mse(p_hat, rare_probs)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          stopifnot(!is.na(ret))

          return(ret)
     }
)

mse <- function(x, y){
     
     stopifnot(is.numeric(x) | is.integer(x))
     stopifnot(is.numeric(y) | is.integer(y))
     stopifnot(all(!is.na(x)))
     stopifnot(all(!is.na(y)))
     n <- length(x)
     stopifnot(n >= 1)
     stopifnot(n == length(y))

     return(sum((x - y)^2)/n)
}

getEvalQuantities <- function(model, out, X){

     checkGetEvalQuantitiesInputs(model, out, X)

     # For methods that rely on ordinalNet, we simply generate predictions
     # using the built-in prediction function. Those methods don't return
     # alpha_hat, so that's how we can tell whether we can use the ordinalNet
     # prediction function.
     if(any(!is.na(out$alpha_hat))){

          # Confirm that paramter estimates needed for calculating probability
          # estimates are available in method output, as expected
          stopifnot("intercepts" %in% names(model@params))
          stopifnot(all(c("alpha_hat", "Beta_hat") %in% names(out)))

          stopifnot(length(out$alpha_hat) == model$K - 1)
          stopifnot(all(dim(out$Beta_hat) == c(model$p, model$K - 1)))
          stopifnot(all(!is.na(out$alpha_hat)))
          stopifnot(all(!is.na(out$Beta_hat)))

          if(model$p > 1){
               p_hat <- 1 - plogis(out$alpha_hat[model$K - 1] +
                    X %*% out$Beta_hat[, model$K - 1])[, 1]
          } else{
               p_hat <- 1 - plogis(out$alpha_hat[model$K - 1] +
                    out$Beta_hat[, model$K - 1]*X)
          }
     } else{
          # Confirm that the method output a valid ordinalNet object
          stopifnot("fit" %in% names(out))

          stopifnot(!is.na(out$fit))
          stopifnot("ordinalNet" %in% class(out$fit))

          # Having confirmed this is an ordinalNet object, use the provided
          # prediction function to get the rare class probailities

          p_hat <- predict(out$fit, newx=X)[, model$K]
     }

     # Verify predicted probabilities either way
     stopifnot(!is.na(p_hat))
     stopifnot(is.numeric(p_hat) | is.integer(p_hat))
     stopifnot(length(p_hat) == model$n)
     stopifnot(all(p_hat >= 0))
     stopifnot(all(p_hat <= 1))

     # Also get true rare probabilities as a reference
     mu <- model$intercepts[model$K - 1] + X %*%
          out$beta_final[, model$K - 1]

     stopifnot(length(mu) == model$n)

     rare_probs <- 1 - plogis(mu)

     stopifnot(all(!is.na(rare_probs)))
     stopifnot(all(rare_probs >= 0))
     stopifnot(all(rare_probs <= 1))
     stopifnot(length(rare_probs) == model$n)

     return(list(mu=mu, p_hat=p_hat, rare_probs=rare_probs))
}

checkGetEvalQuantitiesInputs <- function(model, out, X){
     stopifnot(all(c("n", "K", "p") %in% names(model@params)))
     stopifnot(all(c("X", "y", "beta_final") %in% names(out)))

     if(model$p > 1){
          stopifnot(model$n == nrow(out$X))
          stopifnot(model$p == ncol(out$X))
          stopifnot(model$p == ncol(X))
     } else{
          stopifnot(length(out$X) == model$n)
          stopifnot(length(X) == model$n)
     }
     
     stopifnot(model$K == length(unique(out$y)))

     stopifnot(model$n == length(out$y))
     stopifnot(model$K >= 3)

     stopifnot(all(dim(out$beta_final) == c(model$p, model$K - 1)))
     stopifnot(all(!is.na(out$beta_final)))

     stopifnot(length(model$intercepts) == model$K - 1)
}

# Average proportion of rare observations
prop_rare_obs <- new_metric("prop_rare_obs",
     "Average proportion of rare observations", metric = function(model, out) {

          stopifnot(all(c("n", "K") %in% names(model@params)))
          stopifnot("y" %in% names(out))

          number_rare_obs <- sum(out$y == model$K)

          stopifnot(all(!is.na(number_rare_obs)))
          stopifnot(all(number_rare_obs > 0))
          stopifnot(all(number_rare_obs < model$n))
          stopifnot(all(number_rare_obs == round(number_rare_obs)))

          ret <- number_rare_obs/model$n

          stopifnot(is.numeric(ret) | is.integer(ret))
          stopifnot(!is.na(ret))
          stopifnot(length(ret) == 1)
          stopifnot(ret >= 0)

          return(ret)
     }
)

cal_osce_gen_data_app <- new_metric("cal_osce_gen_data_app",
     "OSCE (Equal Frequency)", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- getOCE(actual=act_pred$actual, predicted=act_pred$predicted,
               n_bins=n_bins, error="squared")

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for cal_osce_gen_data_app")
          }

          return(ret)
     }
)

get_actual_pred_rare <- function(model, out){
     p <- ncol(out$X_ordnet_test)
     n <- nrow(out$X_ordnet_test)
     K <- length(unique(out$test_y))

     stopifnot(n == length(out$test_y))
     stopifnot(K >= 3)

     stopifnot(p > 1)

     if(any(!is.na(out$alpha_hat))){

          stopifnot(length(out$alpha_hat) == K - 1)
          stopifnot(all(dim(out$Beta_hat) == c(p, K - 1)))

          if(p > 1){
               stopifnot(nrow(out$X_ordnet_test) == n)
               stopifnot(ncol(out$X_ordnet_test) == p)
               p_hat <- 1 - plogis(out$alpha_hat[K - 1] +
                    out$X_ordnet_test %*% out$Beta_hat[, K - 1])[, 1]
          } else{
               stopifnot(length(out$X_ordnet_test) == n)
               p_hat <- 1 - plogis(out$alpha_hat[K - 1] +
                    out$Beta_hat[, K - 1]*out$X_ordnet_test)
          }

     } else{
          stopifnot(!is.na(out$fit))
          stopifnot("ordinalNet" %in% class(out$fit))

          p_hat <- predict(out$fit, newx=out$X_ordnet_test)[, K]
     }

     stopifnot(length(p_hat) == n)

     return(list(actual=as.numeric(out$test_y == model$resp_levels[K]),
          predicted=p_hat))
}

# Modified from https://github.com/cran/CalibratR/blob/master/R/getECE.R
getOCE <- function(actual, predicted, n_bins=10, error="squared",
     trim_end_bins=FALSE){

     # Check inputs
     checkOCEInputs(actual, predicted, n_bins, error, trim_end_bins)

     N <- length(actual)

     # Order predicted probabilities and corresponding labels
     idx <- order(predicted)
     pred_actual <- cbind(predicted[idx], actual[idx])

     stopifnot(nrow(pred_actual) == N)
     # Get matrix of indices for subsamples
     inds_mat <- partitionIndicesFreq(N, n_bins)

     stopifnot(nrow(inds_mat) == n_bins)

     if(trim_end_bins){
          inds_to_delete_later <- c(inds_mat[1, 1]:inds_mat[1, 2],
               inds_mat[n_bins, 1]:inds_mat[n_bins, 2])
          inds_mat <- inds_mat[2:(n_bins - 1), ]
     }
     
     # Vector that will contain the estimated squared errors of each observation
     S <- rep(as.numeric(NA), N)

     groups <- list()

     for(i in 1:nrow(inds_mat)){
          first_ind <- inds_mat[i, 1]
          last_ind <- inds_mat[i, 2]
     
          # Predicted probabilities for ith bin
          group_pred <- pred_actual[first_ind:last_ind, 1]
          # Actual outcomes for ith bin
          group_actual <- pred_actual[first_ind:last_ind, 2]

          n_i <- length(group_pred)

          # Observed proportions of succeses in bin i (estimated probability
          # for observations in bin i)
          observed <- mean(group_actual) #true fraction of pos.instances = prevalence in bin b

          # Get estimated squared errors for each observatino
          for(j in first_ind:last_ind){
               # Make sure this space in S is empty (didn't mess up indices)
               stopifnot(is.na(S[j]))
               if(error == "squared"){
                    S[j] <- (predicted[j] - observed)^2
               } else if(error == "absolute"){
                    S[j] <- abs(predicted[j] - observed)
               }
               
          }
     }
     if(trim_end_bins){
          S <- S[setdiff(1:N, inds_to_delete_later)]
     }
     stopifnot(all(!is.na(S)))
     stopifnot(all(S >= 0))
     return(mean(S))
}

checkOCEInputs <- function(actual, predicted, n_bins, error, trim_end_bins){

     stopifnot(!is.na(actual))
     stopifnot(is.numeric(actual) | is.integer(actual))
     N <- length(actual)
     stopifnot(all(actual %in% c(0, 1)))

     stopifnot(!is.na(predicted))
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(length(predicted) == N)
     stopifnot(all(predicted >= 0))
     stopifnot(all(predicted <= 1))

     stopifnot(!is.na(n_bins))
     stopifnot(is.integer(n_bins) | is.numeric(n_bins))
     stopifnot(n_bins > 0)
     stopifnot(n_bins == round(n_bins))

     stopifnot(!is.na(error))
     stopifnot(is.character(error))
     stopifnot(length(error) == 1)
     stopifnot(error %in% c("squared", "absolute"))

     stopifnot(!is.na(trim_end_bins))
     stopifnot(is.logical(trim_end_bins))
     stopifnot(length(trim_end_bins) == 1)

}

partitionIndicesFreq <- function(n, n_bins){
     n_bins <- min(n, n_bins)
     ret <- matrix(as.integer(NA), nrow=n_bins, ncol=2)
     remainder <- n %% n_bins
     quotient <- floor(n/n_bins)

     stopifnot(n == n_bins*quotient + remainder)

     ret[1, 1] <- 1
     if(remainder >= 1){
          ret[1, 2] <- quotient + 1
     } else{
          ret[1, 2] <- quotient
     }

     for(i in 2:n_bins){
          stopifnot(is.na(ret[i, 1]))
          stopifnot(is.na(ret[i, 2]))
          ret[i, 1] <- ret[i - 1, 2] + 1
          if(remainder >= i){
               ret[i, 2] <- ret[i - 1, 2] + quotient + 1
          } else{
               ret[i, 2] <- ret[i - 1, 2] + quotient
          }
     }
     # Check output
     ret <- checkIndicesMatOutput(ret, n)
     return(ret)
}

checkIndicesMatOutput <- function(mat, n){
     n_bins <- nrow(mat)
     stopifnot(ncol(mat) == 2)
     stopifnot(mat[1, 2] >= mat[1, 1])
     for(i in 2:n_bins){
          stopifnot(mat[i, 2] >= mat[i, 1])
          stopifnot(mat[i, 1] == mat[i - 1, 2] + 1)
     }
     if(mat[n_bins, 2] > n){
          mat[n_bins, 2] <- n
          mat[n_bins, 1] <- min(mat[n_bins, 1], n)
          for(i in (n_bins - 1):1){
               mat[i, 2] <- min(mat[i, 2], mat[i + 1, 1] - 1)
               mat[i, 1] <- min(mat[i, 1], mat[i, 2])
          }
     }

     # Check output
     stopifnot(all(!is.na(mat)))
     stopifnot(mat[1, 1] == 1)
     stopifnot(mat[n_bins, 2] == n)
     inds <- integer()
     for(i in 1:n_bins){
          inds <- c(inds, mat[i, 1]:mat[i, 2])
     }
     stopifnot(identical(inds, 1L:n))
     return(mat)
}