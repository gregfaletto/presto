## @knitr metrics

# just fdr (tweak cutoff?)
# area under partial precision-recall AUC curve
# are we predicting small probabilities well? big probabilities well?

# Beta MSE

beta_mse <- new_metric("beta_mse", "Beta MSE", metric = function(model, out) {
          beta_hat <- out$theta_hat[2:(model$p + 1)]
          stopifnot(length(beta_hat) == model$p)

          beta_mses <- (beta_hat - out$beta_final)^2
          stopifnot(all(!is.na(beta_mses)))
          stopifnot(all(beta_mses >= 0))

          return(sum(beta_mses)/model$p)
	}
)

rare_beta_mse <- new_metric("rare_beta_mse", "Rare Beta MSE",
     metric = function(model, out) {

          stopifnot(all(c("p", "K") %in% names(model@params)))
          stopifnot(all(c("Beta_hat", "beta_final") %in% names(out)))

          stopifnot(model$K >= 3)
          stopifnot(model$p > 1)

          stopifnot(all(dim(out$Beta_hat) == c(model$p, model$K - 1)))
          stopifnot(all(dim(out$beta_final) == c(model$p, model$K - 1)))

          stopifnot(all(!is.na(out$Beta_hat)))
          stopifnot(all(!is.na(out$beta_final)))

          beta_mses <- (out$Beta_hat[, model$K - 1] -
               out$beta_final[, model$K - 1])^2

          stopifnot(all(!is.na(beta_mses)))
          stopifnot(all(beta_mses >= 0))

          ret <- sum(beta_mses)/model$p

          stopifnot(is.numeric(ret) | is.integer(ret))
          stopifnot(length(ret) == 1)
          stopifnot(!is.na(ret))
          stopifnot(ret >= 0)

          return(ret)
     }
)


          

# # Common intercept MSE

# com_intcp_mse <- new_metric("com_intcp_mse", "Common intercept MSE",
#      metric = function(model, out) {
#           alpha_hat <- out$theta_hat[1]

#           common_intercept_mse <- (alpha_hat - model$intercepts[1])^2

#           stopifnot(!is.na(common_intercept_mse))
#           stopifnot(common_intercept_mse >= 0)

#           return(common_intercept_mse)
#      }
# )

# Rare intercept MSE

rare_intcp_mse <- new_metric("rare_intcp_mse", "Rare intercept MSE",
     metric = function(model, out) {

          stopifnot(all(c("K", "intercepts") %in% names(model@params)))
          stopifnot("alpha_hat" %in% names(out))

          stopifnot(!is.na(out$alpha_hat[model$K - 1]))
          stopifnot(!is.na(model$intercepts[model$K - 1]))
 
          rare_intercept_mse <- (out$alpha_hat[model$K - 1] -
               model$intercepts[model$K - 1])^2

          stopifnot(is.numeric(rare_intercept_mse) | is.integer(rare_intercept_mse))
          stopifnot(!is.na(rare_intercept_mse))
          stopifnot(length(rare_intercept_mse) == 1)
          stopifnot(rare_intercept_mse >= 0)

          return(rare_intercept_mse)
     }
)

# Average proportion of rare observations

prop_rare_obs <- new_metric("prop_rare_obs", "Average proportion of rare observations",
     metric = function(model, out) {

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

            # if(K >= 4){
            #     number_less_rare_obs[i] <- sum(y == K - 1)
            # }


# Get confusion matrix entries

conf_mat <- function(actual, predicted, threshold=0.5){
     # Verify inputs
     stopifnot(is.numeric(actual) | is.integer(actual))
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(all(actual %in% c(0, 1)))
     n <- length(actual)
     stopifnot(n == length(predicted))
     stopifnot(all(!is.na(actual)))
     stopifnot(all(!is.na(predicted)))

     counts <- integer(4)
     names(counts) <- c("true_neg", "true_pos", "false_neg", "false_pos")

     for(k in 0:1){
          inds_k <- actual == k
          stopifnot(sum(inds_k) >= 0)

          if(k == 0){
               counts["true_neg"]<- sum(predicted[inds_k] <= threshold)
               counts["false_pos"] <- sum(predicted[inds_k] > threshold)
          } else{
               counts["true_pos"] <- sum(predicted[inds_k] > threshold)
               counts["false_neg"] <- sum(predicted[inds_k] <= threshold)
          }
          
     }

     stopifnot(is.integer(counts))
     stopifnot(length(counts) == 4)
     stopifnot(identical(names(counts), c("true_neg", "true_pos", "false_neg",
          "false_pos")))
     stopifnot(sum(counts) == n)

     for(i in 1:4){
          stopifnot(!is.na(counts[i]))
          stopifnot(counts[i] == round(counts[i]))
          stopifnot(counts[i] >= 0)
          stopifnot(counts[i] <= n)
     }
     return(counts)
}

# False positive rate (low is good)--worked okay in simulations, better than logit but 
# worse than proportional odds by a statistically significant margin

fpr <- function(actual, predicted){

     conf_mat <- conf_mat(actual, predicted)

     ret <- conf_mat["false_pos"]/(conf_mat["false_pos"] + conf_mat["true_neg"])

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)
     stopifnot(ret <= 1)

     return(ret)
}

fpr_gen <- new_metric("fpr_gen", "False positive rate", metric = function(model,
     out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- fpr(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# False negative rate (lower is good)--PRESTO was better than proportional odds
# but worse than logit by statistically significant margins

fnr <- function(actual, predicted){

     conf_mat <- conf_mat(actual, predicted)

     ret <- conf_mat["false_neg"]/(conf_mat["false_neg"] + conf_mat["true_pos"])

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)
     stopifnot(ret <= 1)

     return(ret)
}

fnr_gen <- new_metric("fnr_gen", "False Negative Rate",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- fnr(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Cohen's kappa (higher is good) PRESTO better than proportional odds but
# worse than logit by statistially significant margins
ckap <- function(actual, predicted){
     # https://en.wikipedia.org/wiki/Cohen's_kappa#Binary_classification_confusion_matrix
     conf_mat <- conf_mat(actual, predicted)

     tp <- conf_mat["true_pos"]
     fp <- conf_mat["false_pos"]
     tn <- conf_mat["true_neg"]
     fn <- conf_mat["false_neg"]

     numerator <- 2*(tp*tn - fp*fn)
     denom <- (tp + fp)*(fp + tn) + (tp + fn)*(fn + tn)
     ret <- numerator/denom

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret <= 1)

     return(ret)
}

ckap_gen <- new_metric("ckap_gen", "Cohen's Kappa",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- ckap(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

ckap_gen_data_app <- new_metric("ckap_gen_data_app", "Cohen's Kappa",
     metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- ckap(actual=act_pred$actual, predicted=act_pred$predicted)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for ckap_gen_data_app")
          }

          return(ret)
     }
)

# MCC (higher is better)--PRESTO did worse than logit by statistically
# significant margins
mcc <- function(actual, predicted){
     # https://en.wikipedia.org/wiki/Phi_coefficient#Machine_learning
     conf_mat <- conf_mat(actual, predicted)

     tp <- conf_mat["true_pos"]
     fp <- conf_mat["false_pos"]
     tn <- conf_mat["true_neg"]
     fn <- conf_mat["false_neg"]

     numerator <- tp*tn - fp*fn

     denom <- sqrt(tp + fp)
     stopifnot(!is.na(denom))
     denom <- denom*sqrt(fp + tn)
     stopifnot(!is.na(denom))
     denom <- denom*sqrt(tp + fn)
     stopifnot(!is.na(denom))
     denom <- denom*sqrt(fn + tn)
     stopifnot(!is.na(denom))

     if(denom == 0){
          denom <- 1
     }
     ret <- numerator/denom

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret <= 1)

     return(ret)
}

mcc_gen <- new_metric("mcc_gen", "MCC",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- mcc(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

mcc_gen_data_app <- new_metric("mcc_gen_data_app", "MCC",
     metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- mcc(actual=act_pred$actual, predicted=act_pred$predicted)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for mcc_gen_data_app")
          }

          return(ret)
     }
)

# False discovery rate (lower is better)--DID PRETTY WELL IN SIMULATIONS! Mean
# lower for PRESTO than for both other methods, but confidence intervals
# overlapped. Difference between PRESTO and logit had a p-value of 9.140327e-05
# (two sample t test assuming unequal variances), difference between PRESTO and
# proportional odds had a p-value of 0.06302532.

fdr <- function(actual, predicted){

     conf_mat <- conf_mat(actual, predicted)

     if(conf_mat["false_pos"] + conf_mat["true_pos"] == 0){
          return(0)
     }
     ret <- conf_mat["false_pos"]/(conf_mat["false_pos"] + conf_mat["true_pos"])

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)
     stopifnot(ret <= 1)

     return(ret)
}

fdr_gen <- new_metric("fdr_gen", "False Discovery Rate",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- fdr(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

rare_prob_fdr_data_app <- new_metric("rare_prob_fdr_data_app",
     "Rare Event FDR", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          ret <- fdr(actual=act_pred$actual, predicted=act_pred$predicted)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for rare_prob_fdr_data_app")
          }

          return(ret)
     }
)

# Sensitivity/recall/true positive rate (higher is better)--DIDN'T WORK IN
# SIMULATIONS (PRESTO was higher than proportional odds but lower than logistic
# regression)

tpr <- function(actual, predicted){

     conf_mat <- conf_mat(actual, predicted)

     ret <- conf_mat["true_pos"]/(conf_mat["false_neg"] + conf_mat["true_pos"])

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)
     stopifnot(ret <= 1)

     return(ret)
}

tpr_gen <- new_metric("tpr_gen", "Sensitivity", metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- tpr(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Specificity/true negative rate (higher is better)--did pretty well in
# simulations, but a little worse than proportional odds (non-overlapping
# confidence intervals)

tnr <- function(actual, predicted){

     conf_mat <- conf_mat(actual, predicted)

     n <- length(actual)

     ret <- conf_mat["true_neg"]/(conf_mat["false_pos"] + conf_mat["true_neg"])

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)
     stopifnot(ret <= 1)

     return(ret)
}

tnr_gen <- new_metric("tnr_gen", "Specificity",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- tnr(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Precision (higher is better)--PRESTO did better than proportional odds but
# worse than logit. A little surprising because if there is at least one
# predicted positive then precision = 1 - fdr, but I think in an appreciable
# number of simulations at least somet methods might have predicted no positives

precision <- function(actual, predicted){

     conf_mat <- conf_mat(actual, predicted)

     n <- length(actual)

     denom <- max(conf_mat["false_pos"] + conf_mat["true_pos"], 1)

     ret <- conf_mat["true_pos"]/denom

     if(conf_mat["false_pos"] + conf_mat["true_pos"] == 0){
          stopifnot(ret == 0)
     }

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)
     stopifnot(ret <= 1)

     return(ret)
}

precision_gen <- new_metric("precision_gen", "Precision",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- precision(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

rare_prob_precision_data_app <- new_metric("rare_prob_precision_data_app",
     "Rare Event Precision", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          ret <- precision(actual=act_pred$actual, predicted=act_pred$predicted)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for rare_prob_precision_data_app")
          }

          return(ret)
     }
)

# F score (higher is better)--PRESTO did much worse than logit
fscore <- function(actual, predicted){
     precision <- precision(actual, predicted)

     if(precision == 0){
          return(0)
     }

     recall <- tpr(actual, predicted)

     if(recall == 0){
          return(0)
     }

     ret <- 2*precision*recall/(precision + recall)

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)
     stopifnot(ret <= 1)

     return(ret)
}

fscore_gen <- new_metric("fscore_gen", "F Score",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- fscore(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Adjusted GM (higher is better)--didn't work in simulations (worse than logit,
# better than proportional odds)
# https://pubmed.ncbi.nlm.nih.gov/22809416/
adj_gm <- function(actual, predicted){

     conf_mat <- conf_mat(actual, predicted)

     n <- length(actual)

     Q <- (conf_mat["true_neg"] + conf_mat["false_pos"])/n

     rec <- tpr(actual, predicted)

     if(rec == 0){
          return(0)
     }

     spe <- tnr(actual, predicted)

     gm <- sqrt(rec*spe)

     ret <- (gm + spe*Q)/(1 + Q)

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)
     stopifnot(ret <= 1)

     return(ret)
}

adj_gm_gen <- new_metric("adj_gm_gen", "Adjusted GM",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- adj_gm(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

library(CalibratR)

# On test set: for n_bins = 25, PRESTO does better than logit but worse than
# proportional odds. For n_bins = 25, same thing, PRESTO worse than proportional
# odds by a statistically significant margin.

# Uniform differences:
# Calibration Expected Calibration Error (low is good)--with n_bins=25, 
# PRESTO is statistically tied with PO but loses by a statistically significant
# margin to logit. Similar results with n_bins=100 or 250.
cal_gen <- new_metric("cal_gen", "ECE",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- getECE(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat, n_bins=100)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

cal_gen_data_app <- new_metric("cal_gen_data_app", "ECE",
     metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- getECE(actual=act_pred$actual, predicted=act_pred$predicted,
               n_bins=n_bins)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for cal_gen_data_app")
          }

          return(ret)
     }
)

# Modified from https://github.com/cran/CalibratR/blob/master/R/getECE.R
getESCE <- function(actual, predicted, n_bins=10){ #equal frequency bins

  predicted <- predicted
  labels <- actual
  idx <- order(predicted)
  pred_actual <- (cbind(predicted[idx], labels[idx]))

  N <- nrow(pred_actual)
  rest <- N%%n_bins
  S <- 0
  W <- c()
  B <- min(N,n_bins) #if less then n_bins elements in data set, then use that number of bins
  groups <- list()

  for (i in 1:B){ #i von 1 bis B
    if (i <= rest){ #put rest elements into each bin
      group_pred <- (pred_actual[(((i-1)*ceiling(N/n_bins)+1) : (i*ceiling(N/n_bins))),1])
      group_actual <- (pred_actual[(((i-1)*ceiling(N/n_bins)+1) : (i*ceiling(N/n_bins))),2])
    }
    else {
      group_pred <- (pred_actual[((rest+(i-1)*floor(N/n_bins)+1) : (rest+i*floor(N/n_bins))),1])#group size=N/B
      group_actual <- (pred_actual[((rest+(i-1)*floor(N/n_bins)+1) : (rest+i*floor(N/n_bins))),2])
      }

    n_ <- length(group_pred)
    expected <- mean(group_pred) #mean of predictions in bin b
    observed <- mean(group_actual) #true fraction of pos.instances = prevalence in bin b

    S[i] <- (observed-expected)^2 #squared difference of observed value-predicted value in bin
    W[i] <- n_/N #empirical frequence of all instances that fall into bin i, should be equal when using equal freq binning approach
    groups[[i]] <- group_pred

  }

  mean_prediction <- lapply(groups, mean)
  min_group <- lapply(groups, min)
  max_group <- lapply(groups, max)

  res <- t(S)%*%W
  return(as.numeric(res))
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

getCutoffsLinearWidth <- function(probs, n_bins){
     max_prob <- max(probs)
     min_prob <- min(probs)

     bin_width <- (max_prob - min_prob)/n_bins

     bin_cutoffs <- numeric(n_bins + 1)
     bin_cutoffs[1] <- min_prob 
     bin_cutoffs[n_bins + 1] <- max_prob
     for(i in 2:n_bins){
          bin_cutoffs[i] <- bin_cutoffs[i - 1] + bin_width
     }
     stopifnot(identical(bin_cutoffs, sort(bin_cutoffs)))

     return(bin_cutoffs)
}


partitionIndicesWidth <- function(n, n_bins, probs, width="linear"){
     stopifnot(n == length(probs))
     stopifnot(!is.na(width))
     stopifnot(is.character(width))
     stopifnot(length(width) == 1)
     stopifnot(width %in% c("linear", "log"))

     n_bins <- min(n, n_bins)
     ret <- matrix(as.integer(NA), nrow=n_bins, ncol=2)
     
     if(width == "linear"){
          bin_cutoffs <- getCutoffsLinearWidth(probs, n_bins)
     } else if(width == "log"){
          if(any(probs <= 0)){
               # Find smallest probability greater than 0, and set any 
               # probabilities equal to 0 equal to this
               nonzero_probs <- probs[probs > 0]
               probs[probs <= 0] <- min(nonzero_probs)
          }
          if(any(probs >= 1)){
               # Find smallest probability less than 1, and set any 
               # probabilities equal to 1 equal to this
               nonone_probs <- probs[probs < 1]
               probs[probs >= 1] <- max(nonone_probs)
          }
          stopifnot(all(probs > 0))
          stopifnot(all(probs < 1))
          bin_cutoffs <- exp(getCutoffsLinearWidth(log(probs), n_bins))
          # Correct possible issues due to rounding
          bin_cutoffs[n_bins + 1] <- max(bin_cutoffs[n_bins + 1], max(probs))
          bin_cutoffs[1] <- min(bin_cutoffs[1], min(probs))
     }

     ret[1, 1] <- 1
     ret[1, 2] <- max(sum(probs <= bin_cutoffs[2]), 1)

     for(i in 2:n_bins){
          stopifnot(is.na(ret[i, 1]))
          stopifnot(is.na(ret[i, 2]))
          ret[i, 1] <- ret[i - 1, 2] + 1
          ret[i, 2] <- max(sum(probs <= bin_cutoffs[i + 1]), ret[i, 1])
     }
     if(ret[i, 2] == n -1){
          ret[i, 2] <- n
     }
     stopifnot(ret[i, 2] >= n)
     # Make sure last entry is at most n, and correct backwards if needed
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

checkOCEInputs <- function(actual, predicted, n_bins, error, binning,
     trim_end_bins){

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

     stopifnot(!is.na(binning))
     stopifnot(is.character(binning))
     stopifnot(length(binning) == 1)
     stopifnot(binning %in% c("equal_freq", "equal_lin_width",
          "equal_log_width"))

     stopifnot(!is.na(trim_end_bins))
     stopifnot(is.logical(trim_end_bins))
     stopifnot(length(trim_end_bins) == 1)

}

# Modified from https://github.com/cran/CalibratR/blob/master/R/getECE.R
getOCE <- function(actual, predicted, n_bins=10, error="squared",
     binning="equal_freq", trim_end_bins=FALSE){

     # Check inputs
     checkOCEInputs(actual, predicted, n_bins, error, binning, trim_end_bins)

     N <- length(actual)

     # Order predicted probabilities and corresponding labels
     idx <- order(predicted)
     pred_actual <- cbind(predicted[idx], actual[idx])

     stopifnot(nrow(pred_actual) == N)
     # Get matrix of indices for subsamples
     if(binning == "equal_freq"){
          inds_mat <- partitionIndicesFreq(N, n_bins)
     } else if(binning == "equal_lin_width"){
          inds_mat <- partitionIndicesWidth(N, n_bins, pred_actual[, 1],
               width="linear")
     } else if(binning == "equal_log_width"){
          inds_mat <- partitionIndicesWidth(N, n_bins, pred_actual[, 1],
               width="log")
     }

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

# For n_bins = 10, 25, 100, PRESTO does worse than both proportional odds and
# logit
cal_esce_gen <- new_metric("cal_esce_gen", "ESCE",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- getESCE(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat, n_bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

cal_esce_gen_data_app <- new_metric("cal_esce_gen_data_app", "ESCE",
     metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- getESCE(actual=act_pred$actual, predicted=act_pred$predicted,
               n_bins=n_bins)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for cal_esce_gen_data_app")
          }

          return(ret)
     }
)

# With n_bins = 10, PRESTO does significantly better than logit but significantly
# worse than proportional odds. With n_bins = 25, similar with larger gap between
# proportional odds and PRESTO. Similar with n_bins = 100.
cal_osce_gen <- new_metric("cal_osce_gen", "OSCE",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- getOCE(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat, n_bins=25, error="squared",
               binning="equal_freq")

          stopifnot(!is.na(ret))

          return(ret)
     }
)

cal_osce_gen_data_app <- new_metric("cal_osce_gen_data_app",
     "OSCE (Equal Frequency)", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- getOCE(actual=act_pred$actual, predicted=act_pred$predicted,
               n_bins=n_bins, error="squared", binning="equal_freq")

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for cal_osce_gen_data_app")
          }

          return(ret)
     }
)

# PRESTO worse than proportional odds at statistically significant difference;
# better than logit
cal_osce_high_gen <- new_metric("cal_osce_high_gen", "OSCE (High prob)",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          labels <- as.numeric(out$y == as.character(model$K))

          p_rare <- mean(labels)

          # likely_inds <- p_hat >= 0.01
          likely_inds <- p_hat >= p_rare

          ret <- getOCE(actual=labels[likely_inds],
               predicted=p_hat[likely_inds], n_bins=100, error="squared",
               binning="equal_freq")

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# n_bins = 10: PRESTO significantly better that logit but worse than
# proportional odds. Similar for n_bins=25, n_bins = 100.
cal_osce_trim_gen <- new_metric("cal_osce_trim_gen", "OSCE",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- getOCE(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat, n_bins=100, error="squared",
               binning="equal_freq", trim_end_bins=TRUE)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Seems to work worse than squared error
cal_oce_gen <- new_metric("cal_oce_gen", "OCE",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- getOCE(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat, n_bins=25, error="absolute",
               binning="equal_freq")

          stopifnot(!is.na(ret))

          return(ret)
     }
)

cal_oce_gen_data_app <- new_metric("cal_oce_gen_data_app", "OCE",
     metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- getOCE(actual=act_pred$actual, predicted=act_pred$predicted,
               n_bins=n_bins, error="absolute", binning="equal_freq")

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for cal_oce_gen_data_app")
          }

          return(ret)
     }
)

cal_oce_high_gen <- new_metric("cal_oce_high_gen", "OCE (High prob)",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          labels <- as.numeric(out$y == as.character(model$K))

          # p_rare <- mean(labels)

          likely_inds <- p_hat >= 0.01

          ret <- getOCE(actual=labels[likely_inds],
               predicted=p_hat[likely_inds], n_bins=25, error="absolute",
               binning="equal_freq")

          stopifnot(!is.na(ret))

          return(ret)
     }
)

cal_osce_width_gen <- new_metric("cal_osce_width_gen", "OSCE (Equal Width)",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- getOCE(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat, n_bins=25, error="squared",
               binning="equal_lin_width")

          stopifnot(!is.na(ret))

          return(ret)
     }
)

cal_osce_width_gen_data_app <- new_metric("cal_osce_width_gen_data_app",
     "OSCE (Equal Width)", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- getOCE(actual=act_pred$actual, predicted=act_pred$predicted,
               n_bins=n_bins, error="squared", binning="equal_lin_width")

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for cal_osce_width_gen_data_app")
          }

          return(ret)
     }
)

# n_bins = 10: PRESTO does better than logit but worse than proportional odds
# by significant margins (weird because proportional odds is worse than both
# methods when looking at MSE of actual probabilities). n_bins = 25: qualitatively
# similar. also for n_bins = 100.
cal_osce_log_width_gen <- new_metric("cal_osce_log_width_gen",
     "OSCE (Equal Log Width)", metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- getOCE(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat, n_bins=25, error="squared",
               binning="equal_log_width")

          stopifnot(!is.na(ret))

          return(ret)
     }
)

cal_osce_log_width_gen_data_app <- new_metric("cal_osce_log_width_gen_data_app",
     "OSCE (Equal Log Width)", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- getOCE(actual=act_pred$actual, predicted=act_pred$predicted,
               n_bins=n_bins, error="squared", binning="equal_log_width")

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for cal_osce_log_width_gen_data_app")
          }

          return(ret)
     }
)

# Significantly worse for PRESTO than proportional odds
cal_osce_high_log_width_gen  <- new_metric("cal_osce_high_log_width_gen",
     "OSCE (High prob, equal log width)",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          labels <- as.numeric(out$y == as.character(model$K))

          # p_rare <- mean(labels)

          likely_inds <- p_hat >= 0.01

          ret <- getOCE(actual=labels[likely_inds],
               predicted=p_hat[likely_inds], n_bins=100, error="squared",
               binning="equal_log_width")

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# n_bins = 10: PRESTO does better than logit but worse than proportional odds
# by significant margins. Similar for n_bins = 25, 100.
cal_osce_log_width_trim_gen <- new_metric("cal_osce_log_width_trim_gen",
     "OSCE (Equal Log Width)", metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- getOCE(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat, n_bins=100, error="squared",
               binning="equal_log_width", trim_end_bins=TRUE)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# This function is documented in the R package, but it looks like the authors
# forgot to export it. From
# https://github.com/cran/CalibratR/blob/master/R/getECE.R
get_ECE_equal_width <- function(actual, predicted, bins=10){ #equal width bins

  pred_actual <- cbind(predicted, actual)

  if(all(predicted<=1) && all(predicted>=0)){
    hist_x <- hist(pred_actual[,1], breaks=seq(0,1,1/bins), plot=F)
  }
  else{
    hist_x <- hist(pred_actual[,1], breaks=bins, plot=F)
  }

  breaks_y <- hist_x$breaks
  y_true <- hist(subset(pred_actual[,1], pred_actual[,2]=="1"), breaks=breaks_y, plot=F)
  divided <- cut(pred_actual[,1], breaks=c(hist_x$breaks), label = seq(1,length(y_true$mids)), include.lowest = T)
  prediction_in_bin <- list()
  expected <- c()

    for (i in as.numeric(levels(divided))){
      prediction_in_bin[[i]] <- pred_actual[which(divided==i),1]
      expected[i] <- mean(prediction_in_bin[[i]]) #mean prediction in that bin
      #expected[i] <- hist_x$mids[i] #hist mids as mean prediction in that bin
    }


  counts_all <- hist_x$counts
  counts_true <- y_true$counts
  zeros <- which(counts_all==0)

  prevalence <- counts_true/counts_all
  prevalence[zeros] <- 0 #set prevalence to 0 when no observations are in the bin
  expected[zeros] <- hist_x$mids[zeros] #set expectation to the mid bin point, when no elements are in bin

  S_2 <- abs(prevalence-expected)
  W_2 <- counts_all/(length(predicted))


  return(as.numeric(t(S_2)%*%W_2))
}

# PRESTO does worse than both logit and proportional odds for n_bins = 100, 25,
# 10
cal_gen_ew <- new_metric("cal_gen_ew", "ECE (EW)",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- get_ECE_equal_width(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat, bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Modified from https://github.com/cran/CalibratR/blob/master/R/getMCE.R
getMSCE <- function(actual, predicted, n_bins=10){

  predicted <- predicted
  labels <- actual
  idx <- order(predicted)
  pred_actual <- (cbind(predicted[idx], actual[idx]))
  N <- nrow(pred_actual)
  rest <- N%%n_bins
  B <- min(N,n_bins)

  S <- 0
  W <- c()
  for (i in 1:B){ #i von 1 bis B
    if (i <= rest){ #put rest elements into each bin
      group_pred <- (pred_actual[(((i-1)*ceiling(N/n_bins)+1) : (i*ceiling(N/n_bins))),1])
      group_actual <- (pred_actual[(((i-1)*ceiling(N/n_bins)+1) : (i*ceiling(N/n_bins))),2])
    }
    else {
      group_pred <- (pred_actual[((rest+(i-1)*floor(N/n_bins)+1) : (rest+i*floor(N/n_bins))),1])#group size=N/B
      group_actual <- (pred_actual[((rest+(i-1)*floor(N/n_bins)+1) : (rest+i*floor(N/n_bins))),2])
    }

    n <- length(group_pred)
    expected <- mean(group_pred) #mean of predictions in bin b
    observed <- mean(group_actual) #true fraction of pos.instances = prevalence in bin b

    S[i] <- (observed-expected)^2 # squared difference of observed value-predicted value in bin
    W[i] <- n/N #empirical frequence of all instances that fall into bin i, should be pretty much the same among all bins
  }

  res <- max(S*W)
  return(res)
}

# On test set: PRESTO did worse than other methods but not by statistically significant
# margin. (Could try larger test n, or re-writing code so that every method
# uses the same test set in order to reduce variance.) (n_bins = 10)
# For n_bins = 25, PRESTO did a little better than logistic regression, otherwise
# similar.

# Closest results for n_bins = 25 rather than 100 or 10, but PRESTO still worse
# than other methods
cal_msce_gen <- new_metric("cal_msce_gen", "MSCE",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- getMSCE(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat, n_bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)



# From https://github.com/cran/CalibratR/blob/master/R/getMCE.R
get_MCE_equal_width <- function(actual, predicted, bins=10){ #equal width bins

  predicted <- predicted
  labels <- actual
  idx <- order(predicted)
  pred_actual <- (cbind(predicted[idx], labels[idx]))

  hist_x <- hist(pred_actual[,1],breaks=bins, plot=F)
  breaks_y <- hist_x$breaks
  y_true <- hist(subset(pred_actual[,1], pred_actual[,2]=="1"),breaks=breaks_y, plot=F)
  divided <- cut(pred_actual[,1], breaks=c(hist_x$breaks),label = seq(1,length(y_true$mids)),include.lowest = T)
  prediction_in_bin <- list()
  expected <- c()

  for (i in as.numeric(levels(divided))){
    prediction_in_bin[[i]] <- pred_actual[which(divided==i),1]
    #expected[i] <- hist_x$mids[i] #mean prediction in that bin
    expected[i] <- mean(pred_actual[which(divided==i),1]) #mean prediction in that bin
  }

  counts_all <- hist_x$counts
  counts_true <- y_true$counts
  zeros <- which(counts_all==0)

  prevalence <- counts_true/counts_all
  prevalence[zeros] <- 0 #set prevalence to 0 when no observations are in the bin
  expected[zeros] <- hist_x$mids[zeros] #set expectation to the mid bin point, when no elements are in bin

  S_2 <- abs(prevalence-expected)
  W_2 <- counts_all/(length(predicted))
  return(max(S_2*W_2))
}

# On test set: n_bins = 25, PRESTO did worse than proportional odds but better
# than logit.

# Uniform differences:
# PRESTO is statistically tied with PO (p = 0.42) but loses by a statistically
# significant margin to logit.
cal_max_gen <- new_metric("cal_max_gen", "MCE",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- getMCE(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat, n_bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

cal_max_gen_data_app <- new_metric("cal_max_gen_data_app", "MCE",
     metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          n_bins <- max(round(length(act_pred$actual)/100), 10)

          ret <- getMCE(actual=act_pred$actual, predicted=act_pred$predicted,
               n_bins=n_bins)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for cal_max_gen_data_app")
          }

          return(ret)
     }
)

# Uniform differences:
# (Low is good) PRESTO beats logit by statistically significant margins and
# ties statistically with proportional odds
cal_high_prob_gen <- new_metric("cal_high_prob_gen", "ECE (High Prob)",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          labels <- as.numeric(out$y == as.character(model$K))

          p_rare <- mean(labels)

          likely_inds <- p_hat >= p_rare

          ret <- getECE(actual=labels[likely_inds],
               predicted=p_hat[likely_inds], n_bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Uniform differences:
# (Low is good) PRESTO beats logit by a statistically significant margin and
# ties with proportional odds (p = 0.089)
cal_1_prob_gen <- new_metric("cal_1_prob_gen", "ECE (1 Prob)",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          labels <- as.numeric(out$y == as.character(model$K))

          likely_inds <- p_hat >= 0.01

          ret <- getECE(actual=labels[likely_inds],
               predicted=p_hat[likely_inds], n_bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Uniform differences:
# (Low is good) PRESTO loses to logit by a statistically significant margin but
# beats proportional odds by a statistically significant margin
cal_low_prob_gen <- new_metric("cal_low_prob_gen", "ECE (Low Prob)",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          labels <- as.numeric(out$y == as.character(model$K))

          p_rare <- mean(labels)

          likely_inds <- p_hat < p_rare

          ret <- getECE(actual=labels[likely_inds],
               predicted=p_hat[likely_inds], n_bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Uniform differences:
# PRESTO beats logit by a statistically significant margin and statistically
# ties with proportional odds (p-value 0.249)
cal_max_high_prob_gen <- new_metric("cal_max_high_prob_gen", "MCE (High Prob)",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          labels <- as.numeric(out$y == as.character(model$K))

          p_rare <- mean(labels)

          likely_inds <- p_hat >= p_rare

          ret <- getMCE(actual=labels[likely_inds],
               predicted=p_hat[likely_inds], n_bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Uniform differences:
# PRESTO beats logit by a statistically significant margin and ties statistically
# with proportional odds (p = 0.241)
cal_max_1_prob_gen <- new_metric("cal_max_1_prob_gen", "MCE (1 Prob)",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          labels <- as.numeric(out$y == as.character(model$K))

          likely_inds <- p_hat >= .01

          ret <- getMCE(actual=labels[likely_inds],
               predicted=p_hat[likely_inds], n_bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Uniform differences:
# PRESTO loses to logit by a statistially significant margin (p = 0.021) but
# beats proportional odds by a statistically significant margin (p = 2.8e-14)
cal_max_low_prob_gen <- new_metric("cal_max_low_prob_gen", "MCE (Low Prob)",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          labels <- as.numeric(out$y == as.character(model$K))

          p_rare <- mean(labels)

          likely_inds <- p_hat < p_rare

          ret <- getMCE(actual=labels[likely_inds],
               predicted=p_hat[likely_inds], n_bins=25)

          stopifnot(!is.na(ret))

          return(ret)
     }
)



# PR AUC (higher is better, supposedly better when there is a rare class)--
# PRESTO did better than proportional odds but worse than logit (p-value
# 4.059465e-06)

# library(PRROC)

pr_auc_gen <- new_metric("pr_auc_gen", "PR AUC", metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          pos_cases <- out$y == as.character(model$K)

          stopifnot(length(pos_cases) == length(p_hat))

          response <- as.numeric(out$y == as.character(model$K))

          ret <- pr.curve(scores.class0=qlogis(p_hat[pos_cases]),
               scores.class1=qlogis(p_hat[!pos_cases]))$auc.integral

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# AUC (higher is better)--PRESTO was higher than proportional odds but lower
# than logit
# library(pROC)

auc_gen <- new_metric("auc_gen", "AUC",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- suppressMessages(auc(response=as.numeric(out$y == as.character(model$K)),
               predictor=p_hat))

          stopifnot(!is.na(ret))

          return(ret)
     }
)

probs <- new_metric("probs", "Probabilities",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          return(p_hat)
     }
)

labels <- new_metric("labels", "Labels",
     metric = function(model, out) {

          y <- as.numeric(out$y == as.character(model$K))

          return(y)
     }
)

# Partial AUC (higher is better)--PRESTO did worse than logit but better than
# proportional odds with lower bound on sensitivity of 0.5, similar with 0.75,
# similar with 0.9

part_auc_gen <- new_metric("part_auc_gen", "Partial AUC",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- suppressMessages(auc(response=as.numeric(out$y == as.character(model$K)),
               predictor=p_hat, partial.auc=c(0.99, 1),
               partial.auc.focus="sensitivity"))

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Specificity Partial AUC (higher is better)--PRESTO did worse than logit
# but better than proportional odds with lower bound on specificity of 0.5,
# similar with 0.75, similar with 0.9)

spec_part_auc_gen <- new_metric("spec_part_auc_gen", "Specificity Partial AUC",
     metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- suppressMessages(auc(response=as.numeric(out$y == as.character(model$K)),
               predictor=p_hat, partial.auc=c(0.9, 1),
               partial.auc.focus="specificity"))

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Rare class accuracy (higher is better)--presto did better than proportional
# odds but worse than logit in simulations, but not by statistically significant
# margin (p-value: 0.254751)

accuracy <- function(actual, predicted){

     conf_mat <- conf_mat(actual, predicted)

     n <- length(actual)

     ret <- (conf_mat["true_pos"] + conf_mat["true_neg"])/n

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)
     stopifnot(ret <= 1)

     return(ret)

     # stopifnot(is.numeric(actual) | is.integer(actual))
     # stopifnot(is.numeric(predicted) | is.integer(predicted))
     # stopifnot(all(actual %in% c(0, 1)))
     # stopifnot(length(actual) == length(predicted))
     # stopifnot(all(!is.na(actual)))
     # stopifnot(all(!is.na(predicted)))

     # accurate <- 0

     # for(k in 0:1){
     #      inds_k <- actual == k
     #      stopifnot(sum(inds_k) >= 0)
     #      if(sum(inds_k) == 0){
     #           next
     #      }
     #      if(k == 0){
     #           accurate <- accurate + sum(predicted[inds_k] <= 0.5)
     #      } else{
     #           accurate <- accurate + sum(predicted[inds_k] > 0.5)
     #      }
          
     # }

     # stopifnot(!is.na(accurate))
     # stopifnot(is.numeric(accurate) | is.integer(accurate))
     # stopifnot(length(accurate) == 1)
     # stopifnot(accurate >= 0)

     # ret <- accurate/length(actual)

     # stopifnot(!is.na(ret))
     # stopifnot(is.numeric(ret) | is.integer(ret))
     # stopifnot(length(ret) == 1)
     # stopifnot(ret >= 0)
     # stopifnot(ret <= 1)

     # return(ret)
}

generateTestData <- function(model, out, n_tries=500){
     
     num_iters <- 0
     no_X <- TRUE

     while((num_iters < n_tries) & no_X){
          dat <- sim_data(n=model$n, p=model$p, K=model$K,
               intercepts=model$intercepts, beta_mat=out$beta_final)

          # Was an X generated so that every probability for every class
          # is nonnegative?
          no_X <- any(is.na(dat))

          num_iters <- num_iters + 1
     }
     
     if(no_X){
          stop(paste("No X generated with all nonnegative probabilities in",
          n_tries, "attempts."))
     }

     return(list(X=dat$X, y=dat$y, probs=dat$probs))
}

# On test set: PRESTO outperformed other methods but not by statistically significant
# margin. (Could try larger test n, or re-writing code so that every method
# uses the same test set in order to reduce variance.) True probability method
# didn't do significantly better than presto either.
acc_gen <- new_metric("acc_gen", "Accuracy",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- accuracy(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

rare_prob_acc_data_app <- new_metric("rare_prob_acc_data_app",
     "Rare Event Accuracy", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          ret <- accuracy(actual=act_pred$actual, predicted=act_pred$predicted)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for rare_prob_acc_data_app")
          }

          return(ret)
     }
)

# On test set: PRESTO outperformed other methods but not by statistically significant
# margin. (Could try larger test n, or re-writing code so that every method
# uses the same test set in order to reduce variance.) True probability method
# didn't do significantly better than presto either.

# Rare class probability Brier score (low is good)--presto did worse than logit
# in simulations (p-value 0.03841468) but better than proportional odds
brier <- function(actual, predicted){
     stopifnot(is.numeric(actual) | is.integer(actual))
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(all(actual %in% c(0, 1)))
     stopifnot(length(actual) == length(predicted))
     stopifnot(all(!is.na(actual)))
     stopifnot(all(!is.na(predicted)))

     loss <- 0

     for(k in 0:1){
          inds_k <- actual == k
          stopifnot(sum(inds_k) >= 0)
          if(sum(inds_k) == 0){
               next
          }
          loss <- loss + sum((k - predicted[inds_k])^2)
     }

     stopifnot(!is.na(loss))
     stopifnot(is.numeric(loss) | is.integer(loss))
     stopifnot(length(loss) == 1)
     stopifnot(loss >= 0)

     ret <- loss/length(actual)

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)

     return(ret)
}

brier_gen <- new_metric("brier_gen", "Brier Score",
     metric = function(model, out) {

          # Generate test X and y
          test_data <- generateTestData(model, out)

          p_hat <- getEvalQuantities(model, out, test_data$X)$p_hat

          ret <- brier(actual=as.numeric(test_data$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

rare_prob_brier_data_app <- new_metric("rare_prob_brier_data_app",
     "Rare Probability Brier Score", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          ret <- brier(actual=act_pred$actual, predicted=act_pred$predicted)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for rare_prob_brier_data_app")
          }

          return(ret)
     }
)

rare_prob_strat_brier_pos_data_app <- new_metric("rare_prob_strat_brier_pos_data_app",
     "Stratified Brier Score (Positive)", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          pos_inds <- act_pred$actual == 1

          if(sum(pos_inds) > 0){
               ret <- brier(actual=act_pred$actual[pos_inds],
                    predicted=act_pred$predicted[pos_inds])
          } else{
               ret <- 0
          }

          if(is.na(ret)){
               stop("NA detected as output for rare_prob_strat_brier_pos_data_app")
          }

          return(ret)
     }
)

rare_prob_strat_brier_neg_data_app <- new_metric("rare_prob_strat_brier_neg_data_app",
     "Stratified Brier Score (Negative)", metric = function(model, out) {

          act_pred <- get_actual_pred_rare(model, out)

          neg_inds <- act_pred$actual == 0

          if(sum(neg_inds) > 0){
               ret <- brier(actual=act_pred$actual[neg_inds],
                    predicted=act_pred$predicted[neg_inds])
          } else{
               ret <- 0
          }

          if(is.na(ret)){
               stop("NA detected as output for rare_prob_strat_brier_neg_data_app")
          }

          return(ret)
     }
)

# Kolmogorov-Smirnov test (higher is better)--PRESTO did worse than logit,
# better than proportional odds
# https://towardsdatascience.com/evaluating-classification-models-with-kolmogorov-smirnov-ks-test-e211025f5573

ks <- function(actual, predicted){
     stopifnot(is.numeric(actual) | is.integer(actual))
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(all(actual %in% c(0, 1)))
     n <- length(actual)
     stopifnot(n == length(predicted))
     stopifnot(all(!is.na(actual)))
     stopifnot(all(!is.na(predicted)))

     ks <- ks.test(x=predicted[actual == 0], y=predicted[actual == 1])$statistic

     stopifnot(!is.na(ks))
     stopifnot(is.numeric(ks) | is.integer(ks))
     stopifnot(length(ks) == 1)
     stopifnot(ks >= 0)
     stopifnot(ks <= 1)

     return(ks)
}

ks_gen <- new_metric("ks_gen", "KS Statistic", metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- ks(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

empDist <- function(y, data){
     # Calculate the empirical distribution value for y given a data set
     ret <- numeric(length(y))
     for(i in 1:length(ret)){
          ret[i] <- 1/length(data)*sum(data <= y[i])

          stopifnot(!is.na(ret[i]))
          stopifnot(is.numeric(ret[i]) | is.integer(ret[i]))
          stopifnot(length(ret[i]) == 1)
          stopifnot(ret[i] >= 0)
          stopifnot(ret[i] <= 1)
     }
     
     return(ret)
}

# Anderson-Darling Test (higher is better)--not getting finite results for any
# method, using either DescTools implementation or goftest
# library(DescTools)
ad <- function(actual, predicted){
     stopifnot(is.numeric(actual) | is.integer(actual))
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(all(actual %in% c(0, 1)))
     n <- length(actual)
     stopifnot(n == length(predicted))
     stopifnot(all(!is.na(actual)))
     stopifnot(all(!is.na(predicted)))

     func1 <- function(t){
          if(sum(actual == 1) > 0){
               return(empDist(t, predicted[actual == 0]))
          }
          return(1)
     }

     ad <- ad.test(x=predicted[actual == 1], null=func1)$statistic

     stopifnot(!is.na(ad))
     stopifnot(is.numeric(ad) | is.integer(ad))
     stopifnot(length(ad) == 1)
     stopifnot(ad >= 0)

     return(ad)
}

ad_gen <- new_metric("ad_gen", "AD Statistic", metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- ad(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Cramer-von Mises Test (higher is better)--PRESTO does better than proportional
# odds but worse than logit
# library(goftest)
cvm <- function(actual, predicted){
     stopifnot(is.numeric(actual) | is.integer(actual))
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(all(actual %in% c(0, 1)))
     n <- length(actual)
     stopifnot(n == length(predicted))
     stopifnot(all(!is.na(actual)))
     stopifnot(all(!is.na(predicted)))

     func1 <- function(t){
          if(sum(actual == 1) > 0){
               return(empDist(t, predicted[actual == 1]))
          }
          return(1)
     }

     cvm <- cvm.test(x=predicted[actual == 0], null=func1)$statistic

     stopifnot(!is.na(cvm))
     stopifnot(is.numeric(cvm) | is.integer(cvm))
     stopifnot(length(cvm) == 1)
     stopifnot(cvm >= 0)

     return(cvm)
}

cvm_gen <- new_metric("cvm_gen", "CvM Statistic", metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- cvm(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

          stopifnot(!is.na(ret))

          return(ret)
     }
)

# Rare class probability log loss (lower is better)--PRESTO did worse than logit
# by statistically significant margin, better than proportional odds

log_loss <- function(actual, predicted){
     stopifnot(is.numeric(actual) | is.integer(actual))
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(all(actual %in% c(0, 1)))
     stopifnot(length(actual) == length(predicted))
     stopifnot(all(!is.na(actual)))
     stopifnot(all(!is.na(predicted)))

     # score <- -(actual * log(predicted) + (1 - actual) * log(1 - predicted))
     # score[actual == predicted] <- 0
     # score[is.nan(score)] <- Inf
     # return(mean(score))

     loss <- 0

     for(k in 0:1){
          inds_k <- actual == k
          stopifnot(sum(inds_k) >= 0)
          if(sum(inds_k) == 0){
               next
          }
          if(k == 0){
               if(any(1 - predicted[inds_k] <= 0)){
                    return(NA)
               }
               loss <- loss + sum(-log(1 - predicted[inds_k]))
          } else{
               if(any(predicted[inds_k] <= 0)){
                    return(NA)
               }
               loss <- loss + sum(-log(predicted[inds_k]))
          }
     }

     stopifnot(!is.na(loss))
     stopifnot(is.numeric(loss) | is.integer(loss))
     stopifnot(length(loss) == 1)
     stopifnot(loss >= 0)

     ret <- loss/length(actual)

     stopifnot(!is.na(ret))
     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)

     return(ret)
}

# Multinomial log loss

multi_log_loss <- function(actual, predicted){
     stopifnot(is.numeric(actual) | is.integer(actual))
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(all(actual >= 1))
     stopifnot(all(actual == round(actual)))
     K <- max(actual)
     stopifnot(length(unique(actual)) == K)

     n <- length(actual)

     # stopifnot(all(predicted >= 0))
     # stopifnot(all(predicted <= 1))
     stopifnot(nrow(predicted) == n)
     stopifnot(ncol(predicted) == K)

     stopifnot(all(!is.na(actual)))
     stopifnot(all(!is.na(predicted)))

     stopifnot(all(abs(rowSums(predicted) - 1) < 10^(-6)))

     loss <- 0

     for(k in 1:K){
          inds_k <- actual == k
          stopifnot(sum(inds_k) > 0)
          if(any(predicted[inds_k, k] <= 0)){
               return(NA)
          }
          loss <- loss + sum(-log(predicted[inds_k, k]))
     }

     stopifnot(loss >= 0)

     ret <- loss/n

     stopifnot(is.numeric(ret) | is.integer(ret))
     stopifnot(!is.na(ret))
     stopifnot(length(ret) == 1)
     stopifnot(ret >= 0)

     return(ret)
}

log_loss_weighted <- function(actual, predicted, weights){
     stopifnot(is.numeric(actual) | is.integer(actual))
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(is.numeric(weights) | is.integer(weights))
     stopifnot(all(actual >= 0))
     stopifnot(all(actual <= 1))
     stopifnot(all(predicted >= 0))
     stopifnot(all(predicted <= 1))
     stopifnot(all(weights >= 0))
     stopifnot(length(actual) == length(predicted))
     stopifnot(length(actual) == length(weights))
     stopifnot(all(!is.na(actual)))
     stopifnot(all(!is.na(predicted)))
     stopifnot(all(!is.na(weights)))
     stopifnot(sum(weights) > 0)

     score <- -(actual * log(predicted) + (1 - actual) * log(1 - predicted))
     score[actual == predicted] <- 0
     score[is.nan(score)] <- Inf

     stopifnot(length(score) == length(weights))

     weights <- weights/sum(weights)

     return(mean(weights*score))
}

rare_prob_log_loss <- new_metric("rare_prob_log_loss",
     "Rare Probability Log Loss", metric = function(model, out) {
          alpha_hat <- out$theta_hat[1]

          beta_hat <- out$theta_hat[2:(model$p + 1)]
          stopifnot(length(beta_hat) == model$p)

          if(model$p > 1){
               stopifnot(nrow(out$X) == model$n)
               stopifnot(ncol(out$X) == model$p)
               p_hat <- 1 - plogis(alpha_hat + out$X %*% beta_hat)[, 1]
          } else{
               stopifnot(length(out$X) == model$n)
               p_hat <- 1 - plogis(alpha_hat + beta_hat*out$X)
          }
          stopifnot(length(p_hat) == model$n)

          return(log_loss(actual=as.numeric(out$y == model$K),
               predicted=p_hat))
     }
)

rare_prob_log_loss_gen <- new_metric("rare_prob_log_loss_gen",
     "Rare Probability Log Loss", metric = function(model, out) {

          p_hat <- getEvalQuantities(model, out, out$X)$p_hat

          ret <- log_loss(actual=as.numeric(out$y == as.character(model$K)),
               predicted=p_hat)

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

mspe <- function(predicted, actual){
     
     stopifnot(is.numeric(predicted) | is.integer(predicted))
     stopifnot(is.numeric(actual) | is.integer(actual))
     stopifnot(all(!is.na(predicted)))
     stopifnot(all(!is.na(actual)))
     n <- length(predicted)
     stopifnot(n >= 1)
     stopifnot(n == length(actual))

     if(any(actual == 0)){
          stop("MSPE metric provided with an actual value equal to 0; MSPE not defined")
     }

     ret <- sum(((predicted - actual)/actual)^2)/n

     stopifnot(!is.na(ret))
     stopifnot(ret >= 0)

     return(ret)
}

mae <- function(x, y){
     
     stopifnot(is.numeric(x) | is.integer(x))
     stopifnot(is.numeric(y) | is.integer(y))
     stopifnot(all(!is.na(x)))
     stopifnot(all(!is.na(y)))
     n <- length(x)
     stopifnot(n >= 1)
     stopifnot(n == length(y))

     return(sum(abs(x - y))/n)
}

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

rare_prob_mspe_gen <- new_metric("rare_prob_mspe_gen",
     "Rare Probability MSPE", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          ret <- mspe(p_hat, rare_probs)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          stopifnot(!is.na(ret))

          return(ret)
     }
)

rare_prob_mae_gen <- new_metric("rare_prob_mae_gen",
     "Rare Probability MAE", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          ret <- mae(p_hat, rare_probs)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          stopifnot(!is.na(ret))

          return(ret)
     }
)

rare_prob_mse_gen_only_rare <- new_metric("rare_prob_mse_gen_only_rare",
     "Rare Probability MSE (Only Rare)", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          rare_inds <- out$y == as.character(model$K)

          stopifnot(sum(rare_inds) < model$n)
          stopifnot(sum(rare_inds) >= 1)

          actual <- rare_probs[rare_inds]
          pred <- p_hat[rare_inds]

          stopifnot(length(actual) == length(pred))
          stopifnot(length(actual) == sum(rare_inds))

          ret <- mse(actual, pred)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          stopifnot(!is.na(ret))

          return(ret)
     }
)

rare_prob_mae_gen_only_rare <- new_metric("rare_prob_mae_gen_only_rare",
     "Rare Probability MAE (Only Rare)", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          rare_inds <- out$y == as.character(model$K)

          stopifnot(sum(rare_inds) < model$n)
          stopifnot(sum(rare_inds) >= 1)

          actual <- rare_probs[rare_inds]
          pred <- p_hat[rare_inds]

          stopifnot(length(actual) == length(pred))
          stopifnot(length(actual) == sum(rare_inds))

          ret <- mae(actual, pred)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          stopifnot(!is.na(ret))

          return(ret)
     }
)

rare_prob_mse_gen_near_db <- new_metric("rare_prob_mse_gen_near_db",
     "Rare Probability MSE (Near DB)", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          mu <- ret_list$mu
          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          rare_inds <- (mu >= -1) & (mu <= 1)

          pred <- p_hat[rare_inds]
          actual <- rare_probs[rare_inds]

          stopifnot(length(actual) == length(pred))

          if(sum(rare_inds) == 0){
               return(NA)
          } else{

               ret <- mse(pred, actual)

               # Should be possibe to get NA in this case--all predicted probabilities
               # should be nonnegative regardless of model
               stopifnot(!is.na(ret))

               return(ret)
          }
     }
)

rare_prob_mae_gen_near_db <- new_metric("rare_prob_mae_gen_near_db",
     "Rare Probability MAE (Near DB)", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          mu <- ret_list$mu
          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          rare_inds <- (mu >= -1) & (mu <= 1)

          pred <- p_hat[rare_inds]
          actual <- rare_probs[rare_inds]

          stopifnot(length(actual) == length(pred))

          if(sum(rare_inds) == 0){
               return(NA)
          } else{

               ret <- mae(pred, actual)

               # Should be possibe to get NA in this case--all predicted probabilities
               # should be nonnegative regardless of model
               stopifnot(!is.na(ret))

               return(ret)
          }
     }
)

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

rare_prob_mse_gen_high_prob <- new_metric("rare_prob_mse_gen_high_prob",
     "Rare Probability MSE (High prob units)", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)
          
          rare_inds <- rare_probs >= 1/3

          pred <- p_hat[rare_inds]
          actual <- rare_probs[rare_inds]

          stopifnot(length(actual) == length(pred))

          if(sum(rare_inds) == 0){
               return(NA)
          } else{

               ret <- mse(pred, actual)

               # Should be possibe to get NA in this case--all predicted probabilities
               # should be nonnegative regardless of model
               stopifnot(!is.na(ret))

               return(ret)
          }
     }
)

rare_prob_mae_gen_high_prob <- new_metric("rare_prob_mae_gen_high_prob",
     "Rare Probability MAE (High prob units)", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)
          
          rare_inds <- rare_probs >= 1/3

          pred <- p_hat[rare_inds]
          actual <- rare_probs[rare_inds]

          stopifnot(length(actual) == length(pred))

          if(sum(rare_inds) == 0){
               return(NA)
          } else{

               ret <- mae(pred, actual)

               # Should be possibe to get NA in this case--all predicted probabilities
               # should be nonnegative regardless of model
               stopifnot(!is.na(ret))

               return(ret)
          }
     }
)



rare_prob_log_loss_gen_only_rare <- new_metric("rare_prob_log_loss_gen_only_rare",
     "Rare Probability Log Loss (Only Rare)", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          # Only select indices with observed rare outcomes

          rare_inds <- out$y == as.character(model$K)

          actual <- as.numeric(out$y == as.character(model$K))

          actual <- actual[rare_inds]

          pred <- p_hat[rare_inds]

          stopifnot(sum(rare_inds) > 0)
          stopifnot(length(actual) == length(pred))
          stopifnot(length(unique(actual)) == 1)

          ret <- log_loss(actual=actual, predicted=pred)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          stopifnot(!is.na(ret))

          return(ret)
     }
)

rare_prob_log_loss_gen_near_db <- new_metric("rare_prob_log_loss_gen_near_db",
     "Rare Probability Log Loss (Near DB)", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          mu <- ret_list$mu
          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          rare_inds <- (mu >= -1) & (mu <= 1)

          actual <- as.numeric(out$y == as.character(model$K))
          actual <- actual[rare_inds]
          pred <- p_hat[rare_inds]

          stopifnot(length(actual) == length(pred))

          if(sum(rare_inds) == 0){
               return(NA)
          } else{
               ret <- log_loss(actual=actual, predicted=pred)

               # Should be possibe to get NA in this case--all predicted probabilities
               # should be nonnegative regardless of model
               stopifnot(!is.na(ret))

               return(ret)
          }
     }
)

rare_prob_log_loss_gen_high_prob <- new_metric("rare_prob_log_loss_gen_high_prob",
     "Rare Probability Log Loss (High prob units)", metric = function(model, out) {

          ret_list <- getEvalQuantities(model, out, out$X)

          p_hat <- ret_list$p_hat
          rare_probs <- ret_list$rare_probs

          rm(ret_list)

          rare_inds <- rare_probs >= 1/3

          actual <- as.numeric(out$y == as.character(model$K))
          actual <- actual[rare_inds]
          pred <- p_hat[rare_inds]

          stopifnot(length(actual) == length(pred))

          if(sum(rare_inds) == 0){
               return(NA)
          } else{
               ret <- log_loss(actual=actual, predicted=pred)

               # Should be possibe to get NA in this case--all predicted probabilities
               # should be nonnegative regardless of model
               stopifnot(!is.na(ret))

               return(ret)
          }
     }
)

rare_prob_log_loss_data_app <- new_metric("rare_prob_log_loss_data_app",
     "Rare Probability Log Loss", metric = function(model, out) {

          p <- ncol(out$X_ordnet)
          n <- nrow(out$X_ordnet)
          K <- ncol(out$test_ymat)

          stopifnot(n == nrow(out$test_ymat))
          stopifnot(K >= 3)

          stopifnot(p > 1)

          stopifnot(length(out$theta_hat) == p + 1)
          alpha_hat <- out$theta_hat[1]

          beta_hat <- out$theta_hat[2:(p + 1)]
          stopifnot(length(beta_hat) == p)

          actual <- rep(as.integer(NA), 2*n)
          predicted <- rep(as.numeric(NA), 2*n)
          weights <- rep(as.numeric(NA), 2*n)

          for(i in 1:n){
               ind_neg <- (i - 1)*2 + 1
               ind_pos <- 2*i

               actual[ind_neg] <- 0
               actual[ind_pos] <- 1

               predicted[ind_neg] <- 1 - plogis(alpha_hat + out$X_ordnet[i, ] %*% beta_hat)[, 1]
               predicted[ind_pos] <- predicted[ind_neg]

               weights[ind_neg] <- sum(out$test_ymat[i, 1:(K - 1)])
               weights[ind_pos] <- out$test_ymat[i, K]
          }

          return(log_loss_weighted(actual=actual, predicted=predicted,
               weights=weights))
     }
)

rare_prob_log_loss_gen_data_app <- new_metric("rare_prob_log_loss_gen_data_app",
     "Rare Probability Log Loss", metric = function(model, out) {

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

          ret <- log_loss(actual=as.numeric(out$test_y == model$resp_levels[K]),
               predicted=p_hat)

          # Should be possibe to get NA in this case--all predicted probabilities
          # should be nonnegative regardless of model
          if(is.na(ret)){
               stop("NA detected as output for rare_prob_log_loss_gen_data_app")
          }

          return(ret)
     }
)

log_lik_gen_data_app <- new_metric("log_lik_gen_data_app",
     "Overall Log Likeilhood", metric = function(model, out) {

          p <- ncol(out$X_ordnet_test)
          n <- nrow(out$X_ordnet_test)
          K <- length(unique(out$test_y))

          stopifnot(n == length(out$test_y))
          stopifnot(K >= 3)

          stopifnot(p > 1)

          if(any(!is.na(out$alpha_hat))){

               stopifnot(length(out$alpha_hat) == K - 1)
               stopifnot(all(dim(out$Beta_hat) == c(p, K - 1)))

               stopifnot(nrow(out$X_ordnet_test) == n)
               stopifnot(ncol(out$X_ordnet_test) == p)

               p_hat_mat <- matrix(as.numeric(NA), nrow=n, ncol=K)

               p_hat_mat[, 1] <- plogis(out$alpha_hat[1] +
                    out$X_ordnet_test %*% out$Beta_hat[, 1])[, 1]

               for(k in 2:(K - 1)){
                    p_hat_mat[, k] <- plogis(out$alpha_hat[k] +
                         out$X_ordnet_test %*% out$Beta_hat[, k])[, 1] -
                         plogis(out$alpha_hat[k - 1] +
                         out$X_ordnet_test %*% out$Beta_hat[, k - 1])[, 1]
               }

               p_hat_mat[, K] <- 1 - plogis(out$alpha_hat[K - 1] +
                    out$X_ordnet_test %*% out$Beta_hat[, K - 1])[, 1]

               # print("out$alpha_hat:")
               # print(out$alpha_hat)
               # print("out$Beta_hat:")
               # print(out$Beta_hat)
               # print("p_hat_mat[is.na(p_hat_mat:)]")
               # print(p_hat_mat[is.na(p_hat_mat)])
               # print("sum(is.na(p_hat_mat)):")
               # print(sum(is.na(p_hat_mat)))

               # print("head(p_hat_mat):")
               # print(head(p_hat_mat))
               # print("tail(p_hat_mat):")
               # print(tail(p_hat_mat))

               # print("n:")
               # print(n)
               # print("p:")
               # print(p)
               # print("K:")
               # print(K)

               stopifnot(all(!is.na(p_hat_mat)))
               stopifnot(all(abs(rowSums(p_hat_mat) - 1) < 10^(-6)))
          } else{
               stopifnot(!is.na(out$fit))
               stopifnot("ordinalNet" %in% class(out$fit))

               p_hat_mat <- predict(out$fit, newx=out$X_ordnet_test)

               stopifnot(ncol(p_hat_mat) == K)
               stopifnot(nrow(p_hat_mat) == n)
          }

          return(multi_log_loss(actual=as.integer(out$test_y),
               predicted=p_hat_mat))
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



log_lik_gen <- new_metric("log_lik_gen",
     "Overall Log Likeilhood", metric = function(model, out) {


          stopifnot(all(c("n", "K", "p") %in% names(model@params)))
          stopifnot(all(c("X", "y") %in% names(out)))

          if(model$p > 1){
               stopifnot(model$n == nrow(out$X))
               stopifnot(model$p == ncol(out$X))
          } else{
               stopifnot(length(out$X) == model$n)
          }
          
          stopifnot(model$K == length(unique(out$y)))

          stopifnot(model$n == length(out$y))
          stopifnot(model$K >= 3)

          if(any(!is.na(out$alpha_hat))){

               stopifnot("intercepts" %in% names(model@params))
               stopifnot(all(c("alpha_hat", "Beta_hat") %in% names(out)))

               stopifnot(length(out$alpha_hat) == model$K - 1)
               stopifnot(all(dim(out$Beta_hat) == c(model$p, model$K - 1)))
               stopifnot(all(!is.na(out$alpha_hat)))
               stopifnot(all(!is.na(out$Beta_hat)))

               p_hat_mat <- matrix(as.numeric(NA), nrow=model$n, ncol=model$K)

               p_hat_mat[, 1] <- plogis(out$alpha_hat[1] +
                    out$X %*% out$Beta_hat[, 1])[, 1]

               for(k in 2:(model$K - 1)){
                    p_hat_mat[, k] <- plogis(out$alpha_hat[k] +
                         out$X %*% out$Beta_hat[, k])[, 1] -
                         plogis(out$alpha_hat[k - 1] +
                         out$X %*% out$Beta_hat[, k - 1])[, 1]
               }

               p_hat_mat[, model$K] <- 1 - plogis(out$alpha_hat[model$K - 1] +
                    out$X %*% out$Beta_hat[, model$K - 1])[, 1]

               stopifnot(all(!is.na(p_hat_mat)))
               stopifnot(all(abs(rowSums(p_hat_mat) - 1) < 10^(-6)))
          } else{

               stopifnot("fit" %in% names(out))
               stopifnot(!is.na(out$fit))
               stopifnot("ordinalNet" %in% class(out$fit))

               p_hat_mat <- predict(out$fit, newx=out$X)

               stopifnot(ncol(p_hat_mat) == model$K)
               stopifnot(nrow(p_hat_mat) == model$n)
          }

          return(multi_log_loss(actual=as.integer(out$y),
               predicted=p_hat_mat))
     }
)

# if(K >= 4){
#      less_rare_intercept_mses[i] <- (intercepts[K - 2] - unname(model$zeta[K - 2]))^2
# }

# Estimated probabilities


