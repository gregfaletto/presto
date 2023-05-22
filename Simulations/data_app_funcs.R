# setwd("/Users/gregfaletto/Library/CloudStorage/OneDrive-USCMarshallSchoolofBusiness/Proportional odds generalization/Code/Simulations")

prepOrdNet <- function(df, final_y_names, cum_y_names){

     # Check intput

     stopifnot(is.character(cum_y_names))
     stopifnot(all(!is.na(cum_y_names)))
     stopifnot(all(cum_y_names %in% colnames(df)))

     K <- length(cum_y_names)

     stopifnot(K >= 3)

     stopifnot(is.character(final_y_names))
     stopifnot(all(!is.na(final_y_names)))
     stopifnot(length(final_y_names) == K)

     stopifnot(is.data.frame(df))

     n_obs <- nrow(df)
     
     response_mat <- df[, cum_y_names]

     y_mat <- matrix(as.numeric(NA), nrow=n_obs, ncol=K)

     colnames(y_mat) <- final_y_names

     y_mat[, K] <- response_mat[, K]

     for(k in (K - 1):1){
          y_mat[, k] <- response_mat[, k] - response_mat[, k + 1]
     }

     # Prepare ordinalNet predictor matrix

     # # sim <- new_simulation("data_app_basic", "Data application")
     # X_ordnet <- model.matrix(~ xyz_campaign_id + age + gender + interest,
     #      data=df)

     # # sim <- new_simulation("data_app_basic_ints", "Data application")
     # X_ordnet <- model.matrix(~ xyz_campaign_id + age + gender +
     #      xyz_campaign_id*age + xyz_campaign_id*gender, data=df)

     # sim <- new_simulation("data_app", "Data application")
     # X_ordnet <- model.matrix(~ xyz_campaign_id + age + gender, data=df)

     # sim <- new_simulation("data_app_all_ints", "Data application")
     X_ordnet <- model.matrix(~ xyz_campaign_id + age + gender +
          interest + xyz_campaign_id*age + xyz_campaign_id*gender +
          xyz_campaign_id*interest, data=df)

     X_ordnet <- X_ordnet[, colnames(X_ordnet) != "(Intercept)"]

     # Check output

     stopifnot(nrow(X_ordnet) == n_obs)
     stopifnot(all(dim(y_mat) == c(n_obs, K)))
     stopifnot(all(!is.na(y_mat)))
     stopifnot(all(y_mat >= 0))

     return(list(X_ordnet=X_ordnet, y_mat=y_mat))
}

make_ypolr <- function(y_mat){

	# Columns of y_mat must be labeled with categories and must be in
	# increasing order.
	# Check input

	stopifnot(all(!is.na(y_mat)))
	stopifnot(all(y_mat >= 0))

	K <- ncol(y_mat)
	n_obs <- nrow(y_mat)

	stopifnot(!is.null(colnames(y_mat)))

	var_names <- colnames(y_mat)

	stopifnot(length(var_names) == K)
	stopifnot(is.character(var_names))

	# Prepare MASS::polr responses (vector of factors, as well as vector of 
	# weights)

	y_polr <- rep(as.character(NA), n_obs*K)
	weights_polr <- rep(as.numeric(NA), n_obs*K)

	for(k in 1:K){
	  inds <- ((1:n_obs) - 1)*K + k
	  stopifnot(all(is.na(y_polr[inds])))
	  stopifnot(all(is.na(weights_polr[inds])))

	  y_polr[inds] <- var_names[k]
	  weights_polr[inds] <- y_mat[, k]
	}


	# Check output

	stopifnot(all(!is.na(y_polr)))
	stopifnot(all(var_names %in% y_polr))
	stopifnot(length(unique(y_polr)) == K)
	stopifnot(length(y_polr) == n_obs*K)

	stopifnot(sum(weights_polr) == sum(y_mat))
	stopifnot(all(!is.na(weights_polr)))
	stopifnot(all(weights_polr >= 0))
	stopifnot(length(weights_polr) == length(y_polr))

	y_polr <- factor(y_polr, levels=var_names, ordered=TRUE)

	stopifnot(length(weights_polr) == length(y_polr))

	return(list(y_polr=y_polr, weights_polr=weights_polr))

}

pre_df_polr <- function(y_mat, X_ordnet){

	# Check inputs

	stopifnot(all(!is.na(y_mat)))
	stopifnot(all(y_mat >= 0))

	n_obs <- nrow(y_mat)
	K <- ncol(y_mat)

	stopifnot(is.matrix(X_ordnet))
	stopifnot(nrow(X_ordnet) == n_obs)

	p <- ncol(X_ordnet)

	ret_list <- make_ypolr(y_mat)

	y_polr <- ret_list$y_polr
	weights_polr <- ret_list$weights_polr

	stopifnot(length(weights_polr) == length(y_polr))

	rm(ret_list)

	# # Prepare logistic regression responses
	# y_logit <- as.factor(y_polr == my_responses[K])
	# levels(y_logit) <- c("NoPurchase", "Purchase")

	# Prepare MASS::polr and glm predictor matrix

	X_polr <- X_ordnet[rep(1:n_obs , each=K), ]
	stopifnot(nrow(X_polr) == n_obs*K)

	X_polr <- X_polr[weights_polr > 0, ]
	stopifnot(nrow(X_polr) <= length(weights_polr))
	stopifnot(nrow(X_polr) == sum(weights_polr > 0))

	y_polr <- y_polr[weights_polr > 0]
	stopifnot(nrow(X_polr) == length(y_polr))

	# y_logit <- y_logit[weights_polr > 0]
	# stopifnot(length(y_logit) == length(y_polr))

	weights_polr <- weights_polr[weights_polr > 0]
	stopifnot(length(weights_polr) == nrow(X_polr))

	stopifnot(sum(weights_polr) == sum(y_mat))

	df_polr <- data.frame(X_polr, y_polr)
	# df_logit <- data.frame(X_polr, y_logit)

	return(list(df_polr=df_polr, weights_polr=weights_polr))
}

#### Functions for subsampling

subsampCountVec <- function(count_vec, sample_size, names=NA){
     # Takes a subsample from a vector of counts. If count_vec has names,
     # then returned fector will have same names. If names are provided,
     # any names of count_vec will be overwritten internally in the same order
     # as names, and the returned vector will also have these names.

     stopifnot(!is.matrix(count_vec))

     remove_names_later <- FALSE

     if(any(!is.na(names))){
          stopifnot(is.character(names))
          stopifnot(length(names) == length(count_vec))
          names(count_vec) <- names
     }

     # count_vec must be named
     if(is.null(names(count_vec))){
          remove_names_later <- TRUE
          names(count_vec) <- paste("x", 1:length(count_vec), sep="")
     }

     stopifnot(length(unique(names(count_vec))) == length(count_vec))
     stopifnot(all(count_vec >= 0))
     stopifnot(all(count_vec == round(count_vec)))

     stopifnot(is.integer(sample_size) | is.numeric(sample_size))
     stopifnot(length(sample_size) == 1)
     stopifnot(sample_size > 0)
     stopifnot(sample_size <= sum(count_vec))
     stopifnot(sample_size == round(sample_size))

     # first, draw from the possible indexes (does not create the full vector)
     draw <- sample.int(sum(count_vec), sample_size)

     # then assign indices back to original group
     items <- findInterval(draw - 1, c(0, cumsum(count_vec)),
          rightmost.closed=TRUE)

     # extract sample    
     obs <- names(count_vec)[items]
     ret_vec <- as.matrix(table(obs))[, 1]
     stopifnot(all(names(ret_vec) %in% names(count_vec)))
     stopifnot(length(unique(names(count_vec))) >= length(unique(names(ret_vec))))
     if(length(unique(names(count_vec))) > length(unique(names(ret_vec)))){
          missing_names <-setdiff(unique(names(count_vec)), unique(names(ret_vec)))
          n_missing <- length(missing_names)
          ret_vec <- c(ret_vec, rep(0, n_missing))
          names(ret_vec)[(length(ret_vec) - n_missing + 1):
               length(ret_vec)] <- missing_names
     }

     stopifnot(all(names(count_vec) %in% names(ret_vec)))
     stopifnot(length(unique(names(count_vec))) == length(unique(names(ret_vec))))
     
     ret_vec <- ret_vec[order(factor(names(ret_vec), levels=names(count_vec)))]

     # Check output
     stopifnot(all(ret_vec <= count_vec))
     stopifnot(all(ret_vec >= 0))
     stopifnot(length(ret_vec) == length(count_vec))
     stopifnot(all(names(ret_vec) == names(count_vec)))
     stopifnot(sum(ret_vec) == sample_size)

     if(any(!is.na(names))){
          stopifnot(all(names(ret_vec) == names))
     }

     if(remove_names_later){
          ret_vec <- unname(ret_vec)
     }

     return(ret_vec)
     }

# marbleCounts <- c(red = 5, green = 3, blue = 2)

# marbleSample <- subsampCountVec(marbleCounts, sample_size = 4)

# marbleCounts
# marbleSample

# # Unnamed example

# marbleCounts2 <- c(5, 3, 2)

# names <- c("purple", "orange", "yellow")

# marbleSample2 <- subsampCountVec(marbleCounts2, sample_size=6)

# marbleCounts2
# marbleSample2

# marbleSample3 <- subsampCountVec(marbleCounts2, sample_size=6, names=names)

# marbleCounts2
# marbleSample3

subsampCountMat <- function(count_mat, sample_size){
     # Takes a subsample from a matrix of counts. Columns of matrix must be 
     # named.

     stopifnot(is.matrix(count_mat))

     n <- nrow(count_mat)
     p <- ncol(count_mat)

     stopifnot(n > 1)
     stopifnot(p > 1)
     stopifnot(!is.null(colnames(count_mat)))
     stopifnot(length(unique(colnames(count_mat))) == p)
     if(any(colSums(count_mat) == 0)){
          stop("Provided count_mat must have at least one nonzero entry in every column")
     }

     vec <- as.vector(count_mat)

     stopifnot(length(vec) == n*p)

     if(is.null(rownames(count_mat))){
          # Make up names to correspond to rows
          rownames(count_mat) <- paste("x", 1:n, sep="")
     }

     stopifnot(length(unique(rownames(count_mat))) == n)

     stopifnot(is.integer(sample_size) | is.numeric(sample_size))
     stopifnot(length(sample_size) == 1)
     stopifnot(sample_size > 0)
     stopifnot(sample_size <= sum(count_mat))
     stopifnot(sample_size == round(sample_size))

     # Name vec with names according to column name and row names

     names(vec) <- rep(as.character(NA), n*p)

     for(j in 1:p){
          inds <- (j - 1)*n + 1:n
          stopifnot(all(is.na(names(vec)[inds])))
          names(vec)[inds] <- paste(rownames(count_mat),
               colnames(count_mat)[j], sep="_")
     }

     stopifnot(all(!is.na(names(vec))))

     # Subsample these indices

     # Want to ensure that every column has at least one observation
     all_cols_flag <- FALSE
     iter_count <- 0

     while(!all_cols_flag){

          if(iter_count >= 10){
               stop("Unable to create subsampled matrix with at least one observation in each column after 10 attempts; breaking loop here.")
          }

          subsamp_vec <- subsampCountVec(vec, sample_size=sample_size)

          # Return to matrix

          subsamp_mat <- matrix(subsamp_vec, nrow=n, ncol=p)
          colnames(subsamp_mat) <- colnames(count_mat)

          # Check output
          stopifnot(all(subsamp_mat <= count_mat))
          stopifnot(all(subsamp_mat >= 0))
          stopifnot(all(dim(subsamp_mat) == dim(count_mat)))
          stopifnot(all(colnames(subsamp_mat) == colnames(count_mat)))
          stopifnot(sum(subsamp_mat) == sample_size)

          all_cols_flag <- all(colSums(subsamp_mat) > 0)
          iter_count <- iter_count + 1
     }

     return(subsamp_mat)
}

conf_int <- function(est, se, alpha=0.05, round=5){
     stopifnot(length(est) == 1)
     stopifnot(length(se) == 1)
     stopifnot(length(alpha) == 1)
     stopifnot(length(round) == 1)
     stopifnot(se >= 0)
     stopifnot(round >= 0)
     stopifnot(alpha > 0)
     stopifnot(alpha < 1)
     stopifnot(round == round(round))
     lb <- est - qnorm(1 - alpha/2)*se
     ub <- est + qnorm(1 - alpha/2)*se
     paste("(", round(lb, round), ", ", round(ub, round), ")", sep="")
}

# n <- 8
# p <- 3
# stopifnot(p <= length(letters))
# sample_mat <- matrix(sample(1:(n*p), size=n*p), nrow=n, ncol=p)
# colnames(sample_mat) <- letters[1:p]

# samp_size <- n*p*2

# subsamp_mat <- subsampCountMat(sample_mat, samp_size)


# sample_mat

# subsamp_mat










