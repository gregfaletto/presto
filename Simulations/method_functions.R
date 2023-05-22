# qs used for this iteration

### Proportional Odds

prop_odds_meth <- new_method("prop_odds_meth", "PO",
	method = function(model, draw) {

	# theta_hat <- rep(as.numeric(NA), model$p + 1)
	# # Proportional odds model
	# df <- data.frame(y=draw$y, X=draw$X)
	# polr_model <- MASS::polr(y ~., df, method="logistic")

	# # Save estimates
	# theta_hat[1] <- unname(polr_model$zeta[model$K - 1])
	# # Correct for sign flip in MASS:polr model
	# theta_hat[2:(model$p + 1)] <- -unname(coef(polr_model))

	# stopifnot(all(!is.na(theta_hat)))
	# names(theta_hat) <- c("alpha", paste("beta", 1:model$p,
	# 	sep="_"))

	# return(list(theta_hat=theta_hat, X=draw$X, y=draw$y, probs=draw$probs,
	# 	beta_final=draw$beta_final))

		alpha_hat <- rep(as.numeric(NA), model$K - 1)
		Beta_hat <- matrix(rep(as.numeric(NA), model$p*(model$K - 1)),
			nrow=model$p, ncol=model$K - 1)

		# Proportional odds model

		df <- data.frame(y=draw$y, X=draw$X)
		prop_odds_model <- MASS::polr(y ~., df, method="logistic")

		# Save estimates
		alpha_hat <- prop_odds_model$zeta[1:(model$K - 1)]
		names(alpha_hat) <- names(prop_odds_model$zeta[1:(model$K - 1)])
		# Correct for sign flip in MASS:polr model
		for(k in 1:(model$K - 1)){
			Beta_hat[, k] <- -coef(prop_odds_model)
		}

		stopifnot(length(coef(prop_odds_model)) == model$p)

		rownames(Beta_hat) <- names(coef(prop_odds_model))
		colnames(Beta_hat) <- names(prop_odds_model$zeta[1:(model$K - 1)])

		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))
		stopifnot(length(alpha_hat) == model$K - 1)
		stopifnot(all(dim(Beta_hat) == c(model$p, model$K - 1)))

		return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat, X=draw$X, y=draw$y,
			probs=draw$probs, beta_final=draw$beta_final, fit=NA))

	}#,
	# settings = list(cutoff = 0.7, q=10)
)

# True probabilities (reference)

true_probs_meth <- new_method("true_probs_meth", "GT",
	method = function(model, draw) {

		mu <- model$intercepts[model$K - 1] + draw$X %*%
			draw$beta_final[, model$K - 1]

		stopifnot(length(mu) == model$n)

		rare_probs <- 1 - plogis(mu)

		stopifnot(all(!is.na(rare_probs)))
		stopifnot(all(rare_probs >= 0))
		stopifnot(all(rare_probs <= 1))
		stopifnot(length(rare_probs) == model$n)

		return(list(alpha_hat=model$intercepts, Beta_hat=draw$beta_final,
			X=draw$X, y=draw$y, probs=draw$probs, beta_final=draw$beta_final,
			fit=NA))
	}
)

### Proportional Odds for data application (matrix response with weights)

prop_odds_data_analysis <- new_method("prop_odds_data_analysis",
	"PO", method = function(model, draw) {

		K <- length(model$resp_names)

		p <- ncol(draw$train_df_polr) - 1

		theta_hat <- rep(as.numeric(NA), p + 1)

		# Proportional odds model

		prop_odds_model <- MASS::polr(y_polr ~., data=draw$train_df_polr,
			weights=draw$train_weights_polr)

		# Save estimates
		theta_hat[1] <- unname(prop_odds_model$zeta[K - 1])
		# Correct for sign flip in MASS:polr model
		theta_hat[2:(p + 1)] <- -unname(coef(prop_odds_model))

		stopifnot(all(!is.na(theta_hat)))
		names(theta_hat) <- c("alpha", paste("beta", 1:p, sep="_"))

		return(list(theta_hat=theta_hat, X_ordnet=draw$X_ordnet,
			test_weights_polr=draw$weights_polr_test,
			test_ymat=draw$test_ymat))
	}#,
	# settings = list(cutoff = 0.7, q=10)
)

### Proportional Odds for data application (typical vector response)

prop_odds_data_analysis_vec <- new_method("prop_odds_data_analysis_vec",
	"PO", method = function(model, draw) {

		K <- length(model$resp_levels)

		p <- ncol(draw$X_ordnet_train)

		alpha_hat <- rep(as.numeric(NA), K - 1)
		Beta_hat <- matrix(rep(as.numeric(NA), p*(K - 1)), nrow=p, ncol=K - 1)

		# Proportional odds model

		formula <- as.formula(paste(model$resp_name, " ~ .", sep=""))

		prop_odds_model <- MASS::polr(formula, data=draw$train_df_polr)

		# print("length(coef(prop_odds_model)):")
		# print(length(coef(prop_odds_model)))
		# print("p:")
		# print(p)
		# print("colnames(draw$X_ordnet_train:")
		# print(colnames(draw$X_ordnet_train))
		# print("names(coef(prop_odds_model)):")
		# print(names(coef(prop_odds_model)))
		stopifnot(length(coef(prop_odds_model)) == p)

		# Save estimates
		alpha_hat <- prop_odds_model$zeta[1:(K - 1)]
		names(alpha_hat) <- names(prop_odds_model$zeta[1:(K - 1)])
		# Correct for sign flip in MASS:polr model
		for(k in 1:(K - 1)){
			Beta_hat[, k] <- -coef(prop_odds_model)
		}

		stopifnot(length(coef(prop_odds_model)) == p)
		rownames(Beta_hat) <- names(coef(prop_odds_model))
		colnames(Beta_hat) <- names(prop_odds_model$zeta[1:(K - 1)])

		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))

		return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat,
			X_ordnet_test=draw$X_ordnet_test, test_y=draw$test_y, fit=NA))
	}#,
	# settings = list(cutoff = 0.7, q=10)
)

### Logistic regression on rare class

logit_meth <- new_method("logit_meth", "Logit",
	method = function(model, draw) {

	# theta_hat <- rep(as.numeric(NA), model$p + 1)
	# # Logistic regression model
	# df <- data.frame(y=draw$y, X=draw$X)
	# logit_model <- glm(y ~., family=binomial, data=df)

	# # Save estimates
	# stopifnot(length(coef(logit_model)) == model$p + 1)
	# theta_hat[1] <- unname(coef(logit_model)["(Intercept)"])
	# # Correct for sign flip in MASS:polr model
	# theta_hat[2:(model$p + 1)] <- unname(coef(logit_model)[2:(model$p + 1)])

	# stopifnot(all(!is.na(theta_hat)))
	# names(theta_hat) <- c("alpha", paste("beta", 1:model$p,
	# 	sep="_"))

	# return(list(theta_hat=theta_hat, X=draw$X, y=draw$y, probs=draw$probs,
	# 	beta_final=draw$beta_final))

		alpha_hat <- rep(as.numeric(NA), model$K - 1)
		Beta_hat <- matrix(rep(as.numeric(NA), model$p*(model$K - 1)),
			nrow=model$p, ncol=model$K - 1)

		for(k in 1:(model$K - 1)){
			# Logistic regression model
			df_k <- data.frame(X=draw$X)
			df_k$y <- as.integer(as.integer(draw$y) >= k + 1)

			stopifnot(all(df_k$y %in% c(0, 1)))
			stopifnot(sum(df_k$y == 0) > 0)
			stopifnot(sum(df_k$y == 1) > 0)

			df_k$y <- as.factor(df_k$y)

			logit_model <- glm(y ~., family=binomial, data=df_k)

			# Save estimates (need to flip signs to match convention)

			stopifnot(length(coef(logit_model)) == model$p + 1)
			stopifnot(all(!is.na(coef(logit_model))))

			alpha_hat[k] <- -unname(coef(logit_model)["(Intercept)"])
			# Correct for sign flip
			Beta_hat[, k] <- -unname(coef(logit_model)[2:(model$p + 1)])

			rownames(Beta_hat) <- names(coef(logit_model)[2:(model$p + 1)])
	
		}

		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))
		stopifnot(length(alpha_hat) == model$K - 1)
		stopifnot(all(dim(Beta_hat) == c(model$p, model$K - 1)))

		return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat, X=draw$X, y=draw$y,
			probs=draw$probs, beta_final=draw$beta_final, fit=NA))
	}#,
	# settings = list(cutoff = 0.7, q=10)
)

### Our method!!!

fused_polr <- new_method("fused_polr", "PRESTO",
	method = function(model, draw) {

		stopifnot(model$K >= 3)

		# theta_hat <- rep(as.numeric(NA), model$p + 1)

		# # Use part of data to get lambda

		# tune_inds <- 1:round(model$tune_prop*model$n)
		# train_inds <- setdiff(1:model$n, tune_inds)

		# stopifnot(length(tune_inds) + length(train_inds) == model$n)

		# lambda_star <- getLambdaFusedPolr(X=draw$X[tune_inds, ],
		# 	y=draw$y[tune_inds], printProgress=FALSE)

		# # Fit model on rest of data using selected lambda and get parameter estimates

		# # Fit model on rest of data using selected lambda

		# ret_list <- getParamEstsMod(X=draw$X[train_inds, ],
		# 	y=draw$y[train_inds], lambdaVals=lambda_star, printIter=FALSE,
		# 	K=model$K)

		# Use data to get lambda

		lambda_star <- getLambdaFusedPolr(X=draw$X, y=draw$y,
			printProgress=FALSE)

		# Fit model on rest of data using selected lambda and get parameter estimates

		# Fit model on rest of data using selected lambda

		ret_list <- getParamEstsMod(X=draw$X, y=draw$y, lambdaVals=lambda_star,
			printIter=FALSE, K=model$K)

		alpha_hat <- ret_list$alpha_hat

		Beta_hat <- ret_list$Beta_hat

		rm(ret_list)

		stopifnot(length(alpha_hat) == model$K - 1)
		stopifnot(nrow(Beta_hat) == model$p)
		stopifnot(ncol(Beta_hat) == model$K - 1)

		# theta_hat[1] <- alpha_hat[model$K - 1]

		# theta_hat[2:(model$p + 1)] <- Beta_hat[, model$K - 1]

		# stopifnot(all(!is.na(theta_hat)))
		
		# names(theta_hat) <- c("alpha", paste("beta", 1:model$p,
		# 	sep="_"))

		# return(list(theta_hat=theta_hat, X=draw$X, y=draw$y, probs=draw$probs,
		# 	beta_final=draw$beta_final))

		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))
		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))

		return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat, X=draw$X, y=draw$y,
			probs=draw$probs, beta_final=draw$beta_final, fit=NA))
	}
)


# L2 version

fused_polr_l2 <- new_method("fused_polr_l2", "PRESTO_L2",
	method = function(model, draw) {

		stopifnot(model$K >= 3)

		# theta_hat <- rep(as.numeric(NA), model$p + 1)

		# # Use part of data to get lambda

		# tune_inds <- 1:round(model$tune_prop*model$n)
		# train_inds <- setdiff(1:model$n, tune_inds)

		# stopifnot(length(tune_inds) + length(train_inds) == model$n)

		# lambda_star <- getLambdaFusedPolr(X=draw$X[tune_inds, ],
		# 	y=draw$y[tune_inds], printProgress=FALSE)

		# # Fit model on rest of data using selected lambda and get parameter estimates

		# # Fit model on rest of data using selected lambda

		# ret_list <- getParamEstsMod(X=draw$X[train_inds, ],
		# 	y=draw$y[train_inds], lambdaVals=lambda_star, printIter=FALSE,
		# 	K=model$K)

		# Use data to get lambda

		lambda_star <- getLambdaFusedPolr(X=draw$X, y=draw$y,
			printProgress=FALSE, alpha_param=0)

		# Fit model on rest of data using selected lambda and get parameter estimates

		# Fit model on rest of data using selected lambda

		ret_list <- getParamEstsMod(X=draw$X, y=draw$y, lambdaVals=lambda_star,
			printIter=FALSE, K=model$K, alpha_param=0)

		alpha_hat <- ret_list$alpha_hat

		Beta_hat <- ret_list$Beta_hat

		rm(ret_list)

		stopifnot(length(alpha_hat) == model$K - 1)
		stopifnot(nrow(Beta_hat) == model$p)
		stopifnot(ncol(Beta_hat) == model$K - 1)

		# theta_hat[1] <- alpha_hat[model$K - 1]

		# theta_hat[2:(model$p + 1)] <- Beta_hat[, model$K - 1]

		# stopifnot(all(!is.na(theta_hat)))
		
		# names(theta_hat) <- c("alpha", paste("beta", 1:model$p,
		# 	sep="_"))

		# return(list(theta_hat=theta_hat, X=draw$X, y=draw$y, probs=draw$probs,
		# 	beta_final=draw$beta_final))

		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))
		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))

		return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat, X=draw$X, y=draw$y,
			probs=draw$probs, beta_final=draw$beta_final, fit=NA))
	}
)




### Our method for data application (response matrix with weights)

fused_polr_data_analysis <- new_method("fused_polr_data_analysis", "PRESTO",
	method = function(model, draw) {

		K <- length(model$resp_names)
		stopifnot(K >= 3)

		p <- ncol(draw$X_ordnet)

		theta_hat <- rep(as.numeric(NA), p + 1)

		# Use part of data to get lambda
		# Want to ensure that every column has at least one observation (the
	    # function subsampCountMat ensures this already for y_mat_tune, so just
	    # need to confirm for y_mat_train)
	    all_cols_flag <- FALSE
	    iter_count <- 0

	    n_tune <- round(model$tune_prop*sum(draw$train_ymat))

	    while(!all_cols_flag){
	          
	        if(iter_count >= 10){
	            stop("Unable to create subsampled matrix with at least one observation in each column after 10 attempts; breaking loop here.")
	        }

	        y_mat_tune <- subsampCountMat(count_mat=draw$train_ymat,
	        	sample_size=n_tune)

	        y_mat_train <- draw$train_ymat - y_mat_tune

	        stopifnot(all(y_mat_train <= draw$train_ymat))
	        stopifnot(all(y_mat_train >= 0))
	        stopifnot(all(dim(y_mat_train) == dim(draw$train_ymat)))
	        stopifnot(all(colnames(y_mat_train) == colnames(draw$train_ymat)))
	        stopifnot(sum(y_mat_train) == sum(draw$train_ymat) - sum(y_mat_tune))

	        all_cols_flag <- all(colSums(y_mat_train) > 0)
	        iter_count <- iter_count + 1
	    }

	    lambda_star <- getLambdaFusedPolr(X=draw$X_ordnet,
			y=y_mat_tune, printProgress=FALSE)

		# fit_tune <- ordinalNetTuneMod(draw$X_ordnet, y_mat_tune,
		# 	family="cumulative", link="logit", printProgress=FALSE)

		# briers <- rowMeans(fit_tune$brier)

		# lambda_star <- max(fit_tune$fit$lambdaVals[briers == min(briers)])

		# Fit model on rest of data using selected lambda

		ret_list <- getParamEstsMod(X=draw$X_ordnet,
			y=y_mat_train, lambdaVals=lambda_star, printIter=FALSE,
			K=K)

		alpha_hat <- ret_list$alpha_hat

		Beta_hat <- ret_list$Beta_hat

		rm(ret_list)

		stopifnot(length(alpha_hat) == K - 1)
		stopifnot(nrow(Beta_hat) == p)
		stopifnot(ncol(Beta_hat) == K - 1)

		theta_hat[1] <- alpha_hat[K - 1]

		theta_hat[2:(model$p + 1)] <- Beta_hat[, K - 1]

		stopifnot(all(!is.na(theta_hat)))
		
		names(theta_hat) <- c("alpha", paste("beta", 1:model$p,
			sep="_"))

		# fit <- ordinalNetMod(draw$X_ordnet, y_mat_train,
		# 	family="cumulative", link="logit", printIter=FALSE,
		# 	lambdaVals=lambda_star)

		# # Get estimates

		# theta_hat[1] <- fit$coefs[1, K - 1]

		# psi <- fit$coefs[1, K:length(fit$coefs)]

		# stopifnot(length(psi) == p*(K - 1))

		# beta <- psi[1:p]

		# for(k in 2:(K - 1)){
		#      beta <- beta + psi[((k-1)*p + 1):(k*p)]
		# }

		# theta_hat[2:(p + 1)] <- beta

		# stopifnot(all(!is.na(theta_hat)))
		
		# names(theta_hat) <- c("alpha", paste("beta", 1:p,
		# 	sep="_"))

		return(list(theta_hat=theta_hat, X_ordnet=draw$X_ordnet,
			test_weights_polr=draw$weights_polr_test,
			test_ymat=draw$test_ymat))
	}
)

getLambdaFusedPolr <- function(X, y, printProgress=FALSE, nLambda=20,
	lambdaMinRatio=.01, alpha_param=1){
	# Check input

	stopifnot(nrow(X) == length(y))

	fit_tune <- ordinalNetTuneMod(X, y, family="cumulative", link="logit",
		printProgress=printProgress, nLambda=nLambda,
		lambdaMinRatio=lambdaMinRatio, alpha=alpha_param)

	briers <- rowMeans(fit_tune$brier)

	lambda_star <- max(fit_tune$fit$lambdaVals[briers == min(briers)])

	# Check output

	stopifnot(is.numeric(lambda_star) | is.integer(lambda_star))
	stopifnot(!is.na(lambda_star))
	stopifnot(length(lambda_star) == 1)
	stopifnot(lambda_star >= 0)

	return(lambda_star)
}

getLambdaOrdNet <- function(X, y, printProgress=FALSE, nLambda=20,
	lambdaMinRatio=.01){
	# Check input

	stopifnot(nrow(X) == length(y))

	fit_tune <- ordinalNetTune(X, y, family="cumulative", link="logit",
		printProgress=printProgress, nLambda=nLambda,
		lambdaMinRatio=lambdaMinRatio)

	briers <- rowMeans(fit_tune$brier)

	lambda_star <- max(fit_tune$fit$lambdaVals[briers == min(briers)])

	# Check output

	stopifnot(is.numeric(lambda_star) | is.integer(lambda_star))
	stopifnot(!is.na(lambda_star))
	stopifnot(length(lambda_star) == 1)
	stopifnot(lambda_star >= 0)

	return(lambda_star)
}

getParamEstsMod <- function(X, y, lambdaVals, K, printIter=FALSE, alpha_param=1){

	fit <- ordinalNetMod(X, y, family="cumulative", link="logit",
		printIter=printIter, lambdaVals=lambdaVals, alpha=alpha_param)

	p <- ncol(X)

	# Get estimates

	alpha_hat <- fit$coefs[1, 1:(K - 1)]
	names(alpha_hat) <- names(fit$coefs[1, 1:(K - 1)])

	psi <- fit$coefs[1, K:length(fit$coefs)]

	stopifnot(length(psi) == p*(K - 1))

	Beta_hat <- matrix(as.numeric(NA), nrow=p, ncol=K - 1)

	Beta_hat[, 1] <- psi[1:p]

	for(k in 2:(K - 1)){
	     Beta_hat[, k] <- Beta_hat[, k - 1] + psi[((k-1)*p + 1):(k*p)]
	}

	stopifnot(nrow(Beta_hat) == ncol(X))
	rownames(Beta_hat) <- colnames(X)
	colnames(Beta_hat) <- names(fit$coefs[1, 1:(K - 1)])

	stopifnot(all(!is.na(alpha_hat)))
	stopifnot(all(!is.na(Beta_hat)))


	return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat))
}

getParamEstsOrdNet <- function(X, y, lambdaVals, K, printIter=FALSE){

	fit <- ordinalNet::ordinalNet(X, y,
		family="cumulative", parallelTerms=TRUE, nonparallelTerms=TRUE,
		link="logit", printIter=printIter, lambdaVals=lambdaVals)

	p <- ncol(X)

	# Get estimates

	alpha_hat <- fit$coefs[1, 1:(K - 1)]
	names(alpha_hat) <- names(fit$coefs[1, 1:(K - 1)])

	psi <- fit$coefs[1, K:length(fit$coefs)]

	stopifnot(length(psi) == p*(K - 1))

	Beta_hat <- matrix(as.numeric(NA), nrow=p, ncol=K - 1)

	Beta_hat[, 1] <- psi[1:p]

	for(k in 2:(K - 1)){
	     Beta_hat[, k] <- Beta_hat[, k - 1] + psi[((k-1)*p + 1):(k*p)]
	}

	stopifnot(nrow(Beta_hat) == ncol(X))
	rownames(Beta_hat) <- colnames(X)
	colnames(Beta_hat) <- names(fit$coefs[1, 1:(K - 1)])

	stopifnot(all(!is.na(alpha_hat)))
	stopifnot(all(!is.na(Beta_hat)))


	return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat))
}

### Our method for data application (typical vector response)

fused_polr_data_analysis_vec <- new_method("fused_polr_data_analysis_vec",
	"PRESTO", method = function(model, draw) {

		K <- length(model$resp_levels)
		stopifnot(K >= 3)

		p <- ncol(draw$X_ordnet_train)
		n_obs <- nrow(draw$X_ordnet_train)

		alpha_hat <- rep(as.numeric(NA), K - 1)
		Beta_hat <- matrix(rep(as.numeric(NA), p*(K - 1)), nrow=p, ncol=K - 1)


		if(model$tune_prop > 0){
			# Use part of data to get lambda
		
			ret_list <- splitIntoTwoSetsVec(X=draw$X_ordnet_train,
				y=draw$train_y, prop=model$tune_prop)

			tune_X_ordnet <- ret_list$X_prop

			train_X_ordnet <- ret_list$X_comp

			y_tune <- ret_list$y_prop

			y_train <- ret_list$y_comp

			rm(ret_list)

		    lambda_star <- getLambdaFusedPolr(X=tune_X_ordnet, y=y_tune,
		    	printProgress=FALSE, nLambda=100, lambdaMinRatio=.0001)

			# Fit model on rest of data using selected lambda

			ret_list <- getParamEstsMod(X=train_X_ordnet, y=y_train,
				lambdaVals=lambda_star, printIter=FALSE, K=K)

			alpha_hat <- ret_list$alpha_hat

			Beta_hat <- ret_list$Beta_hat

			rm(ret_list)
		} else{

			# Use all of data to get lambda

		    lambda_star <- getLambdaFusedPolr(X=draw$X_ordnet_train,
		    	y=draw$train_y, printProgress=FALSE, nLambda=100,
		    	lambdaMinRatio=.0001)

			# Fit model on rest of data using selected lambda

			ret_list <- getParamEstsMod(X=draw$X_ordnet_train, y=draw$train_y,
				lambdaVals=lambda_star, printIter=FALSE, K=K)

			alpha_hat <- ret_list$alpha_hat

			Beta_hat <- ret_list$Beta_hat

			rm(ret_list)
		}

		return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat,
			X_ordnet_test=draw$X_ordnet_test, test_y=draw$test_y, fit=NA))
	}
)

### Logistic regression on rare class for data application (typical vector
# response)

logit_meth_gen <- new_method("logit_meth_gen", "Logit",
	method = function(model, draw) {

		K <- length(model$resp_levels)

		p <- ncol(draw$X_ordnet_train)

		alpha_hat <- rep(as.numeric(NA), K - 1)
		Beta_hat <- matrix(rep(as.numeric(NA), p*(K - 1)), nrow=p, ncol=K - 1)

		for(k in 1:(K - 1)){
			# Logistic regression model
			df_k <- draw$train_df_polr
			df_k$y <- as.integer(as.integer(df_k[, model$resp_name]) >= k + 1)

			stopifnot(all(df_k$y %in% c(0, 1)))
			stopifnot(sum(df_k$y == 0) > 0)
			stopifnot(sum(df_k$y == 1) > 0)

			df_k <- df_k[, !(colnames(df_k) == model$resp_name)]
			df_k$y <- as.factor(df_k$y)

			logit_model <- glm(y ~., family=binomial, data=df_k)

			# Save estimates (need to flip signs to match convention)

			stopifnot(length(coef(logit_model)) == p + 1)
			stopifnot(all(!is.na(coef(logit_model))))

			alpha_hat[k] <- -unname(coef(logit_model)["(Intercept)"])
			# Correct for sign flip
			Beta_hat[, k] <- -unname(coef(logit_model)[2:(p + 1)])

			rownames(Beta_hat) <- names(coef(logit_model)[2:(p + 1)])
	
		}

		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))

		return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat,
			X_ordnet_test=draw$X_ordnet_test, test_y=draw$test_y, fit=NA))
	}#,
	# settings = list(cutoff = 0.7, q=10)
)


lasso_logit_gen <- new_method("lasso_logit_gen", "L1 Logit",
	method = function(model, draw) {

		K <- length(model$resp_levels)

		p <- ncol(draw$X_ordnet_train)

		# Get lambda (based only on rare category)
		rare_y <- as.integer(as.integer(draw$train_df_polr[,
			model$resp_name]) == model$K)

		stopifnot(all(rare_y %in% c(0, 1)))
		stopifnot(sum(rare_y == 0) > 0)
		stopifnot(sum(rare_y == 1) > 0)
		cv_ret <- cv.glmnet(draw$X_ordnet_train, rare_y, family = "binomial")
		lambda <- cv_ret$lambda.min
	

		lasso_logit <- glmnet(draw$X_ordnet_train, rare_y, family = "binomial")


		return(list(alpha_hat=NA, Beta_hat=NA, X_ordnet_test=draw$X_ordnet_test,
			test_y=draw$test_y, fit=lasso_logit, s=lambda))
	}#,
	# settings = list(cutoff = 0.7, q=10)
)



### Semi-parallel ordinalNet for data application (typical vector response)

semipar_ordnet_vec <- new_method("semipar_ordnet_vec",
	"ON", method = function(model, draw) {

		stopifnot(model$K >= 3)

		# alpha_hat <- rep(as.numeric(NA), model$K - 1)
		# Beta_hat <- matrix(rep(as.numeric(NA), model$p*(model$K - 1)),
		# 	nrow=model$p, ncol=model$K - 1)

		# if(model$tune_prop > 0){
		# 	# Use part of data to get lambda
		
		# 	ret_list <- splitIntoTwoSetsVec(X=draw$X_ordnet_train,
		# 		y=draw$train_y, prop=model$tune_prop)

		# 	tune_X_ordnet <- ret_list$X_prop

		# 	train_X_ordnet <- ret_list$X_comp

		# 	y_tune <- ret_list$y_prop

		# 	y_train <- ret_list$y_comp

		# 	rm(ret_list)

		# 	lambda_star <- getLambdaOrdNet(X=tune_X_ordnet, y=y_tune,
		# 		printProgress=FALSE, nLambda=100, lambdaMinRatio=0.0001)

		# 	# Fit model on rest of data using selected lambda

		# 	fit <- ordinalNet::ordinalNet(train_X_ordnet, y_train,
		# 		family="cumulative", parallelTerms=TRUE, nonparallelTerms=TRUE,
		# 		link="logit", printIter=FALSE, lambdaVals=lambda_star)

		# 	# Don't need estimates--will use ordinalNet predict function later
		# } else{

		# 	lambda_star <- getLambdaOrdNet(X=draw$X_ordnet_train,
		# 		y=draw$train_y, printProgress=FALSE, nLambda=100,
		# 		lambdaMinRatio=0.0001)

		# 	# Fit model on all of data using selected lambda

		# 	fit <- ordinalNet::ordinalNet(draw$X_ordnet_train, draw$train_y,
		# 		family="cumulative", parallelTerms=TRUE, nonparallelTerms=TRUE,
		# 		link="logit", printIter=FALSE, lambdaVals=lambda_star)

		# 	# Don't need estimates--will use ordinalNet predict function later
		# }

		lambda_star <- getLambdaOrdNet(X=draw$X, y=draw$y,
			printProgress=FALSE)

		# Fit model on rest of data using selected lambda

		# ret_list <- getParamEstsOrdNet(X=draw$X, y=draw$y,
		# 	lambdaVals=lambda_star, printIter=FALSE, K=model$K)

		fit <- ordinalNet::ordinalNet(draw$X, draw$y,
				family="cumulative", parallelTerms=TRUE, nonparallelTerms=TRUE,
				link="logit", printIter=FALSE, lambdaVals=lambda_star)

		# alpha_hat <- ret_list$alpha_hat

		# Beta_hat <- ret_list$Beta_hat

		# rm(ret_list)

		# stopifnot(length(alpha_hat) == model$K - 1)
		# stopifnot(nrow(Beta_hat) == model$p)
		# stopifnot(ncol(Beta_hat) == model$K - 1)

		# stopifnot(all(!is.na(alpha_hat)))
		# stopifnot(all(!is.na(Beta_hat)))
		# stopifnot(all(!is.na(alpha_hat)))
		# stopifnot(all(!is.na(Beta_hat)))
		
		# return(list(alpha_hat=NA, Beta_hat=NA, X_ordnet_test=draw$X_ordnet_test,
		# 	test_y=draw$test_y, fit=fit))

		stopifnot(!is.na(fit))
    	stopifnot("ordinalNet" %in% class(fit))

		return(list(alpha_hat=NA, Beta_hat=NA, X=draw$X, y=draw$y,
			probs=draw$probs, beta_final=draw$beta_final, fit=fit))
	}
)


### Non-parallel ordinalNet for data application (typical vector response)

nonpar_ordnet_vec <- new_method("nonpar_ordnet_vec",
	"ordinalNet (non-parallel)", method = function(model, draw) {

		K <- length(model$resp_levels)
		stopifnot(K >= 3)

		p <- ncol(draw$X_ordnet_train)
		n_obs <- nrow(draw$X_ordnet_train)

		alpha_hat <- rep(as.numeric(NA), K - 1)
		Beta_hat <- matrix(rep(as.numeric(NA), p*(K - 1)), nrow=p, ncol=K - 1)

		if(model$tune_prop > 0){
			# Use part of data to get lambda
		
			ret_list <- splitIntoTwoSetsVec(X=draw$X_ordnet_train,
				y=draw$train_y, prop=model$tune_prop)

			tune_X_ordnet <- ret_list$X_prop

			train_X_ordnet <- ret_list$X_comp

			y_tune <- ret_list$y_prop

			y_train <- ret_list$y_comp

			rm(ret_list)

			lambda_star <- getLambdaOrdNet(X=tune_X_ordnet, y=y_tune,
				printProgress=FALSE, nLambda=100, lambdaMinRatio=0.0001)

			# Fit model on rest of data using selected lambda

			fit <- ordinalNet::ordinalNet(train_X_ordnet, y_train,
				family="cumulative", parallelTerms=FALSE, nonparallelTerms=TRUE,
				link="logit", printIter=FALSE, lambdaVals=lambda_star)

			# Don't need estimates--will use ordinalNet predict function later
		} else{

			lambda_star <- getLambdaOrdNet(X=draw$X_ordnet_train,
				y=draw$train_y, printProgress=FALSE, nLambda=100,
				lambdaMinRatio=0.0001)

			# Fit model on all of data using selected lambda

			fit <- ordinalNet::ordinalNet(draw$X_ordnet_train, draw$train_y,
				family="cumulative", parallelTerms=FALSE, nonparallelTerms=TRUE,
				link="logit", printIter=FALSE, lambdaVals=lambda_star)

			# Don't need estimates--will use ordinalNet predict function later
		}
		
		return(list(alpha_hat=NA, Beta_hat=NA, X_ordnet_test=draw$X_ordnet_test,
			test_y=draw$test_y, fit=fit))
	}
)





