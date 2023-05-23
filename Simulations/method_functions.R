### Logistic regression on rare class

logit_meth <- new_method("logit_meth", "Logit",
	method = function(model, draw) {

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
	}
)

### Proportional Odds

prop_odds_meth <- new_method("prop_odds_meth", "PO",
	method = function(model, draw) {

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

	}
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
	}
)

### PRESTO
fused_polr <- new_method("fused_polr", "PRESTO",
	method = function(model, draw) {

		stopifnot(model$K >= 3)

		# Use data to get lambda
		lambda_star <- getLambdaFusedPolr(X=draw$X, y=draw$y,
			printProgress=FALSE)

		# Fit model on rest of data using selected lambda and get parameter estimates
		ret_list <- getParamEstsMod(X=draw$X, y=draw$y, lambdaVals=lambda_star,
			printIter=FALSE, K=model$K)

		alpha_hat <- ret_list$alpha_hat

		Beta_hat <- ret_list$Beta_hat

		rm(ret_list)

		stopifnot(length(alpha_hat) == model$K - 1)
		stopifnot(nrow(Beta_hat) == model$p)
		stopifnot(ncol(Beta_hat) == model$K - 1)
		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))
		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))

		return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat, X=draw$X, y=draw$y,
			probs=draw$probs, beta_final=draw$beta_final, fit=NA))
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

# L2 version
fused_polr_l2 <- new_method("fused_polr_l2", "PRESTO_L2",
	method = function(model, draw) {

		stopifnot(model$K >= 3)

		# Use data to get lambda
		lambda_star <- getLambdaFusedPolr(X=draw$X, y=draw$y,
			printProgress=FALSE, alpha_param=0)

		# Fit model on rest of data using selected lambda and get parameter estimates
		ret_list <- getParamEstsMod(X=draw$X, y=draw$y, lambdaVals=lambda_star,
			printIter=FALSE, K=model$K, alpha_param=0)

		alpha_hat <- ret_list$alpha_hat

		Beta_hat <- ret_list$Beta_hat

		rm(ret_list)

		stopifnot(length(alpha_hat) == model$K - 1)
		stopifnot(nrow(Beta_hat) == model$p)
		stopifnot(ncol(Beta_hat) == model$K - 1)

		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))
		stopifnot(all(!is.na(alpha_hat)))
		stopifnot(all(!is.na(Beta_hat)))

		return(list(alpha_hat=alpha_hat, Beta_hat=Beta_hat, X=draw$X, y=draw$y,
			probs=draw$probs, beta_final=draw$beta_final, fit=NA))
	}
)

### PRESTO for data application (typical vector response)
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