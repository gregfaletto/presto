## @knitr models
relax_prop_odds_model_rand <- function(n, p, K, intercepts, beta, dev_size,
    dev_prob){
    # Number of classes
    stopifnot(is.numeric(K) | is.integer(K))
    stopifnot(K == round(K))
    stopifnot(K >= 3)
    stopifnot(is.numeric(intercepts) | is.integer(intercepts))
    stopifnot(length(intercepts) >= K - 1)
    if(length(intercepts) > K - 1){
        intercepts <- intercepts[1:(K - 1)]
    }
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))
    stopifnot(is.numeric(p) | is.integer(p))
    stopifnot(p == round(p))
    stopifnot(p >= 1)
    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) >= p)
    if(length(beta) > p){
        beta <- beta[1:p]
    }

    stopifnot(is.numeric(dev_prob) | is.integer(dev_prob))
    stopifnot(!is.na(dev_prob))
    stopifnot(length(dev_prob) == 1)
    stopifnot(dev_prob >= 0)
    stopifnot(dev_prob <= 1)
    
    my_model <- new_model(name = "relax_prop_odds_mod", 
        label = sprintf("Relaxed proportional odds model (n = %s,
            p = %s, K = %s, rare intercept = %s, deviation size = %s)",
            n, p, K, intercepts[K-1], dev_size),
        params = list(n = n, p = p, K = K, intercepts=intercepts,
            beta = beta, dev_prob=dev_prob, dev_size = dev_size),
        simulate = sim_relaxed_prop_odds_rand
    )
    return(my_model)
}

sim_relaxed_prop_odds_rand <- function(n, p, K, intercepts, beta, nsim,
    dev_prob, dev_size=0){

    stopifnot(K >= 3)
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))

    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) == p)

    stopifnot(is.numeric(dev_prob) | is.integer(dev_prob))
    stopifnot(!is.na(dev_prob))
    stopifnot(length(dev_prob) == 1)
    stopifnot(dev_prob >= 0)
    stopifnot(dev_prob <= 1)


    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()

    for(i in 1:nsim){

        # If beta_mat that is generated results in negative probabilities,
        # try again with another random draw, up to 10 times before giving up

        ret_i <- NA
        n_tries <- 0

        while(any(is.na(ret_i))){

            if(n_tries >= 100){
                print(paste("p =", p, ", K =", K, "dev_num =", dev_num,
                    "dev_size =", dev_size))
                print("intercepts:")
                print(intercepts)
                print("colMeans(probs):")
                print(colMeans(probs))
                stop("At least one probability was negative")
            }

            # Generate beta_mat

            beta_mat <- matrix(0, nrow=p, ncol=K - 1)

            beta_mat[, 1] <- beta

            # Select random indices for base beta vector to keep
            inds <- as.logical(rbinom(p, size=1, prob=dev_prob))
            reject_inds <- which(!inds)
            
            beta_mat[reject_inds, 1] <- 0
            

            # Generate sparse random deviations
            dev_mat <- matrix(0, nrow=p, ncol=K - 2)

            for(k in 1:(K - 2)){
                # Select random indices for deviations
                inds <- which(as.logical(rbinom(p, size=1, prob=dev_prob)))
                dev_mat[inds, k] <- dev_size
            }

            # Create random sign flips to multiply deviations by
            sign_mat_content <- sample(c(-1, 1), size=p*(K - 2), replace=TRUE)
            sign_mat <- matrix(sign_mat_content, nrow=p, ncol=K - 2)

            # Multiply deviations by random sign flips
            dev_mat <- dev_mat * sign_mat

            # Add sparse random deviations to beta_mat
            for(k in 2:(K-1)){
                beta_mat[, k] <- beta_mat[, k - 1] + dev_mat[, k - 1]
            }

            # Generate X, y, probabilities
            ret_i <- sim_data(n, p, K, intercepts, beta_mat, dev_num,
                dev_size)

            n_tries <- n_tries + 1

        }

        ret_list[[i]] <- ret_i
        
    }
    
    return(ret_list)
}

relax_prop_odds_unif_model <- function(n, p, K, intercepts, beta, dev_size){
    # Number of classes
    stopifnot(is.numeric(K) | is.integer(K))
    stopifnot(K == round(K))
    stopifnot(K >= 3)
    stopifnot(is.numeric(intercepts) | is.integer(intercepts))
    stopifnot(length(intercepts) >= K - 1)
    if(length(intercepts) > K - 1){
        intercepts <- intercepts[1:(K - 1)]
    }
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))
    stopifnot(is.numeric(p) | is.integer(p))
    stopifnot(p == round(p))
    stopifnot(p >= 1)
    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) >= p)
    if(length(beta) > p){
        beta <- beta[1:p]
    }

    stopifnot(is.numeric(dev_size) | is.integer(dev_size))
    stopifnot(!is.na(dev_size))
    stopifnot(length(dev_size) == 1)
    stopifnot(dev_size > 0)
    
    my_model <- new_model(name = "relax_prop_odds_mod_unif", 
        label = sprintf("Relaxed proportional odds model with uniform noise (n = %s,
            p = %s, K = %s, rare intercept = %s, deviation size = %s)",
            n, p, K, intercepts[K-1], dev_size),
        params = list(n = n, p = p, K = K, intercepts=intercepts,
            beta = beta, dev_size=dev_size),
        simulate = sim_relaxed_unif_prop_odds
    )
    return(my_model)
}

sim_relaxed_unif_prop_odds <- function(n, p, K, intercepts, beta, nsim,
    dev_size){

    stopifnot(K >= 3)
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))

    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) == p)

    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()

    # dev_num <- round(sqrt(p))

    for(i in 1:nsim){

        # If beta_mat that is generated results in negative probabilities,
        # try again with another random draw, up to 10 times before giving up
        ret_i <- NA
        n_tries <- 0

        while(any(is.na(ret_i))){

            if(n_tries >= 100){
                print(paste("p =", p, ", K =", K))
                print("intercepts:")
                print(intercepts)
                stop("At least one probability was negative")
            }

            # Generate beta_mat

            beta_mat <- matrix(0, nrow=p, ncol=K - 1)

            # Generate uniform deviations for base beta
            beta_mat[, 1] <- runif(p, min=-dev_size, max=dev_size)

            # Generate uniform deviations
            dev_mat <- matrix(runif(p*(K - 2), min=-dev_size, max=dev_size),
                nrow=p, ncol=K - 2)

            # Add uniform deviations to beta_mat
            for(k in 2:(K-1)){
                beta_mat[, k] <- beta_mat[, k - 1] + dev_mat[, k - 1]
            }

            # Generate X, y, probabilities
            ret_i <- sim_data(n, p, K, intercepts, beta_mat)

            n_tries <- n_tries + 1

            ret_list[[i]] <- ret_i

        }

    }
    
    return(ret_list)
}

gen_data_app_model_prediabetes <- function(df, tune_prop, test_prop, resp_name,
    resp_levels, x_formula, omitted_x, age_cutoff){

    stopifnot(is.data.frame(df))
    stopifnot(all(!is.na(df)))

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop >= 0)
    stopifnot(test_prop < 1)

    stopifnot(is.numeric(tune_prop) | is.integer(tune_prop))
    stopifnot(length(tune_prop) == 1)
    stopifnot(tune_prop >= 0)
    stopifnot(tune_prop < 1)

    stopifnot(is.character(resp_levels))
    stopifnot(all(!is.na(resp_levels)))

    K <- length(resp_levels)

    stopifnot(K >= 3)

    stopifnot(is.character(resp_name))
    stopifnot(length(resp_name) == 1)
    

    if(any(!is.na(omitted_x))){
        stopifnot(is.character(omitted_x))
        stopifnot(omitted_x %in% colnames(df))
    } else{
        stopifnot(all(is.na(omitted_x)))
    }
    
    my_model <- new_model(name = "dat_app_mod", 
                label = "Data Application",
                params = list(df=df, test_prop=test_prop, resp_name=resp_name,
                    resp_levels=resp_levels, tune_prop=tune_prop,
                    x_formula=x_formula, omitted_x=omitted_x,
                    age_cutoff=age_cutoff),
                simulate = gen_sim_data_app_prediabetes
    )

    return(my_model)
}

gen_sim_data_app_prediabetes <- function(df, test_prop, resp_name, resp_levels, x_formula,
    omitted_x, age_cutoff, nsim){

    n <- nrow(df)

    response <- rep("NoDiabetes", n)

    # Takes about an hour for each, so I can run all of these in about 7 hours

    response[(df$Age_PreDiabetes < age_cutoff) &
         (df$Age_Diabetes >= age_cutoff)] <- "PreDiabetes"

    response[df$Age_Diabetes < age_cutoff] <- "Diabetes"

    df$response <- factor(response, levels=RESP_LEVELS, ordered=TRUE)

    # print("category proportions")
    # print(round(summary(df$response)/n, 4))

    df <- df[, !(colnames(df) %in% c("Age_PreDiabetes",
         "Age_Diabetes", "Time_Pre_To_Diabetes",
         "PreDiabetes_Checks_Before_Diabetes"))]

    df$Sex <- as.factor(df$Sex)

    # Check inputs

    stopifnot(is.data.frame(df))

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop > 0)
    stopifnot(test_prop < 1)

    stopifnot(is.character(resp_levels))
    stopifnot(all(!is.na(resp_levels)))

    stopifnot(is.character(resp_name))
    stopifnot(length(resp_name) == 1)
    stopifnot(resp_name %in% colnames(df))

    K <- length(resp_levels)

    stopifnot(K >= 3)

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop >= 0)
    stopifnot(test_prop < 1)

    if(any(!is.na(omitted_x))){
        stopifnot(is.character(omitted_x))
        stopifnot(omitted_x %in% colnames(df))
    } else{
        stopifnot(all(is.na(omitted_x)))
    }
    
    # Create design matrix and response matrix for ordinalNet

    X_ordnet <- model.matrix(x_formula, data=df)

    stopifnot(nrow(X_ordnet) == nrow(df))

    X_ordnet <- X_ordnet[, colnames(X_ordnet) != "(Intercept)"]

    p <- ncol(X_ordnet)

    n_obs <- nrow(X_ordnet)

    stopifnot(is.matrix(X_ordnet))
    stopifnot(all(!is.na(X_ordnet)))
    stopifnot(n_obs == nrow(df))

    ret <- list()

    for(i in 1:nsim){

        # Split y_mat into training and test sets

        ret_list <- splitIntoTwoSetsVec(X=X_ordnet, y=df[, resp_name],
            prop=test_prop, n_attempts=100)

        X_ordnet_test <- ret_list$X_prop
        X_ordnet_train <- ret_list$X_comp
        y_test <- ret_list$y_prop
        y_train <- ret_list$y_comp
        train_inds <- ret_list$comp_inds

        rm(ret_list)

        # Create data.frame for MASS::polr function

        if(all(!is.na(omitted_x))){
            df <- df[, !(colnames(df) %in% omitted_x)]
        }
        
        df_polr_train <- df[train_inds, ]
        # No need for test df--will just extract coefficients

        # Check outputs
        stopifnot(ncol(X_ordnet_train) == ncol(X_ordnet_test))

        stopifnot(nrow(X_ordnet_train) == length(y_train))
        stopifnot(nrow(df_polr_train) == length(y_train))

        stopifnot(nrow(X_ordnet_test) == length(y_test))

        stopifnot(is.factor(y_train))
        stopifnot(is.factor(y_test))

        stopifnot(length(unique(y_train)) == K)
        stopifnot(length(unique(y_test)) == K)

        ret[[i]] <- list(train_df_polr=df_polr_train,
            X_ordnet_train=X_ordnet_train, train_y=y_train,
            test_y=y_test, X_ordnet_test=X_ordnet_test)
    }

    return(ret)
}

gen_data_app_model <- function(df, tune_prop, test_prop, resp_name,
    resp_levels, x_formula, omitted_x){

    stopifnot(is.data.frame(df))
    stopifnot(all(!is.na(df)))

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop >= 0)
    stopifnot(test_prop < 1)

    stopifnot(is.numeric(tune_prop) | is.integer(tune_prop))
    stopifnot(length(tune_prop) == 1)
    stopifnot(tune_prop >= 0)
    stopifnot(tune_prop < 1)

    stopifnot(is.character(resp_levels))
    stopifnot(all(!is.na(resp_levels)))

    K <- length(resp_levels)

    stopifnot(K >= 3)

    stopifnot(is.character(resp_name))
    stopifnot(length(resp_name) == 1)
    stopifnot(resp_name %in% colnames(df))

    if(any(!is.na(omitted_x))){
        stopifnot(is.character(omitted_x))
        stopifnot(omitted_x %in% colnames(df))
    } else{
        stopifnot(all(is.na(omitted_x)))
    }
    
    my_model <- new_model(name = "dat_app_mod", 
                label = "Data Application",
                params = list(df=df, test_prop=test_prop, resp_name=resp_name,
                    resp_levels=resp_levels, tune_prop=tune_prop,
                    x_formula=x_formula, omitted_x=omitted_x),
                simulate = gen_sim_data_app
    )

    return(my_model)
}

gen_sim_data_app <- function(df, test_prop, resp_name, resp_levels, x_formula,
    omitted_x, nsim){

    # Check inputs

    stopifnot(is.data.frame(df))

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop > 0)
    stopifnot(test_prop < 1)

    stopifnot(is.character(resp_levels))
    stopifnot(all(!is.na(resp_levels)))

    stopifnot(is.character(resp_name))
    stopifnot(length(resp_name) == 1)
    stopifnot(resp_name %in% colnames(df))

    K <- length(resp_levels)

    stopifnot(K >= 3)

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop >= 0)
    stopifnot(test_prop < 1)

    if(any(!is.na(omitted_x))){
        stopifnot(is.character(omitted_x))
        stopifnot(omitted_x %in% colnames(df))
    } else{
        stopifnot(all(is.na(omitted_x)))
    }
    

    # Create design matrix and response matrix for ordinalNet
    X_ordnet <- model.matrix(x_formula, data=df)

    stopifnot(nrow(X_ordnet) == nrow(df))

    X_ordnet <- X_ordnet[, colnames(X_ordnet) != "(Intercept)"]

    p <- ncol(X_ordnet)

    n_obs <- nrow(X_ordnet)

    stopifnot(is.matrix(X_ordnet))
    stopifnot(all(!is.na(X_ordnet)))
    stopifnot(n_obs == nrow(df))

    ret <- list()

    for(i in 1:nsim){

        # Split y_mat into training and test sets

        ret_list <- splitIntoTwoSetsVec(X=X_ordnet, y=df[, resp_name],
            prop=test_prop, n_attempts=100)

        X_ordnet_test <- ret_list$X_prop
        X_ordnet_train <- ret_list$X_comp
        y_test <- ret_list$y_prop
        y_train <- ret_list$y_comp
        train_inds <- ret_list$comp_inds

        rm(ret_list)

        # Create data.frame for MASS::polr function

        if(all(!is.na(omitted_x))){
            df <- df[, !(colnames(df) %in% omitted_x)]
        }
        
        df_polr_train <- df[train_inds, ]
        # No need for test df--will just extract coefficients

        # Check outputs

        stopifnot(ncol(X_ordnet_train) == ncol(X_ordnet_test))

        stopifnot(nrow(X_ordnet_train) == length(y_train))
        stopifnot(nrow(df_polr_train) == length(y_train))

        stopifnot(nrow(X_ordnet_test) == length(y_test))

        stopifnot(is.factor(y_train))
        stopifnot(is.factor(y_test))

        stopifnot(length(unique(y_train)) == K)
        stopifnot(length(unique(y_test)) == K)

        ret[[i]] <- list(train_df_polr=df_polr_train,
            X_ordnet_train=X_ordnet_train, train_y=y_train,
            test_y=y_test, X_ordnet_test=X_ordnet_test)
    }

    return(ret)
}

sim_data <- function(n, p, K, intercepts, beta_mat, dev_num=NA, dev_size=NA){

    stopifnot(K >= 3)
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))

    stopifnot(is.numeric(beta_mat) | is.integer(beta_mat))
    stopifnot(nrow(beta_mat) == p)
    stopifnot(ncol(beta_mat) == K - 1)

    # If we don't observe at least one observation from each class, toss
    # out simulation and try again
    all_classes <- FALSE
    counter <- 0
    while(!all_classes){

        if(counter >= 10){
            return(NA)
        }
        
        # Generate uniformly distributed X
        if(p > 1){
            X <- matrix(runif(n*p, min=-1, max=1), nrow=n, ncol=p)
        } else{
            X <- runif(n, min=-1, max=1)
        }

        # Probabilities of each class
        probs <- matrix(as.numeric(NA), nrow=n, ncol=K)
        if(p > 1){
            probs[, 1] <- plogis(intercepts[1] + X %*% beta_mat[, 1])
        } else{
            probs[, 1] <- plogis(intercepts[1] + beta_mat[, 1]*X)
        }

        if(p > 1){
            for(k in 2:(K-1)){
                if(k > 2){
                    probs[, k] <- plogis(intercepts[k] + X %*% beta_mat[, k]) -
                        rowSums(probs[, 1:(k - 1)])
                } else{
                    probs[, k] <- plogis(intercepts[k] + X %*% beta_mat[, k]) -
                        probs[, 1]
                }
                # stopifnot(all(probs[, k] > 0))
            }
        } else{
            for(k in 2:(K-1)){
                if(k > 2){
                    probs[, k] <- plogis(intercepts[k] + beta_mat[, k]*X) -
                        rowSums(probs[, 1:(k - 1)])
                } else{
                    probs[, k] <- plogis(intercepts[k] + beta_mat[, k]*X) -
                        probs[, 1]
                }
                # stopifnot(all(probs[, k] > 0))
            } 
        }

        stopifnot(all(rowSums(probs[, 1:(K-1)]) <= 1))

        probs[, K] <- 1 - rowSums(probs[, 1:(K-1)])

        stopifnot(all(!is.na(probs)))
        stopifnot(all(abs(rowSums(probs) - 1) < 10^(-7)))

        if(any(probs < 0)){
            return(NA)
        }
        if(any(probs > 1)){
            return(NA)
        }

        stopifnot(all(probs >= 0))
        stopifnot(all(probs <= 1))

        # Generate y
        y <- rep(as.numeric(NA), n)
        for(j in 1:n){
            y[j] <- which(rmultinom(n=1, size=1, prob=probs[j,])[, 1] != 0)
        }
        stopifnot(all(!is.na(y)))
        stopifnot(all(y %in% 1:K))

        # Make sure at least one observation from each class; otherwise, toss 
        # out this simulation
        all_classes <- (length(unique(y)) == K)

        stopifnot(length(unique(y)) <= K)

        counter <- counter + 1
    }

    # y <- as.factor(y)
    y <- factor(y, levels=1:K, ordered=TRUE)

    
    return(list(X=X, y=y, probs=probs, beta_final=beta_mat))
}

splitIntoTwoSetsVec <- function(X, y, prop, n_attempts=10){

     # Given a proportion prop, returns a matrix X and vector y that constitute
     # that proportion of the data in (X_prop, y_prop) as well as another
     # (X_comp, y_comp) containing the complement of (X_prop, y_prop) (the 
     # remainder of the data). Ensures at least one outcome from every category
     # in y occur in both y_prop and y_comp (or throws an error if after 10
     # attempts no such set is found).

     # Check inputs

     stopifnot(is.numeric(prop) | is.integer(prop))
     stopifnot(length(prop) == 1)
     stopifnot(prop > 0)
     stopifnot(prop < 1)

     stopifnot(is.matrix(X))
     n_obs <- nrow(X)

     stopifnot(n_obs >= 4)

     stopifnot(length(y) == n_obs)
     stopifnot(is.factor(y))

     K <- length(unique(y))

     stopifnot(K >= 3)
     stopifnot(K < n_obs/2)

     # Want to ensure that in both subsamples there is at least one
     # observation from every class
     all_cols_flag <- FALSE
     iter_count <- 0

     n_prop <- round(prop*n_obs)

     stopifnot(n_prop >= 2)
     stopifnot(n_prop <= n_obs - 2)

     while(!all_cols_flag){
          
          if(iter_count >= n_attempts){
               err_message <- paste("Unable to create subsampled matrix with at least one observation in each column after",
                    n_attempts, "attempts; breaking loop here.")
               stop(err_message)
          }

          prop_inds <- sample(n_obs, n_prop)
          comp_inds <- setdiff(1:n_obs, prop_inds)

          X_prop <- X[prop_inds, ]
          y_prop <- y[prop_inds]

          X_comp <- X[comp_inds, ]
          y_comp <- y[comp_inds]

          stopifnot(length(unique(y_prop)) <= K)
          stopifnot(length(unique(y_comp)) <= K)

          all_cols_flag <- (length(unique(y_prop)) == K) &
               (length(unique(y_comp)) == K)
          iter_count <- iter_count + 1
     }

     colnames(X_prop) <- colnames(X)
     rownames(X_prop) <- rownames(X)[prop_inds]

     colnames(X_comp) <- colnames(X)
     rownames(X_comp) <- rownames(X)[comp_inds]

     names(y_comp) <- names(y)[comp_inds]
     names(y_prop) <- names(y)[prop_inds]

     # Check outputs

     stopifnot(ncol(X_prop) == ncol(X_comp))

     stopifnot(nrow(X_prop) == length(y_prop))
     stopifnot(nrow(X_comp) == length(y_comp))

     stopifnot(is.factor(y_prop))
     stopifnot(is.factor(y_comp))

     stopifnot(length(unique(y_prop)) == K)
     stopifnot(length(unique(y_comp)) == K)

     stopifnot(length(y_prop) + length(y_comp) == n_obs)

     return(list(X_prop=X_prop, X_comp=X_comp, y_prop=y_prop, y_comp=y_comp,
          prop_inds=prop_inds, comp_inds=comp_inds))
}