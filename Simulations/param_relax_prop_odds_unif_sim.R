# setwd("/Users/gregfaletto/OneDrive - USC Marshall School of Business/Proportional odds generalization/Code/Simulations")

# dev.off()
rm(list=ls())

library(simulator)
library(MASS)
library(ordinalNet)
library(parallel)

dir_main <- getwd()
dir_code <- "/Users/gregfaletto/OneDrive - USC Marshall School of Business/Proportional odds generalization/Code/ordinalNet modified/R"

setwd(dir_code)
source("cdIn.R")
source("cdOut.R")
source("links.R")
source("mirlsNet.R")
source("misc.R")
source("ordinalNet-methods.R")
source("ordinalNet.R")
source("ordinalNetCV.R")
source("ordinalNetTune.R")

setwd(dir_main)

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")


t0 <- Sys.time()

set.seed(2390812)

# n = 2000, dev_size=0.5: proportional odds almost as good as ours, logistic
# regression clearly worse

# n = 5000, dev_size=0.5, 50 sims: mixed bag. ours estimates rare probabilities better
# than proportional odds everywhere, as well as beta. our method seems to
# do the worst at estimating the rare intercept. overall, logistic regression
# seems to beat us by a hair on estimating rare probabilties, and at 
# estimating probabilties only among observations in rare category. we do
# beat logistic regression near the rare decision boundary, and for high
# probability units.

# n = 3000, 50 sims: overall pretty promising. our method does better than proportional
# odds. for logistic regression, looks very similar but mean is worse for ours
# due to outliers. on only rare observations we seem to do a little better.
# definitely better than logistic regression on observations
# near decision boundary and high probability units

# n = 2500 seems to work! doing 100 simulations

# probabilities are estimated better by PRESTO, as is beta vector, though intercept
# is better for logit than PRESTO (not by a statistically significan margin,
# p-value 0.4327586), and much worse for PRESTO than proportional odds (t-stat
# over 14)

# Initialize parallel processing--only works on Mac or Unix
n_cores <- detectCores() - 1

intcpt_list <- list(c(0, 3.5, 5.5), c(0, 4, 6), c(0, 4, 6.5))

param_relax_prop_odds_unif_more2 <- new_simulation("param_relax_prop_odds_unif_more2",
     "Relaxed Proportional Odds (Uniform noise)")

param_relax_prop_odds_unif_more2 <- generate_model(param_relax_prop_odds_unif_more2,
     relax_prop_odds_unif_model, n = 2500,
          p = 10, 
          K = 4,
          # K = as.list(c(4,5)),
          # intercepts=c(0, 4, 6, 10),
          intercepts=intcpt_list,
          beta = rep(2, 20),
          dev_size=0.5
          , vary_along=c("intercepts")
          )

param_relax_prop_odds_unif_more2 <- simulate_from_model(param_relax_prop_odds_unif_more2,
     nsim = 100
          , index = 1:n_cores
          )

print("")
print("")
print("")
print("")
print("")
print("")
print("")
print("")
print("Done generating data! Now running methods (in parallel)...")
print("")
print("")
print("")
print("")
print("")
print("")
print("")
print("")

param_relax_prop_odds_unif_more2 <- run_method(param_relax_prop_odds_unif_more2,
     list(logit_meth, prop_odds_meth
          , fused_polr
          # , semipar_ordnet_vec
          )
     , parallel=list(socket_names=n_cores)
     )
     # evaluate(list(beta_mse, rare_intcp_mse, rare_prob_log_loss, prop_rare_obs))
     # rare_beta_mse, rare_intcp_mse, r

param_relax_prop_odds_unif_more2 <- evaluate(param_relax_prop_odds_unif_more2,
     list(prop_rare_obs, rare_prob_mse_gen
          # , rare_prob_log_loss_gen, rare_prob_mse_gen_only_rare,
          # prop_rare_obs, log_lik_gen, rare_prob_log_loss_gen_only_rare,
          # rare_prob_log_loss_gen_near_db, rare_prob_mse_gen_near_db,
          # rare_prob_log_loss_gen_high_prob, rare_prob_mse_gen_high_prob,
          # rare_prob_mae_gen, rare_prob_mae_gen_only_rare,
          # rare_prob_mae_gen_near_db, rare_prob_mae_gen_high_prob
          ))

     ### TODO: add metric that weights rare class and common class equally
     # regardless of number of observations?

# param_relax_prop_odds <- new_simulation("param_relax_prop_odds",
#      "Relaxed Proportional Odds") |>
#      generate_model(relax_prop_odds_model, n = 5000, p = as.list(c(10, 20)),
#           K = as.list(c(4, 5)),
#           intercepts=c(0, 4, 6, 10), beta = rep(1, 20),
#           dev_size=as.list(c(.3, .7)), dev_prob=as.list(c(.2, .3)),
#           vary_along = c("K", "dev_size", "p", "dev_prob")) |>
#      simulate_from_model(nsim = 1) |>
#      run_method(list(logit_meth, prop_odds_meth, fused_polr)) |>
#      # evaluate(list(beta_mse, rare_intcp_mse, rare_prob_log_loss, prop_rare_obs))
#      evaluate(list(rare_beta_mse, rare_intcp_mse, rare_prob_log_loss_gen,
#           prop_rare_obs, log_lik_gen))

# param_relax_prop_odds <- new_simulation("param_relax_prop_odds",
#      "Relaxed Proportional Odds") |>
#      generate_model(relax_prop_odds_model, n = 5000, p = as.list(c(10, 20)),
#           K = as.list(c(4, 5)),
#           intercepts=c(0, 4, 6, 10), beta = rep(1, 20),
#           dev_size=as.list(c(.1, .3, .5, .7)), dev_prob=as.list(c(.1, .2, .3)),
#           vary_along = c("K", "dev_size", "p", "dev_prob")) |>
#      simulate_from_model(nsim = 1) |>
#      run_method(list(logit_meth, prop_odds_meth, fused_polr)) |>
#      # evaluate(list(beta_mse, rare_intcp_mse, rare_prob_log_loss, prop_rare_obs))
#      evaluate(list(rare_beta_mse, rare_intcp_mse, rare_prob_log_loss_gen,
#           prop_rare_obs, log_lik_gen))

# if(K >= 4){
#      less_rare_intercept_mses <- rep(as.numeric(NA), nsims)
#      number_less_rare_obs <- rep(as.numeric(NA), nsims)
# }

save_simulation(param_relax_prop_odds_unif_more2)

# param_relax_prop_odds_unif <- load_simulation("param_relax_prop_odds_unif")

# param_relax_prop_odds_unif <- run_method(param_relax_prop_odds_unif,
#     methods=c(lassoSS_phat_cor)) %>% evaluate(list(phat, labels))

# 

# param_relax_prop_odds_unif <- param_relax_prop_odds_unif |> evaluate(probs)

print("Done! Total time for simulations:")
t1 <- Sys.time()
print(t1 - t0)

# tabulate_eval(param_relax_prop_odds_unif, "rare_prob_mse_gen")

# subset_simulation(sim=param_relax_prop_odds_unif
#                   , methods = c("prop_odds_meth", "fused_polr", "logit_meth")) |> evaluate(rare_prob_mse_gen) |>
#     tabulate_eval("rare_prob_mse_gen")

# param_relax_prop_odds_rand <- load_simulation("param_relax_prop_odds_rand")

# subset_simulation(param_relax_prop_odds_unif, methods=c("prop_odds_meth",
#      "fused_polr", "logit_meth")) |> plot_eval("rare_prob_mse_gen") +
#      ggtitle("Uniform Differences")

# subset_simulation(param_relax_prop_odds_unif, methods=c("prop_odds_meth",
#      "fused_polr", "logit_meth")) |> tabulate_eval("rare_prob_mse_gen")

# subset_simulation(sim=param_relax_prop_odds_unif, index=1) %>%
#      plot_eval("prop_rare_obs")

# param_relax_prop_odds_unif |> evaluate(tpr_gen) |> plot_eval("tpr_gen")
# param_relax_prop_odds_unif |> evaluate(tpr_gen) |> tabulate_eval("tpr_gen")

# subset_simulation(sim=param_relax_prop_odds_unif
#      , methods = c("prop_odds_meth", "fused_polr")) %>%
#      plot_eval("rare_prob_mse_gen")

# subset_simulation(sim=param_relax_prop_odds_unif
#      , methods = c("prop_odds_meth", "fused_polr")) %>%
#      tabulate_eval("rare_prob_log_loss")

# subset_simulation(sim=param_relax_prop_odds_unif
#      # , methods = c("lasso", "lasso_refit")
#      , K == 4, p == 10, dev_size == 0.7, dev_prob == .3) %>%
#      plot_eval("prop_rare_obs")

# subset_simulation(sim=param_relax_prop_odds_unif
#      # , methods = c("lasso", "lasso_refit")
#      , n==1000) %>%
#      plot_eval("rare_prob_log_loss_gen")

# print(plot_eval(param_relax_prop_odds_unif, "prop_rare_obs"))
# print(plot_eval(param_relax_prop_odds_unif, "rare_beta_mse"))
# print(plot_eval(param_relax_prop_odds_unif, "rare_intcp_mse"))
print(plot_eval(param_relax_prop_odds_unif_more2, "rare_prob_mse_gen"))
# print(plot_eval(param_relax_prop_odds_unif, "log_lik_gen"))
# print(plot_eval(param_relax_prop_odds_unif, "rare_prob_mse_gen_only_rare"))
# print(plot_eval(param_relax_prop_odds_unif, "rare_prob_mse_gen_near_db"))
# print(plot_eval(param_relax_prop_odds_unif, "rare_prob_mse_gen_high_prob"))


# # print("Average beta MSE:")
# # print(mean(common_intercept_mses))

# # print("Average common intercept MSE:")
# # print(mean(beta_mses))

# # print("Average rare intercept MSE:")
# # print(mean(rare_intercept_mses))

# # print("Average number of rare class observations:")
# # print(mean(number_rare_obs))


# Plot functions
setwd("/Users/gregfaletto/Dropbox/Jacob and Greg/AISTATS 2023")
source("sim_eval_function.R")
setwd(dir_main)
create_plots(param_relax_prop_odds_unif_more2)



