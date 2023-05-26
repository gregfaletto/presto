# setwd("/Users/gregfaletto/Documents/GitHub/presto")

if(!is.null(dev.list())){
     dev.off()
}
rm(list=ls())

library(simulator)
library(MASS)
library(parallel)
library(cowplot)
library(ggplot2)
library(stargazer)

dir_main <- getwd()
dir_ordnet <- paste(dir_main, "/ordinalNet modified", sep="")
setwd(dir_ordnet)

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
dir_code <- paste(dir_main, "/Simulations", sep="")
setwd(dir_code)

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")
source("sim_eval_function.R")

setwd(dir_main)

t0 <- Sys.time()

set.seed(2390812)

# Initialize parallel processing--only works on Mac or Unix.
n_cores <- 7
stopifnot(n_cores <= detectCores())
if(n_cores == 7){
     nsims <- 100
} else{
     nsims <- 700
}


# About 0.88 minutes per simulation
intcpt_list <- list(c(0, 3, 5), c(0, 3.5, 5.5), c(0, 4, 6))

sparse_sim <- new_simulation("sparse_sim", "Relaxed Proportional Odds")

sparse_sim <- generate_model(sparse_sim, relax_prop_odds_model_rand, n = 2500,
          p = 10, K = 4, intercepts=intcpt_list, beta = rep(.5, 10), dev_size=.5,
          dev_prob=list(1/3, 1/2), vary_along=c("intercepts", "dev_prob"))

sparse_sim <- simulate_from_model(sparse_sim, nsim = nsims, index = 1:n_cores)

print("")
print("")
print("")
print("")
print("")
print("")
print("")
print("")
print("Done generating data! Now running methods (in parallel)...")
print("Time elapsed:")
print(Sys.time() - t0)
print("")
print("")
print("")
print("")
print("")
print("")
print("")
print("")

sparse_sim <- run_method(sparse_sim, list(logit_meth, prop_odds_meth,
     fused_polr, fused_polr_l2), parallel=list(socket_names=n_cores))


print("")
print("")
print("Done running methods! Now evaluating...")
print("Time elapsed:")
print(Sys.time() - t0)
print("")
print("")
print("")

sparse_sim <- evaluate(sparse_sim, list(prop_rare_obs, rare_prob_mse_gen))

save_simulation(sparse_sim)

print("Done! Total time for simulations:")
t1 <- Sys.time()
print(t1 - t0)

sparse_plots_1_2 <- create_sparse_plots(subset_simulation(sparse_sim,
     methods=c("logit_meth", "prop_odds_meth", "fused_polr")))

# Figure 6 (sparsity 1/2)
fig_6 <- sparse_plots_1_2$main_plot

# Figure 7 (sparsity 1/2)
fig_7 <- sparse_plots_1_2$supp_plot

sparse_plots_1_3 <- create_sparse_plots(subset_simulation(sparse_sim,
     methods=c("logit_meth", "prop_odds_meth", "fused_polr")), plots=c(2, 4, 6))

# Figure 1 (sparsity 1/3)
fig_1 <- sparse_plots_1_3$main_plot

# Figure 5 (sparsity 1/3)
fig_5 <- sparse_plots_1_3$supp_plot

ret <- df_sim_stats(subset_simulation(sparse_sim, methods=c("logit_meth",
     "prop_odds_meth", "fused_polr")), methods_to_compare=c("logit_meth",
     "prop_odds_meth" ))

# Table 3
stargazer(ret$t_d_df, summary=FALSE)

# Table 4
stargazer(ret$summary_df, summary=FALSE)


# Ridge presto results

# Figure 10
fig_10 <- create_sparse_plot2_ridge(sparse_sim, plots=1)

# Figure 11
fig_11 <- create_sparse_plot2_ridge(sparse_sim, plots=2)

ret_ridge <- df_sim_stats(subset_simulation(sparse_sim, methods=c("logit_meth",
     "prop_odds_meth", "fused_polr", "fused_polr_l2")),
     methods_to_compare=c("logit_meth", "prop_odds_meth", "fused_polr_l2"))

# Table 3
stargazer(ret_ridge$t_d_df, summary=FALSE)

# Table 4
stargazer(ret_ridge$summary_df, summary=FALSE)



