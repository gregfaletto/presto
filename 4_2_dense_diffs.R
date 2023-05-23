# setwd("/Users/gregfaletto/Documents/GitHub/presto")

if(!is.null(dev.list())){
     dev.off()
}
rm(list=ls())

library(simulator)
library(MASS)
library(ordinalNet)
library(parallel)
library(cowplot)
library(ggplot2)

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

# Initialize parallel processing--only works on Mac or Unix
n_cores <- detectCores() - 1

intcpt_list <- list(c(0, 2.5, 4.5), c(0, 3, 5), c(0, 3.5, 5.5), c(0, 4, 6))

dense_sim <- new_simulation("dense_sim",
     "Relaxed Proportional Odds (Uniform noise)")

dense_sim <- generate_model(dense_sim, relax_prop_odds_unif_model, n = 2500,
     p = 10, K = 4, intercepts=intcpt_list, beta = rep(2, 10), dev_size=0.5,
     vary_along=c("intercepts"))

dense_sim <- simulate_from_model(dense_sim, nsim = 1, index = 1:n_cores)

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

dense_sim <- run_method(dense_sim, list(logit_meth, prop_odds_meth, fused_polr),
     parallel=list(socket_names=n_cores))

dense_sim <- evaluate(dense_sim, list(prop_rare_obs, rare_prob_mse_gen))

# save_simulation(dense_sim)

print("Done! Total time for simulations:")
t1 <- Sys.time()
print(t1 - t0)

print(plot_eval(dense_sim, "rare_prob_mse_gen"))

# Plot functions
create_plots(dense_sim)

df_sim_stats(dense_sim)

