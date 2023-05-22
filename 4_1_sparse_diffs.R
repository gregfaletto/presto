# setwd("/Users/gregfaletto/Documents/GitHub/presto")

# dev.off()
rm(list=ls())

library(simulator)
library(MASS)
library(ordinalNet)
library(parallel)

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

setwd(dir_main)

t0 <- Sys.time()

set.seed(2390812)

# Initialize parallel processing--only works on Mac or Unix
n_cores <- detectCores() - 1

# About 0.88 minutes per simulation
intcpt_list <- list(c(0, 3, 5), c(0, 3.5, 5.5), c(0, 4, 6))

sparse_sim <- new_simulation("sparse_sim", "Relaxed Proportional Odds")

sparse_sim <- generate_model(sparse_sim, relax_prop_odds_model_rand, n = 2500,
          p = 10, K = 4, intercepts=intcpt_list, beta = rep(1, 10), dev_size=.5,
          dev_prob=1/3, vary_along=c("intercepts"))

sparse_sim <- simulate_from_model(sparse_sim, nsim = 100, index = 1:n_cores)

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

sparse_sim <- run_method(sparse_sim, list(logit_meth, prop_odds_meth,
     fused_polr), parallel=list(socket_names=n_cores))

sparse_sim <- evaluate(sparse_sim, list(rare_prob_mse_gen))

save_simulation(sparse_sim)

print("Done! Total time for simulations:")
t1 <- Sys.time()
print(t1 - t0)

print(plot_eval(sparse_sim, "rare_prob_mse_gen"))







