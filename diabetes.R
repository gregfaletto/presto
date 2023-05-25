# setwd("/Users/gregfaletto/Documents/GitHub/presto")

if(!is.null(dev.list())){
     dev.off()
}
rm(list=ls())

library(simulator)
library(MASS)
# MLDataR package contains data set
library(MLDataR)
library(parallel)
library(ggplot2)
library(cowplot)

# Initialize parallel processing--only works on Mac or Unix
n_cores <- 7
stopifnot(n_cores <= detectCores())
if(n_cores == 7){
     nsims <- 7
} else{
     nsims <- 49
}

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

# Load data

data(PreDiabetes)

PreDiabetes <- as.data.frame(PreDiabetes)

# Prepare data

print("Preparing data...")

n <- nrow(PreDiabetes)

age_cutoff <- 65

# A little over 8 hours for nsim = 50 across 7 cores on laptop
# About 45 minutes for nsim = 5 across 7 cores on laptop
# Run experiments

print("Running simulations...")

RESP_LEVELS <- c("NoDiabetes", "PreDiabetes", "Diabetes")

t0 <- Sys.time()

set.seed(123571)

sim <- new_simulation("prediabetes_data_app", "Data application") |>
     generate_model(gen_data_app_model_prediabetes, df=PreDiabetes,
          test_prop=0.1, resp_name="response", resp_levels=RESP_LEVELS,
          tune_prop=0, x_formula= ~ Age + Sex + IMD_Decile + BMI + HbA1C,
          omitted_x = NA, age_cutoff=as.list(5*(6:13)),
          vary_along="age_cutoff") |>
     simulate_from_model(nsim = nsims, index = 1:n_cores) |>
     run_method(list(prop_odds_data_analysis_vec, logit_meth_gen,
          fused_polr_data_analysis_vec), parallel=list(socket_names=n_cores))

sim <- sim |> evaluate(list(cal_osce_gen_data_app))

save_simulation(sim)

print("Done! Total time for simulations:")
t1 <- Sys.time()
print(t1 - t0)

plot_eval_by(sim, "cal_osce_gen_data_app", varying = "age_cutoff") +
     xlab("Age cutoff") + ggtitle(NULL) + ylab("Estimated Rare Probability MSE")

tabulate_eval(sim, "cal_osce_gen_data_app", se_format="None",
     format_args=list(digits=3))

