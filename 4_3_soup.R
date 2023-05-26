# setwd("/Users/gregfaletto/Documents/GitHub/presto")

if(!is.null(dev.list())){
     dev.off()
}
rm(list=ls())

library(simulator)
library(MASS)
# ordinal package contains data set
library(ordinal)
library(parallel)
library(cowplot)
library(ggplot2)

# Initialize parallel processing--only works on Mac or Unix.
n_cores <- 7
stopifnot(n_cores <= detectCores())
if(n_cores == 7){
     nsims <- 50
} else{
     nsims <- 350
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

data(soup)

print("Preparing data...")

resp_levels <- paste("sure", 1:6, sep="_")

soup$SURENESS <- factor(soup$SURENESS, ordered=TRUE)

levels(soup$SURENESS) <- resp_levels

# Reverse ordering since sure_1 is the rare category
resp_levels <- paste("sure", 6:1, sep="_")

qualities <- soup$SURENESS

soup$SURENESS <- factor(qualities, levels=resp_levels, ordered=TRUE)

stopifnot(all(levels(soup$SURENESS) == resp_levels))
stopifnot(all(!is.na(soup$SURENESS)))

print("Category probabilities:")
print(summary(soup$SURENESS)/length(soup$SURENESS))

print("n:")
print(nrow(soup))

print("Running simulations...")

t0 <- Sys.time()

set.seed(123571)

sim <- new_simulation("soup_data_app", "Data application") |>
     generate_model(gen_data_app_model, df=soup, test_prop=0.1,
          resp_name="SURENESS", resp_levels=resp_levels, tune_prop=0,
          x_formula= ~  PROD + DAY + SOUPTYPE + SOUPFREQ + COLD + EASY +
          GENDER + AGEGROUP + LOCATION, omitted_x = c("RESP", "PRODID")) |>
     simulate_from_model(nsim = nsims, index = 1:n_cores) |>
     run_method(list(prop_odds_data_analysis_vec, logit_meth_gen,
          fused_polr_data_analysis_vec), parallel=list(socket_names=n_cores)) 

sim <- sim |> evaluate(list(cal_osce_gen_data_app))

save_simulation(sim)

print("Done! Total time for simulations:")
t1 <- Sys.time()
print(t1 - t0)

# Create plot for paper (Figure 3)
fig_3 <- create_data_app_plots(sim)

# Generate sample means for paper (as reported in Section 4.3)
df_data_app_stats(sim)

