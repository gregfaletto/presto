# setwd("/Users/gregfaletto/OneDrive - USC Marshall School of Business/Proportional odds generalization/Code/Data Application")

if(!is.null(dev.list())){
     dev.off()
}
rm(list=ls())

library(simulator)
library(MASS)
library(ordinalNet)
# ORION package contains data set
library(ordinal)
# library(doParallel)
library(parallel)

# Initialize parallel processing--only works on Mac or Unix
n_cores <- detectCores() - 1
# registerDoParallel(cores = n_cores)

dir_main <- getwd()
dir_code <- "/Users/gregfaletto/OneDrive - USC Marshall School of Business/Proportional odds generalization/Code/ordinalNet modified/R"
# dir_data <- "/Users/gregfaletto/Library/CloudStorage/OneDrive-USCMarshallSchoolofBusiness/Proportional odds generalization/Possible Data Sets/Sales conversion optimization"
dir_sims <- "/Users/gregfaletto/OneDrive - USC Marshall School of Business/Proportional odds generalization/Code/Simulations"


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

setwd(dir_sims)

source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

# Load functions
source("data_app_funcs.R")

setwd(dir_main)

# Load data

data(soup)

# Prepare data

# Results on all outcomes (49 simulations):
#
# cal_osce_log_width_gen_data_app: PRESTO outperforms other methods. T-stat for
# beating proportional odds: 1.153897 T-stat for beating logit: 2.333834

# cal_osce_width_gen_data_app: PRESTO outperforms other methods. T-stat for
# beating proportional odds: 0.4663336. T-stat for beating logit: 3.047715.
#
# cal_osce_gen_data_app: PRESTO outperforms other methods. T-stat for beating
# proportional odds: 1.140587. T-stat for beating logit: 3.73311.
#
# cal_gen_data_app: PRESTO performs worse than other methods.
#
# cal_max_gen_data_app: PRESTO performs worse than other methods.
#
# cal_esce_gen_data_app: PRESTO performs worse than other methods.
#
# Rare event accuracy basically identical for all
# Rare probability Brier score worse for us






# Results on collapsed outcomes (875 simulations) (took 11.01346 hours),
# categories 2 - 5 merged together, features "RESP" and "PRODID" omitted: 
#
# cal_osce_log_width_gen_data_app: PRESTO outperforms other methods. T-stat for
# beating proportional odds: 1.454514. T-stat for beating logit: 4.945855.
#
# cal_osce_width_gen_data_app: PRESTO outperforms other methods. T-stat for
# beating proportional odds: 0.9012153. T-stat for beating logit: 10.02963.
#
# cal_osce_gen_data_app: PRESTO outperforms other methods. T-stat for beating
# proportional odds: 1.791337. T-stat for beating logit: 10.87275.
#
# cal_gen_data_app: PRESTO performs worse than other methods.
#
# cal_max_gen_data_app: PRESTO performs worse than other methods.
#
# cal_esce_gen_data_app: PRESTO performs worse than other methods.
#
# Rare event accuracy basically identical for all
# Rare probability Brier score worse for us

# To do:
# * PRODID omitted because not feasible to have observations at every level on
# every subsample, so sometimes coefficients are omitted altogether. (RESP has
# missing observations)
# (could also omit cold and easy as paper does)

print("Preparing data...")

resp_levels <- paste("sure", 1:6, sep="_")

soup$SURENESS <- factor(soup$SURENESS, ordered=TRUE)

levels(soup$SURENESS) <- resp_levels

# Reverse ordering since sure_1 is the rare category
resp_levels <- paste("sure", 6:1, sep="_")

qualities <- soup$SURENESS

# qualities <- character(nrow(soup))

# qualities[soup$SURENESS %in% c("sure_1")] <- "SureRef"
# qualities[soup$SURENESS %in% c("sure_2", "sure_3", "sure_4",
#      "sure_5")] <- "NotSure"
# qualities[soup$SURENESS %in% c("sure_6")] <- "SureNotRef"

# stopifnot(all(!is.na(qualities)))
# resp_levels <- c("SureNotRef", "NotSure", "SureRef")
# stopifnot(all(qualities %in% resp_levels))

# print("qualities:")
# print(summary(qualities))
# print(str(qualities))

soup$SURENESS <- factor(qualities, levels=resp_levels, ordered=TRUE)

stopifnot(all(levels(soup$SURENESS) == resp_levels))

stopifnot(all(!is.na(soup$SURENESS)))




print("Category probabilities:")
print(summary(soup$SURENESS)/length(soup$SURENESS))

print("n:")
print(nrow(soup))

# Took a little over 6 minutes per simulation with original categories, a little
# over 2 with shortened categories.

# Parallel: 125*7 simulations (took 11.01346 hours)
# 25*7 simulations took a little under 2.5 hours.

# All outcomes, parallel: 7*7 simulations took 2.230848 hours.

# Run experiments

print("Running simulations...")

t0 <- Sys.time()

set.seed(123571)

sim <- new_simulation("soup_data_app", "Data application") |>
     generate_model(gen_data_app_model, df=soup, test_prop=0.1,
          resp_name="SURENESS", resp_levels=resp_levels,
          tune_prop=0,
          x_formula= ~  PROD + DAY + SOUPTYPE + SOUPFREQ + COLD + EASY + GENDER + AGEGROUP + LOCATION ,
          omitted_x = c("RESP", "PRODID")) |>
     simulate_from_model(nsim = 50 #40
          , index = 1:n_cores
          ) |>
     run_method(list(prop_odds_data_analysis_vec
          , logit_meth_gen
          , fused_polr_data_analysis_vec
          # , semipar_ordnet_vec
          )
     , parallel=list(socket_names=n_cores)
     ) 

sim <- sim |> evaluate(list(
     # rare_prob_brier_data_app,
          # rare_prob_acc_data_app
          # , cal_osce_log_width_gen_data_app
          , cal_osce_gen_data_app
          , cal_osce_width_gen_data_app
          , prop_rare_obs))

save_simulation(sim)

print("Done! Total time for simulations:")
t1 <- Sys.time()
print(t1 - t0)


print(plot_eval(sim, "cal_osce_gen_data_app"))
print(plot_eval(sim, "cal_osce_width_gen_data_app"))
# print(plot_eval(sim, "rare_prob_acc_data_app"))
# print(plot_eval(sim, "rare_prob_brier_data_app"))
# print(plot_eval(sim, "cal_osce_log_width_gen_data_app"))

tabulate_eval(sim, "cal_osce_gen_data_app")


setwd("/Users/gregfaletto/Dropbox/Jacob and Greg/AISTATS 2023")
source("sim_eval_function.R")
setwd(dir_main)
create_data_app_plots(sim)



