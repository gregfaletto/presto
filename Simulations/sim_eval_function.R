

create_plots <- function(sim){
  require(cowplot)
  require(ggplot2)

  rare_probs <- round(100*rare_probs(sim), 2)
  
  titles <- paste("Rare Proportion: ", rare_probs, "%", sep="")
  
  # main text plot
  
  boxplot_2 <- sim |> subset_simulation(subset=2) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[2]) + xlab(NULL) +
    scale_y_log10()
  
  e_df <- as.data.frame(evals(sim))
  stopifnot("rare_prob_mse_gen" %in% colnames(e_df))
  stopifnot("prop_rare_obs" %in% colnames(e_df))
  
  ratio_plot_1 <- sim |> subset_simulation(subset=1) |> evals() |>
    as.data.frame() |> plot_ratios(titles[1]) + xlab(NULL) + scale_y_log10()
  
  ratio_plot_2 <- sim |> subset_simulation(subset=2) |> evals() |>
    as.data.frame() |> plot_ratios(titles[2]) + xlab(NULL) + scale_y_log10()
  
  ratio_plot_3 <- sim |> subset_simulation(subset=3) |> evals() |>
    as.data.frame() |> plot_ratios(titles[3]) + xlab(NULL) + scale_y_log10()
  
  plot_1 <- cowplot::plot_grid(boxplot_2, ratio_plot_2, ratio_plot_3,
    ratio_plot_1, ncol = 2, nrow = 2)
  
  # supplement plot
  boxplot_1 <- sim |> subset_simulation(subset=1) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[1]) + xlab(NULL) +
    scale_y_log10()
  
  
  boxplot_3 <- sim |> subset_simulation(subset=3) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[3]) + xlab(NULL) +
    scale_y_log10()
  
  
  supp_plot <- cowplot::plot_grid(boxplot_1, boxplot_3, ncol = 2, nrow = 1)
  
  return(list(main_plot=plot_1, supp_plot=supp_plot))
}

plot_ratios <- function(df, title, loss_name="rare_prob_mse_gen",
                        ylab="MSE Ratio (PRESTO/Other)",
                        meth_names=c("fused_polr", "logit_meth",
                                     "prop_odds_meth")){
  require(ggplot2)
  presto_losses <- df[df$Method == meth_names[1], loss_name]
  logit_losses <- df[df$Method == meth_names[2], loss_name]
  po_losses <- df[df$Method == meth_names[3], loss_name]

  logit_ratio <- presto_losses/logit_losses
  po_ratio <- presto_losses/po_losses
  
  n <- length(logit_ratio)
  stopifnot(n == length(po_ratio))
  
  labels <- c(rep("Logit", n), rep("PO", n))

  df_gg <- data.frame(Labels=labels, MSE=c(logit_ratio, po_ratio))

  plot <- ggplot(df_gg, aes(x=Labels, y=MSE)) + geom_boxplot() +
    ggtitle(title) + geom_hline(yintercept=1, color="red",
                                      linetype="dashed") + ylab(ylab) +
    scale_y_log10()

  return(plot)
}

create_data_app_plots <- function(sim){
  require(cowplot)
  require(ggplot2)
  
  # sim <- sim |> evaluate(list(cal_osce_gen_data_app))
  
  boxplot <- sim |> plot_eval("cal_osce_gen_data_app") + ggtitle(NULL) +
    xlab(NULL) + ylab("Estimated Rare Probability MSE")
  
  e_df <- as.data.frame(evals(sim))
  stopifnot("cal_osce_gen_data_app" %in% colnames(e_df))
  
  ratio_plot <- sim |> evals() |> as.data.frame() |>
    plot_ratios(NULL, loss_name="cal_osce_gen_data_app",
                meth_names=c("fused_polr_data_analysis_vec", "logit_meth_gen",
                             "prop_odds_data_analysis_vec"),
                ylab="Estimated MSE Ratio (PRESTO/Other)") + xlab(NULL)
  
  plot_1 <- cowplot::plot_grid(boxplot, ratio_plot, ncol = 2)
  
  return(plot_1)
}

df_sim_stats <- function(sim, presto_meth="fused_polr",
  methods_to_compare=c("logit_meth", "prop_odds_meth")){
  
  # sim <- sim |> evaluate(list(rare_prob_mse_gen, prop_rare_obs))

  rare_probs <- round(100*rare_probs(sim), 2)
  rare_prob_titles <- paste(rare_probs, "%", sep="")

  e_df <- evals(sim) |> as.data.frame()

  metric <- "rare_prob_mse_gen"
  stopifnot(metric %in% colnames(e_df))
  models <- model(sim)

  model_names <- unique(e_df$Model)
  n_models <- length(model_names)
  
  stopifnot(length(rare_prob_titles) == n_models)
  
  # model_labels <- rep(as.character(NA), n_models)
  method_names <- unique(e_df$Method)

  n_comps <- length(methods_to_compare)

  stopifnot(presto_meth %in% method_names)
  stopifnot(all(methods_to_compare %in% method_names))
  stopifnot(!(presto_meth %in% methods_to_compare))
  
  for(i in 1:n_models){
    ind_i <- models[[i]]@name == model_names
    stopifnot(sum(ind_i) == 1)
    stopifnot(which(ind_i) == i)
    # model_labels[ind_i] <- models[[i]]@label
  }



  presto_meth_mse_means <- rep(as.numeric(NA), n_models)
    presto_meth_mse_ses <- rep(as.numeric(NA), n_models)

     n_comps <- length(methods_to_compare)

    comp_means <- matrix(as.numeric(NA), n_models, n_comps)
    comp_ses <- matrix(as.numeric(NA), n_models, n_comps)
    p_values <- matrix(as.character(NA), n_models, n_comps)

    colnames(p_values) <- methods_to_compare
    colnames(comp_means) <- methods_to_compare
    colnames(comp_ses) <- methods_to_compare


    # e_df_presto_meth <- e_df[e_df$Method == presto_meth, ]
    for(i in 1:n_models){
      e_df_i <- e_df[e_df$Model == model_names[i], ]

      presto_mses <- e_df_i[e_df_i$Method == presto_meth, metric]
      sample_size <- length(presto_mses)
      
      print("sample size:")
      print(sample_size)


      presto_meth_mse_means[i] <- mean(presto_mses)
      presto_meth_mse_ses[i] <- sd(presto_mses)/sqrt(sample_size)

      for(j in 1:n_comps){

        comp_mses <- e_df_i[e_df_i$Method == methods_to_compare[j], metric]
        stopifnot(length(comp_mses) == sample_size)
        

        comp_means[i, j] <- mean(comp_mses)
        comp_ses[i, j] <- sd(comp_mses)/sqrt(sample_size)

        p_values_i_j <- t.test(x=presto_mses, y=comp_mses, alternative="less",
                              paired=TRUE, var.equal=FALSE)$p.value

        p_values[i, j] <- signif(p_values_i_j, digits=3)

      } 
    }
  
  stopifnot(all(!is.na(rare_prob_titles)))

  t_df <- data.frame(rare_prob_titles, p_values)
  stopifnot(ncol(t_df) == n_comps + 1)
  stopifnot(nrow(t_df) == n_models)
  colnames(t_df) <- c("Rare Proportion", methods_to_compare)

  presto_meth_mse_means <- signif(presto_meth_mse_means, digits=3)
  presto_meth_mse_means <- formatC(presto_meth_mse_means, format="e", digits=2)
  presto_meth_mse_ses <- signif(presto_meth_mse_ses, digits=2)

  comp_means <- signif(comp_means, digits=3)
  comp_means <- formatC(comp_means, format="e", digits=2)
  comp_ses <- signif(comp_ses, digits=2)

  presto_stats <- paste(presto_meth_mse_means, " (", presto_meth_mse_ses, ")",
    sep="")

  print("str(comp_means):")
  print(str(comp_means))
  comp_stats <- matrix(character(), n_models, n_comps)
  for(j in 1:n_comps){
    comp_stats[, j] <- paste(comp_means[, j], " (", comp_ses[, j], ")", sep="")
  }


  print("str(comp_stats):")
  print(str(comp_stats))

  print("dim(comp_stats):")
  print(dim(comp_stats))

  stopifnot(nrow(comp_stats) == n_models)
  stopifnot(ncol(comp_stats) == n_comps)

  mean_se_df <- data.frame(rare_prob_titles, presto_stats, comp_stats)
  colnames(mean_se_df) <- c("Rare Proportion", "PRESTO", methods_to_compare)
  
  return(list(t_d_df=t_df, summary_df=mean_se_df))
    
}

rare_probs <- function(sim){
  rare_prob_df <- sim |> evals() |> as.data.frame()
  
  mod_names <- unique(rare_prob_df$Model)

  print("mod_names:")
  print(mod_names)
  
  n_models <- length(mod_names)
  
  stopifnot(n_models >= 1)
  
  rare_probs <- rep(as.numeric(NA), n_models)
  
  for(i in 1:n_models){
    rare_probs[i] <- mean(rare_prob_df[rare_prob_df$Model == mod_names[i],
                                       "prop_rare_obs"])
  }
  if(any(is.na(rare_probs))){
    print("total na:")
    print(sum(is.na(rare_probs)))
    print("total rare_probs:")
    print(length(rare_probs))
  }
  stopifnot(all(!is.na(rare_probs)))
  stopifnot(all(rare_probs >= 0))
  stopifnot(all(rare_probs <= 1))
  
  return(rare_probs)
}

df_data_app_stats <- function(sim){
  
  e_df <- evals(sim) |> as.data.frame()
  stopifnot("cal_osce_gen_data_app" %in% colnames(e_df))
  models <- model(sim)
  
  model_names <- unique(e_df$Model)
  n_models <- length(model_names)
  
  # model_labels <- rep(as.character(NA), n_models)
  method_names <- unique(e_df$Method)

  stopifnot(n_models == 1)
 
  presto_mses <- e_df[e_df$Method == "fused_polr_data_analysis_vec",
                        "cal_osce_gen_data_app"]
  sample_size <- length(presto_mses)
  
  print("sample size:")
  print(sample_size)
  
  logit_mses <- e_df[e_df$Method == "logit_meth_gen",
                       "cal_osce_gen_data_app"]
  stopifnot(length(logit_mses) == sample_size)
  
  po_mses <- e_df[e_df$Method == "prop_odds_data_analysis_vec",
                    "cal_osce_gen_data_app"]
  stopifnot(length(po_mses) == sample_size)
  # Summary statistics
  
  presto_mean <- mean(presto_mses)
  logit_mean <- mean(logit_mses)
  po_mean <- mean(po_mses)
  
  # Sample means and standard errors
  presto_mean <- signif(presto_mean, digits=3)
  
  logit_mean <- signif(logit_mean, digits=3)
  
  po_mean <- signif(po_mean, digits=3)

  mean_df <- c(presto_mean, logit_mean, po_mean)
  names(mean_df) <- c("PRESTO", "Logit", "PO")
  
  return(mean_df)
  
}
