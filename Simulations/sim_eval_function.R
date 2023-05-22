library(cowplot)
library(ggplot2)
# sim <- load_simulation("param_relax_prop_odds_unif_more2")
# sim <- param_relax_prop_odds_unif_more2
# e <- evals(sim)
# e_df <- as.data.frame(e)
# models <- model(sim)

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



        # e_df_presto_meth_j <- e_df_presto_meth[e_df_presto_meth$Model ==
        # model_names[i], ]


        # presto_mses <- e_df_i[e_df_i$Method == presto_meth, metric]
        # sample_size <- length(presto_mses)
        
        # print("sample size:")
        # print(sample_size)


        # presto_means[i] <- mean(presto_mses)
        # presto_ses[i] <- sd(presto_mses)/sqrt(sample_size)



        #   model_size_res <- e_df_presto_meth |> dplyr::slice(j)
        #   stopifnot(nrow(model_size_res) == n_sims)
        #   presto_meth_mses[j, ] <- model_size_res[, metric][[metric]]
        #   presto_meth_mse_means[j] <- mean(model_size_res[, metric][[metric]], na.rm=TRUE)
        #   presto_meth_mse_means[j] <- signif(presto_meth_mse_means[j], digits=3)
        #   sample_size_j <- sum(!is.na(model_size_res[, metric][[metric]]))
        #   if(sample_size_j > 0){
        #       presto_meth_mse_ses[j] <- sd(model_size_res[, metric][[metric]],
        #           na.rm=TRUE)/sqrt(sample_size_j)
        #       presto_meth_mse_ses[j] <- signif(presto_meth_mse_ses[j], digits=3)
        #   }

      }



      
        
    }

   









  
  # logit_ps <- rep(as.numeric(NA), n_models)
  # po_ps <- rep(as.numeric(NA), n_models)
  
  # logit_t_ps <- rep(as.numeric(NA), n_models)
  # po_t_ps <- rep(as.numeric(NA), n_models)
  
  # presto_means <- rep(as.numeric(NA), n_models)
  # presto_ses <- rep(as.numeric(NA), n_models)
  
  # logit_means <- rep(as.numeric(NA), n_models)
  # logit_ses <- rep(as.numeric(NA), n_models)
  
  # po_means <- rep(as.numeric(NA), n_models)
  # po_ses <- rep(as.numeric(NA), n_models)


  
  # for(i in 1:n_models){
  #   e_df_i <- e_df[e_df$Model == model_names[i], ]
    
  #   presto_mses <- e_df_i[e_df_i$Method == "fused_polr", metric]
  #   sample_size <- length(presto_mses)
    
  #   print("sample size:")
  #   print(sample_size)


  #   presto_means[i] <- mean(presto_mses)
  #   presto_ses[i] <- sd(presto_mses)/sqrt(sample_size)

  #   for(j in 1:n_comps){
  #     logit_mses <- e_df_i[e_df_i$Method == "logit_meth", metric]
  #     stopifnot(length(logit_mses) == sample_size)
      
  #     po_mses <- e_df_i[e_df_i$Method == "prop_odds_meth", metric]
  #     stopifnot(length(po_mses) == sample_size)
      

  #     # Logit t
  #     logit_t_ps[i] <- t.test(x=presto_mses, y=logit_mses, alternative="less",
  #                             paired=TRUE, var.equal=FALSE)$p.value

  #     # Proportional odds t
  #     po_t_ps[i] <- t.test(x=presto_mses, y=po_mses, alternative="less",
  #                          paired=TRUE, var.equal=FALSE)$p.value
      
  #     # Summary statistics
      
      
        
  #     logit_means[i] <- mean(logit_mses)
  #     logit_ses[i] <- sd(logit_mses)/sqrt(sample_size)
      
  #     po_means[i] <- mean(po_mses)
  #     po_ses[i] <- sd(po_mses)/sqrt(sample_size)
  #   }
    
    

  # }
  
  stopifnot(all(!is.na(rare_prob_titles)))
  
  # logit_t_ps <- signif(logit_t_ps, digits=3)
  # logit_t_ps <- as.character(formatC(logit_t_ps, format="e", digits=2))
  # po_t_ps <- signif(po_t_ps, digits=3)
  # po_t_ps <- as.character(formatC(po_t_ps, format="e", digits=2))

  t_df <- data.frame(rare_prob_titles, p_values)
  stopifnot(ncol(t_df) == n_comps + 1)
  stopifnot(nrow(t_df) == n_models)
  colnames(t_df) <- c("Rare Proportion", methods_to_compare)
  
  # t_df <- data.frame(rare_prob_titles, logit_t_ps, po_t_ps)
  # colnames(t_df) <- c("Rare Proportion", methods_to_compare)
  
  # Sample means and standard errors

  presto_meth_mse_means <- signif(presto_meth_mse_means, digits=3)
  presto_meth_mse_means <- formatC(presto_meth_mse_means, format="e", digits=2)
  presto_meth_mse_ses <- signif(presto_meth_mse_ses, digits=2)

  comp_means <- signif(comp_means, digits=3)
  comp_means <- formatC(comp_means, format="e", digits=2)
  comp_ses <- signif(comp_ses, digits=2)




  # presto_means <- signif(presto_means, digits=3)
  # presto_means <- formatC(presto_means, format="e", digits=2)
  # presto_ses <- signif(presto_ses, digits=2)
  
  # logit_means <- signif(logit_means, digits=3)
  # logit_means <- formatC(logit_means, format="e", digits=2)
  # logit_ses <- signif(logit_ses, digits=2)
  
  # po_means <- signif(po_means, digits=3)
  # po_means <- formatC(po_means, format="e", digits=2)
  # po_ses <- signif(po_ses, digits=2)

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

  
  # presto_stats <- paste(presto_means, " (", presto_ses, ")", sep="")
  # logit_stats <- paste(logit_means, " (", logit_ses, ")", sep="")
  # po_stats <- paste(po_means, " (", po_ses, ")", sep="")
  
  # mean_se_df <- data.frame(rare_prob_titles, presto_stats, logit_stats,
  #                          po_stats)

  mean_se_df <- data.frame(rare_prob_titles, presto_stats, comp_stats)
  colnames(mean_se_df) <- c("Rare Proportion", "PRESTO", methods_to_compare)
  
  return(list(t_d_df=t_df, summary_df=mean_se_df))
    
}

# print_sim_stats <- function(e_df, models, rare_props, dist){
# 
#   model_names <- unique(e_df$Model)
#   n_models <- length(model_names)
#   model_labels <- rep(as.character(NA), n_models)
#   method_names <- unique(e_df$Method)
#   for(i in 1:n_models){
#     ind_i <- models[[i]]@name == model_names
#     stopifnot(sum(ind_i) == 1)
#     stopifnot(which(ind_i) == i)
#     model_labels[ind_i] <- models[[i]]@label
#   }
#   
#   # Pairwise t_statistics
#   
#   for(i in 1:n_models){
#     e_df_i <- e_df[e_df$Model == model_names[i], ]
#     
#     # Logit t statistic
#     stopifnot("fused_polr" %in% method_names)
#     stopifnot("logit_meth" %in% method_names)
    # logit_results <- t.test(x=e_df_i[e_df_i$Method == "fused_polr", "rare_prob_mse_gen"],
    #                         y=e_df_i[e_df_i$Method == "logit_meth", "rare_prob_mse_gen"],
    #                         alternative="less", paired=TRUE, var.equal=FALSE)
#     
#     print("")
#     print("")
#     print("")
#     print("Model:")
#     print(model_labels[i])
#     print("T-statistic for logistic regression vs. PRESTO (one-tailed test of paired differences):")
#     print(logit_results$statistic)
#     print("p-value:")
#     print(logit_results$p.value)
#     
#     # Logit Wilcoxon
#     stopifnot("logit_meth" %in% method_names)
#     logit_results_wilcox <- wilcox.test(x=e_df_i[e_df_i$Method == "fused_polr", "rare_prob_mse_gen"],
#                                      y=e_df_i[e_df_i$Method == "logit_meth", "rare_prob_mse_gen"],
#                                      alternative="less", paired=TRUE)
#     
#     print("")
#     # print("T-statistic for PO vs. PRESTO (one-tailed test of paired differences):")
#     # print(po_results$statistic)
#     print("p-value for Wilcoxon test (logistic regression):")
#     print(logit_results_wilcox$p.value)
#     
#     
#     
#     
#     
#     # Proportional odds t statistic
#     stopifnot("prop_odds_meth" %in% method_names)
#     po_results <- t.test(x=e_df_i[e_df_i$Method == "fused_polr", "rare_prob_mse_gen"],
#                             y=e_df_i[e_df_i$Method == "prop_odds_meth", "rare_prob_mse_gen"],
#                             alternative="less", paired=TRUE, var.equal=FALSE)
#     
#     
#     print("")
#     print("T-statistic for PO vs. PRESTO (one-tailed test of paired differences):")
#     print(po_results$statistic)
#     print("p-value:")
#     print(po_results$p.value)
#     
#     # Proportional odds Wilcoxon
#     stopifnot("prop_odds_meth" %in% method_names)
#     po_results_wilcox <- wilcox.test(x=e_df_i[e_df_i$Method == "fused_polr", "rare_prob_mse_gen"],
#                          y=e_df_i[e_df_i$Method == "prop_odds_meth", "rare_prob_mse_gen"],
#                          alternative="less", paired=TRUE)
#     
#     print("")
#     # print("T-statistic for PO vs. PRESTO (one-tailed test of paired differences):")
#     # print(po_results$statistic)
#     print("p-value for Wilcoxon test (proportional odds):")
#     print(po_results_wilcox$p.value)
# 
#     ratio_plot_i <- plot_ratios(e_df_i, paste(dist,
#                                               " Differences (Rare Proportion: ",
#                                               rare_props[i], ")", sep=""))
#     print(ratio_plot_i)
#     
#     diff_plot_i <- plot_diffs(e_df_i, paste(dist,
#                                             " Differences (Rare Proportion: ",
#                                             rare_props[i], ")", sep=""))
#     print(diff_plot_i)
#   }
# }

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

# plot_diffs <- function(df, title){
#   require(ggplot2)
#   presto_losses <- df[df$Method == "fused_polr", "rare_prob_mse_gen"]
#   logit_losses <- df[df$Method == "logit_meth", "rare_prob_mse_gen"]
#   po_losses <- df[df$Method == "prop_odds_meth", "rare_prob_mse_gen"]
#   
#   logit_ratio <- presto_losses - logit_losses
#   po_ratio <- presto_losses - po_losses
#   labels <- c(rep("Logit", nrow(df)), rep("PO", nrow(df)))
#   
#   df_gg <- data.frame(Labels=labels, MSE=c(logit_ratio, po_ratio))
#   
#   plot <- ggplot(df_gg, aes(x=Labels, y=MSE)) + geom_boxplot() +
#     ggtitle(title) + geom_hline(yintercept=0, color="red",
#                                       linetype="dashed") +
#     ylab("MSE Difference (PRESTO - Other)")
#   
#   return(plot)
# }

# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

# Make four plots for this simulation: one boxplot of main results, then three
# plots of pairwise ratios

rare_probs <- function(sim){
  rare_prob_df <- sim |> evals() |> as.data.frame()
  
  mod_names <- unique(rare_prob_df$Model)

  print("mod_names:")
  print(mod_names)
  
  n_models <- length(mod_names)
  
  stopifnot(n_models >= 1)
  
  rare_probs <- rep(as.numeric(NA), n_models)
  
  # print("unique(rare_prob_df$Model):")
  # print(unique(rare_prob_df$Model))
  # print("colnames(rare_prob_df):")
  # print(colnames(rare_prob_df))
  for(i in 1:n_models){
    # print("mod_names[i]:")
    # print(mod_names[i])
    # print("mod_names[i] %in% rare_prob_df$Model:")
    # print(mod_names[i] %in% rare_prob_df$Model)
    # print("sum(rare_prob_df$Model == mod_names[i]):")
    # print(sum(rare_prob_df$Model == mod_names[i]))
    # print("str(df):")
    # print(str(rare_prob_df[rare_prob_df$Model == mod_names[i],
                                       # "prop_rare_obs"]))
    # rare_probs[i] <- mean(rare_prob_df[rare_prob_df$Model == mod_names[i],
    #                                    "prop_rare_obs"])
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

create_plots <- function(sim){
  require(cowplot)
  require(ggplot2)
  
  # sim <- sim |> evaluate(list(rare_prob_mse_gen, prop_rare_obs))
  
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
    as.data.frame() |> plot_ratios(titles[1]) + xlab(NULL)
  
  ratio_plot_2 <- sim |> subset_simulation(subset=2) |> evals() |>
    as.data.frame() |> plot_ratios(titles[2]) + xlab(NULL)
  
  ratio_plot_3 <- sim |> subset_simulation(subset=3) |> evals() |>
    as.data.frame() |> plot_ratios(titles[3]) + xlab(NULL)
  
  plot_1 <- plot_grid(boxplot_2, ratio_plot_2, ratio_plot_3, ratio_plot_1,
                      ncol = 2, nrow = 2)
  
  # supplement plot
  
  boxplot_1 <- sim |> subset_simulation(subset=1) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[1]) + xlab(NULL) +
    scale_y_log10()
  
  
  boxplot_3 <- sim |> subset_simulation(subset=3) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[3]) + xlab(NULL) +
    scale_y_log10()
  
  
  supp_plot <- plot_grid(boxplot_1, boxplot_3, ncol = 2, nrow = 1)
  
  return(list(main_plot=plot_1, supp_plot=supp_plot))
}

create_plot2 <- function(sim){
  require(cowplot)
  require(ggplot2)
  
  # sim <- sim |> evaluate(list(rare_prob_mse_gen, prop_rare_obs))
  
  rare_probs <- round(100*rare_probs(sim), 3)
  
  titles <- paste("Rare Proportion: ", rare_probs, "%", sep="")

  print("titles:")
  print(titles)
  
  # main text plot
  
  boxplot_2 <- sim |> subset_simulation(subset=2) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[2]) + xlab(NULL) +
    scale_y_log10()
  
  e_df <- as.data.frame(evals(sim))
  stopifnot("rare_prob_mse_gen" %in% colnames(e_df))
  stopifnot("prop_rare_obs" %in% colnames(e_df))
  
  ratio_plot_1 <- sim |> subset_simulation(subset=1) |> evals() |>
    as.data.frame() |> plot_ratios(titles[1]) + xlab(NULL)
  
  ratio_plot_2 <- sim |> subset_simulation(subset=2) |> evals() |>
    as.data.frame() |> plot_ratios(titles[2]) + xlab(NULL)
  

   boxplot_1 <- sim |> subset_simulation(subset=1) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[1]) + xlab(NULL) +
    scale_y_log10()
  # ratio_plot_3 <- sim |> subset_simulation(subset=3) |> evals() |>
  #   as.data.frame() |> plot_ratios(titles[3]) + xlab(NULL)
  
  plot_1 <- plot_grid(boxplot_2, boxplot_1, ratio_plot_2, ratio_plot_1,
                      ncol = 2, nrow = 2)
  
  # supplement plot
  
 
  
  
  # boxplot_3 <- sim |> subset_simulation(subset=3) |>
  #   plot_eval("rare_prob_mse_gen") + ggtitle(titles[3]) + xlab(NULL) +
  #   scale_y_log10()
  
  
  # supp_plot <- plot_grid(boxplot_1, boxplot_3, ncol = 2, nrow = 1)
  
  return(plot_1)
}

create_sparse_plot2 <- function(sim, plots=c(1,3)){
  require(cowplot)
  require(ggplot2)

  n_plots <- length(plots)

  stopifnot(n_plots %in% c(1, 2))
  
  # sim <- sim |> evaluate(list(rare_prob_mse_gen, prop_rare_obs))
  
  rare_probs <- round(100*rare_probs(sim), 3)
  
  titles <- paste("Rare Proportion: ", rare_probs, "%", sep="")

  print("titles:")
  print(titles)
  
  # main text plot
  
  if(n_plots == 2){
    boxplot_2 <- sim |> subset_simulation(subset=plots[2]) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[plots[2]]) + xlab(NULL) +
    scale_y_log10()
  }
  
  
  e_df <- as.data.frame(evals(sim))
  stopifnot("rare_prob_mse_gen" %in% colnames(e_df))
  stopifnot("prop_rare_obs" %in% colnames(e_df))
  
  ratio_plot_1 <- sim |> subset_simulation(subset=plots[1]) |> evals() |>
    as.data.frame() |> plot_ratios(titles[plots[1]]) + xlab(NULL)
  
  if(n_plots == 2){
    ratio_plot_2 <- sim |> subset_simulation(subset=plots[2]) |> evals() |>
    as.data.frame() |> plot_ratios(titles[plots[2]]) + xlab(NULL)
  }
  
  
  boxplot_1 <- sim |> subset_simulation(subset=plots[1]) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[plots[1]]) + xlab(NULL) +
    scale_y_log10()
  
  if(n_plots == 2){
    plot_1 <- plot_grid(boxplot_2, boxplot_1, ratio_plot_2, ratio_plot_1,
                      ncol = 2, nrow = 2)
    } else if(n_plots == 1){
      plot_1 <- plot_grid(boxplot_1,ratio_plot_1,
                      ncol = 2, nrow = 1)
    }
  

  return(plot_1)
}

create_sparse_plots <- function(sim, plots=c(1,3,5)){
  require(cowplot)
  require(ggplot2)

  stopifnot(length(plots) == 3)
  
  # sim <- sim |> evaluate(list(rare_prob_mse_gen, prop_rare_obs))
  
  rare_probs <- round(100*rare_probs(sim), plots[2])
  
  titles <- paste("Rare Proportion: ", rare_probs, "%", sep="")
  
  # main text plot
  
  boxplot_2 <- sim |> subset_simulation(subset=plots[2]) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[plots[2]]) + xlab(NULL) +
    scale_y_log10()
  
  e_df <- as.data.frame(evals(sim))
  stopifnot("rare_prob_mse_gen" %in% colnames(e_df))
  stopifnot("prop_rare_obs" %in% colnames(e_df))
  
  ratio_plot_1 <- sim |> subset_simulation(subset=plots[1]) |> evals() |>
    as.data.frame() |> plot_ratios(titles[plots[1]]) + xlab(NULL)
  
  ratio_plot_2 <- sim |> subset_simulation(subset=plots[2]) |> evals() |>
    as.data.frame() |> plot_ratios(titles[plots[2]]) + xlab(NULL)
  
  ratio_plot_3 <- sim |> subset_simulation(subset=plots[3]) |> evals() |>
    as.data.frame() |> plot_ratios(titles[plots[3]]) + xlab(NULL)
  
  plot_1 <- plot_grid(boxplot_2, ratio_plot_2, ratio_plot_3, ratio_plot_1,
                      ncol = 2, nrow = 2)
  
  # supplement plot
  
  boxplot_1 <- sim |> subset_simulation(subset=plots[1]) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[plots[1]]) + xlab(NULL) +
    scale_y_log10()
  
  
  boxplot_3 <- sim |> subset_simulation(subset=plots[3]) |>
    plot_eval("rare_prob_mse_gen") + ggtitle(titles[plots[3]]) + xlab(NULL) +
    scale_y_log10()
  
  
  supp_plot <- plot_grid(boxplot_1, boxplot_3, ncol = 2, nrow = 1)
  
  return(list(main_plot=plot_1, supp_plot=supp_plot))
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
  
  plot_1 <- plot_grid(boxplot, ratio_plot, ncol = 2)
  
  return(plot_1)
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
  # presto_mean <- formatC(presto_mean, format="e", digits=2)
  
  logit_mean <- signif(logit_mean, digits=3)
  # logit_mean <- formatC(logit_mean, format="e", digits=2)
  
  po_mean <- signif(po_mean, digits=3)
  # po_mean <- formatC(po_mean, format="e", digits=2)
  
  # presto_stats <- paste(presto_means, " (", presto_ses, ")", sep="")
  # logit_stats <- paste(logit_means, " (", logit_ses, ")", sep="")
  # po_stats <- paste(po_means, " (", po_ses, ")", sep="")
  
  mean_df <- c(presto_mean, logit_mean, po_mean)
  names(mean_df) <- c("PRESTO", "Logit", "PO")
  
  return(mean_df)
  
}


# sim <- sim |> evaluate(list(rare_prob_mse_gen, prop_rare_obs))
# dfs <- df_sim_stats(sim)
# plots <- create_plots(sim)
# stargazer(dfs[[3]], summary=FALSE)
