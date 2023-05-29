# presto
Code for ["Predicting Rare Events by Shrinking Towards Proportional Odds"](https://gregoryfaletto.com/research/) (Faletto and Bien, 2023). Forthcoming at Fortieth International Conference on Machine Learning ([ICML 2023](https://icml.cc/)).

In all data experiments, we implemented \textsc{PRESTO} by slightly modifying the code for \texttt{ordinalNet} (version 2.12); see Section H in the appendix of the paper for details, and see the modified code in this repo. In all cases, executing the code in this repo generates the exact plots that appear in the published version of the paper.

# Computation details: synthetic data experiments

The synthetic data experiments from Sections 4.1 ("4_1_sparse_diffs.R") and 4.2 ("4_2_dense_diffs.R") of the paper, as well as Simulation Studies A ("e_sim_study_a.R") and B ("e_sim_study_b.R") in Appendix E, were conducted in R Version 4.2.2 running on macOS 10.15.7 on an iMac with a 3.5 GHz Quad-Core Intel Core i7 processor and 32 GB or RAM. We used the R packages \texttt{MASS} (version 7.3.58.1), \texttt{simulator} (version 0.2.4), \texttt{ggplot2} (version 3.3.6), \texttt{cowplot} (version 1.1.1), and \texttt{stargazer} (version 5.2.3), all available for download on CRAN, as well as the base \texttt{parallel} package (version 4.2.2). We implemented \textsc{PRESTO} by slightly modifying the code for \texttt{ordinalNet} (version 2.12); see Section \ref{sec.est.presto} in the appendix for details, and see the modified code in this repo.

# Computation details: real data experiments

The synthetic data experiments from Sections 4.3 ("4_3_soup.R") and Appendix B ("diabetes.R") of the paper were conducted in R Version 4.3.0 running on macOS Ventura 13.3.1 on a MacBook Pro with a 2.3 GHz Quad-Core Intel Core i5 processor and 16 GB or RAM. We used the R packages \texttt{MASS} (version 7.3.58.4), \texttt{simulator} (version 0.2.4), \texttt{ggplot2} (version 3.4.2), \texttt{cowplot} (version 1.1.1), and \texttt{stargazer} (version 5.2.3), all available for download on CRAN, as well as the base \texttt{parallel} package (version 4.3.0). The data for the real data epxeriment from Section 4.3 is the \texttt{soup} data set from the R \texttt{ordinal} package (version 2022.11.16). The data for the real data experiment from Appendix B is the \texttt{PreDiabetes} data set from the R \texttt{MLDataR} package (version 1.0.1).
