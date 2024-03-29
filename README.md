# presto
Code for ["Predicting Rare Events by Shrinking Towards Proportional Odds"](https://proceedings.mlr.press/v202/faletto23a.html) (Faletto and Bien, 2023). Citation:

Gregory Faletto and Jacob Bien (2023). Predicting Rare Events by Shrinking Towards Proportional Odds. Proceedings of the 40th International Conference on Machine Learning, in Proceedings of Machine Learning Research 202:9547-9602.

# Code implementing PRESTO

To get an implementation of PRESTO, save this repo, use `setwd()` to set your working directory in R to the directory for this repo, and load the file `presto_func.R`:

```
source("presto_func.R")
```
After this, two functions will be available for you in implementing PRESTO:
* The function `presto` fits the PRESTO model. Given a model matrix X and a response (ordered categorical variable) y, `presto` automatically selects a tuning parameter, then estimates the PRESTO model, returning the estimated coefficients.
* The function `predictPrestoProbs` estimates class probabilities on a test data set, given the output of `presto`.

Read the comments of `presto_func.R` for more documentation details. **See the file `example_code.R` for code illustrating how to use these functions:**

```
source("example_code.R")
```

# Details on data experiments from paper

In all data experiments, we implemented PRESTO by slightly modifying the code for `ordinalNet` (version 2.12); see Appendix H in the paper for details, and see the modified code in this repo. In all cases, executing the code in this repo generates the exact plots that appear in the published version of the paper.

## Computation details: synthetic data experiments

The synthetic data experiments from Sections 4.1 (`4_1_sparse_diffs.R`) and 4.2 (`4_2_dense_diffs.R`) of the paper, as well as Simulation Studies A (`e_sim_study_a.R`) and B (`e_sim_study_b.R`) in Appendix E, were conducted in R Version 4.2.2 running on macOS 10.15.7 on an iMac with a 3.5 GHz Quad-Core Intel Core i7 processor and 32 GB or RAM. We used the R packages `MASS` (version 7.3.58.1), `simulator` (version 0.2.4), `ggplot2` (version 3.3.6), `cowplot` (version 1.1.1), and `stargazer` (version 5.2.3), all available for download on CRAN, as well as the base `parallel` package (version 4.2.2).

## Computation details: real data experiments

The real data experiments from Section 4.3 (`4_3_soup.R`) and Appendix B (`diabetes.R`) of the paper were conducted in R Version 4.3.0 running on macOS Ventura 13.3.1 on a MacBook Pro with a 2.3 GHz Quad-Core Intel Core i5 processor and 16 GB or RAM. We used the R packages `MASS` (version 7.3.58.4), `simulator` (version 0.2.4), `ggplot2` (version 3.4.2), `cowplot` (version 1.1.1), and `stargazer` (version 5.2.3), all available for download on CRAN, as well as the base `parallel` package (version 4.3.0). The data for the real data experiment from Section 4.3 is the `soup` data set from the R `ordinal` package (version 2022.11.16). The data for the real data experiment from Appendix B is the `PreDiabetes` data set from the R `MLDataR` package (version 1.0.1).
