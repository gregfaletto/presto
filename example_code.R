set.seed(4512)

rm(list=ls())

# Set your working directory to the folder containing presto_func.R
source("presto_func.R")
# Now the functions presto, which fits the model, and predictPrestoProbs, which
# generates class probabilties for unlabeled data, are loaded.

# Run the below line of code if you haven't already installed the ordinal
# package
# install.packages(“ordinal”)
library(ordinal)

# Load the soup data set
data(soup)

# Take a look at the data
print("Data set:")
str(soup)

print("Response (soup$SURENESS):")

y <- soup$SURENESS

print(str(y))

# Get the proportion of observations that lie in each class
print("Response class proportions:")
n <- length(y)
print(summary(y)/n)

dat <- model.matrix( ~ PROD + DAY + SOUPTYPE + SOUPFREQ + COLD + EASY +
	GENDER + AGEGROUP + LOCATION, soup)

X <- dat[, colnames(dat) != "(Intercept)"]

stopifnot(nrow(X) == n)
stopifnot(ncol(X) == 22)

# Fit model on rest of data using selected lambda and get parameter estimates
# (takes about 3.5 minutes on my laptop)

print("Estimating PRESTO...")
presto_results <- presto(X, y, printIter=TRUE)

print("Estimated coefficients:")
print(presto_results)

print("Estimated probabilities:")
probs <- predictPrestoProbs(presto_results, X)
print(probs)

