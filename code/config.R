library(ggplot2)
library(data.table)
library(magrittr)

cbPalette <- c("#0072B2", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")


bootstrap <- function(y_true, y_pred, n=999){
  # resample with replacement
}

## -------------------------
# Variance explained
# 1 - SS_res / SS_tot
# The value can be negative, because the model can be arbitrary worse.
# sklearn implementation: https://github.com/scikit-learn/scikit-learn/blob/55bf5d93e5674f13a1134d93a11fd0cd11aabcd1/sklearn/metrics/regression.py#L355
## -------------------------
variance_explained <- function(y_true, y_pred){
  y_diff_avg = mean(y_true - y_pred)
  numerator = mean((y_true - y_pred - y_diff_avg) ^ 2)
  y_true_avg = mean(y_true)
  denominator = mean((y_true - y_true_avg) ^ 2)
  output_score = 1 - numerator/denominator
  return(output_score)
}

R_2 <- function(y_true, y_pred){
  cor(y_true, y_pred) ^ 2
}


# metric <- variance_explained
metric <- R_2


## With affine transformation, R^2 is different from variance explained
n <- 100
x <- rnorm(n)

# affine transform + noise
y <- 2*x + 10 + rnorm(n, 0, 0.5)

plot(x, y)

abline(coef = c(10, 2))

variance_explained(y, x)
R_2(y, x)
