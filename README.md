# non-parametric SATT estimation

This is an entry (with [Leo Liu](https://github.com/leonidliu)) in the [causal inference challenge](https://docs.google.com/document/d/1p5xdeJVY5GdBC2ar_3wVjaboph0PemXulnMD5OojOCI/edit) for the 
2016 atlantic causal inference conference.

Our method is intended to calculate the SATT in an automatic and
nonparametric fashion.

For each data set, we do the following:

Step 1: We fit two ten-fold cross-validated lasso regressions,
the first predicting assignment to the treatment using all available
covariates, the second predicting the outcome using all available
covariates (excluding the treatment). We store a vector of the
selected variables for each.

Step 2: Using the covariates that predict treatment, we run
Sekhon's GenMatch algorithm to find 1 analagous control observation
that's as similar as possible to each treatment observation. This
creates a matched sample in which the probability of assignment
to treatment can be thought of as random, and in which
the covariate distributions of treatment and control groups
are balanced across variables that predict treatment assignment.

Step 3: Fit a generalized additive model of the outcome on the
treatment and on variables selected by either of the two lasso regressions.
Variables that have fewer than ten unique values are treated as
categorical and estimated with a series of dummies; otherwise,
variables are assumed to be continuous and fitted using a smoothed
spline with default values.

Step 4: We store the coefficient on the treatment variable as the SATT.
Assuming that the SATT is drawn from a normal distribution, we
calculate the 95% confidence interval as SATT +/- (1.96 * SE).
