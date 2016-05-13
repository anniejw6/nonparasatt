match_samp <- function(x, treat, y,
	nfold = 10, parallel = F,
	alpha = 1,
	cv_cutoff = 0.5, 
	pop.size = 1000, max.generations = 100,
	wait.generations = 10,
	cluster = NULL){

	# x should be a numeric matrix
	# treat should be a numeric vector
	# y should be a numeric vector


	# run lasso to find covariates that predict treatment
	cat('Variables Predicting Treatment\n')
	trt_vars <- pick_glmnet_var(
		x = x, y = treat, 
		nfold = nfold, parallel = parallel,
		alpha = alpha, cv_cutoff = cv_cutoff)

	# run lasso to find covariates that predict outcome
	cat('Variables Predicting Outcome\n')
	out_vars <- pick_glmnet_var(
		x = x, y = y, 
		nfold = nfold, parallel = parallel,
		alpha = alpha, cv_cutoff = cv_cutoff)

	x_small <- x[, trt_vars, drop = F]

	# generate covariate weights using genmatch
	z <- Matching::GenMatch(Tr=treat, X=x_small, 
		pop.size = pop.size,
		max.generations = max.generations, 
		wait.generations = wait.generations,
		cluster = cluster)

	# run matching using genmatch weights
	mgen1 <- Match(Y = y, Tr = treat, X = x_small, Weight.matrix = z)
	
	# matched sample
	df <- data.frame(y = y,
	                 treat = treat,
	                 x = x)
	names(df) <- gsub("x\\.", "", names(df))  # Remove covariate prefix
	df <- df[unique(c(mgen1$index.treated, mgen1$index.control)), ]

	# return matched sample
	list(trt_vars = trt_vars,
		out_vars = out_vars,
		matched = z$matches,
		mout = mgen1,
		df = df)

}

check_balance <- function(matched, nboots = 500){
  
  # make formula
  formula <- paste("treat ~", paste0(matched$trt_vars, collapse = " + "))
  formula <- as.formula(formula)
  
  # get balance stats
  MatchBalance(formula, 
               data = matched$df, 
               match.out = matched$mout, 
               nboots = nboots)
}

satt_est <- function(matched){
  
  # matched: return object from match_samp()
  
  # variables selected: intersection of predictors of treat and outcome
  covariates <- unique(c(matched$trt_vars, matched$out_vars))
  xs <- matched$df[, covariates]
  
  # remove columns with zero variance
  covariates <- covariates[!(covariates %in% zero_variance(xs))]
  xs <- xs[, covariates]

  # classify covariates as factors or continuous
  factors <- get_factor_vars(xs, threshold = 10)
  continuous <- covariates[!(covariates %in% factors)]
  
  # flags
  no_factors <- length(factors) == 0
  no_continuous <- length(continuous) == 0
  no_vars <- length(factors) == 0 & length(continuous) == 0
  
  # convert factors
  if (length(factors) > 0) {
    matched$df[, factors] <- apply(matched$df[, factors], 2, 
                                   function(x) as.factor(x))
  }
  
  # create formula
  if (no_vars == T) {
    formula <- "y ~ treat"
  } else if (no_factors == T) {
    formula_p1 <- "y ~ treat"
    formula_p2 <- paste0("s(", paste0(continuous, collapse = ") + s("), ")")
    formula <- paste(formula_p1, formula_p2, sep = " + ")
  } else if (no_continuous == T) {
    formula_p1 <- "y ~ treat"
    formula_p2 <- paste0(factors, collapse = " + ")
    formula <- paste(formula_p1, formula_p2, sep = " + ")
  } else {
    formula_p1 <- "y ~ treat"
    formula_p2 <- paste0(factors, collapse = " + ")
    formula_p3 <- paste0("s(", paste0(continuous, collapse = ") + s("), ")")
    formula <- paste(formula_p1, formula_p2, formula_p3, sep = " + ")
  }
  formula <- as.formula(formula)
  
  # fit gam
  model <- mgcv::gam(formula,
                     data = matched$df,
                     family = ifelse(length(unique(matched$df$y)) == 2,
                                     "binomial",
                                     "gaussian"))
  
  # get satt and 9% confidence interval
  satt <- model$coefficients[["treat"]]
  se <- summary(model)$se[["treat"]]
  ci.95 <- satt + se * c(-1.96, 1.96)
  
	# return in requested format
  data.frame(est = satt,
             ci_lower = ci.95[1],
             ci_upper = ci.95[2])
}