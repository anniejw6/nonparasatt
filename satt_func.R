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
	df <- df[, unique(c(mgen1$index.treated, mgen1$index.control))]

	# return matched sample
	list(trt_vars = trt_vars,
		out_vars = out_vars,
		matched = z$matches,
		mout = mgen1,
		df = df)

}

check_balance <- function(matched_sample, nboots = 500){
  
  # make formula
  formula <- paste("treat ~", paste0(matched_sample$trt_vars, collapse = " + "))
  formula <- as.formula(formula)
  
  # get balance stats
  MatchBalance(formula, 
               data = matched_sample$df, 
               match.out = matched_sample$mout, 
               nboots = nboots)
}

satt_est <- function(matched_sample){
  
  # classify covariates as factors or continuous
  covariates <- unique(c(matched_sample$trt_vars, matched_sample$out_vars))
  factors <- get_factor_vars(covariates)
  continuous <- covariates[!(factors %in% covariates)]
  
  # create formula
  formula_p1 <- "y ~ treat +"
  formula_p2 <- paste0(factors, collapse = " + ")
  formula_p3 <- paste0("s(", paste0(continuous, collapse = ") + s("), ")")
  formula <- paste(formula_p1, formula_p2, formula_p3)
  
  # fit gam
  model <- mgcv::gam(formula,
                     data = matched_sample$df,
                     family = ifelse(length(unique(matched_sample$df$y)) == 2,
                                     "binomial",
                                     "gaussian"))
  
  # get satt and se
  
	# return satt, confidence interval
  list(satt_nocontrol = matched_sample$mout$est,
       se_nocontrol = matched_sample$mout$se)
}