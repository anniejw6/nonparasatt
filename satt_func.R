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

	# run genmatch
	z <- Matching::GenMatch(Tr=treat, X=x_small, 
		pop.size = pop.size,
		max.generations = max.generations, 
		wait.generations = wait.generations,
		cluster = cluster)

	# check balance
	mgen1 <- Match(Y = y, Tr = treat, X = x_small, Weight.matrix = z)

	# return matched sample
	list(trt_vars = trt_vars,
		out_vars = out_vars,
		matched = z$matches,
		bal = mgen1)

}

satt_est <- function(matched_sample){

	# run gam here ... 

	# return satt, confidence interval
}