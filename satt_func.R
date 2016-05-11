match_samp <- function(x, treat, y,
	nfold = 10, parallel = F,
	alpha = 1,
	cv_cutoff = 0.5){

	# x should be a numeric matrix
	# treat should be a numeric vector
	# y should be a numeric vector


	# run lasso to find covariates that predict treatment
	vars_selected <- pick_glmnet_var(
		x = x, y = treat, 
		nfold = nfold, parallel = parallel,
		alpha = alpha, cv_cutoff = cv_cutoff)
	x_small <- x[, vars_selected, drop = F]

	# run genmatch
	genout <- Matching::GenMatch(Tr=treat, X=x_small)


	# return matched sample

}

satt_est <- function(matched_sample){

	# run gam here ... 

	# return satt, confidence interval
}