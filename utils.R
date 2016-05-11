# Helper function to run LASSO
pick_glmnet_var <- function(x, y, 
	nfold = 10, parallel = F,
	alpha = 1,
	cv_cutoff = 0.5){

	model <- glmnet::cv.glmnet(x = x, y = y, 
		family = "binomial", 
		nfolds = nfold, 
		type.logistic = "modified.Newton", 
		parallel = parallel, 
		alpha = alpha)

  # Find optimal penalty parameter (lambda)
	min.cvm <- model$cvm[ which(model$lambda == model$lambda.min) ]
	if(cv_cutoff != 0){
		new.ixs <- which( model$cvm < (min.cvm + model$cvsd * cv_cutoff) )
	} else {
		new.ixs <- which( model$cvm == min.cvm )
	}

  # Subset to variables selected
	beta.selected <- model$glmnet.fit$beta[, new.ixs[1]]
	vars.selected <- which(abs(beta.selected) > 1e-5)
	vars_selected <- names(vars.selected)

	cat(sprintf('Variables Selected: %s\n', 
		paste(vars_selected, collapse = ', ')
		))
	return(vars_selected)
}