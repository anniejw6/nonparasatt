library(glmnet)
library(Matching)
library(parallel)
library(doParallel)
registerDoParallel(15)

wd <- '~/projects/nonparasatt/'
source(paste0(wd, 'utils.R'))
source(paste0(wd, 'satt_func.R'))

# Parallel cores
cl <- parallel::makeCluster(15)

# Clean covariates
xs <- read.csv(paste0(wd, 'data/x.csv'))
for(i in 1:length(xs)){
	if(length(unique(xs[[i]])) <= 7){
		print(i)
		xs[[i]] <- factor(xs[[i]])
	}
}
xs <- model.matrix(~., xs)
xs <- xs[, -1]

# Initialize empty set
satt_estimates <- data.frame(est = numeric(), 
                             ci_lower = numeric(), 
                             ci_upper = numeric())

# Populate satt_estimates
for (i in 1:20) {
  set.seed(2345)
  df <- read.csv(paste0(wd, 'data/zy_', as.character(i), '.csv'))
  ms <- match_samp(x = xs, 
                   treat = df$z, 
                   y = df$y,
                   cluster = cl, 
                   parallel = T, 
                   pop.size = 5)
  satt_estimates <- rbind(satt_estimates, satt_est(ms))
  write.csv(satt_estimates,
            file = paste0(wd, 'submission/est.csv'),
            row.names = F)
}