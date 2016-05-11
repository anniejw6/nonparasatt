# Setup
library(glmnet)
library(Matching)
library(parallel)

wd <- '~/projects/nonparasatt/'
source(paste0(wd, 'utils.R'))
source(paste0(wd, 'satt_func.R'))

# Parallel cores
cl <- parallel::makeCluster(detectCores() - 1)

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

# Read in dataset 1
df1 <- read.csv(paste0(wd, 'data/zy_1.csv'))

z <- match_samp(x = xs, treat = df1$z, y = df1$y,
	cluster = cl, parallel = T)