# Test script for the internal functions

# Load required libraries
library(MASS)
library(glmnet)
library(randomForest)
library(BiocParallel)

# Load the functions
source('utils.R')
source('JRF_internal.R')
source('ZILGM_internal.R')
source('PCzinb_internal.R')

# Test with small synthetic data
set.seed(123)
n <- 50
p <- 5

# Create simple count data
X <- matrix(rpois(n*p, lambda=3), nrow=n, ncol=p)
colnames(X) <- paste0('Gene', 1:p)
rownames(X) <- paste0('Cell', 1:n)

cat('Testing PCzinb with Poisson method...\n')
# Test PCzinb with poi method (should work)
adj_poi <- try(PCzinb_internal(X, method='poi', maxcard=1, alpha=0.05, extend=TRUE))
if (!inherits(adj_poi, 'try-error')) {
    cat('PCzinb Poisson method: SUCCESS\n')
    cat('Adjacency matrix dimensions:', dim(adj_poi), '\n')
} else {
    cat('PCzinb Poisson method: FAILED\n')
    print(adj_poi)
}

cat('Testing PCzinb with NB method...\n')
# Test PCzinb with nb method (should work)
adj_nb <- try(PCzinb_internal(X, method='nb', maxcard=1, alpha=0.05, extend=TRUE))
if (!inherits(adj_nb, 'try-error')) {
    cat('PCzinb NB method: SUCCESS\n')
    cat('Adjacency matrix dimensions:', dim(adj_nb), '\n')
} else {
    cat('PCzinb NB method: FAILED\n')
    print(adj_nb)
}

# Test JRF with multi-class data
cat('Testing JRF method...\n')
# Create multi-class data for JRF
X_list <- list(
    class1 = matrix(rpois(n*p, lambda=2), nrow=p, ncol=n),
    class2 = matrix(rpois(n*p, lambda=4), nrow=p, ncol=n)
)
rownames(X_list[[1]]) <- paste0('Gene', 1:p)
rownames(X_list[[2]]) <- paste0('Gene', 1:p)

adj_jrf <- try(JRF_internal(X_list, ntree=50, mtry=2, genes.name=paste0('Gene', 1:p), nCores=2))
if (!inherits(adj_jrf, 'try-error')) {
    cat('JRF method: SUCCESS\n')
    cat('Output columns:', colnames(adj_jrf), '\n')
} else {
    cat('JRF method: FAILED\n')
    print(adj_jrf)
}

# Test ZILGM method
cat('Testing ZILGM method...\n')
adj_zilgm <- try(zilgm_internal(t(X), family='Poisson', nlambda=10, nCores=2, verbose=0))
if (!inherits(adj_zilgm, 'try-error')) {
    cat('ZILGM method: SUCCESS\n')
    cat('Number of lambda values:', length(adj_zilgm$lambda), '\n')
} else {
    cat('ZILGM method: FAILED\n')
    print(adj_zilgm)
}

# All done - BiocParallel handles cleanup automatically