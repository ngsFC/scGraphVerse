# Test ZINB methods specifically
library(MASS)
library(glmnet)
library(randomForest)
library(BiocParallel)

# Load the functions
source('utils.R')
source('JRF_internal.R')
source('ZILGM_internal.R')
source('PCzinb_internal.R')

# Create test data
set.seed(123)
n <- 20  # Small sample for testing
p <- 3   # Small number of genes
X <- matrix(rpois(n*p, lambda=3), nrow=n, ncol=p)
colnames(X) <- paste0('Gene', 1:p)
rownames(X) <- paste0('Cell', 1:n)

cat('Testing PCzinb with ZINB0 method...\n')
adj_zinb0 <- try(PCzinb_internal(X, method='zinb0', maxcard=1, alpha=0.1, extend=TRUE, nCores=2))
if (!inherits(adj_zinb0, 'try-error')) {
    cat('PCzinb ZINB0 method: SUCCESS\n')
    cat('Adjacency matrix dimensions:', dim(adj_zinb0), '\n')
} else {
    cat('PCzinb ZINB0 method: FAILED\n')
    print(adj_zinb0)
}

cat('Testing PCzinb with ZINB1 method...\n')
adj_zinb1 <- try(PCzinb_internal(X, method='zinb1', maxcard=1, alpha=0.1, extend=TRUE, nCores=2))
if (!inherits(adj_zinb1, 'try-error')) {
    cat('PCzinb ZINB1 method: SUCCESS\n')
    cat('Adjacency matrix dimensions:', dim(adj_zinb1), '\n')
} else {
    cat('PCzinb ZINB1 method: FAILED\n')
    print(adj_zinb1)
}