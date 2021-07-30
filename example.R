#############################################################################################
##                                                                                         ##
## Supplementary Code                                                                      ##
##                                                                                         ##
## Cong Mu, Angelo Mele, Lingxin Hao, Joshua Cape, Avanti Athreya, and Carey E. Priebe     ##
##                                                                                         ##
#############################################################################################




#### Change the path or working directory if needed
source("algorithm.R")




#### Example 1 - 2-block rank one model with one 5-categorical covariate (Section 5.1)


## Hyperparameters
seed <- 2020
addCovariates <- TRUE
betahatmethod <- "SA"
dmax <- 10

## Number of blocks
K <- 2

## Dimension of latent positions
d <- 1    

## Number of vertices
n <- 1000 

## Latent positions and vertex covariate effect
p <- 0.3
q <- 0.4
latent <- cbind(p, q)
beta <- 0.15
cov <- 5

## Balanced case
pi <- rep(1/K, K)
pi_z <- rep(1/(K*cov), K*cov)

block_size <- round(pi*n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## Generate covariates
block_size_z <- round(pi_z*n)
covariates <- matrix(c(rep(1,block_size_z[1]), rep(2,block_size_z[2]), rep(3,block_size_z[3]), rep(4,block_size_z[4]), rep(5,block_size_z[5]), 
                       rep(1,block_size_z[6]), rep(2,block_size_z[7]), rep(3,block_size_z[8]), rep(4,block_size_z[9]), rep(5,block_size_z[10])))

## Generate adjacency matrix
B <- generateB(latent, K, d, addCovariates, cov, beta)
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
A <- generateA(n, P, seed = seed)

## Apply two algorithms
dhat <- rankMatrix(B)
G <- K * cov
result1 <- Algo1(A, cov, dhat, dmax, G) 
result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)

## Evalute the performance via ARI
ARI1 <- adjustedRandIndex(blocks, result1$tauhat)
ARI2 <- adjustedRandIndex(blocks, result2$tautilde)
data.frame(Algorithms = c("Algo1", "Algo2"), ARI = c(ARI1, ARI2))




#### Example 2 - 2-block homogeneous model with one 5-categorical covariate (Section 5.2)

## Hyperparameters
seed <- 2020
addCovariates <- TRUE
betahatmethod <- "SA"
dmax <- 10

## Number of blocks
K <- 2

## Dimension of latent positions
d <- 2  

## Number of vertices
n <- 1000 

## Latent positions and vertex covariate effect
a <- 0.3
b <- 0.25
latent <- cbind(c(sqrt(a), 0), c(b/sqrt(a), sqrt((a+b)*(a-b)/a)))
beta <- -0.2
cov <- 5

## Balanced case
pi <- rep(1/K, K)
pi_z <- rep(1/(K*cov), K*cov)

block_size <- round(pi*n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## Generate covariates
block_size_z <- round(pi_z*n)
covariates <- matrix(c(rep(1,block_size_z[1]), rep(2,block_size_z[2]), rep(3,block_size_z[3]), rep(4,block_size_z[4]), rep(5,block_size_z[5]), 
                       rep(1,block_size_z[6]), rep(2,block_size_z[7]), rep(3,block_size_z[8]), rep(4,block_size_z[9]), rep(5,block_size_z[10])))

## Generate adjacency matrix
B <- generateB(latent, K, d, addCovariates, cov, beta)
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
A <- generateA(n, P, seed = seed)

## Apply two algorithms
dhat <- rankMatrix(B)
G <- K * cov
result1 <- Algo1(A, cov, dhat, dmax, G) 
result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)

## Evalute the performance via ARI
ARI1 <- adjustedRandIndex(blocks, result1$tauhat)
ARI2 <- adjustedRandIndex(blocks, result2$tautilde)
data.frame(Algorithms = c("Algo1", "Algo2"), ARI = c(ARI1, ARI2))



