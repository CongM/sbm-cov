##################################################
## codeFunc_approx_chernoff_optimize_v2.R ########
##################################################
##################################################
##
## [ ---------- VERSION 2 ---------- ]
##
## README:
## 
## The master function in this code is: approx_chernoff_opt_v2(matrix.B, vec.pi).
## The function outputs a vector of length seven with entries given by
## rho^star_ASE   argmin_k_ASE   argmin_l_ASE   rho^star_LSE   argmin_k_LSE   argmin_l_LSE rho^star_ratio
## for rho^star squantities as in our Network Science ("Chernoff") paper.
## Note that entries two and three form a tuple of (block) indicies corresponding to minimum Chernoff under ASE,
## and similarly, entries five and six form a tuple of (block) indicies corresponding to minimum Chernoff under LSE
##
## Example code with output:
##
##    > source(".../codeFunc_approx_chernoff_optimize.R")
##    > matrix.B <- rbind(c(0.5, 0.4, 0.3), c(0.4, 0.5, 0.3), c(0.3, 0.3, 0.5))
##    > vec.pi <- c(3/10, 5/10, 2/10)
##    > approx_chernoff_opt_v2(matrix.B, vec.pi)
##    rho^star_ASE   argmin_k_ASE   argmin_l_ASE   rho^star_LSE   argmin_k_LSE   argmin_l_LSE rho^star_ratio 
##    0.008163318    1.000000000    2.000000000    0.008152563    1.000000000    2.000000000    1.001319270  


###########################################################
## [BEGIN] Required libraries #############################
###########################################################

require(MASS)
require(Matrix)
require(methods)
require(stats)
require(graphics)
require(grDevices)
require(stringr)
require(utils)

###########################################################
## [BEGIN] Master function: approx_chernoff_opt ###########
###########################################################
## This function computes
## rho^star_ASE, rho^star_LSE, and rho^star_ratio
## as in the Chernoff paper 
###########################################################

approx_chernoff_opt_v2 <- function(matrix.B, vec.pi){
  
  LB <- eigen(matrix.B)
  index.keep <- which(abs(LB$values) > sqrt(.Machine$double.eps))
  vals.B <- LB$values[index.keep]
  vecs.B <- LB$vectors[,index.keep]
  matrix.lp <- vecs.B %*% diag(sqrt(abs(vals.B)), nrow=length(vals.B))
  dim.K <- dim(matrix.lp)[1]
  dim.lp <- dim(matrix.lp)[2]
  Ipq = diag(sign(vals.B))
  
  vec.mu <- colSums(diag(vec.pi) %*% matrix.lp)
  
  delta <- matrix(0, nrow = dim.lp, ncol = dim.lp)
  for(i in 1:dim.K){ delta = delta + vec.pi[i] * (matrix.lp[i,] %o% matrix.lp[i,])}
  delta.inv <- solve(delta)
  
  delta.tilde <- matrix(0, nrow = dim.lp, ncol = dim.lp)
  for(i in 1:dim.K){ delta.tilde = delta.tilde + vec.pi[i] *
    (1/(drop(matrix.lp[i,] %*% Ipq %*% vec.mu))) * (matrix.lp[i,] %o% matrix.lp[i,])}
  delta.tilde.inv <- solve(delta.tilde)
  
  covar_ASE <- function(lp){
    inner.covar <- matrix(0, nrow = dim.lp, ncol = dim.lp)
    for(i in 1:dim.K){ inner.covar = inner.covar + vec.pi[i] *
      drop((lp %*% Ipq %*% matrix.lp[i,]) * (1 - (lp %*% Ipq %*% matrix.lp[i,]))) * 
      matrix.lp[i,] %o% matrix.lp[i,]}
    return(Ipq %*% delta.inv %*% inner.covar %*% delta.inv %*% Ipq)
  }
  
  covar_LSE <- function(lp){
    temp.covar <- matrix(0, nrow = dim.lp, ncol = dim.lp)
    for(i in 1:dim.K){ temp.covar = temp.covar + vec.pi[i] *
      drop(((lp %*% Ipq %*% matrix.lp[i,])*(1 - lp %*% Ipq %*% matrix.lp[i,]))/(lp %*% Ipq %*% vec.mu)) *
      (as.vector(matrix.lp[i,] %*% delta.tilde.inv %*% Ipq / drop(matrix.lp[i,] %*% Ipq %*% vec.mu)) - (lp / drop(2*lp %*% Ipq %*% vec.mu))) %o%
      (as.vector(matrix.lp[i,] %*% delta.tilde.inv %*% Ipq / drop(matrix.lp[i,] %*% Ipq %*% vec.mu)) - (lp / drop(2*lp %*% Ipq %*% vec.mu))) }
    return(temp.covar)
  }
  
  chernoff_ASE <- function(t, lp1, lp2){
    (t*(1-t))*drop(as.vector(lp1 - lp2)%*%solve(t*covar_ASE(lp1) + (1-t)*covar_ASE(lp2))%*%as.vector(lp1 - lp2))
  }
  
  chernoff_LSE <- function(t, lp1, lp2){
    (t*(1-t))*
      drop(
        ((lp1/sqrt(drop(lp1 %*% Ipq %*% vec.mu))) - (lp2/sqrt(drop(lp2 %*% Ipq %*% vec.mu)))) %*%
          solve(t*covar_LSE(lp1) + (1-t)*covar_LSE(lp2)) %*%
          ((lp1/sqrt(drop(lp1 %*% Ipq %*% vec.mu))) - (lp2/sqrt(drop(lp2 %*% Ipq %*% vec.mu))))
      )
  }
  
  all.combs <- t(combn(seq(1:dim.K), 2))
  vec.rhoA <- as.vector(rep(0, length = dim(all.combs)[1]))
  vec.rhoL <- as.vector(rep(0, length = dim(all.combs)[1]))
  
  for(i in 1:dim(all.combs)[1]){
    temp_ASE <- function(t){ chernoff_ASE(t, matrix.lp[all.combs[i,1],], matrix.lp[all.combs[i,2],]) }
    vec.rhoA[i] <- optimize(temp_ASE, interval=c(0,1), maximum=TRUE)$objective
    
    temp_LSE <- function(t){ chernoff_LSE(t, matrix.lp[all.combs[i,1],], matrix.lp[all.combs[i,2],]) }
    vec.rhoL[i] <- optimize(temp_LSE, interval=c(0,1), maximum=TRUE)$objective
  }
  
  argmin.A <- all.combs[which(vec.rhoA == min(vec.rhoA)),]
  argmin.L <- all.combs[which(vec.rhoL == min(vec.rhoL)),]
  
  if (length(argmin.A) > 2 | length(argmin.L) > 2) {
    my.output <- list()
    my.output[[1]] <- min(vec.rhoA)
    my.output[[2]] <- all.combs[which(vec.rhoA == min(vec.rhoA)),]
    my.output[[3]] <- min(vec.rhoL)
    my.output[[4]] <- all.combs[which(vec.rhoL == min(vec.rhoL)),]
    my.output[[5]] <- min(vec.rhoA)/min(vec.rhoL)
    
    names(my.output) <- c("rho^star_ASE",
                          "argmin_pairs_ASE",
                          "rho^star_LSE",
                          "argmin_pairs_LSE",
                          "rho^star_ratio")
  } else {
    my.output <- c(min(vec.rhoA),
                   all.combs[which(vec.rhoA == min(vec.rhoA))[1],],
                   min(vec.rhoL),
                   all.combs[which(vec.rhoL == min(vec.rhoL))[1],],
                   min(vec.rhoA)/min(vec.rhoL))
    names(my.output) <- c("rho^star_ASE",
                          "argmin.k_ASE",
                          "argmin.l_ASE",
                          "rho^star_LSE",
                          "argmin.k_LSE",
                          "argmin.l_LSE",
                          "rho^star_ratio")
  }
  return(my.output)
  
}  


##########################################################
## [END] Master function: approx_chernoff_opt_v2 #########
##########################################################
## [END OF CODE] #########################################
##########################################################