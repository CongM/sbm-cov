#############################################################################################
##                                                                                         ##
## Supplementary Code                                                                      ##
##                                                                                         ##
## Cong Mu, Angelo Mele, Lingxin Hao, Joshua Cape, Avanti Athreya, and Carey E. Priebe     ##
##                                                                                         ##
#############################################################################################




#### Required packages
# install.packages("Matrix")
# install.packages("dplyr")
# install.packages("mclust")
# install.packages("jsonlite")
# install.packages("irlba")
# install.packages("ggplot2")
# install.packages("ggthemes")
# install.packages("ggExtra")
# install.packages("igraph")
# devtools::install_github("neurodata/graphstats")
# devtools::install_github("meleangelo/grdpg")

library(Matrix)
library(dplyr)
library(mclust)
library(jsonlite)
library(irlba)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(igraph)
library(graphstats)
library(grdpg)




#### Algorithm 1 - (Section 3)
## Estimation of induced block assignment including the vertex covariate effect
##
## Input -
## A: An adjacency matrix.
## cov: Possible value that the covariate can take. `2` by default for one binary covariate.
## dhat: Embedding dimension. `NULL` by default to choose by profile likelihood (Remark 2).
## dmax: `nv` for `irlba::irlba`. `10` by default. See `help(irlba, package = "irlba")`.
## G: `G` for `mclust::Mclust`. `1:10` by default. See `help(Mclust, package = "mclust")`.
##
## Output -
## xihat: Block assignments including the vertex covariate effect.
## tauhat: Induced block assignments after accounting for the vertex covariate effect.
##
####

Algo1 <- function(A, cov = 2, dhat = NULL, dmax = 10, G = 1:10) {
  
  output <- list()
  
  ## Step 1
  embedding <- SpectralEmbedding(A, dmax, maxit = 10000, work = 50)
  s <- embedding$D
  dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1]+1, dhat)
  Yhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  
  ## Step 2
  model <- Mclust(Yhat, G, verbose = FALSE)
  xihat <- model$classification
  
  ## Step 3
  Ipq <- getIpq(A, dhat)
  muhats <- matrix(model$parameters$mean, nrow = dhat)
  BZhat <- t(muhats) %*% Ipq %*% muhats
  
  ## Step 4
  model2 <- Mclust(diag(BZhat), ncol(BZhat)/prod(cov), verbose = FALSE)
  phihat <- model2$classification
  
  ## Step 5
  tauhat <- xihat
  for (c in unique(phihat)) {
    ind1 <- which(phihat == c)
    ind2 <- which(xihat %in% ind1)
    tauhat[ind2] <- c
  }
  
  output$xihat <- xihat
  output$tauhat <- tauhat
  
  return(output)
}



#### Algorithm 2 - (Section 3) 
## Estimation of induced block assignment after accounting for the vertex covariate effect
##
## Input -
## A: An adjacency matrix.
## covariates: A matrix of observed covariates.
## betahatmethod: Procedure to estimate vertex covariate effect. "SA" (Step 2a) by default or "WA" (Step 2b).
## cov: Possible value that the covariate can take. `2` by default for one binary covariate.
## dhat: Embedding dimension. `NULL` by default to choose by profile likelihood (Remark 2).
## dmax: `nv` for `irlba::irlba`. `10` by default. See `help(irlba, package = "irlba")`.
## G: `G` for `mclust::Mclust`. `1:10` by default. See `help(Mclust, package = "mclust")`.
##
## Output -
## xihat: Block assignments including the vertex covariate effect.
## tautilde: Induced block assignments after accounting for the vertex covariate effect.
## betahat: Vertex covariate effect.
##
####


Algo2 <- function(A, covariates, betahatmethod = "SA", cov = 2, dhat = NULL, dmax = 10, G = 1:10) {
  
  output <- list()
  
  ## Step 1 (1-4 in Algo1)
  embedding <- SpectralEmbedding(A, dmax, maxit = 10000, work = 50)
  s <- embedding$D
  dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1]+1, dhat)
  Yhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  
  model <- Mclust(Yhat, G, verbose = FALSE)
  xihat <- model$classification
  
  Ipq <- getIpq(A, dhat)
  muhats <- matrix(model$parameters$mean, nrow = dhat)
  BZhat <- t(muhats) %*% Ipq %*% muhats
  
  model2 <- Mclust(diag(BZhat), ncol(BZhat)/prod(cov), verbose = FALSE)
  phihat <- model2$classification
  
  ## Step 2
  if (betahatmethod == "SA") {
    ## Step 2a
    covariates_block <- getBlockCovariates(covariates, xihat)
    result <- estimatebeta(Yhat, muhats, Ipq, cov, covariates_block, xihat, sd = FALSE)
    betahat <- sapply(result$betahats, mean)
  } else {
    ## Step 2b
    result <- estimatebeta2(Yhat, muhats, Ipq, cov, covariates, xihat, sd = FALSE)
    betahat <- sapply(Map('*',result$betahats,result$pis), sum)
  }
  
  ## Step 3
  A_tilde <- getAwithoutCovariates(A, betahat, covariates)
  
  ## Step 4
  embedding_tilde <- SpectralEmbedding(A_tilde, dmax, maxit = 10000, work = 50)
  s_tilde <- embedding_tilde$D
  dhat_tilde <- dimselect(s_tilde)$elbow[1] + 1
  Yhat_tilde <- embedding_tilde$X[,1:dhat_tilde] %*% sqrt(diag(s_tilde[1:dhat_tilde], nrow=dhat_tilde, ncol=dhat_tilde))
  
  ## Step 5
  model_tilde <- Mclust(Yhat_tilde, ncol(BZhat)/prod(cov), verbose = FALSE)
  tautilde <- model_tilde$classification
  
  output$xihat <- xihat
  output$tautilde <- tautilde
  output$betahat <- betahat
  
  return(output)
}



