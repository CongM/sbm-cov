#############################################################################################
##                                                                                         ##
## Supplementary Code                                                                      ##
##                                                                                         ##
## Cong Mu, Angelo Mele, Lingxin Hao, Joshua Cape, Avanti Athreya, and Carey E. Priebe     ##
##                                                                                         ##
#############################################################################################




#### Change the path or working directory if needed
source("algorithm.R")




#### Simulations
col <- c("#F8766D", "#00BA38", "#619CFF")
sha <- c(16, 17, 15)


## Fig. 4 upper right panel
allresults <- data.frame()

addCovariates <- TRUE
betahatmethod <- "SA"
dmax <- 10
m <- 100

K <- 2
d <- 1    

p <- 0.3
q <- 0.668
latent <- cbind(p, q)

beta <- 0.49
cov <- 2

pi <- rep(1/K, K)
pi_z <- rep(1/(K*cov), K*cov)

ns <- seq(from = 100, to = 260, by = 40)
for (n in ns) {
  block_size <- round(pi*n)
  blocks <- c()
  for (k in 1:length(block_size)) {
    blocks <- c(blocks, rep(k, block_size[k]))
  }
  
  block_size_z <- round(pi_z*n)
  covariates <- matrix(c(rep(1,block_size_z[1]), rep(2,block_size_z[2]), rep(1,block_size_z[3]), rep(2,block_size_z[4])))
  
  B <- generateB(latent, K, d, addCovariates, cov, beta)
  P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
  
  dhat <- rankMatrix(B)
  G <- K * cov
  
  for (i in 1:m) {
    cat(i, "\n")
    
    A <- generateA(n, P, seed = i)
    
    start_time1 <- Sys.time()
    result1 <- Algo1(A, cov, dhat, dmax, G)
    end_time1 <- Sys.time()
    Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))
    
    start_time2 <- Sys.time()
    result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
    end_time2 <- Sys.time()
    Elapsed_Algo2 <- as.numeric(difftime(end_time2,start_time2,units="secs"))
    
    A_tilde <- getAwithoutCovariates(A, beta, covariates)
    embedding_tilde <- SpectralEmbedding(A_tilde, dmax, maxit = 10000, work = 50)
    s_tilde <- embedding_tilde$D
    dhat_tilde <- dimselect(s_tilde)$elbow[1] + 1
    Yhat_tilde <- embedding_tilde$X[,1:dhat_tilde] %*% sqrt(diag(s_tilde[1:dhat_tilde], nrow=dhat_tilde, ncol=dhat_tilde))
    model_tilde <- Mclust(Yhat_tilde, G/prod(cov), verbose = FALSE)
    tautilde <- model_tilde$classification
    
    ARI_Algo1 <- adjustedRandIndex(blocks, result1$tauhat)
    ARI_Algo2_betahat <- adjustedRandIndex(blocks, result2$tautilde)
    ARI_Algo2_beta <- adjustedRandIndex(blocks, tautilde)
    
    allresults <- rbind(allresults, data.frame(n, p, q, beta, betahat = result2$betahat, 
                                               ARI_Algo1, ARI_Algo2_beta, ARI_Algo2_betahat, 
                                               Elapsed_Algo1, Elapsed_Algo2, seed = i))
  }
}

save(allresults, file = "figure4ur.rdata")

datA <- allresults %>%
  group_by(n) %>%
  summarise(p = mean(p),
            q = mean(q),
            beta = mean(beta),
            betahat_mean = mean(betahat),
            ARI_Algo1_mean = mean(ARI_Algo1),
            ARI_Algo1_se = sd(ARI_Algo1) / sqrt(m),
            ARI_Algo2_beta_mean = mean(ARI_Algo2_beta),
            ARI_Algo2_beta_se = sd(ARI_Algo2_beta) / sqrt(m),
            ARI_Algo2_betahat_mean = mean(ARI_Algo2_betahat),
            ARI_Algo2_betahat_se = sd(ARI_Algo2_betahat) / sqrt(m),
            Elapsed_Algo1_mean = mean(Elapsed_Algo1),
            Elapsed_Algo1_se = sd(Elapsed_Algo1) / sqrt(m),
            Elapsed_Algo2_mean = mean(Elapsed_Algo2),
            Elapsed_Algo2_se = sd(Elapsed_Algo2) / sqrt(m))


## Fig. 4 lower right panel
allresults <- data.frame()

addCovariates <- TRUE
betahatmethod <- "SA"
dmax <- 10
m <- 100

K <- 2
d <- 1    

p <- 0.3
q <- 0.564
latent <- cbind(p, q)

beta <- 0.49
cov <- 2

pi <- rep(1/K, K)
pi_z <- rep(1/(K*cov), K*cov)

ns <- seq(from = 100, to = 260, by = 40)
for (n in ns) {
  block_size <- round(pi*n)
  blocks <- c()
  for (k in 1:length(block_size)) {
    blocks <- c(blocks, rep(k, block_size[k]))
  }
  
  block_size_z <- round(pi_z*n)
  covariates <- matrix(c(rep(1,block_size_z[1]), rep(2,block_size_z[2]), rep(1,block_size_z[3]), rep(2,block_size_z[4])))
  
  B <- generateB(latent, K, d, addCovariates, cov, beta)
  P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
  
  dhat <- rankMatrix(B)
  G <- K * cov
  
  for (i in 1:m) {
    cat(i, "\n")
    
    A <- generateA(n, P, seed = i)
    
    start_time1 <- Sys.time()
    result1 <- Algo1(A, cov, dhat, dmax, G)
    end_time1 <- Sys.time()
    Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))
    
    start_time2 <- Sys.time()
    result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
    end_time2 <- Sys.time()
    Elapsed_Algo2 <- as.numeric(difftime(end_time2,start_time2,units="secs"))
    
    A_tilde <- getAwithoutCovariates(A, beta, covariates)
    embedding_tilde <- SpectralEmbedding(A_tilde, dmax, maxit = 10000, work = 50)
    s_tilde <- embedding_tilde$D
    dhat_tilde <- dimselect(s_tilde)$elbow[1] + 1
    Yhat_tilde <- embedding_tilde$X[,1:dhat_tilde] %*% sqrt(diag(s_tilde[1:dhat_tilde], nrow=dhat_tilde, ncol=dhat_tilde))
    model_tilde <- Mclust(Yhat_tilde, G/prod(cov), verbose = FALSE)
    tautilde <- model_tilde$classification
    
    ARI_Algo1 <- adjustedRandIndex(blocks, result1$tauhat)
    ARI_Algo2_betahat <- adjustedRandIndex(blocks, result2$tautilde)
    ARI_Algo2_beta <- adjustedRandIndex(blocks, tautilde)
    
    allresults <- rbind(allresults, data.frame(n, p, q, beta, betahat = result2$betahat, 
                                               ARI_Algo1, ARI_Algo2_beta, ARI_Algo2_betahat, 
                                               Elapsed_Algo1, Elapsed_Algo2, seed = i))
  }
}

save(allresults, file = "figure4lr.rdata")

datB <- allresults %>%
  group_by(n) %>%
  summarise(p = mean(p),
            q = mean(q),
            beta = mean(beta),
            betahat_mean = mean(betahat),
            ARI_Algo1_mean = mean(ARI_Algo1),
            ARI_Algo1_se = sd(ARI_Algo1) / sqrt(m),
            ARI_Algo2_beta_mean = mean(ARI_Algo2_beta),
            ARI_Algo2_beta_se = sd(ARI_Algo2_beta) / sqrt(m),
            ARI_Algo2_betahat_mean = mean(ARI_Algo2_betahat),
            ARI_Algo2_betahat_se = sd(ARI_Algo2_betahat) / sqrt(m),
            Elapsed_Algo1_mean = mean(Elapsed_Algo1),
            Elapsed_Algo1_se = sd(Elapsed_Algo1) / sqrt(m),
            Elapsed_Algo2_mean = mean(Elapsed_Algo2),
            Elapsed_Algo2_se = sd(Elapsed_Algo2) / sqrt(m))

dat <- rbind(datA, datB) %>%
  mutate(Chernoff = rep(c("A","B"), each = 5))

pp <- ggplot(dat) + 
  geom_line(aes(x=n, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 3) + 
  geom_point(aes(x=n, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 10) + 
  geom_errorbar(aes(x=n, ymin=ARI_Algo1_mean-ARI_Algo1_se, ymax=ARI_Algo1_mean+ARI_Algo1_se, color = 'Algo 1'), width = 4, size = 3) + 
  geom_line(aes(x=n, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 3) + 
  geom_point(aes(x=n, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 10) + 
  geom_errorbar(aes(x=n, ymin=ARI_Algo2_beta_mean-ARI_Algo2_beta_se, ymax=ARI_Algo2_beta_mean+ARI_Algo2_beta_se, color = 'Algo 2 with beta'), width = 4, size = 3) + 
  geom_line(aes(x=n, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 3) + 
  geom_point(aes(x=n, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 10) + 
  geom_errorbar(aes(x=n, ymin=ARI_Algo2_betahat_mean-ARI_Algo2_betahat_se, ymax=ARI_Algo2_betahat_mean+ARI_Algo2_betahat_se, color = 'Algo 2 with betahat'), width = 4, size = 3) + 
  geom_hline(yintercept = 1, color = 'black', size = 4) + 
  scale_x_continuous(breaks = dat$n) + facet_grid(Chernoff ~ .) + 
  labs(x = 'n', y = 'ARI', color = ' ', shape = ' ') + 
  theme_bw() + 
  scale_color_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = col) + 
  scale_shape_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = sha) + 
  theme(axis.title.x = element_text(size = 60),
        axis.text.x = element_text(size = 58),
        axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 58),
        legend.title = element_text(size = 60),
        legend.text = element_text(size = 58),
        strip.text.y = element_text(size = 58),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = 'top')
pp


## Fig. 5(a)
allresults <- data.frame()

addCovariates <- TRUE
betahatmethod <- "SA"
dmax <- 10
m <- 100

K <- 2
d <- 1    

p <- 0.3
qs <- seq(from = 0.35, to = 0.45, by = 0.025)

beta <- 0.4
cov <- 5

pi <- rep(1/K, K)
pi_z <- rep(1/(K*cov), K*cov)

n <- 2000
block_size <- round(pi*n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

block_size_z <- round(pi_z*n)
covariates <- matrix(c(rep(1,block_size_z[1]), rep(2,block_size_z[2]), rep(3,block_size_z[3]), rep(4,block_size_z[4]), rep(5,block_size_z[5]), 
                       rep(1,block_size_z[6]), rep(2,block_size_z[7]), rep(3,block_size_z[8]), rep(4,block_size_z[9]), rep(5,block_size_z[10])))

for (q in qs) {
  latent <- cbind(p, q)
  
  B <- generateB(latent, K, d, addCovariates, cov, beta)
  P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
  
  dhat <- rankMatrix(B)
  G <- K * cov
  
  for (i in 1:m) {
    cat(i, "\n")
    
    A <- generateA(n, P, seed = i)
    
    start_time1 <- Sys.time()
    result1 <- Algo1(A, cov, dhat, dmax, G)
    end_time1 <- Sys.time()
    Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))
    
    start_time2 <- Sys.time()
    result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
    end_time2 <- Sys.time()
    Elapsed_Algo2 <- as.numeric(difftime(end_time2,start_time2,units="secs"))
    
    A_tilde <- getAwithoutCovariates(A, beta, covariates)
    embedding_tilde <- SpectralEmbedding(A_tilde, dmax, maxit = 10000, work = 50)
    s_tilde <- embedding_tilde$D
    dhat_tilde <- dimselect(s_tilde)$elbow[1] + 1
    Yhat_tilde <- embedding_tilde$X[,1:dhat_tilde] %*% sqrt(diag(s_tilde[1:dhat_tilde], nrow=dhat_tilde, ncol=dhat_tilde))
    model_tilde <- Mclust(Yhat_tilde, G/prod(cov), verbose = FALSE)
    tautilde <- model_tilde$classification
    
    ARI_Algo1 <- adjustedRandIndex(blocks, result1$tauhat)
    ARI_Algo2_betahat <- adjustedRandIndex(blocks, result2$tautilde)
    ARI_Algo2_beta <- adjustedRandIndex(blocks, tautilde)
    
    allresults <- rbind(allresults, data.frame(n, p, q, beta, betahat = result2$betahat, 
                                               ARI_Algo1, ARI_Algo2_beta, ARI_Algo2_betahat, 
                                               Elapsed_Algo1, Elapsed_Algo2, seed = i))
  }
}

save(allresults, file = "figure5a.rdata")

dat <- allresults %>%
  mutate(delta_latent = abs(q-p)) %>%
  group_by(q) %>%
  summarise(p = mean(p),
            n = mean(n),
            delta_latent = mean(delta_latent),
            beta = mean(beta),
            betahat_mean = mean(betahat),
            ARI_Algo1_mean = mean(ARI_Algo1),
            ARI_Algo1_se = sd(ARI_Algo1) / sqrt(m),
            ARI_Algo2_beta_mean = mean(ARI_Algo2_beta),
            ARI_Algo2_beta_se = sd(ARI_Algo2_beta) / sqrt(m),
            ARI_Algo2_betahat_mean = mean(ARI_Algo2_betahat),
            ARI_Algo2_betahat_se = sd(ARI_Algo2_betahat) / sqrt(m),
            Elapsed_Algo1_mean = mean(Elapsed_Algo1),
            Elapsed_Algo1_se = sd(Elapsed_Algo1) / sqrt(m),
            Elapsed_Algo2_mean = mean(Elapsed_Algo2),
            Elapsed_Algo2_se = sd(Elapsed_Algo2) / sqrt(m))

pp <- ggplot(dat) + 
  geom_line(aes(x=delta_latent, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 3) + 
  geom_point(aes(x=delta_latent, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 10) + 
  geom_errorbar(aes(x=delta_latent, ymin=ARI_Algo1_mean-ARI_Algo1_se, ymax=ARI_Algo1_mean+ARI_Algo1_se, color = 'Algo 1'), width = 0.0025, size = 3) + 
  geom_line(aes(x=delta_latent, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 3) + 
  geom_point(aes(x=delta_latent, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 10) + 
  geom_errorbar(aes(x=delta_latent, ymin=ARI_Algo2_beta_mean-ARI_Algo2_beta_se, ymax=ARI_Algo2_beta_mean+ARI_Algo2_beta_se, color = 'Algo 2 with beta'), width = 0.0025, size = 3) + 
  geom_line(aes(x=delta_latent, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 3) + 
  geom_point(aes(x=delta_latent, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 10) + 
  geom_errorbar(aes(x=delta_latent, ymin=ARI_Algo2_betahat_mean-ARI_Algo2_betahat_se, ymax=ARI_Algo2_betahat_mean+ARI_Algo2_betahat_se, color = 'Algo 2 with betahat'), width = 0.0025, size = 3) + 
  geom_hline(yintercept = 1, color = 'black', size = 4) + 
  labs(x = 'q-p', y = 'ARI', color = ' ', shape = ' ') + 
  theme_bw() + 
  scale_color_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = col) + 
  scale_shape_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = sha) + 
  theme(axis.title.x = element_text(size = 60),
        axis.text.x = element_text(size = 58),
        axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 58),
        legend.title = element_text(size = 60),
        legend.text = element_text(size = 58),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = 'top')
pp


## Fig. 5(b)
allresults <- data.frame()

addCovariates <- TRUE
betahatmethod <- "SA"
dmax <- 10
m <- 100

K <- 2
d <- 1    

p <- 0.3
q <- 0.375
latent <- cbind(p, q)

betas <- seq(from = 0.1, to = 0.3, by = 0.05)
cov <- 5

pi <- rep(1/K, K)
pi_z <- rep(1/(K*cov), K*cov)

n <- 2000
block_size <- round(pi*n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

block_size_z <- round(pi_z*n)
covariates <- matrix(c(rep(1,block_size_z[1]), rep(2,block_size_z[2]), rep(3,block_size_z[3]), rep(4,block_size_z[4]), rep(5,block_size_z[5]), 
                       rep(1,block_size_z[6]), rep(2,block_size_z[7]), rep(3,block_size_z[8]), rep(4,block_size_z[9]), rep(5,block_size_z[10])))

for (beta in betas) {
  B <- generateB(latent, K, d, addCovariates, cov, beta)
  P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
  
  dhat <- rankMatrix(B)
  G <- K * cov
  
  for (i in 1:m) {
    cat(i, "\n")
    
    A <- generateA(n, P, seed = i)
    
    start_time1 <- Sys.time()
    result1 <- Algo1(A, cov, dhat, dmax, G)
    end_time1 <- Sys.time()
    Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))
    
    start_time2 <- Sys.time()
    result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
    end_time2 <- Sys.time()
    Elapsed_Algo2 <- as.numeric(difftime(end_time2,start_time2,units="secs"))
    
    A_tilde <- getAwithoutCovariates(A, beta, covariates)
    embedding_tilde <- SpectralEmbedding(A_tilde, dmax, maxit = 10000, work = 50)
    s_tilde <- embedding_tilde$D
    dhat_tilde <- dimselect(s_tilde)$elbow[1] + 1
    Yhat_tilde <- embedding_tilde$X[,1:dhat_tilde] %*% sqrt(diag(s_tilde[1:dhat_tilde], nrow=dhat_tilde, ncol=dhat_tilde))
    model_tilde <- Mclust(Yhat_tilde, G/prod(cov), verbose = FALSE)
    tautilde <- model_tilde$classification
    
    ARI_Algo1 <- adjustedRandIndex(blocks, result1$tauhat)
    ARI_Algo2_betahat <- adjustedRandIndex(blocks, result2$tautilde)
    ARI_Algo2_beta <- adjustedRandIndex(blocks, tautilde)
    
    allresults <- rbind(allresults, data.frame(n, p, q, beta, betahat = result2$betahat, 
                                               ARI_Algo1, ARI_Algo2_beta, ARI_Algo2_betahat, 
                                               Elapsed_Algo1, Elapsed_Algo2, seed = i))
  }
}

save(allresults, file = "figure5b.rdata")

dat <- allresults %>%
  group_by(beta) %>%
  summarise(n = mean(n),
            p = mean(p),
            q = mean(q),
            beta = mean(beta),
            betahat_mean = mean(betahat),
            ARI_Algo1_mean = mean(ARI_Algo1),
            ARI_Algo1_se = sd(ARI_Algo1) / sqrt(m),
            ARI_Algo2_beta_mean = mean(ARI_Algo2_beta),
            ARI_Algo2_beta_se = sd(ARI_Algo2_beta) / sqrt(m),
            ARI_Algo2_betahat_mean = mean(ARI_Algo2_betahat),
            ARI_Algo2_betahat_se = sd(ARI_Algo2_betahat) / sqrt(m),
            Elapsed_Algo1_mean = mean(Elapsed_Algo1),
            Elapsed_Algo1_se = sd(Elapsed_Algo1) / sqrt(m),
            Elapsed_Algo2_mean = mean(Elapsed_Algo2),
            Elapsed_Algo2_se = sd(Elapsed_Algo2) / sqrt(m))

pp <- ggplot(dat) + 
  geom_line(aes(x=beta, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 3) + 
  geom_point(aes(x=beta, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 10) + 
  geom_errorbar(aes(x=beta, ymin=ARI_Algo1_mean-ARI_Algo1_se, ymax=ARI_Algo1_mean+ARI_Algo1_se, color = 'Algo 1'), width = 0.005, size = 3) + 
  geom_line(aes(x=beta, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 3) + 
  geom_point(aes(x=beta, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 10) + 
  geom_errorbar(aes(x=beta, ymin=ARI_Algo2_beta_mean-ARI_Algo2_beta_se, ymax=ARI_Algo2_beta_mean+ARI_Algo2_beta_se, color = 'Algo 2 with beta'), width = 0.005, size = 3) + 
  geom_line(aes(x=beta, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 3) + 
  geom_point(aes(x=beta, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 10) + 
  geom_errorbar(aes(x=beta, ymin=ARI_Algo2_betahat_mean-ARI_Algo2_betahat_se, ymax=ARI_Algo2_betahat_mean+ARI_Algo2_betahat_se, color = 'Algo 2 with betahat'), width = 0.005, size = 3) + 
  geom_hline(yintercept = 1, color = 'black', size = 4) + 
  labs(x = expression(beta), y = 'ARI', color = ' ', shape = ' ') + 
  theme_bw() + 
  scale_color_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = col) + 
  scale_shape_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = sha) + 
  theme(axis.title.x = element_text(size = 60),
        axis.text.x = element_text(size = 58),
        axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 58),
        legend.title = element_text(size = 60),
        legend.text = element_text(size = 58),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = 'top')
pp


## Fig. 6(a)
allresults <- data.frame()

addCovariates <- TRUE
betahatmethod <- "SA"
dmax <- 10
m <- 100

K <- 2
d <- 2  

as <- seq(from = 0.12, to = 0.14, by = 0.005)
b <- 0.1

beta <- 0.2
cov <- 5

pi <- rep(1/K, K)
pi_z <- rep(1/(K*cov), K*cov)

n <- 2000
block_size <- round(pi*n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

block_size_z <- round(pi_z*n)
covariates <- matrix(c(rep(1,block_size_z[1]), rep(2,block_size_z[2]), rep(3,block_size_z[3]), rep(4,block_size_z[4]), rep(5,block_size_z[5]), 
                       rep(1,block_size_z[6]), rep(2,block_size_z[7]), rep(3,block_size_z[8]), rep(4,block_size_z[9]), rep(5,block_size_z[10])))

for (a in as) {
  latent <- cbind(c(sqrt(a), 0), c(b/sqrt(a), sqrt((a+b)*(a-b)/a)))
  
  B <- generateB(latent, K, d, addCovariates, cov, beta)
  P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
  
  dhat <- rankMatrix(B)
  G <- K * cov
  
  for (i in 1:m) {
    cat(i, "\n")
    
    A <- generateA(n, P, seed = i)
    
    start_time1 <- Sys.time()
    result1 <- Algo1(A, cov, dhat, dmax, G)
    end_time1 <- Sys.time()
    Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))
    
    start_time2 <- Sys.time()
    result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
    end_time2 <- Sys.time()
    Elapsed_Algo2 <- as.numeric(difftime(end_time2,start_time2,units="secs"))
    
    A_tilde <- getAwithoutCovariates(A, beta, covariates)
    embedding_tilde <- SpectralEmbedding(A_tilde, dmax, maxit = 10000, work = 50)
    s_tilde <- embedding_tilde$D
    dhat_tilde <- dimselect(s_tilde)$elbow[1] + 1
    Yhat_tilde <- embedding_tilde$X[,1:dhat_tilde] %*% sqrt(diag(s_tilde[1:dhat_tilde], nrow=dhat_tilde, ncol=dhat_tilde))
    model_tilde <- Mclust(Yhat_tilde, G/prod(cov), verbose = FALSE)
    tautilde <- model_tilde$classification
    
    ARI_Algo1 <- adjustedRandIndex(blocks, result1$tauhat)
    ARI_Algo2_betahat <- adjustedRandIndex(blocks, result2$tautilde)
    ARI_Algo2_beta <- adjustedRandIndex(blocks, tautilde)
    
    allresults <- rbind(allresults, data.frame(n, a, b, beta, betahat = result2$betahat, 
                                               ARI_Algo1, ARI_Algo2_beta, ARI_Algo2_betahat, 
                                               Elapsed_Algo1, Elapsed_Algo2, seed = i))
  }
}

save(allresults, file = "figure6a.rdata")

dat <- allresults %>%
  mutate(delta_latent = 2*abs(a-b)) %>%
  group_by(a) %>%
  summarise(b = mean(b),
            n = mean(n),
            delta_latent = mean(delta_latent),
            beta = mean(beta),
            betahat_mean = mean(betahat),
            ARI_Algo1_mean = mean(ARI_Algo1),
            ARI_Algo1_se = sd(ARI_Algo1) / sqrt(m),
            ARI_Algo2_beta_mean = mean(ARI_Algo2_beta),
            ARI_Algo2_beta_se = sd(ARI_Algo2_beta) / sqrt(m),
            ARI_Algo2_betahat_mean = mean(ARI_Algo2_betahat),
            ARI_Algo2_betahat_se = sd(ARI_Algo2_betahat) / sqrt(m),
            Elapsed_Algo1_mean = mean(Elapsed_Algo1),
            Elapsed_Algo1_se = sd(Elapsed_Algo1) / sqrt(m),
            Elapsed_Algo2_mean = mean(Elapsed_Algo2),
            Elapsed_Algo2_se = sd(Elapsed_Algo2) / sqrt(m))

pp <- ggplot(dat) + 
  geom_line(aes(x=delta_latent, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 3) + 
  geom_point(aes(x=delta_latent, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 10) + 
  geom_errorbar(aes(x=delta_latent, ymin=ARI_Algo1_mean-ARI_Algo1_se, ymax=ARI_Algo1_mean+ARI_Algo1_se, color = 'Algo 1'), width = 0.001, size = 3) + 
  geom_line(aes(x=delta_latent, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 3) + 
  geom_point(aes(x=delta_latent, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 10) + 
  geom_errorbar(aes(x=delta_latent, ymin=ARI_Algo2_beta_mean-ARI_Algo2_beta_se, ymax=ARI_Algo2_beta_mean+ARI_Algo2_beta_se, color = 'Algo 2 with beta'), width = 0.001, size = 3) + 
  geom_line(aes(x=delta_latent, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 3) + 
  geom_point(aes(x=delta_latent, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 10) + 
  geom_errorbar(aes(x=delta_latent, ymin=ARI_Algo2_betahat_mean-ARI_Algo2_betahat_se, ymax=ARI_Algo2_betahat_mean+ARI_Algo2_betahat_se, color = 'Algo 2 with betahat'), width = 0.001, size = 3) + 
  geom_hline(yintercept = 1, color = 'black', size = 4) + 
  labs(x = '2(a-b)', y = 'ARI', color = ' ', shape = ' ') + 
  theme_bw() + 
  scale_color_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = col) + 
  scale_shape_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = sha) + 
  theme(axis.title.x = element_text(size = 60),
        axis.text.x = element_text(size = 58),
        axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 58),
        legend.title = element_text(size = 60),
        legend.text = element_text(size = 58),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = 'top')
pp


## Fig. 6(b)
allresults <- data.frame()

addCovariates <- TRUE
betahatmethod <- "SA"
dmax <- 10
m <- 100

K <- 2
d <- 2  

a <- 0.135
b <- 0.1
latent <- cbind(c(sqrt(a), 0), c(b/sqrt(a), sqrt((a+b)*(a-b)/a)))

betas <- seq(from = -0.09, to = -0.05, by = 0.01)
cov <- 5

pi <- rep(1/K, K)
pi_z <- rep(1/(K*cov), K*cov)

n <- 2000
block_size <- round(pi*n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

block_size_z <- round(pi_z*n)
covariates <- matrix(c(rep(1,block_size_z[1]), rep(2,block_size_z[2]), rep(3,block_size_z[3]), rep(4,block_size_z[4]), rep(5,block_size_z[5]), 
                       rep(1,block_size_z[6]), rep(2,block_size_z[7]), rep(3,block_size_z[8]), rep(4,block_size_z[9]), rep(5,block_size_z[10])))

for (beta in betas) {
  B <- generateB(latent, K, d, addCovariates, cov, beta)
  P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
  
  dhat <- rankMatrix(B)
  G <- K * cov
  
  for (i in 1:m) {
    cat(i, "\n")
    
    A <- generateA(n, P, seed = i)
    
    start_time1 <- Sys.time()
    result1 <- Algo1(A, cov, dhat, dmax, G)
    end_time1 <- Sys.time()
    Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))
    
    start_time2 <- Sys.time()
    result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
    end_time2 <- Sys.time()
    Elapsed_Algo2 <- as.numeric(difftime(end_time2,start_time2,units="secs"))
    
    A_tilde <- getAwithoutCovariates(A, beta, covariates)
    embedding_tilde <- SpectralEmbedding(A_tilde, dmax, maxit = 10000, work = 50)
    s_tilde <- embedding_tilde$D
    dhat_tilde <- dimselect(s_tilde)$elbow[1] + 1
    Yhat_tilde <- embedding_tilde$X[,1:dhat_tilde] %*% sqrt(diag(s_tilde[1:dhat_tilde], nrow=dhat_tilde, ncol=dhat_tilde))
    model_tilde <- Mclust(Yhat_tilde, G/prod(cov), verbose = FALSE)
    tautilde <- model_tilde$classification
    
    ARI_Algo1 <- adjustedRandIndex(blocks, result1$tauhat)
    ARI_Algo2_betahat <- adjustedRandIndex(blocks, result2$tautilde)
    ARI_Algo2_beta <- adjustedRandIndex(blocks, tautilde)
    
    allresults <- rbind(allresults, data.frame(n, a, b, beta, betahat = result2$betahat, 
                                               ARI_Algo1, ARI_Algo2_beta, ARI_Algo2_betahat, 
                                               Elapsed_Algo1, Elapsed_Algo2, seed = i))
  }
}

save(allresults, file = "figure6b.rdata")

dat <- allresults %>%
  group_by(beta) %>%
  summarise(n = mean(n),
            a = mean(a),
            b = mean(b),
            beta = mean(beta),
            betahat_mean = mean(betahat),
            ARI_Algo1_mean = mean(ARI_Algo1),
            ARI_Algo1_se = sd(ARI_Algo1) / sqrt(m),
            ARI_Algo2_beta_mean = mean(ARI_Algo2_beta),
            ARI_Algo2_beta_se = sd(ARI_Algo2_beta) / sqrt(m),
            ARI_Algo2_betahat_mean = mean(ARI_Algo2_betahat),
            ARI_Algo2_betahat_se = sd(ARI_Algo2_betahat) / sqrt(m),
            Elapsed_Algo1_mean = mean(Elapsed_Algo1),
            Elapsed_Algo1_se = sd(Elapsed_Algo1) / sqrt(m),
            Elapsed_Algo2_mean = mean(Elapsed_Algo2),
            Elapsed_Algo2_se = sd(Elapsed_Algo2) / sqrt(m))

pp <- ggplot(dat) + 
  geom_line(aes(x=beta, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 3) + 
  geom_point(aes(x=beta, y=ARI_Algo1_mean, color = 'Algo 1', shape = 'Algo 1'), size = 10) + 
  geom_errorbar(aes(x=beta, ymin=ARI_Algo1_mean-ARI_Algo1_se, ymax=ARI_Algo1_mean+ARI_Algo1_se, color = 'Algo 1'), width = 0.001, size = 3) + 
  geom_line(aes(x=beta, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 3) + 
  geom_point(aes(x=beta, y=ARI_Algo2_beta_mean, color = 'Algo 2 with beta', shape = 'Algo 2 with beta'), size = 10) + 
  geom_errorbar(aes(x=beta, ymin=ARI_Algo2_beta_mean-ARI_Algo2_beta_se, ymax=ARI_Algo2_beta_mean+ARI_Algo2_beta_se, color = 'Algo 2 with beta'), width = 0.001, size = 3) + 
  geom_line(aes(x=beta, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 3) + 
  geom_point(aes(x=beta, y=ARI_Algo2_betahat_mean, color = 'Algo 2 with betahat', shape = 'Algo 2 with betahat'), size = 10) + 
  geom_errorbar(aes(x=beta, ymin=ARI_Algo2_betahat_mean-ARI_Algo2_betahat_se, ymax=ARI_Algo2_betahat_mean+ARI_Algo2_betahat_se, color = 'Algo 2 with betahat'), width = 0.001, size = 3) + 
  geom_hline(yintercept = 1, color = 'black', size = 4) + 
  labs(x = expression(beta), y = 'ARI', color = ' ', shape = ' ') + 
  theme_bw() + 
  scale_color_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = col) + 
  scale_shape_manual(labels = c('Algo 1', expression(paste('Algo 2 with ', beta)), expression(paste('Algo 2 with ', widehat(beta)))), values = sha) + 
  theme(axis.title.x = element_text(size = 60),
        axis.text.x = element_text(size = 58),
        axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 58),
        legend.title = element_text(size = 60),
        legend.text = element_text(size = 58),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = 'top')
pp



