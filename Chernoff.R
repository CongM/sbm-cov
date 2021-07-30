#############################################################################################
##                                                                                         ##
## Supplementary Code                                                                      ##
##                                                                                         ##
## Cong Mu, Angelo Mele, Lingxin Hao, Joshua Cape, Avanti Athreya, and Carey E. Priebe     ##
##                                                                                         ##
#############################################################################################




#### Change the path or working directory if needed
source("algorithm.R")
source("codeFunc_approx_chernoff_optimize_v2.R")




#### Chernoff ratio


## Fig. 1
addCovariates <- TRUE
K <- 2
d <- 1

p <- 0.3
qs <- seq(0.31, 0.7, 0.01)
betas <- seq(0.1, 0.5, 0.01)
cov <- 2

pi <- rep(1/K, K)
pi_Z <- rep(1/(K*cov), K*cov)

allresults <- data.frame()
for (q in qs) {
  for (beta in betas) {
    latent <- cbind(p, q)
    B <- t(latent) %*% latent
    B_Z <- generateB(latent, K, d, addCovariates, cov, beta)
    temp1 <- approx_chernoff_opt_v2(B_Z, pi_Z)
    temp2 <- approx_chernoff_opt_v2(B, pi)
    rho1 <- unlist(temp1[1])
    rho2 <- unlist(temp2[1])
    rhostar <- rho1 / rho2
    allresults <- bind_rows(allresults, data.frame(q, beta, rho1, rho2, rhostar))
  }
}

save(allresults, file = "figure1.rdata")

g <- ggplot(allresults, aes(x = q, y = beta, z = rhostar)) +
  geom_tile(aes(fill = rhostar)) +
  stat_contour(aes(color = ..level..), size = 2) +
  scale_fill_distiller(palette = "Spectral", breaks = seq(0.3, 1.1, 0.2), guide = guide_colourbar(ticks.colour = "black", ticks.linewidth = 3)) +
  stat_contour(aes(color = ..level..), breaks = 1, color = "black", size = 3) + 
  labs(y = expression(beta), fill = "Chernoff ratio") + coord_fixed() + theme_bw() + 
  theme(plot.title = element_text(size = 40),
        plot.subtitle = element_text(size = 38),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 38),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 38),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38),
        legend.key.width = unit(3, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "top")

g <- direct.label(g, list("top.pieces", cex=3))


## Fig. 2
K <- 2

b <- 0.1
as <- seq(b+0.01, 0.5, 0.01)
betas <- seq(0.1, 0.5, 0.01)
cov <- 2

pi <- rep(1/K, K)
pi_Z <- rep(1/(K*cov), K*cov)

allresults <- data.frame()
for (a in as) {
  for (beta in betas) {
    B <- matrix(c(a,b,b,a), nrow = 2)
    B_Z <- matrix(c(a+beta,a,b+beta,b,a,a+beta,b,b+beta,b+beta,b,a+beta,a,b,b+beta,a,a+beta), nrow = 4, byrow = TRUE)
    temp1 <- approx_chernoff_opt_v2(B_Z, pi_Z)
    temp2 <- approx_chernoff_opt_v2(B, pi)
    rho1 <- unlist(temp1[1])
    rho2 <- unlist(temp2[1])
    rhostar <- rho1 / rho2
    allresults <- bind_rows(allresults, data.frame(a, beta, rho1, rho2, rhostar))
  }
}

save(allresults, file = "figure2.rdata")

g <- ggplot(allresults, aes(x = a, y = beta, z = rhostar)) +
  geom_tile(aes(fill = rhostar)) +
  stat_contour(aes(color = ..level..), size = 2) +
  scale_fill_distiller(palette = "Spectral", breaks = seq(0.3, 1.1, 0.2), guide = guide_colourbar(ticks.colour = "black", ticks.linewidth = 3)) + 
  stat_contour(aes(color = ..level..), breaks = 1, color = "black", size = 3) + 
  labs(y = expression(beta), fill = "Chernoff ratio") + coord_fixed() + theme_bw() + 
  theme(plot.title = element_text(size = 40),
        plot.subtitle = element_text(size = 38),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 38),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 38),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38),
        legend.key.width = unit(3, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "top")

g <- direct.label(g, list("top.pieces", cex=3))


## Fig. 3
K <- 4

b <- 0.1
as <- seq(b+0.01, 0.5, 0.01)
betas <- seq(0.1, 0.5, 0.01)
cov <- 2

pi <- rep(1/K, K)
pi_Z <- rep(1/(K*cov), K*cov)

allresults <- data.frame()
for (a in as) {
  for (beta in betas) {
    B <- matrix(c(a,b,b,b,
                  b,a,b,b,
                  b,b,a,b,
                  b,b,b,a),
                nrow = 4, byrow = TRUE)
    B_Z <- matrix(c(a+beta,a,b+beta,b,b+beta,b,b+beta,b,
                    a,a+beta,b,b+beta,b,b+beta,b,b+beta,
                    b+beta,b,a+beta,a,b+beta,b,b+beta,b,
                    b,b+beta,a,a+beta,b,b+beta,b,b+beta,
                    b+beta,b,b+beta,b,a+beta,a,b+beta,b,
                    b,b+beta,b,b+beta,a,a+beta,b,b+beta,
                    b+beta,b,b+beta,b,b+beta,b,a+beta,a,
                    b,b+beta,b,b+beta,b,b+beta,a,a+beta),
                  nrow = 8, byrow = TRUE)    
    temp1 <- approx_chernoff_opt_v2(B_Z, pi_Z)
    temp2 <- approx_chernoff_opt_v2(B, pi)
    rho1 <- unlist(temp1[1])
    rho2 <- unlist(temp2[1])
    rhostar <- rho1 / rho2
    allresults <- bind_rows(allresults, data.frame(a, beta, rho1, rho2, rhostar))
  }
}

save(allresults, file = "figure3.rdata")

g <- ggplot(allresults, aes(x = a, y = beta, z = rhostar)) +
  geom_tile(aes(fill = rhostar)) +
  stat_contour(aes(color = ..level..), size = 2) +
  scale_fill_distiller(palette = "Spectral", breaks = seq(0.3, 1.1, 0.2), guide = guide_colourbar(ticks.colour = "black", ticks.linewidth = 3)) + 
  stat_contour(aes(color = ..level..), breaks = 1, color = "black", size = 3) + 
  labs(y = expression(beta), fill = "Chernoff ratio") + coord_fixed() + theme_bw() + 
  theme(plot.title = element_text(size = 40),
        plot.subtitle = element_text(size = 38),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 38),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 38),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38),
        legend.key.width = unit(3, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "top")

g <- direct.label(g, list("top.pieces", cex=3))


## Fig. 4 left panel
load("figure1.rdata")

g <- ggplot(allresults, aes(x = q, y = beta, z = rhostar)) +
  geom_tile(aes(fill = rhostar)) +
  stat_contour(aes(color = ..level..), size = 2) +
  scale_fill_distiller(palette = "Spectral", breaks = seq(0.3, 1.1, 0.2), guide = guide_colourbar(ticks.colour = "black", ticks.linewidth = 3)) +
  stat_contour(aes(color = ..level..), breaks = 1, color = "black", size = 3) + 
  stat_contour(aes(color = ..level..), breaks = 0.91, size = 2) + 
  annotate("text", x = 0.58, y = 0.505, label = "0.91", color = "#3399CC", size = 12.5) +
  geom_point(aes(x=0.564,y=0.49), color = "black", size = 5) +
  annotate("text", x = 0.574, y = 0.49, label = "B", size = 15) + 
  geom_point(aes(x=0.668,y=0.49), color = "black", size = 5) +
  annotate("text", x = 0.678, y = 0.49, label = "A", size = 15) +
  labs(y = expression(beta), fill = "Chernoff ratio") +
  coord_fixed() + theme_bw() + 
  theme(plot.title = element_text(size = 40),
        plot.subtitle = element_text(size = 38),
        axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 38),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 38),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38),
        legend.key.width = unit(2.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "top")

g <- direct.label(g, list("top.pieces", cex=3))



