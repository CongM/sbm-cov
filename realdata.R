#############################################################################################
##                                                                                         ##
## Supplementary Code                                                                      ##
##                                                                                         ##
## Cong Mu, Angelo Mele, Lingxin Hao, Joshua Cape, Avanti Athreya, and Carey E. Priebe     ##
##                                                                                         ##
#############################################################################################




#### Change the path or working directory if needed
source("algorithm.R")




#### Real data experiments


## Connectome data
G <- seq(4, 40, 2)
dmax <- 30
dhat <- NULL
cov <- 2
betahatmethod <- "WA"
files <- dir("./DS72784-weighted-lcc/")

TTresults <- data.frame()
for (i in 1:length(files)) {
  cat(i, "\n")
  
  g <- read_graph(paste0("./DS72784-weighted-lcc/", files[i]), format = "graphml")
  g <- delete_edge_attr(g, "weight")
  A <- g[]
  LR <- V(g)$hemisphere
  GW <- V(g)$tissue
  blocksLR <- ifelse(LR == "left", 1, 2)
  blocksGW <- ifelse(GW == "gray", 1, 2)
  
  start_time1 <- Sys.time()
  result1 <- Algo1(A, cov, dhat, dmax, G)
  end_time1 <- Sys.time()
  Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))
  
  ARI_Algo1_LR <- adjustedRandIndex(blocksLR, result1$tauhat)
  ARI_Algo1_GW <- adjustedRandIndex(blocksGW, result1$tauhat)
  
  ARI_Algo2 <- c()
  Elapsed_Algo2 <- c()
  for (block in c("LR", "GW")) {
    if (block == "LR") {
      blocks <- ifelse(LR == "left", 1, 2)
      covariates <- matrix(ifelse(GW == "gray", 1, 2))
    } else {
      blocks <- ifelse(GW == "gray", 1, 2)
      covariates <- matrix(ifelse(LR == "left", 1, 2))
    }
    start_time2 <- Sys.time()
    result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
    end_time2 <- Sys.time()
    Elapsed_Algo2 <- c(Elapsed_Algo2, as.numeric(difftime(end_time2,start_time2,units="secs")))
    
    ARI_Algo2 <- c(ARI_Algo2, adjustedRandIndex(blocks, result2$tautilde))
  }
  
  ARI_Algo2_LR <- ARI_Algo2[1]
  ARI_Algo2_GW <- ARI_Algo2[2]
  Elapsed_Algo2_LR <- Elapsed_Algo2[1]
  Elapsed_Algo2_GW <- Elapsed_Algo2[2]
  
  TTresults <- rbind(TTresults, data.frame(ARI_Algo1_LR, ARI_Algo2_LR, ARI_Algo1_GW, ARI_Algo2_GW, Elapsed_Algo1, Elapsed_Algo2_LR, Elapsed_Algo2_GW))
}

save(TTresults, file = "TTresults.rdata")

dat <- TTresults %>%
  mutate(ARI12LR = ARI_Algo2_LR - ARI_Algo1_LR, ARI12GW = ARI_Algo2_GW - ARI_Algo1_GW)

pp <- ggplot(dat) + 
  geom_point(aes(x=ARI12LR, y=ARI12GW), alpha = 0.5, size = 5) + 
  geom_hline(yintercept = 0, color = 'red', size = 3) + geom_vline(xintercept = 0, color = 'red', size = 3) + 
  labs(x = 'ARI(Algo2,LR) - ARI(Algo1,LR)', y = 'ARI(Algo2,GW) - ARI(Algo1,GW)') + 
  xlim(-0.12, 0.12) + 
  ylim(-0.12, 0.12) + 
  coord_fixed() + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = 40),
        axis.text.x = element_text(size = 38),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 38),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))
pp <- ggMarginal(pp, type = 'histogram', alpha = 0.5)
pp


## Social network data
SNDresults <- data.frame()

# LastFM data
Dataset <- "LastFM"
raw <- fromJSON("./lasftm_asia/lastfm_asia_features.json")
covariates <- cbind(as.numeric(names(raw))+1,sapply(raw, length)) %>% as.data.frame()
names(covariates) <- c("id", "cov_original")
covariates <- covariates %>%
  arrange(id) %>%
  mutate(cov_discrete = ifelse(cov_original < 200, 1, 
                               ifelse(cov_original < 400, 2,
                                      ifelse(cov_original < 600, 3, 4)))) %>%
  dplyr::select(cov_discrete) %>%
  as.matrix()

el <- read.csv("./lasftm_asia/lastfm_asia_edges.csv") %>% 
  mutate(node_1 = node_1 + 1, node_2 = node_2 + 1) %>%
  as.matrix()
g <- graph_from_edgelist(el, directed = FALSE)
A <- g[]

blocks <- read.csv("./lasftm_asia/lastfm_asia_target.csv") %>%
  mutate(id = id + 1, target = target + 1)

cov <- length(unique(covariates[,1]))
dhat <- NULL
dmax <- 30
G <- prod(cov) * (2:20)
betahatmethod <- "WA"
set.seed(2020)

start_time1 <- Sys.time()
result1 <- Algo1(A, cov, dhat, dmax, G)
end_time1 <- Sys.time()
Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))

start_time2 <- Sys.time()
result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
end_time2 <- Sys.time()
Elapsed_Algo2 <- as.numeric(difftime(end_time2,start_time2,units="secs"))

ARI_Algo_1 <- adjustedRandIndex(blocks$target, result1$tauhat)
ARI_Algo_2 <- adjustedRandIndex(blocks$target, result2$tautilde)

SNDresults <- rbind(SNDresults, data.frame(Dataset, ARI_Algo_1, ARI_Algo_2, Elapsed_Algo1, Elapsed_Algo2))
save(SNDresults, file = "SNDresults.rdata")

# Facebook data
Dataset <- "Facebook"
raw <- fromJSON("./facebook_large/musae_facebook_features.json")
covariates <- cbind(as.numeric(names(raw))+1,sapply(raw, length)) %>% as.data.frame()
names(covariates) <- c("id", "cov_original")
covariates <- covariates %>%
  arrange(id) %>%
  mutate(cov_discrete = ifelse(cov_original < 15, 1, 2)) %>%
  dplyr::select(cov_discrete) %>%
  as.matrix()

el <- read.csv("./facebook_large/musae_facebook_edges.csv") %>% 
  mutate(id_1 = id_1 + 1, id_2 = id_2 + 1) %>%
  as.matrix()
g <- graph_from_edgelist(el, directed = FALSE)
A <- g[]

blocks <- read.csv("./facebook_large/musae_facebook_target.csv") %>%
  mutate(id = id + 1,
         block = ifelse(page_type == "company", 1, 
                        ifelse(page_type == "government", 2, 
                               ifelse(page_type == "politician", 3, 4))))

cov <- length(unique(covariates[,1]))
dhat <- NULL
dmax <- 30
G <- prod(cov) * (2:10)
betahatmethod <- "WA"
set.seed(2020)

start_time1 <- Sys.time()
result1 <- Algo1(A, cov, dhat, dmax, G)
end_time1 <- Sys.time()
Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))

start_time2 <- Sys.time()
result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
end_time2 <- Sys.time()
Elapsed_Algo2 <- as.numeric(difftime(end_time2,start_time2,units="secs"))

ARI_Algo_1 <- adjustedRandIndex(blocks$block, result1$tauhat)
ARI_Algo_2 <- adjustedRandIndex(blocks$block, result2$tautilde)

SNDresults <- rbind(SNDresults, data.frame(Dataset, ARI_Algo_1, ARI_Algo_2, Elapsed_Algo1, Elapsed_Algo2))
save(SNDresults, file = "SNDresults.rdata")

# GitHub data
Dataset <- "GitHub"
raw <- fromJSON("./git_web_ml/musae_git_features.json")
covariates <- cbind(as.numeric(names(raw))+1,sapply(raw, length)) %>% as.data.frame()
names(covariates) <- c("id", "cov_original")
covariates <- covariates %>%
  arrange(id) %>%
  mutate(cov_discrete = ifelse(cov_original < 18, 1, 2)) %>%
  dplyr::select(cov_discrete) %>%
  as.matrix()

el <- read.csv("./git_web_ml/musae_git_edges.csv") %>% 
  mutate(id_1 = id_1 + 1, id_2 = id_2 + 1) %>%
  as.matrix()
g <- graph_from_edgelist(el, directed = FALSE)
A <- g[]

blocks <- read.csv("./git_web_ml/musae_git_target.csv") %>%
  mutate(id = id + 1, block = ml_target + 1)

cov <- length(unique(covariates[,1]))
dhat <- NULL
dmax <- 30
G <- prod(cov) * (2:10)
betahatmethod <- "WA"
set.seed(2020)

start_time1 <- Sys.time()
result1 <- Algo1(A, cov, dhat, dmax, G)
end_time1 <- Sys.time()
Elapsed_Algo1 <- as.numeric(difftime(end_time1,start_time1,units="secs"))

start_time2 <- Sys.time()
result2 <- Algo2(A, covariates, betahatmethod, cov, dhat, dmax, G)
end_time2 <- Sys.time()
Elapsed_Algo2 <- as.numeric(difftime(end_time2,start_time2,units="secs"))

ARI_Algo_1 <- adjustedRandIndex(blocks$block, result1$tauhat)
ARI_Algo_2 <- adjustedRandIndex(blocks$block, result2$tautilde)

SNDresults <- rbind(SNDresults, data.frame(Dataset, ARI_Algo_1, ARI_Algo_2, Elapsed_Algo1, Elapsed_Algo2))
save(SNDresults, file = "SNDresults.rdata")



