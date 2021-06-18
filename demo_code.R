
library(topicmodels)
library(tidytext)
library(tidyverse)
library(quanteda)

source("scripts/setup.R")
source("scripts/helper_functions.R")
source("scripts/gibbs_functions.R")

### generate topic distributions ### 

nT <- 10
#install package with devtools: devtools::install_github("quanteda/quanteda.corpora")
corp_news <- quanteda.corpora::download('data_corpus_guardian')
dfmat_news <- dfm(corp_news, remove_punct = TRUE, remove = stopwords('en')) %>% 
  dfm_remove(c('*-time', '*-timeUpdated', 'GMT', 'BST')) %>% 
  dfm_trim(min_termfreq = 0.95, termfreq_type = "quantile", 
           max_docfreq = 0.1, docfreq_type = "prop")

dfmat_news <- dfmat_news[ntoken(dfmat_news) > 0,]
dtm <- convert(dfmat_news, to = "topicmodels")
lda <- LDA(dtm, k = 10,control = list(seed = 196))

nT <- lda@k

topics <- tidy(lda,matrix = "beta")

topics_tw <- topics %>% 
  group_by(topic) %>% 
  spread(key = topic,value = beta) 

words <- topics_tw$term
tw_true <- topics_tw[,2:(nT+1)] %>% t


### Generate clusters ###
nC <- 4
nUC <- 20
nU <- nUC*nC
nW <- ncol(tw_true)

set.seed(196) 

alpha_true <- matrix(c(1,1,1,1,1,1,1,1,1,1,
                       1,1,1,1,1,0,0,0,0,0,
                       0,0,0,0,0,1,1,1,1,1,
                       0,1,1,1,1,0,0,0,0,0
),nrow = nC,byrow = T)

alpha_true <- alpha_true %>% 
  {./rowSums(.)*100}

#to generate data with a single cluster, uncomment
#alpha_true <- matrix(1,nrow = nC,ncol = nT)

ca_true <- rep(1:nC,times = nUC) %>% sort
ut_true <- sapply(ca_true,function(c) rgamma(n = nT,shape = alpha_true[c,],rate = 1) %>% {./sum(.)}) %>% t

nDperU <- 40
users <- lapply(1:nU,function(u) rep(u,nDperU)) %>% unlist
ta_true <- lapply(users,function(u) sample(1:nT,size=1,prob = ut_true[u,])) %>% unlist

dw <- sapply(ta_true,function(t) rmultinom(n=1,size = 13,prob = tw_true[t,])) %>% t
ut_true_counts <- sapply(1:nU,function(u) sapply(1:nT,function(t) sum(ta_true == t & users ==u))) %>% t 

#sample stLDA-C, see gibbs_functions.R for documentation and parameter descriptions
groundtruth_estimate <- collapsed_gibbs_1topic_clusters(alpha = 1,eta = .1,nu = 1,
                                                        users = users,dw = dw,
                                                        nT = nT,nC = nC,
                                                        niter = 50,
                                                        seed = 196,mcmc_update = T,
                                                        nClusterIter = 100,
                                                        mu_scale = 0,sigma_scale = 100,
                                                        prop_scale_center = 100,alphag_sample_method = "componentwise",
                                                        print_clusters = T)

#sample stLDA, see gibbs_functions.R for documentation and parameter descriptions
# groundtruth_estimate_nocluster <- collapsed_gibbs_1topic(alpha = 1,eta = .1,
#                                                          users = users,dw = dw,
#                                                          nT = nT,
#                                                          niter = 100,
#                                                          seed = 555)

#save resultts
#save(groundtruth_estimate,users,dw,ta_true,ca_true,tw_true,words,file = "output/clda_sims/set1_cldac.Rdata")
#save(groundtruth_estimate_nocluster,users,dw,ta_true,ca_true,tw_true,words,file = "output/clda_sims/set1_clda_100runs.Rdata")

#######################
### Visualizations ####
#######################

#print top 5 words from each topic
groundtruth_estimate[["tw"]] %>% 
  top_topic_words(words = words,n=10) %>% 
  t

#print cluster means with user-level topic estimates overlayed
#grey bars are cluster-level expected values, colored lines are each user's topic distribution
#note that clusters with 1 user do not visualize well

ca_est <- groundtruth_estimate[["ca"]] %>% results_freq_table() %>% apply(1,which.max)
table(ca_est,ca_true)

plot_clusters <- function(ut_mat,cluster_assignment,cluster_alphas,yRange = c(0,.5)){
  cluster_means <- cluster_alphas %>% {./rowSums(.)}
  ut_mat <- ut_mat %>% {./rowSums(.)}
  
  lapply(unique(cluster_assignment),function(c){
    ut_mat %>% 
    {.[cluster_assignment == c,]} %>% 
      t %>% 
      data.frame(Topic = 1:ncol(ut_mat),.) %>% 
      reshape2::melt(id.vars = "Topic") %>% 
      ggplot(aes(x=Topic,y=value)) + 
      geom_line(aes(color=variable)) + 
      guides(color = "none") + 
      geom_bar(data = data.frame(x=1:ncol(ut_mat),y=cluster_means[c,]),aes(x=x,y=y),alpha=.5,stat = "identity") + 
      labs(title = str_c("Cluster ",c," (n=",sum(cluster_assignment == c),")"),y="Probability") + 
      ylim(yRange)
  })
}

clusterPlots <- plot_clusters(ut_mat = groundtruth_estimate[["ut"]] %>% results_array_mean(),
                                 cluster_assignment = groundtruth_estimate[["ca"]] %>% results_freq_table() %>% apply(1,which.max),
                                 cluster_alphas = groundtruth_estimate[["alphag"]] %>% results_array_mean())

clusterPlots %>% gridExtra::grid.arrange(grobs = .)

