
#' Function to run the collapsed Gibbs sampler on twitter data
#' 
#' @param alpha the prior for the user-topic probability vector theta (default set to 1)
#' @param eta the prior for the topic-word probability vector beta (default set to 1)
#' @param users a vector of the users (authors) of each document
#' @param dw a document-word matrix of counts
#' @param nT the desired number of topics, defaults to 2
#' @param niter the number of Gibbs itterations
#' @param nthin save Gibbs output every nthin iterations (defaults to 1)
#' @param nburnin number of iterations to go through before saving results (defaults to 0)
#' @param nupdate print iteration step every nupdate iterations (defaults to 10 prints)
#' @param seed Seed for randomization (optional)

collapsed_gibbs_1topic <- function(alpha=1,eta=1,
                                   users,dw,nT,niter,
                                   nthin=1,nburnin=0,nupdate=max(1,floor(niter/10)),seed=NULL){
  
  
  set.seed(seed)
  
  #extract dimensions
  nU <- users %>% unique %>% length
  nD <- nrow(dw)
  nW <- dim(dw)[2]
  
  #turn users into integers
  users <- users %>% match(unique(users))
  
  #initialize topic assignments
  ta <- sample(1:nT,replace = T,size = nD)
  
  # compute user-topic and topic-word counts
  ut <- matrix(NA,nrow=nU,ncol=nT)
  for(u in 1:nU){
    for(t in 1:nT){
      ut[u,t] <- sum(ta == t & users == u)
    }
  }
  
  tw <- matrix(NA,nrow=nT,ncol=nW)
  for(t in 1:nT){
    tw[t,] <- colSums(dw[ta == t,])
  }
  
  #identify which documents can be resampled faster (no duplicated words)
  no_dups <- apply(dw,1,function(d) max(d) == 1 & sum(d) > 1)
  
  #initilize results matrices
  nsaved <- floor(niter/nthin)
  results_ut <- array(NA,dim = c(nU,nT,nsaved))
  results_tw <- array(NA,dim = c(nT,nW,nsaved))
  results_ta <- array(NA,dim = c(nD,nsaved))
  
  iter_start <- Sys.time()
  for(iter in 1:(niter+nburnin)){
    
    #loop through tweets
    for(d in 1:nD){
      
      #get current status of doc d
      u <- users[d]
      t <- ta[d]
      
      #remove from current counts
      ut[u,t] <- ut[u,t] - 1
      tw[t,] <- tw[t,] - dw[d,]
      
      #compute probabilities
      u_to_t <- log((alpha + ut[u,])/(nT*alpha + sum(ut[u,])))
      if(no_dups[d]){
        w_to_t <- rowSums(log(tw[,dw[d,]>0] + eta)) + #assuming occurence only matters, not number of occurances
          lgamma(rowSums(tw + eta)) - lgamma(rowSums(tw + eta) + sum(dw[d,]))
      } else{
        w_to_t <- lgamma(rowSums(tw + eta)) - rowSums(lgamma(tw+eta)) +
          rowSums(lgamma(tw + eta + dw[rep(d,nT),])) - lgamma(rowSums(tw + eta) + sum(dw[d,]))
      }
      #w_to_t <- log(apply(X=(eta + tw)^rbind(dw[d,],dw[d,]),MARGIN=1,prod)/(rowSums(tw) + nW*eta)^sum(dw[d,]))
      probs <- (u_to_t + w_to_t) %>% {. - max(.)} %>% exp %>% {./sum(.)}
      
      #sample new topic
      t_new <- sample(1:nT,prob=probs,size=1)
      
      #add back to counts
      ta[d] <- t_new
      ut[u,t_new] <- ut[u,t_new] + 1
      tw[t_new,] <- tw[t_new,] + dw[d,]
      
    }
    
    #save results after burnin and every nthin iterations
    if(iter > nburnin && iter %% nthin == 0){
      results_ut[,,iter] <- ut/rowSums(ut)
      results_tw[,,iter] <- tw/rowSums(tw)
      results_ta[,iter] <- ta
    }    
    
    #print update
    if(iter %% nupdate == 0){
      print(str_c("Finished iteration ", iter," at ",Sys.time()))
      print(Sys.time() - iter_start)
      iter_start <- Sys.time()
    }
  }
  
  return(list(ut=results_ut,
              tw=results_tw,
              ta=results_ta))
}

#collapsed_gibbs_1topic(alpha = 1,eta = 1,users = du$user_id,dw = dw,nT = 2,niter = 10,seed=515)

#' Function to run the collapsed Gibbs sampler on twitter data with user-level clusters
#' 
#' @param alpha the prior for the user-topic probability vector theta (default set to 1)
#' @param eta the prior for the topic-word probability vector beta (default set to 1)
#' @param users a vector of the users (authors) of each document
#' @param dw a document-word matrix of counts
#' @param nT the desired number of topics, defaults to 2
#' @param nC the desired number of user clusters, defaults to 2
#' @param niter the number of Gibbs itterations
#' @param nthin save Gibbs output every nthin iterations (defaults to 1)
#' @param nburnin number of iterations to go through before saving results (defaults to 0)
#' @param nupdate print iteration step every nupdate iterations (defaults to 10 prints)
#' @param seed Seed for randomization (optional)
#' @param mcmc_update Option to update alphag parameters with an MCMC step (T) or with the MAP estimate (F)
#' @param nalphag_steps Number of updates to cluster-level paramters to compute between iterations through the corpus
#' @param mu_scale mean for half-normal prior on the scale/precision of each alphag Dirichlet parameters
#' @param sigma_scale standard deviation for half-normal prior on the scale/precision of each alphag Dirichelt parameters  
#' @param prop_scale_center for map_guess and current_center alphag_sample_method, the concnetration of the proposal Dirichlet distribution for the mcmc update
#' @param alphag_sample_method methods for mcmc updates to alphag
#' @param print_cluster option to print the number of users in each cluster after each print itteration
#' @param topic_presets option to pre-specify the topics of certain documents. Must be a vector of length nD with NA for unassigned documents 

collapsed_gibbs_1topic_clusters <- function(alpha=1,eta=1,nu=1,
                                            users,dw,nT,nC,niter,
                                            nthin=1,nburnin=0,nupdate=max(1,floor(niter/10)),seed=NULL,
                                            mcmc_update = T,
                                            nClusterIter = 100,
                                            mu_scale = 100,sigma_scale = 50,
                                            prop_scale_center = 100,
                                            alphag_sample_method = c("map_guess","current_center","componentwise")[3],
                                            print_clusters = F,
                                            topic_presets = NULL,cluster_presets = NULL,
                                            ta_true=NULL){
  
  
  set.seed(seed)
  #set prior for cluster centers
  alpha_center <- rep(alpha,nT)
  
  #extract dimensions
  nU <- users %>% unique %>% length
  nD <- nrow(dw)
  nW <- dim(dw)[2]
  
  #turn users into integers and store list of doc ids
  users <- users %>% match(unique(users))
  
  #initialize topic assignments
  ta <- sample(1:nT,replace = T,size = nD)
  if(!is.null(topic_presets)){
    ta[!is.na(topic_presets)] <- topic_presets[!is.na(topic_presets)]
  }
  
  #initialize user cluster assignments
  ca <- sample(1:nC,replace = T,size = nU)
  if(!is.null(cluster_presets)){
    ca[!is.na(cluster_presets)] <- cluster_presets[!is.na(cluster_presets)]
  }
  #ca <- rep(1:nC,times = rep(nU/nC,nC)) #
  #ca <- rep(1,nU)
  
  # compute user-topic count matrix
  ut <- matrix(NA,nrow=nU,ncol=nT)
  for(u in 1:nU){
    for(t in 1:nT){
      ut[u,t] <- sum(ta == t & users == u)
    }
  }
  
  # compute topic-word count matrix
  tw <- sapply(1:nT,function(t){
    if(sum(ta==t) == 0){
      return(rep(0,nW))
    }else if(sum(ta==t) == 1){
      return(dw[ta==t,])
    }else{
      return(colSums(dw[ta==t,]))  
    }
  }) %>% t
  
  #initialize alphag at within cluster topic counts + 1
  alphag <- 
    sapply(1:nC,function(c){
      if(sum(ca == c) == 1){
        ut[ca == c,] + 1
      }else{
        colSums(ut[ca == c,]) + 1
      }
    }) %>% t 
  
  #scale alphag to sum to a fixed size in each cluster
  alphag <- alphag/rowSums(alphag)*100
  
  #initilze phi at uniform
  phi <- rep(1,times = nC)/nC
  
  #identify which documents can be resampled faster (no duplicated words)
  no_dups <- apply(dw,1,function(d) max(d) == 1 & sum(d) > 1)
  
  #prepare results for saving
  nsaved <- floor(niter/nthin)
  results_ut <- array(NA,dim = c(nU,nT,nsaved))
  results_tw <- array(NA,dim = c(nT,nW,nsaved))
  results_ta <- array(NA,dim = c(nD,nsaved))
  
  results_ca <- array(NA,dim = c(nU,nsaved*nClusterIter))
  results_alphag <- array(NA,dim = c(dim(alphag),nsaved*nClusterIter))
  results_phi <- array(NA,dim = c(nC,nsaved*nClusterIter))
  
  results_iter <- 0
  results_iter_cparams <- 0
  
  #id documents to update 
  if(!is.null(topic_presets)){
    iter_docs <- which(is.na(topic_presets))
  } else{
    iter_docs <- 1:nD
  }
  
  iter_start <- Sys.time()
  for(iter in 1:(niter+nburnin)){
    
    #loop through tweets for new topics
    for(d in sample(iter_docs,size = length(iter_docs),replace = F)){
      
      #get current status of doc d
      u <- users[d]
      t <- ta[d]
      c <- ca[u]
      c_params <- alphag[c,]
      
      #remove from current counts
      ut[u,t] <- ut[u,t] - 1
      if(is.dfm(dw)){ #is.dfm(dw)){ #big dw matrices are passed as dfm objects
        dw_d <- convert(dw[d,],"matrix")
      } else{
        dw_d <- dw[d,]
      }
      tw[t,] <- tw[t,] - dw_d
      
      #compute probabilities
      u_to_t <- log((c_params + ut[u,])/(nT*c_params + sum(ut[u,])))
      if(no_dups[d]){
        w_to_t <- rowSums(log(tw[,dw[d,]>0] + eta)) + #assuming occurence only matters, not number of occurances
          lgamma(rowSums(tw + eta)) - lgamma(rowSums(tw + eta) + sum(dw[d,]))
      } else{
        w_to_t <- lgamma(rowSums(tw + eta)) - rowSums(lgamma(tw+eta)) +
          rowSums(lgamma(tw + eta + dw[rep(d,nT),])) - lgamma(rowSums(tw + eta) + sum(dw[d,]))
      }
      
      probs <- (u_to_t + w_to_t) %>% {. - max(.)} %>% exp %>% {./sum(.)}
      
      #sample new topic
      t_new <- sample(1:nT,prob=probs,size=1)
      
      #add back to counts
      ta[d] <- t_new
      ut[u,t_new] <- ut[u,t_new] + 1
      tw[t_new,] <- tw[t_new,] + dw_d
      
    }
    
    #loop through users for new clusters 
    ca <- sample_ca(nU,ut,nC,alphag,phi,cluster_presets)

    #update phi
    phi <- sample_phi(ca,nC,nu)
    
    if(mcmc_update == F){
      #set new alphag to within cluster topic counts + 1 (should make its own function)
      alphag <- sapply(1:nC,function(c){
        if(sum(ca == c) == 1){
          ut[ca == c,] + 1 #colSums breaks when given a 1 dimensional vector
        }else{
          colSums(ut[ca == c,]) + 1
        }
      }) %>% t
    }
    
    if(mcmc_update == T){    
      #compute metropolis step within the gibbs sampling
      #this step is fast, so repeat it many times to 
      
      for(i in 1:nClusterIter){
        alphag <- sample_alphag(ca,ut,alphag,
                                alpha_center,mu_scale,sigma_scale,prop_scale_center,
                                alphag_sample_method,
                                nC)
        
        #update ca and phi after each step as well
        ca <- sample_ca(nU,ut,nC,alphag,phi,cluster_presets)
        phi <- sample_phi(ca,nC,nu)
        
        if(iter > nburnin && iter %% nthin == 0){
          results_iter_cparams <- results_iter_cparams + 1
          
          results_ca[,results_iter_cparams] <- ca
          results_alphag[,,results_iter_cparams] <- alphag
          results_phi[,results_iter_cparams] <- phi
        }
      }
    }
    
    #save results after burnin and every nthin iterations
    if(iter > nburnin && iter %% nthin == 0){
      results_iter <- results_iter + 1
      
      results_ut[,,results_iter] <- ut/rowSums(ut)
      results_tw[,,results_iter] <- (tw+eta)/rowSums(tw+eta) 
      results_ta[,results_iter] <- ta
      
      if(mcmc_update == F){
        results_ca[,results_iter] <- ca
        results_alphag[,,results_iter] <- alphag
        results_phi[,results_iter] <- phi
      }
    }    
    
    #print update
    if(iter %% nupdate == 0){
      print(str_c("Finished iteration ", iter," at ",Sys.time()))
      print(Sys.time() - iter_start)
      iter_start <- Sys.time()
      if(print_clusters == T) {
        print(table(ca))
        if(!is.null(ta_true)){
          tab <- table(ta,ta_true)
          purity <- apply(tab,1,max)
          print(str_c("topic purity:" ,round(sum(purity)/nD,2)))
        }
      }
    }
  }
  
  return(list(ut=results_ut,
              tw=results_tw,
              ta=results_ta,
              ca=results_ca,
              alphag=results_alphag,
              phi=results_phi,
              words=colnames(dw)))
}




#' Function to run the collapsed Gibbs sampler on twitter data
#' 
#' @param alpha the prior for the user-topic probability vector theta (default set to 1)
#' @param eta the prior for the topic-word probability vector beta (default set to 1)
#' @param users a vector of the users (authors) of each document
#' @param dw a document-word matrix of counts
#' @param nT the desired number of topics, defaults to 2
#' @param niter the number of Gibbs itterations
#' @param nthin save Gibbs output every nthin iterations (defaults to 1)
#' @param nburnin number of iterations to go through before saving results (defaults to 0)
#' @param nupdate print iteration step every nupdate iterations (defaults to 10 prints)
#' @param seed Seed for randomization (optional)

full_gibbs_1topic <- function(alpha=1,eta=1,nu=1,
                              users,dw,nT,nC,niter,
                              nthin=1,nburnin=0,nClusterIter=1,
                              nupdate=max(1,floor(niter/10)),seed=NULL,print_clusters = T,
                              mu_scale = 100,sigma_scale = 50,
                              alphag_sample_method = c("map_guess","current_center","componentwise")[3]){
  
  
  set.seed(seed)
  #set prior for cluster centers
  alpha_center <- rep(alpha,nT)
  
  #extract dimensions
  nU <- users %>% unique %>% length
  nD <- nrow(dw)
  nW <- dim(dw)[2]
  
  #turn users into integers
  users <- users %>% match(unique(users))
  
  #initialize topic assignments (start at truth)
  ta <- dat$ta_true #sample(1:nT,replace = T,size = nD)
  
  #initialize cluster assignment
  ca <- sample(1:nC,size = nU,replace = T)
  
  # compute user-topic and topic-word counts
  ut <- matrix(NA,nrow=nU,ncol=nT)
  for(u in 1:nU){
    for(t in 1:nT){
      ut[u,t] <- sum(ta == t & users == u)
    }
  }
  
  tw <- tw <- sapply(1:nT,function(t){
    if(sum(ta==t) == 0){
      return(rep(0,nW))
    }else if(sum(ta==t) == 1){
      return(dw[ta==t,])
    }else{
      return(colSums(dw[ta==t,]))  
    }
  }) %>% t
  
  #initialize alphag at within cluster topic counts + 1
  alphag <- 
    sapply(1:nC,function(c){
      if(sum(ca == c) == 1){
        ut[ca == c,] + 1
      }else{
        colSums(ut[ca == c,]) + 1
      }
    }) %>% t 
  
  #scale alphag to sum to a fixed size in each cluster
  #alphag <- alphag/rowSums(alphag)*100
  
  #initilze phi at uniform
  phi <- rep(1,times = nC)/nC
  
  #initialize log-beta and log-theta values
  lbeta <- matrix(1/nW,nrow = nT,ncol = nW) %>% log
  ltheta <- matrix(1/nT,nrow = nU,ncol = nT) %>% log
  
  #identify which documents can be resampled faster (no duplicated words)
  no_dups <- apply(dw,1,function(d) max(d) == 1 & sum(d) > 1)
  
  #initilize results matrices
  nsaved <- floor(niter/nthin)
  
  #discrete parameters
  results_ut <- array(NA,dim = c(nU,nT,nsaved))
  results_tw <- array(NA,dim = c(nT,nW,nsaved))
  results_ta <- array(NA,dim = c(nD,nsaved))
  
  #probability parameters
  results_ltheta <- array(NA,dim = c(nU,nT,nsaved))
  results_lbeta <- array(NA,dim = c(nT,nW,nsaved))
  
  #cluster params
  results_ca <- array(NA,dim = c(nU,nsaved))
  results_alphag <- array(NA,dim = c(dim(alphag),nsaved))
  results_phi <- array(NA,dim = c(nC,nsaved))
  
  iter_start <- Sys.time()
  results_iter <- 0
  for(iter in 1:(niter+nburnin)){
    
    #update theta|ta
    ltheta <- apply(ut + alphag[ca,],1,function(r){
      rlogdir(1,alpha = r)
    }) %>% t
    
    #update beta|ta
    lbeta <- apply(tw + eta,1,function(t){
      rlogdir(1,alpha = t)
    }) %>% t
    
    #update ta|theta,beta
    # for (d in 1:nD) {
    #   lprobs <- ltheta[users[d],] + rowSums(t(t(lbeta)*dw[d,]))
    #   #probs <- (lprobs - max(lprobs)) %>% exp %>% {./sum(.)} #can let sample() do this for you
    #   ta[d] <- sample(1:nT,size = 1,prob = exp(lprobs-max(lprobs)))
    # }
    
    ta <- parallel::mclapply(1:nD, function(d){
      lprobs <- ltheta[users[d],] + rowSums(t(t(lbeta)*dw[d,]))
      #probs <- (lprobs - max(lprobs)) %>% exp %>% {./sum(.)} #can let sample() do this for you
      sample(1:nT,size = 1,prob = exp(lprobs-max(lprobs)))
    },mc.cores = parallel::detectCores(),mc.set.seed = seed) %>% unlist
    
    #recompute tw, ut
    ut <- table(c(users,1:nT + nU),c(ta,1:nT)) %>% {.[1:nU,]} 
    tw <- sapply(1:nT,function(t){
      if(sum(ta==t) == 0){
        return(rep(0,nW))
      }else if(sum(ta==t) == 1){
        return(dw[ta==t,])
      }else{
        return(colSums(dw[ta==t,]))  
      }
    }) %>% t

    #update ca,alphag,phi|theta nClusterIter times
    for(j in 1:nClusterIter){
      #ca
      ca <- sapply(1:nU,function(u){
        lprobs <- rep(NA,nC)
        for(c in 1:nC){
          lprobs[c] <- log(phi[c]) + lddir(exp(ltheta[u,]),alpha = alphag[c,])
        }
        if(sum(is.nan(lprobs %>% {.-max(.)} %>% exp))>0){
          print(lprobs)
          print(ltheta[u,])
          print(alphag)
          
          print("Numerical issue (small log probability)")
          
          lprobs[lprobs == Inf] <- 0
          lprobs[lprobs == -Inf] <- 10^-100
        }
        sample(1:nC,size = 1,prob = lprobs %>% {.-max(.)} %>% exp)
      })
      
      #phi
      phi <- sample_phi(ca,nC,nu)
      
      #alphag (hardest one)
      alphag <- sample_alphag.theta(ca,ltheta,alphag,
                                    alpha_center,mu_scale,sigma_scale,
                                    alphag_sample_method = "componentwise",
                                    nC)
    }
    
    #save results after burnin and every nthin iterations
    if(iter > nburnin && iter %% nthin == 0){
      results_iter <- results_iter + 1
      
      results_ut[,,results_iter] <- ut/rowSums(ut)
      results_tw[,,results_iter] <- (tw+eta)/rowSums(tw+eta) 
      results_ta[,results_iter] <- ta
      
      results_ltheta[,,results_iter] <- ltheta
      results_lbeta[,,results_iter] <- lbeta
      
      results_ca[,results_iter] <- ca
      results_alphag[,,results_iter] <- alphag
      results_phi[,results_iter] <- phi
      
    }
    
    # if(iter %% nupdate == 0){
    #   cat("\r",iter)
    # }
    #print update
    if(iter %% nupdate == 0){
      print(str_c("Finished iteration ", iter," at ",Sys.time()))
      print(Sys.time() - iter_start)
      iter_start <- Sys.time()
      if(print_clusters == T) {
        print(table(ca))
      }
    }
  }
  
  return(list(ut=results_ut,
              tw=results_tw,
              ta=results_ta,
              ca=results_ca,
              alphag=results_alphag,
              phi=results_phi,
              ltheta = results_ltheta,
              lbeta = results_lbeta,
              words=colnames(dw)))
}

