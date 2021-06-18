## Sample data
# 
# nT <- 4
# nC <- 2
# nU <- 10
# nDU <- sample(10:20,size = nU)
# 
# alpha <- c(1,1,1)
# 
# alphag <- matrix(NA,nrow=nC,ncol=nT)

#subset a vector to values that appear > some number of times
above_cutoff <- function(x,cutoff){
  counts <- table(x)
  above_cutoff <- names(counts[counts>cutoff])
  x[x %in% above_cutoff]
}

# ddirmult function
lddirmult <- function(x,alpha){
  n <- sum(x)
  lc <- lgamma(n+1) + lgamma(sum(alpha)) - lgamma(n + sum(alpha))
  lkernel <- sum(lgamma(x + alpha) - lgamma(x + 1) - lgamma(alpha))
  lc + lkernel
}

rdir <- function(n,c=NULL,s=NULL,alpha=NULL){
  if(is.null(alpha)){
    g <- replicate(n,rgamma(length(c),shape = c*s,rate = 1)) %>% t
  } else{
    if(length(alpha)==1) return(matrix(1,nrow=n,ncol=1)) #if 1D Dirichlet, return 1
    g <- replicate(n,rgamma(length(alpha),shape = alpha,rate = 1)) %>% t
  }
  d <- g/rowSums(g)
  if(min(d) == 0){
    d[d==0] <- 10^-30
    d <- d/rowSums(d)
  } 
  return(d)
}

rlogdir <- function(n,c=NULL,s=NULL,alpha=NULL){
  if(is.null(alpha)){
    g <- replicate(n,rgamma(length(c),shape = c*s,rate = 1)) %>% t
  } else{
    if(length(alpha)==1) return(matrix(1,nrow=n,ncol=1)) #if 1D Dirichlet, return 1
    g <- replicate(n,rgamma(length(alpha),shape = alpha,rate = 1)) %>% t
  }
  d <- log(g) - log(rowSums(g))

  return(d)
}


lddir <- function(x,c=NULL,s=NULL,alpha=NULL){
  #for null arguments
  if(sum(x) == 0){
    return(0)
  }
  
  if(is.null(alpha)){
    alpha <- c*s
  }
  
  #to ensure lddir(0) = -Inf, pipe below to max(.,-Inf,na.rm=T)
  sum((alpha-1)*log(x)) + lgamma(sum(alpha)) - sum(lgamma(alpha))
}


#mean of a 3d array
results_array_mean <- function(a){
  dims <- dim(a)
  mean_array <- matrix(NA,nrow=dims[1],ncol=dims[2])
  
  for(d1 in 1:dims[1]){
    for(d2 in 1:dims[2]){
      mean_array[d1,d2] <- mean(a[d1,d2,])
    }
  }
  mean_array
}

#compute frequency table for a 1d results parameter
results_freq_table <- function(r){
  levels <- r %>% as.vector() %>% unique %>% sort
  sapply(1:nrow(r),function(id){
    sapply(levels,function(l){
      sum(r[id,] == l)
    })
  }) %>% t %>% {./ncol(r)}
}

#extract top words for each topic from a results vector
top_topic_words <- function(tresults,words,n=10){
  if(length(dim(tresults))==3){
    topic_means <- results_array_mean(tresults)
  } else{
    topic_means <- tresults
  }
  
  nT <- nrow(topic_means)
  lapply(1:nT, function(t){
    cutoff <- sort(topic_means[t,],decreasing = T)[n]
    words[topic_means[t,] >= cutoff]
  }) %>% do.call(what = rbind,args = .)
}

topn_words <- function(topic_dist,words,n=10){
  topic_dist %>% apply(1,function(tdist){
    cutoff <- quantile(tdist,probs = 1-n/length(tdist))
    topwords <- words[tdist >= cutoff]
    topwords <- topwords[order(tdist[tdist >= cutoff],decreasing = T)]
  })
}

write_topn_words <- function(topic_dist,words,TperRow,file,nW=5){
  top_words <- topn_words(topic_dist,words,n = nW)
  nT <- ncol(top_words)
  nrows <- ceiling(nT/TperRow)
  colnames(top_words) <- str_c("Topic ",1:nT)

  if(nrows == 1){
    xtable::xtable(top_words) %>% 
      xtable::print.xtable(floating = F,file = file,include.rownames = F)
  }
  
  if(nrows > 1){
    write(str_c("\\begin{tabular}{",str_c(rep("l",TperRow),collapse=""),"} \n"),file = file,append = F)
    write("\\hline ",file = file,append = T)
    write.table(top_words[,1:TperRow] ,file= file,
                quote = F,sep = " & ",eol = " \\\\ \n",row.names = F,append = T,col.names = T)
    write("\\hline ",file = file,append = T)
    #write remaining rows
    for(r in 2:nrows){
      next_table <- top_words[,((r-1)*TperRow + 1):min(nT,((r)*TperRow ))]
      write.table(next_table,file= file,
                  quote = F,sep = " & ",eol = " \\\\ \n",row.names = F,append = T,col.names = T)
      write("\\hline ",file = file,append = T)
    }
    
    write("\\end{tabular}",file=file,append = T)
  }
}

### clean and prepare twitter data for LDA ###
clean_twitter <- function(twitter_data,cutoff=0,stopwords = stopwords::stopwords()){
  clean_text <- twitter_data$text %>% 
    str_replace_all("\\n"," ") %>% 
    str_replace_all("[^[:alnum:] \\#\\@]","") %>% #keep hashtags and @ symbols, remove all other non-alphanumerics
    str_to_lower() %>% 
    str_squish() 
  
  stopwords_clean <- stopwords %>% str_remove_all("[^[:alnum:]]")
  
  word_list <-  str_c(clean_text,collapse = " ") %>% 
    str_split(" ") %>% {.[[1]]} %>% 
    above_cutoff(cutoff) %>% 
    unique %>%
    {.[!(. %in% stopwords_clean)]}
  
  dfm <- quanteda::dfm(clean_text)
  
  #drop words that did not parse into dfm (generally special characters or non-latin characters)
  word_list <- word_list[word_list %in% dfm@Dimnames$features]
  
  dfm_words <- dfm[,word_list]
  
  dw_mat <- map(1:nrow(dfm_words),~dfm_words[.,] %>% as.vector) %>% do.call(what="rbind",args=.) 
  
  users <- twitter_data$screen_name
  
  #drop documents with 0 words
  zero_word_docs <- rowSums(dw_mat) == 0
  dw_mat <- dw_mat[!zero_word_docs,]
  users <- users[!zero_word_docs]
  original_ids <- (1:nrow(twitter_data))[!zero_word_docs]
  
  #du <- data.frame(doc_id = 1:nrow(twitter_data),user_id = twitter_data$screen_name)
  return(list(dw=dw_mat,words=word_list,users=users,original_ids = original_ids))
}

#compute products for non-gamma versions of sampler
numerator_product <- function(count,posterior){
  if(count == 0) return(1)
  posterior:(posterior + count - 1)
}

numerator_product_vec <- function(count_vec,posterior_vec){
  lapply(1:length(count_vec),function(w) numerator_product(count_vec[w],posterior_vec[w]))
}

#resample cluster assignments

sample_ca <- function(nU,ut,nC,alphag,phi,cluster_presets){
  if(all(!is.na(cluster_presets)) & !is.null(cluster_presets)){
    return(cluster_presets)
  } else{
    ca_new <- sapply(1:nU, function(u){
      if(all(!is.na(cluster_presets[u])) & !is.null(cluster_presets)){
        return(cluster_presets[u])
      } else{
        z <- ut[u,]
        
        lprobs <- sapply(1:nC,function(c) lddirmult(z,alphag[c,]) + log(phi[c]))
        probs <- lprobs %>% {. - max(.)} %>% exp %>% normalize()
        
        c_new <- sample(1:nC,size = 1,prob = probs)
        return(c_new)
      }
    })
    return(ca_new)
  }

}

sample_phi <- function(ca,nC,nu){
  num_c <- sapply(1:nC,function(c) sum(ca == c))
  rgamma(nC,shape = num_c + nu,rate = 1) %>% normalize()
}

sample_alphag <- function(ca,ut,alphag,
                          alpha_center,mu_scale,sigma_scale,prop_scale_center,
                          alphag_sample_method,
                          nC){  
  #parallel::mclapply(1:nC,function(c){
  sapply(1:nC,function(c){
    
    #collect current topics for users in cluster c
    if(sum(ca == c) == 0){
      cluster_topics <- matrix(rep(0,ncol(ut)),nrow=1)
      center_for_prop <- rep(1,ncol(ut))/ncol(ut)
    } else if(sum(ca == c) == 1){
      cluster_topics <- matrix(ut[ca == c,],nrow = 1)
      center_for_prop <- (cluster_topics+1)/sum(cluster_topics+1)
    } else{
      cluster_topics <- ut[ca == c,]
      center_for_prop <- (cluster_topics + 1) %>% {./rowSums(.)} %>% colMeans()
    }
    
    #collect current parameters for cluster c
    alphag_old <- alphag[c,]
    
    ll_current <- apply(cluster_topics,1,FUN = lddirmult,alpha=alphag_old) %>% sum
    
    center_current <- alphag_old %>% {./sum(.)}
    scale_current <- alphag_old %>% sum
    prior_center <- lddir(center_current,alpha = alpha_center)
    prior_scale <- dnorm(scale_current,mean = mu_scale,sd=sigma_scale,log = T)
    ll_current <- ll_current + prior_center + prior_scale
    #sum(dnorm(alphag_old,mean = 0,sd = 100,log = T))
    
    if(alphag_sample_method %in% c("map_guess","current_center")){
      #propose new center from close to the current center
      if(alphag_sample_method == "map_guess"){
        center_prop <- rdir(1,alpha=prop_scale_center*center_for_prop)
        ll_prop <- apply(cluster_topics,1,FUN = lddirmult,alpha=center_prop*scale_current) %>% sum
        prior_center_prop <- lddir(center_prop,alpha = alpha_center)
        
        ll_prop_full <- ll_prop + prior_center_prop + prior_scale
        
        accept_prob <- ll_prop_full - ll_current + 
          lddir(center_current,c=center_for_prop,s=prop_scale_center) - 
          lddir(center_prop,c=center_for_prop,s=prop_scale_center)
      } else if(alphag_sample_method == "current_center"){
        center_prop <- rdir(1,alpha=prop_scale_center*center_current)
        ll_prop <- apply(cluster_topics,1,FUN = lddirmult,alpha=center_prop*scale_current) %>% sum
        prior_center_prop <- lddir(center_prop,alpha = alpha_center)
        
        ll_prop_full <- ll_prop + prior_center_prop + prior_scale
        
        accept_prob <- ll_prop_full - ll_current + 
          lddir(center_current,c=center_prop,s=prop_scale_center) - 
          lddir(center_prop,c=center_current,s=prop_scale_center)
      }
      if(runif(1)<exp(accept_prob)){
        center_current <- center_prop
        ll_current <- ll_prop_full
        prior_center <- prior_center_prop
      }
      
      #next update the scale
      scale_prop <- (scale_current + runif(1,-2,2)) %>% abs
      ll_prop <- apply(cluster_topics,1,FUN = lddirmult,alpha=center_current*scale_prop) %>% sum
      prior_scale_prop <- dnorm(scale_prop,mu_scale,sigma_scale,log = T)
      ll_prop_full <- ll_prop + prior_scale_prop + prior_center
      
      accept_prob <- ll_prop_full - ll_current
      
      if(runif(1)<exp(accept_prob)){
        scale_current <- scale_prop
      }
      
      return(center_current*scale_current)
      
    } else if(alphag_sample_method == "componentwise"){

      alpha_final <- alphag_old
      i <- 1
      for(t in sample(ncol(ut))){
        alpha_prop <- alpha_final
        alpha_prop[t] <-  abs(alpha_final[t] + sample(runif(3,c(-2,-10,-100),c(2,10,100)),prob = c(.60,.30,.10),size = 1))
        
        scale_prop <- sum(alpha_prop)
        center_prop <- alpha_prop/scale_prop
        
        ll_prop <- apply(cluster_topics,1,FUN = lddirmult,alpha=alpha_prop) %>% sum 
        ll_prop_full <- ll_prop + dnorm(scale_prop,mu_scale,sigma_scale,log=T) + lddir(center_prop,alpha = alpha_center)
        
        accept_prob <- ll_prop_full - ll_current
        if(runif(1) < exp(accept_prob)){
          alpha_final <- alpha_prop
          ll_current <- ll_prop_full
        }
      }
    
      return(alpha_final)
  }
    
  }) %>% t #do.call(what = "rbind",args = .)  
}

sample_alphag.theta <- function(ca,ltheta,alphag,
                                alpha_center,mu_scale,sigma_scale,
                                alphag_sample_method = "componentwise",
                                nC){  
  
  #parallel::mclapply(1:nC,function(c){
  sapply(1:nC,function(c){
    
    #collect current topics for users in cluster c in matrix
    if(sum(ca == c) == 0){
      cluster_topics <- matrix(rep(0,ncol(ltheta)),nrow=1)
    } else if(sum(ca == c) == 1){
      cluster_topics <- matrix(ltheta[ca == c,],nrow = 1) %>% exp
    } else{
      cluster_topics <- ltheta[ca == c,] %>% exp
    }
    
    #collect current parameters for cluster c
    alphag_old <- alphag[c,]
    
    ll_current <- apply(cluster_topics,1,FUN = lddir,alpha=alphag_old) %>% sum
    
    center_current <- alphag_old %>% {./sum(.)}
    scale_current <- alphag_old %>% sum
    prior_center <- lddir(center_current,alpha = alpha_center)
    prior_scale <- dnorm(scale_current,mean = mu_scale,sd=sigma_scale,log = T)
    ll_current <- ll_current + prior_center + prior_scale

    #only doing component-wise updates for now, seems to have best mixing
    alpha_final <- alphag_old
    for(t in sample(ncol(ltheta))){
      alpha_prop <- alpha_final
      alpha_prop[t] <-  abs(alpha_final[t] + sample(runif(3,c(-2,-10,-100),c(2,10,100)),prob = c(.60,.30,.10),size = 1))
      
      scale_prop <- sum(alpha_prop)
      center_prop <- alpha_prop/scale_prop
      
      ll_prop <- apply(cluster_topics,1,FUN = lddir,alpha=alpha_prop) %>% sum 
      ll_prop_full <- ll_prop + dnorm(scale_prop,mu_scale,sigma_scale,log=T) + lddir(center_prop,alpha = alpha_center)
      
      accept_prob <- ll_prop_full - ll_current
      if(runif(1) < exp(accept_prob)){
        alpha_final <- alpha_prop
        ll_current <- ll_prop_full
      }
    }
    
    return(alpha_final)
  
    
  }) %>% t
  
}
##########################
### Cluster Evaluation ###
##########################

nmi <- function(tab){
  N <- sum(tab)
  rowProb <- rowSums(tab)/N
  colProb <- colSums(tab)/N
  
  mutInf <- (tab/N)*log(diag(1/rowProb) %*% (tab/N) %*% diag(1/colProb))
  mutInf <- sum(mutInf,na.rm = T)
  
  
  entRow <- sum(-rowProb*log(rowProb))
  entCol <- sum(-colProb*log(colProb))
  
  mutInf/(entRow + entCol)*2
}


################################
### Simulation Study Helpers ###
################################

#simulate data from stlda-C with fixed topic distribuions

make_data <- function(nW,nT,topic_low = 0,topic_high = 10,
                      nC = 3,alpha_true = alphag,
                      nUC=10,nDperU = 50,seed = 515){
  set.seed(seed) 
  
  # simulate topic distributions
  tw_true <- matrix(topic_low,nrow=nT,ncol = nW)
  for(t in 1:nT){
    tw_true[t,((nW/nT)*(t-1) + 1):((nW/nT)*(t))] <- topic_high
  }
  tw_true <- tw_true/rowSums(tw_true)
  
  ### Generate clusters ###
  nU <- nUC*nC
  nW <- ncol(tw_true)
  
  alpha_true <- alpha_true[1:nC,] #for testing, sometimes the hard cluster is dropped
  
  alpha_true <- alpha_true %>% 
    {./rowSums(.)*100} %>% {.+1}
  
  #to generate data with a single cluster, uncomment
  #alpha_true <- matrix(1,nrow = nC,ncol = nT)
  
  #make users
  ca_true <- rep(1:nC,times = nUC) %>% sort
  ut_true <- sapply(ca_true,function(c) rgamma(n = nT,shape = alpha_true[c,],rate = 1) %>% {./sum(.)}) %>% t
  
  #make posts
  #nDperU <- 50
  users <- lapply(1:nU,function(u) rep(u,nDperU)) %>% unlist
  ta_true <- lapply(users,function(u) sample(1:nT,size=1,prob = ut_true[u,])) %>% unlist
  dw <- sapply(ta_true,function(t) rmultinom(n=1,size = 15,prob = tw_true[t,])) %>% t
  ut_true_counts <- sapply(1:nU,function(u) sapply(1:nT,function(t) sum(ta_true == t & users ==u))) %>% t 
  
  return(list(tw_true = tw_true,alpha_true = alpha_true,ca_true = ca_true,
              ut_true = ut_true,users = users,ta_true = ta_true,dw = dw,ut_true_counts = ut_true_counts,
              nT = nT,nC = nC,nW = nW,nU = nUC*nC))
}
