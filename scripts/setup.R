### set up ###

################
### Packages ###
################

require(tidyverse)

#################
### Functions ###
#################

standardize <- function(x){
  (x-mean(x))/sd(x)
}

normalize <- function(x){
  x/sum(x)
}

scale_01 <- function(x){
  y <- x - min(x)
  y/max(y)
}

has_na <- function(x){
  sum(is.na(x))>0
}

#return TRUE for every row that has a duplicate (including the first one)
duplicated_plus_original <- function(data){
  dups_topdown <- duplicated(data)
  dups_bottomup <- duplicated(data[nrow(data):1,])[nrow(data):1]
  pmax(dups_topdown,dups_bottomup) == 1
}

#return TRUE if the element in string matches any of the elements in pattern
str_detect_any <- function(string,pattern){
  results <- map(pattern,function(x) str_detect(string,x)) 
  results <- do.call(rbind,results) %>% apply(MARGIN = 2,max)
  results == 1
}
