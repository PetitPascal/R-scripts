#------------------------------------------------------------------
#### Configurations  ####
#------------------------------------------------------------------

## Installing packages if needed
pack_needed<-c("tidyverse","blockrand")

for (i in 1:length(pack_needed)){
  if(pack_needed[i]%in%.packages(all.available=TRUE)){
    #ok
  }else{
    install.packages(pack_needed[i],repos = "http://cran.us.r-project.org") #installation package(s) manquant(s)
  }
}


## Loading packages
library(tidyverse)
library(blockrand)

## Setting work directory
here::here("Randomization scheme")

#------------------------------------------------------------------
#### Creating a randomized  ####
#------------------------------------------------------------------

## Setting a reproducible random seed generator
set.seed(1234)

## Randomization example: stratification by a binary factor, 15 in stratum, 3 treatments

# for the first level/category of the binary factor (i.e., factor_1)

random_tab1 <- blockrand(n=15,                   # minimum number of subjects to randomize
                         num.levels = 3,         # number of treatments or factor levels to randomize between
                         id.prefix='Factor1_',   # identifier (optional)
                         block.prefix='Factor1_',# block/group name (optional)
                         stratum='Factor1',      # character string specifying the stratum being generated (optional)
                         block.sizes = 5,        # sizes of blocks/groups to use (i.e., number of subjects per group/block)
                         levels=c("A","B", "C"), # character vector of labels for the different treatments or factor levels
                         uneq.beg=F,             # should an unequal block be used at the beginning of the randomization
                         uneq.mid=F,             # should an unequal block be used in the middle
                         uneq.min=0,             # minimum difference between the most and least common levels in an unequal block/group
                         uneq.maxit=10)          # maximum number of tries to get uneq.min difference

# for the second level/category of the binary factor (i.e., factor_2)

random_tab2 <- blockrand(n=15,                   # minimum number of subjects to randomize
                         num.levels = 3,         # number of treatments or factor levels to randomize between
                         id.prefix='Factor2_',   # identifier (optional)
                         block.prefix='Factor2_',# block/group name (optional)
                         stratum='Factor2',      # character string specifying the stratum being generated (optional)
                         block.sizes = 5,        # sizes of blocks/groups to use (i.e., number of subjects per group/block)
                         levels=c("A","B", "C"), # character vector of labels for the different treatments or factor levels
                         uneq.beg=F,             # should an unequal block be used at the beginning of the randomization
                         uneq.mid=F,             # should an unequal block be used in the middle
                         uneq.min=0,             # minimum difference between the most and least common levels in an unequal block/group
                         uneq.maxit=10)          # maximum number of tries to get uneq.min difference

## Final randomization table
as_tibble(bind_rows(random_tab1,random_tab2)) %>% select(-c(block.id,block.size))
