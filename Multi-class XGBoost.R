#-------------------------------------------------------------------------------
## Reproducible & generalizable example: multi-class XGBoost
# - example 1: Predicting epidemic cluster (based on epidemic curves)
# - example 2: Predicting chemical exposure cluster (based on age, BMI, physical activity)
# - example 3: Predicting cardiovascular risk cluster (based on age, BMI, physical activity, BP, cholesterol)
# use this script as a template for your own dataset by replacing the simulated data

#-------------------------------------------------------------------------------

#----------------------------------------------------------------
#### Configurations ####
#----------------------------------------------------------------

## Disabling memory torture
gctorture(FALSE)

## Installing packages if needed
pack_needed<-c("tidyverse","xgboost","reticulate","here","epimdr")
for (i in 1:length(pack_needed)){
  if(pack_needed[i]%in%.packages(all.available=TRUE)){
  }else{
    install.packages(pack_needed[i])
  }
}

## Package loading
library(tidyverse)
library(xgboost)
library(here)
library(reticulate)

## Setting the working directory
here::here("XGBoost - multiclass")

## Creating the conda 'r-reticulate' environment if it does not exist (to allow Python use with R)
if (!"r-reticulate" %in% conda_list()$name) {
  conda_create("r-reticulate")
}

# ## Installing the Python packages if they are not already installed
# conda_install("r-reticulate",
#               packages = c("scikit-learn", "xgboost", "numpy", "pandas", "joblib",
#                            "transformers", "datasets", "torch",
#                            "shap","matplotlib"),
#               pip = TRUE)

## Using the conda 'r-reticulate' environment
use_condaenv("r-reticulate", required = TRUE)

## Importing the Python packages
sklearn <- import("sklearn.multioutput")
xgb <- import("xgboost")
XGBClassifier <- import("xgboost")$XGBClassifier
np <- import("numpy")
pd <- import("pandas")
joblib <- import("joblib")
sklearn_kernels <- import("sklearn.gaussian_process.kernels")
sklearn_gp <- import("sklearn.gaussian_process")

#----------------------------------------------------------------
#### Creating two smiluated datasets ####
#----------------------------------------------------------------

#- - - - 
### Example 1 - epidemics

## Creating a dataset for training the model

n_obs <- 1000    # number of observations
t_max <- 30      # number of time points
time_cols <- paste0("t_", 0:t_max)

set.seed(123)

# Simulating epidemic curves
curve_matrix <- replicate(t_max + 1,
                          round(dnorm(1:n_obs,
                                mean=runif(1, 20, 80),
                                sd=runif(1, 8, 20))*runif(1, 100, 500)))

data_train_epidemic <- as_tibble(curve_matrix)
names(data_train_epidemic) <- time_cols

# Assigning cluster labels (e.g., 3 epidemic types)
data_train_epidemic$cluster <- sample(0:2, n_obs, replace=TRUE)  # 3 clusters: 0,1,2

## Test dataset
set.seed(123)

N_obs <- 10          # number of observations
t_max <- 30      # number of time points

data_test_SIR <- map_dfr(1:N_obs, function(i) {
  
  # Random initial conditions
  S <- numeric(t_max + 1)
  I <- numeric(t_max + 1)
  R <- numeric(t_max + 1)
  
  S[1] <- sample(800:1200, 1)
  I[1] <- sample(5:20, 1)
  R[1] <- 0
  Npop <- S[1] + I[1] + R[1]
  
  # Random parameters
  beta <- runif(1, 0.2, 1.0)
  mu   <- runif(1, 0.05, 0.3)
  
  # Stochastic SIR simulation
  for (t in 1:t_max) {
    # new infections cannot exceed S[t]
    new_inf <- rbinom(1, size = S[t], prob = 1 - exp(-beta * I[t] / Npop))
    # new recoveries cannot exceed I[t]
    new_rec <- rbinom(1, size = I[t], prob = 1 - exp(-mu))
    
    S[t+1] <- S[t] - new_inf
    I[t+1] <- I[t] + new_inf - new_rec
    R[t+1] <- R[t] + new_rec
  }
  
  tibble(
    id   = paste0("Obs", i),
    time = 0:t_max,
    S = S,
    I = I,
    R = R
  )
})

# ensuring no negative values
data_test_SIR_save<-data_test_SIR %>% mutate(time=if_else(time<0,0,time),
                                             S=if_else(S<0,0,S),
                                             I=if_else(I<0,0,I),
                                             R=if_else(R<0,0,R))

data_test_SIR<-data_test_SIR_save %>% select(id,I,time) %>% pivot_wider(names_from=time,values_from=I,names_prefix = "t_")

data_test_epidemic <- data_test_SIR %>% select(-id)

#- - - - 
### Example 2 - Chemical exposure

set.seed(123)
n_expo <- 1000 # number of observations to simulate

data_expo_train <- tibble(age = sample(20:70,n_expo,replace=TRUE),
                          BMI = round(runif(n_expo,18,35),1),
                          activity = sample(0:10,n_expo,replace=TRUE))

data_expo_train$cluster <- sample(0:2, n_expo, replace=TRUE) # 3 exposure clusters

data_expo_test <- tibble(age = sample(20:70,n_test,replace=TRUE),
                         BMI = round(runif(n_test,18,35),1),
                         activity = sample(0:10,n_test,replace=TRUE))

#- - - - 
### Example 3 - Cardiovascular risk

set.seed(123)
n_bp <- 1000  # number of observations to simulate

data_bp_train <- tibble(age = sample(20:70,n_bp,replace=TRUE),
                        BMI = round(runif(n_bp,18,35),1),
                        activity = sample(0:10,n_bp,replace=TRUE),
                        systolic_BP = round(runif(n_bp,100,160)),
                        diastolic_BP = round(runif(n_bp,60,100)),
                        cholesterol = round(runif(n_bp,150,250)))

data_bp_train$cluster <- sample(0:2, n_bp, replace=TRUE) # 3 risk clusters

data_bp_test <- tibble(age = sample(20:70,n_test,replace=TRUE),
                       BMI = round(runif(n_test,18,35),1),
                       activity = sample(0:10,n_test,replace=TRUE),
                       systolic_BP = round(runif(n_test,100,160)),
                       diastolic_BP = round(runif(n_test,60,100)),
                       cholesterol = round(runif(n_test,150,250)))

#----------------------------------------------------------------
#### Creating the models ####
#----------------------------------------------------------------

## Preparing Python arrays

# Example 1 - Epidemics
X_epidemic <- np$array(as.matrix(data_train_epidemic %>% select(-cluster)))
y_epidemic <- np$array(as.integer(data_train_epidemic$cluster))
X_test_epidemic <- np$array(as.matrix(data_test_epidemic))
n_class_epidemic <- as.integer(length(unique(data_train_epidemic$cluster)))

# Example 2 - Chemical exposure
X_expo <- np$array(as.matrix(data_expo_train %>% select(-cluster)))
y_expo <- np$array(as.integer(data_expo_train$cluster))
X_test_expo <- np$array(as.matrix(data_expo_test))
n_class_expo <- as.integer(length(unique(data_expo_train$cluster)))

# Example 3 - Cardiovascular risk
X_bp <- np$array(as.matrix(data_bp_train %>% select(-cluster)))
y_bp <- np$array(as.integer(data_bp_train$cluster))
X_test_bp <- np$array(as.matrix(data_bp_test))
n_class_bp <- as.integer(length(unique(data_bp_train$cluster)))

## Defining XGBClassifier models
py_run_string(sprintf("
from xgboost import XGBClassifier

model1 = XGBClassifier(objective='multi:softprob',
    num_class=%d,
    n_estimators=500,
    max_depth=5,
    learning_rate=0.1,
    min_child_weight=10,
    colsample_bytree=0.8,
    subsample=0.8)

model2 = XGBClassifier(objective='multi:softprob',
    num_class=%d,
    n_estimators=500,
    max_depth=5,
    learning_rate=0.1,
    min_child_weight=10,
    colsample_bytree=0.8,
    subsample=0.8)

model3 = XGBClassifier(objective='multi:softprob',
    num_class=%d,
    n_estimators=500,
    max_depth=5,
    learning_rate=0.1,
    min_child_weight=10,
    colsample_bytree=0.8,
    subsample=0.8)
", n_class_epidemic, n_class_expo, n_class_bp))

## Fitting the models
mod1 <- py$model1$fit(X_epidemic, y_epidemic) # model 1
mod2 <- py$model2$fit(X_expo, y_expo) # model 2
mod3 <- py$model3$fit(X_bp, y_bp) # model 3

## Predicting probabilities

# Model 1
pred_epidemic <- mod1$predict_proba(X_test_epidemic) 
colnames(pred_epidemic)<-paste("Cluster",1:ncol(pred_epidemic),sep="")
pred_epidemic<-pred_epidemic %>% as_tibble() %>% mutate(id=paste("obs ",1:nrow(pred_epidemic),sep=""))

cluster_prob<-pred_epidemic %>%
  pivot_longer(cols=-c(id)) %>% 
  group_by(id) %>% 
  filter(value==max(value,na.rm=T)) %>%
  ungroup %>%
  select(id,name) %>%
  rename(Most_likely_cluster=name)

left_join(pred_epidemic,cluster_prob,by="id") %>% 
  mutate(age=data_expo_test %>% select(age) %>% pull,
         BMI=data_expo_test %>% select(BMI) %>% pull,
         activity=data_expo_test %>% select(activity) %>% pull)

# Model 2
pred_expo<- mod2$predict_proba(X_test_expo)
colnames(pred_expo)<-paste("Cluster",1:ncol(pred_expo),sep="")
pred_expo<-pred_expo %>% as_tibble() %>% mutate(id=paste("obs ",1:nrow(pred_expo),sep=""))

cluster_prob<-pred_expo %>%
              pivot_longer(cols=-c(id)) %>% 
              group_by(id) %>% 
              filter(value==max(value,na.rm=T)) %>%
              ungroup %>%
              select(id,name) %>%
              rename(Most_likely_cluster=name)

as_tibble(data.frame(left_join(pred_expo,cluster_prob,by="id"),I=data_test_SIR %>% select(-id)))

# Model 3
pred_bp<- mod3$predict_proba(X_test_bp)
colnames(pred_bp)<-paste("Cluster",1:ncol(pred_bp),sep="")
pred_bp<-pred_bp %>% as_tibble() %>% mutate(id=paste("obs ",1:nrow(pred_bp),sep=""))

cluster_prob<-pred_bp %>%
  pivot_longer(cols=-c(id)) %>% 
  group_by(id) %>% 
  filter(value==max(value,na.rm=T)) %>%
  ungroup %>%
  select(id,name) %>%
  rename(Most_likely_cluster=name)

left_join(pred_bp,cluster_prob,by="id") %>% 
  mutate(age=data_bp_test %>% select(age) %>% pull,
         BMI=data_bp_test %>% select(BMI) %>% pull,
         activity=data_bp_test %>% select(activity) %>% pull)
