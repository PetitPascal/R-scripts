#-------------------------------------------------------------------------------
## Reproducible & generalizable example: continuous multioutput XGBoost
# - example 1: Predicting epidemic parameters (beta - transmission rate, mu - recovery rate, and S0 - initial susceptible population) from simulated epidemic curves
# - example 2: Predicting chemical concentration/exposure based on age, BMI, and physical activity level
# - example 3: Predicting blood pressure and cholesterol based on age, BMI, and physical activity level
# - use this script as a template for your own dataset by replacing the simulated data
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
here::here("XGBoost - multi-output")

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
np <- import("numpy")
pd <- import("pandas")
joblib <- import("joblib")
sklearn_kernels <- import("sklearn.gaussian_process.kernels")
sklearn_gp <- import("sklearn.gaussian_process")

#----------------------------------------------------------------
#### Creating two smiluated datasets ####
#----------------------------------------------------------------

#- - - - 
### Example 1 - epidemics (SIR simulation)

## Creating a dataset for training the model

n_obs <- 1000    # number of observations
t_max <- 30      # number of time points
time_cols <- paste0("t_", 0:t_max)

set.seed(123)

# Parameter grid
data_mod <- tibble(id   = paste0("St", rep(1:10, each = 10), "_", 1:n_obs),
                   beta = runif(n_obs, 0.5, 2.0),
                   mu   = runif(n_obs, 0.0, 0.3),
                   S0   = seq(80, 120, length.out = n_obs))

# Generating simple epidemic curves
curve_matrix <- replicate(t_max + 1,
                          round(dnorm(1:n_obs,
                                mean=runif(1, 20, 80),
                                sd=runif(1, 8, 20))*runif(1, 100, 500)))

curve_df <- as_tibble(curve_matrix)
names(curve_df) <- time_cols

data_train_SIR <- bind_cols(data_mod, curve_df)

## Creating a dataset for testing the model
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

#- - - - 
### Example 2 - Chemical exposure

set.seed(123)
n_expo <- 1000 # number of observations to simulate

data_expo_train <- tibble(id = paste0("P",1:n_expo),
                          age = sample(20:70,n_expo,replace=TRUE),
                          BMI = round(runif(n_expo,18,35),1),
                          activity = sample(0:10,n_expo,replace=TRUE))

data_expo_train <- data_expo_train %>%
  mutate(exposure_A = 0.5*age + 1.2*BMI - 0.8*activity + rnorm(n_expo,0,5),
         exposure_B = 2*BMI + 0.3*age + 0.5*activity + rnorm(n_expo,0,3))

data_expo_test <- tibble(id = paste0("P_test",1:10),
                         age = sample(20:70,10,replace=TRUE),
                         BMI = round(runif(10,18,35),1),
                         activity = sample(0:10,10,replace=TRUE))

#- - - - 
### Example 3 - Blood pressure and cholesterol

set.seed(123)
n_bp <- 1000  # number of observations to simulate

data_bp_train <- tibble(id = paste0("B",1:n_bp),
                        age = sample(20:70,n_bp,replace=TRUE),
                        BMI = round(runif(n_bp,18,35),1),
                        activity = sample(0:10,n_bp,replace=TRUE)) %>%
  mutate(systolic_BP  = 0.7*age + 1.5*BMI - 0.5*activity + rnorm(n_bp,0,8),
         diastolic_BP = 0.4*age + 1.0*BMI - 0.3*activity + rnorm(n_bp,0,5),
         cholesterol  = 0.3*age + 2*BMI - 0.2*activity + rnorm(n_bp,0,10))

data_bp_test <- tibble(id = paste0("B_test",1:10),
                       age = sample(20:70,10,replace=TRUE),
                       BMI = round(runif(10,18,35),1),
                       activity = sample(0:10,10,replace=TRUE))

#----------------------------------------------------------------
#### Creating the models ####
#----------------------------------------------------------------

## Preprocessing

# model 1
X_mod1 <- as.matrix(data_train_SIR %>% select(-c(id,beta,mu,S0)))
Y_mod1 <- as.matrix(data_train_SIR %>% select(beta,S0,mu))

# model 2
X_mod2 <- as.matrix(data_expo_train %>% select(age,BMI,activity))
Y_mod2 <- as.matrix(data_expo_train %>% select(exposure_A,exposure_B))

# model 3
X_mod3 <- as.matrix(data_bp_train %>% select(age,BMI,activity))
Y_mod3 <- as.matrix(data_bp_train %>% select(systolic_BP,diastolic_BP,cholesterol))

## Converting data to Python

# model 1
X_py_mod1 <- np$array(X_mod1)
Y_py_mod1 <- np$array(Y_mod1)

# model 2
X_py_mod2 <- np$array(X_mod2)
Y_py_mod2 <- np$array(Y_mod2)

# model 3
X_py_mod3 <- np$array(X_mod3)
Y_py_mod3 <- np$array(Y_mod3)

## Defining the models and fitting them with Python

py_run_string("
from sklearn.multioutput import MultiOutputRegressor
from xgboost import XGBRegressor

# Creating 3 models
model1 = MultiOutputRegressor(XGBRegressor(objective='reg:squarederror', n_estimators=500, 
             max_depth=5, 
             learning_rate=0.1,
             min_child_weight=10,
             colsample_bytree=0.8,
             subsample=0.8))
             
model2 = MultiOutputRegressor(XGBRegressor(objective='reg:squarederror', n_estimators=500, 
             max_depth=5, 
             learning_rate=0.1,
             min_child_weight=10,
             colsample_bytree=0.8,
             subsample=0.8))
                          
model3 = MultiOutputRegressor(XGBRegressor(objective='reg:squarederror', n_estimators=500, 
             max_depth=5, 
             learning_rate=0.1,
             min_child_weight=10,
             colsample_bytree=0.8,
             subsample=0.8))
")

## Training the models
mod1 <- py$model1$fit(X_py_mod1, Y_py_mod1) # model 1
mod2 <- py$model2$fit(X_py_mod2, Y_py_mod2) # model 2
mod3 <- py$model3$fit(X_py_mod3, Y_py_mod3) # model 3

#----------------------------------------------------------------
#### Testing the models ####
#----------------------------------------------------------------

## Converting test data to Python
test_mod1 <- r_to_py(as.matrix(data_test_SIR %>% select(-id))) # model 1
test_mod2 <- r_to_py(as.matrix(data_expo_test %>% select(-id))) # model 2
test_mod3 <- r_to_py(as.matrix(data_bp_test %>% select(-id))) # model 3

## Predicting

# model 1
prediction_mod1 <- mod1$predict(test_mod1)
colnames(prediction_mod1)<-c("beta","S0","mu")
results_mod1<-as_tibble(data.frame(prediction_mod1,
                                   I=data_test_SIR %>% select(-id)))
results_mod1

# model 2
prediction_mod2 <- mod2$predict(test_mod2)
colnames(prediction_mod2)<-c("exposure_A","exposure_B")

results_mod2<-as_tibble(data.frame(age=data_expo_test %>% select(age) %>% pull,
                                   BMI=data_expo_test %>% select(BMI) %>% pull,
                                   activity=data_expo_test %>% select(activity) %>% pull,
                                   prediction_mod2))
results_mod2

# model 3
prediction_mod3 <- mod3$predict(test_mod3)
colnames(prediction_mod3)<-c("systolic_BP","diastolic_BP","cholesterol")

results_mod3<-as_tibble(data.frame(age=data_bp_test %>% select(age) %>% pull,
                                   BMI=data_bp_test %>% select(BMI) %>% pull,
                                   activity=data_bp_test %>% select(activity) %>% pull,
                                   prediction_mod3))
results_mod3
