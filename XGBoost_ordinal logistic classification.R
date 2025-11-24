#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing packages if needed
pack_needed<-c("data.table","tidyverse","mllrnrs","broom","doParallel","foreach",
               "splitTools","conflicted","grid","gridExtra","RColorBrewer","mlbench",
               "mlexperiments","caret","MLmetrics","patchwork",
               "xgboost","parallel","here","irr")
for (i in 1:length(pack_needed)){
  if(pack_needed[i]%in%.packages(all.available=TRUE)){
  }else{
    install.packages(pack_needed[i])
  }
}

## Package loading
library(tidyverse)
library(data.table)
library(broom)
library(doParallel)
library(foreach)
library(pROC)
library(grid)
library(gridExtra)
library(conflicted)
library(RColorBrewer)
library(mlexperiments)
library(mllrnrs)
library(mlbench)
library(caret)
library(MLmetrics)
library(xgboost)
library(patchwork)
library(here)

## Preventing package conflicts
conflict_prefer("select", "dplyr") 
conflict_prefer("filter", "dplyr") 
conflict_prefer("slice", "dplyr")
conflict_prefer("alpha", "scales")
conflict_prefer("auc", "pROC")
conflict_prefer("roc", "pROC")

## Setting the working directory
here::here("XGBoost - Ordinal logistic classification - nested CV")

## Setting seed
seed <- 123

## Setting cores
if (isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_")))) {
  ncores <- 2L
} else {
  ncores <- ifelse(
    test = parallel::detectCores() > 4,
    yes = 4L,
    no = ifelse(
      test = parallel::detectCores() < 2L,
      yes = 1L,
      no = parallel::detectCores()
    )
  )
}

## Setting mlexperiments package options
options("mlexperiments.bayesian.max_init" = 10L)
options("mlexperiments.optim.xgb.nrounds" = 100L)
options("mlexperiments.optim.xgb.early_stopping_rounds" = 10L)

#----------------------------------------------------------------
#### Creating functions ####

## Function for plotting SHAP values
plot.shap.summary <- function(data_long){
  x_bound <- max(abs(data_long$value))
  require('ggforce')
  plot1 <- ggplot(data = data_long)+
    coord_flip() + 
    geom_sina(aes(x = variable, y = value, color = stdfvalue)) +
    geom_text(data = unique(data_long[, c("variable", "mean_value"), with = F]),
              aes(x = variable, y=-Inf, label = sprintf("%.3f", mean_value)),
              size = 3, alpha = 0.7,
              hjust = -0.2, 
              fontface = "bold") + # bold
    scale_color_gradient(low="#FFCC33", high="#6600CC", 
                         breaks=c(0,1), labels=c("Low","High")) +
    theme_bw() + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position="bottom") + 
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-x_bound, x_bound)) +
    scale_x_discrete(limits = rev(levels(data_long$variable)) 
    ) + 
    labs(y = "SHAP value (impact on model output)", x = "", color = "Feature value") 
  return(plot1)
}

## Function for standardizing feature values into the same range
std1 <- function(x){
  return ((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}

## Function for formatting the display of summary statistics
test_format<-function(x){
  x<-as.numeric(x)
  sign_x<-if_else(x<0,"neg","pos")
  x<-abs(x)
  x_raw<-x
  
  if(is.na(x)|x==""|x=="NA"|is.infinite(x)){
    x_raw<-0
  }
  
  if(x_raw>=100){
    virg_pos<-str_locate(as.character(x_raw),"[.]")[1]
    
    if(!is.na(virg_pos)&as.numeric(substr(x_raw,virg_pos+1,virg_pos+1))>=5){
      x<-x+1
      x<-as.numeric(substr(x,1,virg_pos-1))
    }
    
  }
  
  if(is.na(x)|x==""|x=="NA"|is.infinite(x)){
    x<-""
  }else{
    
    if(x<0.01|x>=10000){
      x<-format(signif(x,3),scientific = T)
      
      if(nchar(x)==8&substr(x,4,4)!="e"){
        if(as.numeric(substr(x,4,4))>=5&x<0.01){
          
          x1<-as.numeric(substr(x,1,4))+0.1
          x2<-substr(x,5,8)
          x<-paste(substr(x1,1,4),x2,sep="")
        }else{
          x<-paste(substr(x,1,4),substr(x,5,8),sep="")
        }
      }else{
        x<-x
      }
      
    }else{
      
      x_save<-x
      x<-signif(x,3)
      
      if(nchar(x)==6){
        x<-as.numeric(substr(x,1,5))
      }
      
      if(nchar(x)==5){
        if(as.numeric(substr(x,5,5))>=5){
          x<-x+0.01
          x<-substr(x,1,4)
        }else{
          x<-substr(x,1,4)
        }
      }else{
        if(x>=1000){
          x<-as.character(signif(x_save,4))
        }else{
          x<-as.character(x)
        }
      }
    }
    
  }
  
  if(sign_x=="neg"&x_raw!=0){
    x<-paste("-",x,sep="")
  }
  
  if(x=="0e+00"){
    x<-"0"
  }
  
  return(x)
}

# Function for determining the direction (i.e.: how a feature impacts the model)
determine_direction <- function(df, feature_col, shap_col,
                                majority_threshold = 0.55,
                                n_bins = 200) {
  
  x <- df[[feature_col]]
  s <- df[[shap_col]]
  
  valid <- complete.cases(x, s)
  x <- x[valid]; s <- s[valid]
  
  if (length(unique(x)) <= 1 || length(unique(s)) <= 1)
    return("neutral")
  
  #---------------------------
  # BINARY / SPARSE FEATURES
  #---------------------------
  if (length(unique(x)) <= 2 || quantile(x, 0.75) == 0) {
    
    group0 <- s[x == min(x)]
    group1 <- s[x == max(x)]
    
    n0 <- length(group0); n1 <- length(group1)
    if (n0 == 0 || n1 == 0) return("neutral")
    
    group0_sorted <- sort(group0)
    
    pos_count <- sum(findInterval(group1, group0_sorted, left.open = TRUE))
    n_total <- as.double(n0) * as.double(n1)
    
    prop <- pos_count / n_total
    if (prop >= majority_threshold) return("promoting predictor")
    if (prop <= (1 - majority_threshold)) return("mitigating predictor")
    return("neutral")
  }
  
  #---------------------------
  # CONTINUOUS FEATURES
  #---------------------------
  
  bins <- cut(x, breaks = n_bins, include.lowest = TRUE)
  
  bx <- tapply(x, bins, mean)
  bs <- tapply(s, bins, mean)
  bc <- as.numeric(tapply(s, bins, length))
  
  keep <- !is.na(bx)
  bx <- bx[keep]
  bs <- bs[keep]
  bc <- bc[keep]
  
  B <- length(bx)
  if (B < 2) return("neutral")
  
  # Vectorized pairwise comparisons
  bx_mat_i <- matrix(bx, B, B)
  bx_mat_j <- t(bx_mat_i)
  
  bs_mat_i <- matrix(bs, B, B)
  bs_mat_j <- t(bs_mat_i)
  
  # Determining which bin has higher feature value
  higher_x <- bx_mat_i > bx_mat_j
  
  # Counting sample pairs
  bc_mat_i <- matrix(as.numeric(bc), B, B)
  bc_mat_j <- t(bc_mat_i)
  n_pairs  <- bc_mat_i * bc_mat_j
  
  total_pairs <- sum(n_pairs[higher_x])
  
  if (total_pairs == 0) return("neutral")
  
  pos_pairs <- sum(n_pairs[higher_x & (bs_mat_i > bs_mat_j)])
  neg_pairs <- sum(n_pairs[higher_x & (bs_mat_i < bs_mat_j)])
  
  prop_pos <- pos_pairs / total_pairs
  prop_neg <- neg_pairs / total_pairs
  
  if (prop_pos >= majority_threshold) return("promoting predictor")
  if (prop_neg >= majority_threshold) return("mitigating predictor")
  return("neutral")
}

## not including function (opposite function of %in%)
`%ni%`<-Negate('%in%')

#----------------------------------------------------------------
#### Data import and preprocessing ####

# Simulating a dataset
n <- 500
df <- data.frame(
  Age = rnorm(n, 50, 10),
  BMI = rnorm(n, 25, 5),
  Gender = factor(sample(c(0, 1), n, replace = TRUE)),
  Treatment = factor(sample(c("A", "B", "C"), n, replace = TRUE)))

# Ordinal outcome: 0,1,2
df$outcome <- factor(sample(0:2, n, replace = TRUE), ordered = TRUE)

## Creating the clean dataset
clean_data<-df
clean_data_save<-clean_data

## Creating a table with the variable type
var_type_tab<-as_tibble(data.frame(feature=colnames(clean_data),
                                   Var_type=c("continuous","continuous","binary","categorical","ordered")))

## Reorganising the dataset to ensure the outcome is the first column
clean_data<-clean_data %>% select(outcome,colnames(clean_data)[which(colnames(clean_data) %ni% c("outcome"))])

#--------------------------
## Transforming the clean dataset into data.table
dataset <- clean_data |> data.table::as.data.table()

## Creating vectors with the independent and dependent variable names
feature_cols <- colnames(dataset)[2:ncol(dataset)]
target_col <- "outcome"

## Identifying unordered factors
categorical_var<-var_type_tab %>% filter(feature %in% colnames(clean_data)) %>% 
  filter(Var_type %in% c("categorical")) %>% select(feature) %>% pull

## One-hot encoding of unordered factors
if(length(categorical_var) > 0) {
  dummies <- model.matrix(~ . -1, data = dataset[, ..categorical_var])
  dataset <- cbind(dataset[, !..categorical_var, with = FALSE], dummies)
}

## Updating feature columns after encoding
feature_cols <- setdiff(colnames(dataset), "outcome")

# Converting outcome to 0-based integers
dataset$outcome <- as.numeric(as.character(dataset$outcome))
num_class <- length(unique(dataset$outcome))

## Creating the stratified 70/30 train-test split
set.seed(seed)
data_split <- splitTools::partition(
  y = dataset[[target_col]],
  p = c(train = 0.7, test = 0.3),
  type = "stratified",
  seed = seed)

## Preparing training and test datasets

# train dataset
X_train <- as.matrix(dataset[data_split$train, ..feature_cols])
X_train <- apply(X_train, 2, as.numeric)  # convert each column to numeric
y_train <- as.integer(dataset[data_split$train, get(target_col)])

# test dataset
X_test  <- as.matrix(dataset[data_split$test, ..feature_cols])
X_test <- apply(X_test, 2, as.numeric)
y_test  <- as.integer(dataset[data_split$test, get(target_col)])

#----------------------------------------------------------------
#### 10-fold nested cross-validation (CV) ####

#--------------------------
## Performing the outer CV on training set

# Creating the outer folds
outer_folds <- splitTools::create_folds(y_train, k = 10, type = "stratified", seed = seed)

# Creating empty lists to store results
outer_results <- list()
best_params_all <- list()

for (outer_idx in seq_along(outer_folds)) { # for each outer fold
  
  # selecting an outer fold
  val_idx <- outer_folds[[outer_idx]]
  train_idx_cv <- setdiff(seq_len(nrow(X_train)), val_idx)
  
  # creating training and test data sets for a given outer fold
  X_tr <- X_train[train_idx_cv, , drop = FALSE]
  y_tr <- y_train[train_idx_cv]
  X_val <- X_train[val_idx, , drop = FALSE]
  y_val <- y_train[val_idx]

  #--------------------------
  ## Performing the inner CV for tuning
  
  # Creating the inner folds
  inner_folds <- splitTools::create_folds(y_tr, k = 10, type = "stratified", seed = seed)
  
  # Creating random parameter space for hyperparameter tuning
  param_list <- expand.grid(
    subsample = seq(0.5, 1, 0.25),
    colsample_bytree = seq(0.5, 1, 0.25),
    min_child_weight = c(1, 5, 10),
    learning_rate = c(0.05, 0.1, 0.3),
    max_depth = c(3, 5, 7)) %>% dplyr::slice_sample(n = 30, replace = TRUE) # Limiting space to 30 combinations for computational efficiency and environmental sustainability considerations
  
  # Initializing best parameters
  best_auc <- -Inf
  best_params <- NULL
  
  for (i in 1:nrow(param_list)) {

    # Initializing cross validation
    xgb_cv <- mlexperiments::MLCrossValidation$new(
      learner = mllrnrs::LearnerXgboost$new(metric_optimization_higher_better = TRUE),
      fold_list = inner_folds,
      ncores = 2,
      seed = 0)

    # Defining the learner arguments
    xgb_cv$learner_args <- c(
      as.list(param_list[i, ]),
      list(objective = "multi:softprob", eval_metric = "mlogloss", nrounds = 100L,num_class=num_class))
    xgb_cv$performance_metric <- mlexperiments::metric("mae")
    xgb_cv$set_data(x = X_tr, y = y_tr)
    
    # Executing CV tuning
    res_cv <- xgb_cv$execute()
    
    # Evaluating using MAE (ordinal-aware)
    fold_mae <- sapply(res_cv$results$folds, function(f) {
      pred <- max.col(f[[5]]$predictions) - 1
      mean(abs(pred - f[[5]]$labels))
    })
    mean_mae <- mean(fold_mae)
    
    if (!is.na(mean_mae) && mean_mae > best_score) {  # minimizing MAE
      best_score <- mean_mae
      best_params <- params
    }
  }
  
  # If no params are selected, fallback to default
  if (is.null(best_params)) {
    best_params <- data.frame(
      subsample = 1,
      colsample_bytree = 1,
      min_child_weight = 1,
      learning_rate = 0.1,
      max_depth = 3)
  }
  
  best_params_all[[outer_idx]] <- best_params
  
  ## Evaluating on outer fold
  
  # Building the model
  dtrain <- xgboost::xgb.DMatrix(data = X_tr, label = y_tr)
  dval   <- xgboost::xgb.DMatrix(data = X_val, label = y_val)
  
  xgb_outer <- xgboost::xgb.train(
    params = as.list(best_params),
    data = dtrain,
    nrounds = 100,
    objective = "multi:softprob",
    eval_metric = "mlogloss",
    num_class=num_class,
    verbose = 0)
  
  # Assessing the performance
  val_pred <- predict(xgb_outer, newdata = dval)
  val_pred <- matrix(val_pred, ncol = num_class, byrow = TRUE)
  pred_class <- max.col(val_pred) - 1
  MAE <- mean(abs(pred_class - y_val))
  
  ## Saving CV results
  outer_results[[outer_idx]] <- data.frame(Fold = outer_idx, MAE = MAE)
}

# Combining results from outer CV into a summary table
outer_summary <- dplyr::bind_rows(outer_results)
outer_summary

#--------------------------
## Retraining final model on full training set

# Picking the best params (e.g., from highest outer AUROC)
best_idx <- which.min(outer_summary$MAE)
final_params <- best_params_all[[best_idx]]

# Retraining final model
dtrain_full <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
xgb_final <- xgboost::xgb.train(
  params = as.list(final_params),
  data = dtrain_full,
  nrounds = 100,
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class=num_class,
  verbose = 1)

# Saving the final model
saveRDS(xgb_final, file = "Patient journey_XGBoost_model.rds")

#--------------------------
## Final evaluation on holdout test dataset

# Predicted probabilities
pred_prob <- predict(xgb_final, newdata = X_test)
pred_prob <- matrix(pred_prob, ncol = num_class, byrow = TRUE)

# Predicted classes
pred_class <- max.col(pred_prob) - 1
pred_f <- factor(pred_class, levels = sort(unique(y_test)))
true_f <- factor(y_test, levels = sort(unique(y_test)))

# Calculating performance metrics
multi_acc<-yardstick::accuracy_vec(truth = true_f, estimate = pred_f)                   # Accuracy
multi_F1<-yardstick::f_meas_vec(truth = true_f, estimate = pred_f, estimator = "macro") # Macro F1 score
multi_bal_acc<-yardstick::bal_accuracy_vec(truth = true_f, estimate = pred_f)           # Balanced accuracy
precision_macro <- yardstick::precision_vec(true_f, pred_f, estimator = "macro")        # Macro precision
recall_macro <- yardstick::recall_vec(true_f, pred_f, estimator = "macro")              # Macro recall
sens_macro <- yardstick::sens_vec(true_f, pred_f, estimator = "macro")                  # Macro sensitivity
spec_macro <- yardstick::spec_vec(true_f, pred_f, estimator = "macro")                  # Macro specificity
multi_kappa<-caret::confusionMatrix(pred_f, true_f)$overall["Kappa"]                    # Cohen's kappa unweighted
QWK <- irr::kappa2(cbind(pred_class, y_test), weight = "squared")$value                 # Quadratic weighted kappa
multi_logloss <- MLmetrics::MultiLogLoss(y_pred = pred_prob, y_true = y_test)           # Multi-class Logloss
Multi_MAE<-yardstick::mae_vec(truth = y_test, estimate = pred_class)                    # MAE

colnames(pred_prob) <- levels(true_f)
multi_AUC<-as.vector(pROC::multiclass.roc(true_f, pred_prob)$auc)                       # Multi-class AUC

model_eval<-data.frame(MAE=Multi_MAE,
                       Accuracy=multi_acc,
                       Macro_F1_score=multi_F1,
                       Balanced_accuracy=multi_bal_acc,
                       Macro_precision=precision_macro,
                       Macro_recall=recall_macro,
                       Macro_sensitivity=sens_macro,
                       Macro_specificity=spec_macro,
                       Cohen_kappa_unweighted=multi_kappa,
                       QWK=QWK,
                       Multi_class_logloss=multi_logloss,
                       Multi_class_AUC=multi_AUC)

model_eval<-model_eval %>% rowwise %>% mutate_if(is.numeric,test_format) %>% ungroup %>%
            pivot_longer(cols=everything(),names_to="Metric",values_to="Value")

#--------------------------
## Saving CV summary table, confusion matrix and best hyperparameters
write.table(outer_summary,paste("CV summary table_",Sys.Date(),".csv"),sep=";",col.names=T,row.names=F)
write.table(model_eval,paste("Model eval_",Sys.Date(),".csv"),sep=";",col.names=T,row.names=F)
write.table(as_tibble(final_params),paste("Best hyperparameters_",Sys.Date(),".csv"),sep=";",col.names=T,row.names=F)

#----------------------------------------------------------------
#### Univariate association analysis ####

var_test<-colnames(X_test)

Association_tab_res<-c()

for(n_var in 1:length(var_test)){
  
  test_tmpo<-as_tibble(X_test) %>% select(which(colnames(X_test)==var_test[n_var]))
  colnames(test_tmpo)<-"tmp"
  test_tmpo<-bind_cols(y_test,test_tmpo)
  colnames(test_tmpo)[1]<-"outcome"
  test_tmpo$outcome<-factor(test_tmpo$outcome,ordered = T,levels=sort(unique(test_tmpo$outcome)))
  
  clm_model<-ordinal::clm(data=test_tmpo, formula=outcome~.)
  pred_class <- predict(clm_model, type = "class")$fit
  pred_class <- as.numeric(pred_class)
  true_class <- as.numeric(test_tmpo$outcome)
  
  mod_res<-broom::tidy(clm_model,exponentiate=T,conf.int=T) %>%
    filter(coef.type!="intercept") %>%
    mutate(logloss=-mean(rowSums(model.matrix(~ outcome - 1, data=test_tmpo) * log(pmax(pmin(predict(clm_model, type = "prob")$fit, 1 - 1e-15), 1e-15)))),
           nom_test=as_tibble(data.frame(term=rownames(ordinal::nominal_test(clm_model)),
                                         ordinal::nominal_test(clm_model))) %>% filter(term=="tmp") %>% pull(Pr..Chi.),
           scale_test=as_tibble(data.frame(term=rownames(ordinal::scale_test(clm_model)),
                                           ordinal::scale_test(clm_model))) %>% filter(term=="tmp") %>% pull(Pr..Chi.),
           accuracy=mean(pred_class == true_class),
           MAE=mean(abs(pred_class - true_class)),
           kappa=irr::kappa2(data.frame(
             pred = as.numeric(pred_class),
             obs  = as.numeric(test_tmpo$outcome)),
             weight = "squared")$value) %>%
    mutate(term=var_test[n_var],
           direction=case_when(conf.low>=1~"positive association",
                               conf.high<1~"negative association",
                               T~"uncertain association")) %>% 
    bind_cols(performance::performance(clm_model)) %>%
    rowwise %>% 
    mutate_if(is.numeric,test_format) %>%
    ungroup %>%
    mutate(OR=paste(estimate," [",conf.low,"; ",conf.high,"]",sep="")) %>%
    mutate(OR=if_else(OR==" [; ]","",OR)) %>%
    select(-coef.type)

  Association_tab_res<-bind_rows(Association_tab_res,mod_res)
  
}

#----------------------------------------------------------------
#### SHAP analysis ####

#- - - - - - - - - -
## Calculating SHAP values
contr <- predict(xgb_final, as.matrix(X_test[,xgb_final$feature_names]), predcontrib = TRUE)

# Reshaping predicted probabilities into n_samples x num_class
pred_prob_mat <- matrix(pred_prob, ncol = num_class, byrow = TRUE)

# Aggregating SHAP across classes (keeping SHAP distributions), weighted by their probability
if(is.list(contr)==T){
  shap<-list()
  for(i in 1:length(contr)){

    # Get SHAP for class i
    shap_i <- as_tibble(contr[[i]])
    
    # Multiply each row by predicted probability of that class
    shap_i_weighted <- shap_i * pred_prob_mat[, i]
    
    shap[[i]] <- shap_i_weighted
    
  }
  shap<-bind_rows(shap)
  shap_contrib<-as.data.table(shap)
}else{
  shap<-as_tibble(contr)
  shap_contrib <- as.data.table(contr)
}

#- - - - - - - - - -
## Removing BIAS term
shap_contrib <- shap_contrib[, !grepl("bias", names(shap_contrib), ignore.case = TRUE), with=FALSE]

#- - - - - - - - - -
## Computing mean SHAP score
mean_shap_score <- colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing = T)]

#- - - - - - - - - -
## Reshaping SHAP values to long format
shap_score<-shap_contrib
shap_score_sub <- as.data.table(shap_score)
shap_score_sub <- shap_score_sub[, names(mean_shap_score), with = F]
shap_score_long <- melt.data.table(shap_score_sub, measure.vars = colnames(shap_score_sub))

#- - - - - - - - - -
## Matching feature values
fv_sub <- as.data.table(X_test)[, names(mean_shap_score), with = F]
fv_sub_long <- melt.data.table(fv_sub, measure.vars = colnames(fv_sub))
fv_sub_long[, stdfvalue := std1(value), by = "variable"]
names(fv_sub_long) <- c("variable", "rfvalue", "stdfvalue" )

#- - - - - - - - - -
## Merging  SHAP and feature values
shap_long2 <- cbind(shap_score_long, fv_sub_long[,c('rfvalue','stdfvalue')])
shap_long2[, mean_value := mean(abs(value)), by = variable]
setkey(shap_long2, variable)

#- - - - - - - - - -
## Computting summary statistics per feature
Shap_val<-as_tibble(shap_long2) %>% 
  group_by(variable) %>% 
  mutate(IC_2.5=quantile(abs(value),0.025),
         IC_97.5=quantile(abs(value),0.975),
         mean_val=mean(abs(value),na.rm=T)) %>%
  select(variable,mean_val,IC_2.5,IC_97.5) %>% 
  distinct %>% 
  arrange(desc(mean_val)) %>%
  mutate(mean_SHAP=paste(signif(mean_val,3)," (",signif(IC_2.5,3),"; ",signif(IC_97.5,3),")",sep="")) %>% ungroup

colnames(Shap_val)[1]<-"feature"

#- - - - - - - - - -
## Formmating SHAP values
Shap_val<-Shap_val %>% rowwise() %>% mutate(mean_val=test_format(mean_val),
                                            IC_2.5=test_format(IC_2.5),
                                            IC_97.5=test_format(IC_97.5)) %>%
  ungroup %>% mutate(mean_SHAP=paste(mean_val," (",IC_2.5,"; ",IC_97.5,")",sep="")) %>% 
  mutate(mean_SHAP=if_else(mean_SHAP=="0 (0; 0)"|mean_SHAP=="0e+00 (0e+00; 0e+00)","0",mean_SHAP))

#- - - - - - - - - -
## Pivoting SHAP and features long format
shap_tempo <- shap %>% 
  rowid_to_column("id") %>%
  pivot_longer(-id, names_to = "feature", values_to = "shap_value")

tempopo <- X_test[, xgb_final$feature_names] %>%
  as_tibble() %>%
  rowid_to_column("id") %>%
  pivot_longer(-id, names_to = "feature", values_to = "feature_value")

tempopo <- left_join(tempopo, shap_tempo, by = c("id", "feature"))

#- - - - - - - - - -
## SHAP summary and merging labels
tempopo<-left_join(tempopo,Shap_val,by=base::intersect(colnames(tempopo),colnames(Shap_val))) %>% distinct

test3<-Shap_val %>% arrange(desc(as.numeric(mean_val)))

#- - - - - - - - - -
## SHAP directionality
direction_impact<-tempopo %>% group_by(feature) %>%
  group_map(~tibble(feature=.y$feature,
                    direction=determine_direction(.x,"feature_value","shap_value")
  )) %>%
  bind_rows

test<-left_join(tempopo,direction_impact,by="feature")

test<-left_join(test,Shap_val,by=base::intersect(colnames(test),colnames(Shap_val)))

test_tmp<-test %>% mutate(impact=direction) %>%
  select(feature, mean_val,IC_2.5,IC_97.5,impact,mean_SHAP) %>% 
  distinct() %>%
  mutate(mean_val=as.numeric(mean_val),
         IC_2.5=as.numeric(IC_2.5),
         IC_97.5=as.numeric(IC_97.5)) %>%
  mutate(mean_val=if_else(impact=="mitigating predictor",-mean_val,mean_val),
         IC_2.5=if_else(impact=="mitigating predictor",-IC_2.5,IC_2.5),
         IC_97.5=if_else(impact=="mitigating predictor",-IC_97.5,IC_97.5))

test_tmp<-test_tmp %>% arrange(desc(mean_val)) %>% filter(mean_val!=0) %>% filter(!is.na(mean_val))

#- - - - - - - - - -
## Plotting SHAP

# reordoring factors
ordre_def1<-test_tmp %>% filter(mean_val>=0) %>% arrange(desc(IC_2.5))
ordre_def2<-test_tmp %>% filter(mean_val<0) %>% arrange(desc(IC_97.5))
ordre_def<-bind_rows(ordre_def1,ordre_def2)

test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def$feature)))
test_tmp$impact<-factor(test_tmp$impact,levels=c("mitigating predictor","neutral","promoting predictor"),
                        labels=c("mitigating predictor","neutral","promoting predictor"))

# Plot 1: error bar
plot1<-ggplot(test_tmp,aes(y=feature,x=mean_val,col=impact))+
  geom_point()+
  geom_errorbar(aes(xmin = IC_2.5, xmax = IC_97.5))+
  scale_x_continuous("mean |SHAP value|",label=function(x) abs(x))+
  geom_text(aes(x=IC_97.5,y=feature,label=abs(IC_97.5)),size=5,position=position_nudge(x = if_else(test_tmp$IC_97.5>=0,
                                                                                                   round(max(test_tmp$IC_97.5,na.rm=T)/10),
                                                                                                   -round(max(test_tmp$IC_97.5,na.rm=T)/10))))+
  scale_y_discrete("")+
  theme_bw()+
  scale_color_manual("Direction:",
                     na.value="white",
                     values=rev(c("neutral"="black","promoting predictor"="#C35C33","mitigating predictor"="#40B696")))+
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.position="bottom",
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

# Exporting plot 1
ggsave(plot1,
       file=paste("SHAP_errorbar_",
                  Sys.Date(),".pdf",sep=""),
       dpi=600,width=60,height=30,units = "cm",limitsize=F)

# Plot 2: barplot
ordre_def<-test_tmp %>% arrange(desc(mean_val))
test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def$feature)))

plot2<-ggplot(test_tmp,aes(y=feature,x=mean_val,fill=impact))+
  geom_bar(stat = "identity",col='black')+
  geom_text(aes(x=mean_val,y=feature,label=mean_SHAP),size=6,position=position_nudge(x = if_else(test_tmp$mean_val>=0,
                                                                                                 round(max(test_tmp$mean_val,na.rm=T)/10),
                                                                                                 -round(max(test_tmp$mean_val,na.rm=T)/10))))+
  scale_x_continuous("mean |SHAP value|",label=function(x) abs(x))+
  scale_y_discrete("")+
  theme_bw()+
  scale_fill_manual("Direction:",
                    na.value="white",
                    values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2","neutral"="white"))+
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text = element_text(size=16,color="black", face = "bold"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.position="bottom",
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

# Exporting plot 2
ggsave(plot2,
       file=paste("SHAP_barplot_",
                  Sys.Date(),".pdf",sep=""),
       dpi=600,width=60,height=30,units = "cm",limitsize=F)

#- - - - - - - - - -
## SHAP and CLM comparison

tmp_CLM<-Association_tab_res %>% select(term,OR,direction,p.value,R2_Nagelkerke)
colnames(tmp_CLM)<-c("feature","OR [95% CI]","CLM - association direction","CLM - p-value","CLM - Pseudo-R2")

test_SHAP_tmp<-test %>%  
  select(feature,mean_SHAP,direction,mean_val) %>% 
  distinct %>%
  arrange(desc(abs(as.numeric(mean_val)))) %>%
  rowwise %>%
  mutate_if(is.numeric,test_format) %>%
  ungroup %>%
  select(-c(mean_val))

colnames(test_SHAP_tmp)<-c("feature","mean_abs_SHAP [95% CI]","SHAP direction")

Compa_CLM_SHAP<-left_join(test_SHAP_tmp,tmp_CLM,by="feature")

#----------------------------------------------------------------
#### Saving SHAP results

# SHAP values
fwrite(Shap_val,paste("Shap_val_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
fwrite(test_tmp,paste("test_tmp_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
fwrite(test,paste("test_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)
fwrite(shap,paste("shap_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)
fwrite(shap_score_sub,paste("shap_score_sub_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)
fwrite(Association_tab_res,paste("Univariate CLM regression_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
fwrite(Compa_CLM_SHAP,paste("Compa_CLM_SHAP_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
