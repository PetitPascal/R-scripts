#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing packages if needed
pack_needed<-c("data.table","tidyverse","mllrnrs","broom","doParallel","foreach",
               "splitTools","conflicted","grid","gridExtra","RColorBrewer","mlbench",
               "mlexperiments","caret","MLmetrics","patchwork","CardioDataSets","performance",
               "xgboost","parallel","here")
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
here::here("XGBoost - Linear regression - nested CV")

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
    
    if(x<0.05|x>=10000){
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
    return("undefined")
  
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

## Importing the clean dataset
clean_data<-as_tibble(na.omit(CardioDataSets::cardioRiskFactors_df))

## Renaming the outcome variable
colnames(clean_data)[which(colnames(clean_data)=="Uric")]<-"outcome"
clean_data_save<-clean_data

## Creating a table with the variable type
var_type_tab<-as_tibble(data.frame(feature=colnames(clean_data),
                                   Var_type=c("continuous","continuous","continuous","dichotomous","continuous","continuous","continuous","continuous",
                                              "continuous","continuous","continuous","dichotomous","continuous","continuous")))

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
if (length(categorical_var) > 0) {
  for (col in categorical_var) {
    dataset[[col]] <- addNA(dataset[[col]])
  }
  dummies <- model.matrix(~ -1 + ., data = dataset[, ..categorical_var])
  dataset <- cbind(dataset[, !..categorical_var], dummies)
}

## Creating the stratified 70/30 train-test split
set.seed(seed)
data_split <- splitTools::partition(
  y = dataset[[target_col]],
  p = c(train = 0.7, test = 0.3),
  type = "stratified",
  seed = seed
)

## Preparing training and test datasets

# train dataset
X_train <- as.matrix(dataset[data_split$train, setdiff(colnames(dataset), target_col), with = FALSE])
y_train <- as.integer(dataset[data_split$train, get(target_col)])

# test dataset
X_test  <- as.matrix(dataset[data_split$test, setdiff(colnames(dataset), target_col), with = FALSE])
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
  
  # Handling imbalance data
  Weight <- max(table(y_tr)) / min(table(y_tr))
  
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
    max_depth = c(3, 5, 7),
    scale_pos_weight = Weight
  ) %>% dplyr::slice_sample(n = 30, replace = TRUE) # Limiting space to 30 combinations for computational efficiency and environmental sustainability considerations
  
  # Initializing best parameters
  best_auc <- -Inf
  best_params <- NULL
  
  for (i in 1:nrow(param_list)) {
    
    # Initializing cross validation
    xgb_cv <- mlexperiments::MLCrossValidation$new(
      learner = mllrnrs::LearnerXgboost$new(metric_optimization_higher_better = TRUE),
      fold_list = inner_folds,
      ncores = 2,
      seed = 0
    )
    
    # Defining the learner arguments
    xgb_cv$learner_args <- c(
      as.list(param_list[i, ]),
      list(objective = "reg:squarederror", eval_metric = "rmse", nrounds = 100L)
    )
    xgb_cv$performance_metric <- mlexperiments::metric("rmse")
    xgb_cv$set_data(x = X_tr, y = y_tr)
    
    # Executing CV tuning
    res_cv <- xgb_cv$execute()
    
    # Calculating the mean AUROC
    mean_auc <- mean(sapply(res_cv$results$folds, function(f) f[[5]]$performance))
    
    if (!is.na(mean_auc) && mean_auc > best_auc) {
      best_auc <- mean_auc
      best_params <- param_list[i, ]
    }
  }
  
  # If no params are selected, fallback to default
  if (is.null(best_params)) {
    best_params <- data.frame(
      subsample = 1,
      colsample_bytree = 1,
      min_child_weight = 1,
      learning_rate = 0.1,
      max_depth = 3,
      scale_pos_weight = Weight
    )
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
    objective = "reg:squarederror",
    eval_metric = "rmse",
    verbose = 0
  )
  
  # Assessing the performance
  val_pred <- predict(xgb_outer, newdata = dval)
  R2   <- 1 - sum((y_val - val_pred)^2) / sum((y_val - mean(y_val))^2)
  MAE  <- Metrics::mae(y_val, val_pred)
  RMSE <- Metrics::rmse(y_val, val_pred)
  r_pearson<-cor(y_val, val_pred)
  CCC_val<-ccc_value <-(2 * cov(y_val, val_pred, use = "complete.obs"))/(var(y_val, na.rm = TRUE) + var(val_pred, na.rm = TRUE) + (mean(y_val, na.rm = TRUE) - mean(val_pred, na.rm = TRUE))^2)
  
  ## Saving CV results
  outer_results[[outer_idx]] <- data.frame(Fold = outer_idx, R2 = R2, MAE = MAE,RMSE=RMSE,r_pearson=r_pearson,CCC=CCC_val)
}

# Combining results from outer CV into a summary table
outer_summary <- dplyr::bind_rows(outer_results)
outer_summary

#--------------------------
## Retraining final model on full training set

# Picking the best params (e.g., from lowest outer RMSE)
best_idx <- which.min(outer_summary$RMSE)
final_params <- best_params_all[[best_idx]]

# Retraining final model
dtrain_full <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
xgb_final <- xgboost::xgb.train(
  params = as.list(final_params),
  data = dtrain_full,
  nrounds = 100,
  objective = "reg:squarederror",
  eval_metric = "rmse",
  verbose = 1
)

# Saving the final model
saveRDS(xgb_final, file = "XGBoost_model.rds")

#--------------------------
## Final evaluation on holdout test dataset

# Predicting outcome in holdout test dataset
val_pred <- predict(xgb_final, newdata = X_test)
R2   <- 1 - sum((y_test - val_pred)^2) / sum((y_test - mean(y_test))^2)
MAE  <- Metrics::mae(y_test, val_pred)
RMSE <- Metrics::rmse(y_test, val_pred)
r_pearson<-cor(y_test, val_pred)
CCC_val<-ccc_value <-(2 * cov(y_test, val_pred, use = "complete.obs"))/(var(y_test, na.rm = TRUE) + var(val_pred, na.rm = TRUE) + (mean(y_test, na.rm = TRUE) - mean(val_pred, na.rm = TRUE))^2)

# Saving the model characteristics into a table
Mod_charact<-as_tibble(data.frame(R2=R2,
                                  MAE=MAE,
                                  RMSE=RMSE,
                                  r_pearson=r_pearson,
                                  CCC=CCC_val))

#--------------------------
## Saving CV summary table, confusion matrix and best hyperparameters
write.table(outer_summary,paste("CV summary table_",Sys.Date(),".csv"),sep=";",col.names=T,row.names=F)
write.table(Mod_charact,paste("Model characteristics_",Sys.Date(),".csv"),sep=";",col.names=T,row.names=F)
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
  
  glm_model<-glm(outcome ~ tmp, data = test_tmpo,family="gaussian")
  mod_res<-broom::tidy(glm_model,exponentiate=F,conf.int=T) %>%
    filter(term!="(Intercept)") %>%
    mutate(term=var_test[n_var],
           direction=case_when(conf.low>=0~"positive association",
                               conf.high<0~"negative association",
                               T~"uncertain association")) %>% 
    bind_cols(performance::performance(glm_model)) %>%
    rowwise %>% 
    mutate_if(is.numeric,test_format) %>%
    ungroup %>%
    mutate(beta=paste(estimate," [",conf.low,"; ",conf.high,"]",sep="")) %>%
    mutate(beta=if_else(beta==" [; ]","",beta))
  
  Association_tab_res<-bind_rows(Association_tab_res,mod_res)
  
}

#----------------------------------------------------------------
#### SHAP analysis ####

#- - - - - - - - - -
## Calculating SHAP values
contr <- predict(xgb_final, as.matrix(X_test[,xgb_final$feature_names]), predcontrib = TRUE)
shap<-as_tibble(contr)
shap_contrib <- as.data.table(contr)

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
test_save<-xgb.ggplot.shap.summary(as.matrix(X_test[,xgb_final$feature_names]), contr, model = xgb_final,top_n=25)
test_save<-as_tibble(test_save$data)

tempopo<-left_join(tempopo,Shap_val,by=base::intersect(colnames(tempopo),colnames(Shap_val)))

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
                     values=rev(c("neutral"="white","promoting predictor"="#C35C33","mitigating predictor"="#40B696")))+
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
                    values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2"))+
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

# Plot 3: dispersion

bornes<-test %>% select(feature_value) %>% summarize(min=floor(min(feature_value,na.rm=T)),
                                                     max=ceiling(max(feature_value,na.rm=T)))
bornes_shap<-test %>% select(shap_value) %>% summarize(min=min(shap_value,na.rm=T),max=max(shap_value,na.rm=T))

test<-test %>% group_by(feature) %>% 
  mutate(feature_value2=(feature_value-min(feature_value,na.rm=T))/(max(feature_value,na.rm=T)-min(feature_value,na.rm=T))) %>%
  ungroup

test3<-test %>% filter(mean_SHAP!="0")
ordre_Def<-test3 %>% select(feature,mean_val) %>% distinct %>% mutate(mean_val=as.numeric(mean_val)) %>% arrange(mean_val)
test3$feature<-factor(test3$feature,level=ordre_Def$feature)
test2<-test3 %>% select(feature,mean_SHAP) %>% distinct

plot3<-ggplot2::ggplot(test3, ggplot2::aes(x = feature, y = shap_value,colour = feature_value2)) + 
  ggplot2::geom_jitter(width = 0.1) + 
  scale_color_gradientn("Feature value",labels=c("Low", "High"),
                        breaks=c(0,1),
                        limits=c(0,1),
                        colours=c("blue",scales::alpha("blue",0.75),
                                  scales::alpha("blue",0.5),
                                  scales::alpha("blue",0.35),
                                  scales::alpha("#FFA5001A",0.25),
                                  scales::alpha("red",0.35),
                                  scales::alpha("red",0.5),
                                  scales::alpha("red",0.75),
                                  "red"),
                        na.value ="lightgrey",
                        guide = guide_colorbar(barheight = 50,
                                               title.theme = element_text(angle = 90, hjust = 0.5, vjust = 0,size = 16, face = "bold"),
                                               title.position = "left"
                        ))+
  ggplot2::geom_abline(slope = 0,intercept = 0, colour = "darkgrey") + 
  scale_x_discrete("Feature")+
  scale_y_continuous("SHAP value (impact on model output)",limits=c(bornes_shap$min,bornes_shap$max*1.5))+
  ggplot2::coord_flip()+
  geom_text(data=test2,aes(x = feature, y = bornes_shap$max*1.35,label=mean_SHAP),size=5,col="black")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill=scales::alpha('#009E73', 0.2), colour="black", size=1),
        axis.text = element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

# Exporting plot 3
ggsave(plot3,
       file=paste("SHAP_dispersion_",
                  Sys.Date(),".pdf",sep=""),
       dpi=600,width=60,height=30,units = "cm",limitsize=F)

#- - - - - - - - - -
## SHAP and GLM comparison

tmp_GLM<-Association_tab_res %>% select(term,beta,direction,R2,p.value)
colnames(tmp_GLM)<-c("feature","Beta [95% CI]","GLM - association direction","GLM - R2","GLM -p-value")

test_SHAP_tmp<-test %>%  
  select(feature,mean_SHAP,direction,mean_val) %>% 
  distinct %>%
  arrange(desc(abs(as.numeric(mean_val)))) %>%
  rowwise %>%
  mutate_if(is.numeric,test_format) %>%
  ungroup %>%
  select(-c(mean_val))

colnames(test_SHAP_tmp)<-c("feature","mean_abs_SHAP [95% CI]","SHAP direction")

Compa_GLM_SHAP<-left_join(test_SHAP_tmp,tmp_GLM,by="feature")

#----------------------------------------------------------------
#### Saving SHAP results

# SHAP values
fwrite(Shap_val,paste("Shap_val_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
fwrite(test_tmp,paste("test_tmp_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
fwrite(test,paste("test_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)
fwrite(shap,paste("shap_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)
fwrite(shap_score_sub,paste("shap_score_sub_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)
fwrite(Association_tab_res,paste("Univariate GLM regression_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
fwrite(Compa_GLM_SHAP,paste("Compa_GLM_SHAP_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
