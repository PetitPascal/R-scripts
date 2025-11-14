#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing packages if needed
pack_needed<-c("data.table","tidyverse","mllrnrs","mlsurvlrnrs","survival","splitTools","conflicted","mlexperiments","kdry","survminer","timeROC","cluster","pec",
               "factoextra","R6","xgboost","mgcv","quantreg","parallel","here")
for (i in 1:length(pack_needed)){
  if(pack_needed[i]%in%.packages(all.available=TRUE)){
  }else{
    install.packages(pack_needed[i])
  }
}

## Package loading
library(mgcv)
library(quantreg)
library(mllrnrs)
library(mlsurvlrnrs)
library(mlexperiments)
library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(splitTools)
library(conflicted)
library(timeROC)
library(cluster)
library(factoextra)
library(pec)
library(xgboost)
library(here)

## Preventing package conflicts
conflict_prefer("select", "dplyr") 
conflict_prefer("filter", "dplyr") 
conflict_prefer("slice", "dplyr")
conflict_prefer("alpha", "scales")

## Setting the working directory
here::here("Survival XGBoost Cox analysis")

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

## Function for standardizing feature values into the same range
std1 <- function(x, clip = TRUE) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    return(rep(0.5, length(x)))
  }
  scaled <- (x - rng[1]) / (rng[2] - rng[1])
  if (clip) scaled <- pmin(pmax(scaled, 0), 1)
  return(scaled)
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

## Normalization function
normalize_custom <- function(x, new_min, new_max) {
  scaled <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  new_min + scaled * (new_max - new_min)
}

## Function for determining the feature direction effect (i.e.: how a feature impacts the model)
determine_direction<-function(df,feature_col,shap_col){
  
  # extracting vectors
  x<-df[[feature_col]]
  s<-df[[shap_col]]
  
  # removing NA
  valid<-complete.cases(x,s)
  x<-x[valid]
  s<-s[valid]
  
  # for cases when all values are identical
  if(length(unique(x))<=1||length(unique(s))<=1){
    return("undefined")
  }
  
  # for cases when the feature is binary or sparse
  if(length(unique(x))<=2||quantile(x,0.75)==0){
    mean_diff<-mean(s[x>0],na.rm=T)-mean(s[x<=0],na.rm=T)
    
    return(ifelse(mean_diff>0,'risk factor',ifelse(mean_diff<0,"protective factor","neutral")))
  }
  
  # for cases when the feature is continuous
  tryCatch({
    gam_model<-gam(s~s(x), method='REML')
    derivs<-derivatives(gam_model, term='s(x)')
    mean_deriv<-mean(derivs$derivative,na.rm=T)
    
    return(ifelse(mean_deriv>0,'risk factor',ifelse(mean_deriv<0,"protective factor","neutral")))
  }, error=function(e){
    
    # if GAM fails, Spearman correlation is computed
    rho<-suppressWarnings(cor(x,s,method="spearman"))
    if(is.na(rho)) return("undefined")
    if(abs(rho)<0.05) return("neutral")
    return(ifelse(rho>0,"risk factor","protective factor"))
  }
  )
}

#----------------------------------------------------------------
#### Data import and preprocessing ####

## Importing the clean dataset
clean_data<-as_tibble(na.omit(survival::colon) %>% filter(etype == 2) %>% select(-id) %>% mutate(rx=case_when(rx=="Obs"~0,rx=="Lev"~1,T~2)))

## Renaming outcome variables
colnames(clean_data)[which(colnames(clean_data)=="status")]<-"Disease_status"
colnames(clean_data)[which(colnames(clean_data)=="time")]<-"time_to_diagnosis"

## Transforming the clean dataset into data.table
dataset <- as.data.table(clean_data)

## Creating vectors with the independent and dependent variable names
surv_cols <- c("Disease_status", "time_to_diagnosis")
feature_cols <- setdiff(colnames(dataset), c(surv_cols,"ID"))

## Creating split vector and data partition
split_vector <- splitTools::multi_strata(df = dataset[, ..surv_cols], strategy = "kmeans", k = 4)
data_split <- splitTools::partition(y = split_vector, p = c(train = 0.7, test = 0.3), type = "stratified", seed = seed)

## Ensuring no individuals are present in both training and test datasets
good_split<-length(intersect(data_split$train, data_split$test))==0
print(good_split)

#- - - - - - - - - -
## Preparing training and test datasets

## Creating the training set
train_x <- model.matrix(~ -1 + ., dataset[data_split$train, ..feature_cols])
train_y <- Surv(time = dataset[data_split$train, time_to_diagnosis], event = dataset[data_split$train, Disease_status])

## Creating the test set
test_x <- model.matrix(~ -1 + ., dataset[data_split$test, ..feature_cols])
test_y <- Surv(time = dataset[data_split$test, time_to_diagnosis], event = dataset[data_split$test, Disease_status])

## Computing sample weights to account for imbalance data
event_counts <- dataset[data_split$train, .N, by = Disease_status][, weight := 1 / N]
train_data <- dataset[data_split$train]
train_weights <- event_counts[train_data, on = "Disease_status", weight]

## Creating cross-validation (CV) folds
split_vector_train <- splitTools::multi_strata(dataset[data_split$train, ..surv_cols], strategy = "kmeans", k = 4)
fold_list <- splitTools::create_folds(y = split_vector_train, 
                                      k = 5,  # 5-fold CV
                                      type = "stratified",
                                      seed = seq(1,5,1)*123) # different random seed for each fold

stopifnot(length(intersect(data_split$train, data_split$test)) == 0) # Ensuring no overlap

#----------------------------------------------------------------
#### Customing the survival XGBoost Cox learner with sample weights support ####

LearnerSurvXgboostCoxWeighted <- R6::R6Class(
  classname = "LearnerSurvXgboostCoxWeighted",
  inherit = mlsurvlrnrs::LearnerSurvXgboostCox,
  public = list(
    sample_weight = NULL,
    
    set_sample_weight = function(w) {
      self$sample_weight <- w
    },
    
    train = function(x, y, ...) {
      
      if (!is.null(self$sample_weight)) {
        super$train(x = x, y = y, sample_weight = self$sample_weight, ...)
      } else {
        super$train(x = x, y = y, ...)
      }
    }
  )
)

#----------------------------------------------------------------
#### Creating and executing the model ####

#- - - - - - - - - -
## Defining the learner arguments
learner_args <- list(
  objective = "survival:cox",
  eval_metric = "cox-nloglik",
  nthread = ncores
)

#- - - - - - - - - -
## Setting arguments for the predict function and performance metric
predict_args <- NULL
performance_metric <- c_index
performance_metric_args <- NULL
return_models <- T

#- - - - - - - - - -
## Creating a random-like grid for hyperparameter tuning (random search is not currently possible with mlexperiments::MLCrossValidation$new)

# Grid search
parameter_grid <- expand.grid(
  learning_rate = c(0.05, 0.1),
  max_depth = c(3, 4, 5),
  min_child_weight = c(5, 10, 20), # to avoid overfitting because small values cause instability
  subsample = c(0.7, 0.9, 1.0), 
  colsample_bytree = c(0.7, 0.9, 1.0),
  nrounds = c(50, 100, 200)
)

# Limiting grid search to a maximum of 10 rows for computational efficiency and environmental sustainability considerations
if (nrow(parameter_grid) > 10) {
  set.seed(seed)
  parameter_grid <- kdry::mlh_subset(parameter_grid, sample(seq_len(nrow(parameter_grid)), 10))
}

#- - - - - - - - - -
## Hyperparameter tuning

# Creating the learner and setting weights
learner <- LearnerSurvXgboostCoxWeighted$new(metric_optimization_higher_better = FALSE)
learner$set_sample_weight(train_weights)

# Initializing hyperparameter tuning
tuner <- mlexperiments::MLTuneParameters$new(
  learner = learner,
  strategy = "grid",
  ncores = ncores,
  seed = seed
)

# Setting tuning arguments
tuner$parameter_grid <- parameter_grid
tuner$learner_args <- learner_args
tuner$split_type <- "stratified"
tuner$split_vector <- split_vector_train
tuner$learner$set_sample_weight(train_weights)
tuner$set_data(x = train_x, y = train_y)

# Executing tuning
tuner_results_grid <- tuner$execute(k = 5)

# Identifying best tuning
best_setting <- tuner$results$best.setting

#- - - - - - - - - -
## Initializing cross validation
validator <- mlexperiments::MLCrossValidation$new(
  learner = learner,
  fold_list = fold_list,
  ncores = ncores,
  seed = seed
)

#- - - - - - - - - -
## Setting validator arguments
validator$learner_args <- tuner$results$best.setting[-1]
validator$learner$set_sample_weight(train_weights) # taking weight into account in case of imbalance dataset
validator$predict_args <- NULL
validator$performance_metric <- c_index
validator$performance_metric_args <- NULL
validator$return_models <- T
validator$set_data(x = train_x,y = train_y) # Setting data for validator (training set)

#- - - - - - - - - -
## Executing training
validator_results <- validator$execute()

#- - - - - - - - - -
## Predicting outcome in holdout test dataset
preds_xgboost <- mlexperiments::predictions(object = validator, newdata = test_x)

#- - - - - - - - - -
## Evaluating performance on holdout test dataset
perf_xgboost <- mlexperiments::performance(object=validator,
                                           prediction_results = preds_xgboost,
                                           y_ground_truth = test_y)

#- - - - - - - - - -
## Identifying the best model (highest C index)
best_mod<-perf_xgboost %>% as_tibble() %>% filter(performance==max(performance)) %>% slice(1) %>% pull(model)
best_model <- validator$results$fold[[best_mod]]$model

#----------------------------------------------------------------
#### Retraining final model on full training set

# Extracting best hyperparameters (excluding model ID)
best_params <- tuner$results$best.setting[-1]

# Initializing new model
final_model <- LearnerSurvXgboostCoxWeighted$new(metric_optimization_higher_better = FALSE)
final_model$set_sample_weight(train_weights)

# Training on full training set
final_mod<-final_model$fit(
  x = train_x,
  y = train_y,
  seed = seed,
  ncores = ncores,
  objective = best_params$objective,
  eval_metric = best_params$eval_metric,
  subsample = best_params$subsample,
  colsample_bytree = best_params$colsample_bytree,
  min_child_weight = best_params$min_child_weight,
  learning_rate = best_params$learning_rate,
  nrounds = best_params$nrounds,
  max_depth = best_params$max_depth
)

# Predicting on training and test sets
preds_train_final <- predict(final_mod, xgboost::xgb.DMatrix(data = as.matrix(train_x)))
preds_test_final <- predict(final_mod, xgboost::xgb.DMatrix(data = as.matrix(test_x)))

# Evaluating performance
c_index_train_final <- c_index(predictions = preds_train_final, ground_truth = train_y)
c_index_test_final  <- c_index(predictions = preds_test_final,  ground_truth = test_y)

#----------------------------------------------------------------
#### Cox model evaluation ####

#- - - - - - - - - -
## Grouping individuals into groups

# Creating a data frame with predictions and outcomes
risk_df <- data.frame(risk_score = preds_xgboost[["mean"]],
                      time = test_y[, "time"],
                      status = test_y[, "status"])

# Creating sex groups
risk_df$sex <- as.factor(unlist(dataset[data_split$test, "sex"]))

# Creating age groups
risk_df$age <- dataset[data_split$test, "age"]
risk_df$age_group <- ifelse(risk_df$age < 60, "<60", "≥60")

# Combining age group and sex
risk_df$sex_age <- interaction(risk_df$age_group, risk_df$sex, sep = "_")
risk_df$sex_age <- as.factor(risk_df$sex_age)

# Creating risk tertile groups
risk_df <- risk_df %>%
           mutate(risk_group = ntile(risk_score, 3)) %>%
           mutate(risk_group = factor(risk_group, labels = c("Low", "Medium", "High")))

#- - - - - - - - - -
## Computing hazard ratios

# Fiting Cox model on test data using sex
cox_model_sex <- coxph(Surv(time, status) ~ sex, data = risk_df)
summary(cox_model_sex)

# Fiting Cox model on test data using age group
cox_model_age <- coxph(Surv(time, status) ~ age_group, data = risk_df)
summary(cox_model_age)

# Fiting Cox model on test data using age and sex groups combined
cox_model_sex_age <- coxph(Surv(time, status) ~ sex_age, data = risk_df)
summary(cox_model_sex_age)

# Fiting Cox model on test data using the tertile groups -> hazard ratios for medium and high risk vs. low risk group
cox_model_tertiles <- coxph(Surv(time, status) ~ risk_group, data = risk_df)
summary(cox_model_tertiles)

#- - - - - - - - - -
## Kaplan-Meier plot by sex

# KM curve by predicted sex
km_sex <- survfit(Surv(time, status) ~ sex, data = risk_df)

# Plot
survminer::ggsurvplot(
  km_sex,
  data = risk_df,
  risk.table = TRUE,
  pval = TRUE,
  palette = c("#40B696", "#F1CB0E"),
  legend.title = "Sex",
  legend.labs = c("Male", "Female"),
  xlab = "Time",
  ylab = "Survival Probability")

#- - - - - - - - - -
## Kaplan-Meier plot by age groups

# KM curve by predicted age
km_age <- survfit(Surv(time, status) ~ age_group, data = risk_df)

# Plot
survminer::ggsurvplot(
  km_age,
  data = risk_df,
  risk.table = TRUE,
  pval = TRUE,
  palette = c("#40B696", "#F1CB0E"),
  legend.title = "Age Group",
  legend.labs = c("<60", "≥60"))

#- - - - - - - - - -
## Age and sex interaction

# KM curve - fitting survival curve
surv_fit_sex_age <- survfit(Surv(time, status) ~ sex_age, data = risk_df)

# Plot
survminer::ggsurvplot(
  fit = surv_fit_sex_age,
  data = risk_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = F,
  palette = "Set1",
  legend.title = "Predicted Risk × Age",
  xlab = "Time",
  ylab = "Survival probability",
  ggtheme = theme_minimal())

#- - - - - - - - - -
## Kaplan-Meier plot by risk tertile group

# KM curve by predicted risk tertile
km_fit <- survfit(Surv(time, status) ~ risk_group, data = risk_df)

# Plot
survminer::ggsurvplot(
  km_fit,
  data = risk_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  palette = c("#40B696", "#F1CB0E", "red"),
  legend.title = "Predicted Risk Group",
  legend.labs = c("Low", "Medium", "High"),
  xlab = "Time",
  ylab = "Survival probability")

#----------------------------------------------------------------
#### Feature importance ####

#- - - - - - - - - -
## Extracting features considered in the best model
names <- final_mod$feature_names

## Creating the importance matrix
importance_matrix <- xgb.importance(names, model = final_mod)
importance_matrix<-as_tibble(importance_matrix)

#- - - - - - - - - -
## Selecting top 50 predictors
importance_matrix2<-as.data.table(importance_matrix %>% slice(1:50))

#- - - - - - - - - -
## Plotting the feature importance
plot1<-xgb.ggplot.importance(importance_matrix2)+
  theme_bw()+
  scale_y_continuous("Importance")+
  ggtitle("")+
  theme(strip.text = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill='#009E73', colour="black", size=1),
        axis.text = element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title=element_text(size = 16, face = "bold"),
        legend.position = c(.87, .8),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid")) 

## Exporting the graph
ggsave(plot1,
       file=paste("Importance_Top 50_",
                  Sys.Date(),".pdf",sep=""),
       dpi=600,width=40,height=30,units = "cm",limitsize=F)

#----------------------------------------------------------------
#### SHAP analysis ####

#- - - - - - - - - -
## Calculating SHAP values
contr <- predict(final_mod, as.matrix(test_x[,final_mod$feature_names]), predcontrib = TRUE)
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
fv_sub <- as.data.table(test_x)[, names(mean_shap_score), with = F]
fv_sub_long <- melt.data.table(fv_sub, measure.vars = colnames(fv_sub))
fv_sub_long[, stdfvalue := std1(value), by = "variable"]
names(fv_sub_long) <- c("variable", "rfvalue", "stdfvalue")

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
  select(variable,mean_value,IC_2.5,IC_97.5) %>% 
  distinct %>% 
  arrange(desc(mean_value)) %>%
  mutate(mean_SHAP=paste(signif(mean_value,3)," (",signif(IC_2.5,3),"; ",signif(IC_97.5,3),")",sep="")) %>% ungroup

colnames(Shap_val)[1]<-"feature"

#- - - - - - - - - -
## Formmating SHAP values
Shap_val<-Shap_val %>% rowwise() %>% mutate(mean_value=test_format(mean_value),
                                            IC_2.5=test_format(IC_2.5),
                                            IC_97.5=test_format(IC_97.5)) %>%
  ungroup %>% mutate(mean_SHAP=paste(mean_value," (",IC_2.5,"; ",IC_97.5,")",sep="")) %>% 
  mutate(mean_SHAP=if_else(mean_SHAP=="0 (0; 0)"|mean_SHAP=="0e+00 (0e+00; 0e+00)","0",mean_SHAP))

#- - - - - - - - - -
## Pivoting SHAP and features long format
shap_tempo <- shap %>% 
  rowid_to_column("id") %>%
  pivot_longer(-id, names_to = "feature", values_to = "shap_value")

tempopo <- test_x[, final_mod$feature_names] %>%
  as_tibble() %>%
  rowid_to_column("id") %>%
  pivot_longer(-id, names_to = "feature", values_to = "feature_value")

tempopo <- left_join(tempopo, shap_tempo, by = c("id", "feature"))

#- - - - - - - - - -
## SHAP summary and merging labels
test_save<-xgb.ggplot.shap.summary(as.matrix(test_x[,final_mod$feature_names]), contr, model = final_mod,top_n=50)
test_save<-as_tibble(test_save$data)

tempopo<-left_join(tempopo,Shap_val,by=base::intersect(colnames(tempopo),colnames(Shap_val)))

test3<-Shap_val %>% arrange(desc(as.numeric(mean_value)))

#- - - - - - - - - -
## SHAP directionality summary
direction_impact<-test_save %>% group_by(feature) %>%
  group_map(~tibble(feature=.y$feature,
                    direction=determine_direction(.x,"feature_value","shap_value")
  )) %>%
  bind_rows

test<-left_join(tempopo,direction_impact,by="feature")

test_tmp<-test %>% mutate(impact=direction) %>%
  select(feature, mean_value,IC_2.5,IC_97.5,impact) %>% 
  distinct() %>%
  mutate(mean_value=as.numeric(mean_value),
         IC_2.5=as.numeric(IC_2.5),
         IC_97.5=as.numeric(IC_97.5)) %>%
  mutate(mean_value=if_else(impact=="protective factor",-mean_value,mean_value),
         IC_2.5=if_else(impact=="protective factor",-IC_2.5,IC_2.5),
         IC_97.5=if_else(impact=="protective factor",-IC_97.5,IC_97.5))

test_tmp<-test_tmp %>% arrange(desc(mean_value)) %>% filter(mean_value!=0)

#- - - - - - - - - -
## Plotting SHAP

# reordoring factors
ordre_def1<-test_tmp %>% filter(mean_value>=0) %>% arrange(desc(IC_2.5))
ordre_def2<-test_tmp %>% filter(mean_value<0) %>% arrange(desc(IC_97.5))
ordre_def<-bind_rows(ordre_def1,ordre_def2)

test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def$feature)))
test_tmp$impact<-factor(test_tmp$impact,levels=c("protective factor","neutral","risk factor"),labels=c("protective factor","no impact","risk factor"))

# Plot 1: error bar
plot1<-ggplot(test_tmp,aes(y=feature,x=mean_value,col=impact))+
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
                     values=rev(c("white","#C35C33","#40B696")))+
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
       file=paste("SHAP_top 50_v1_",
                  Sys.Date(),".pdf",sep=""),
       dpi=600,width=60,height=30,units = "cm",limitsize=F)

# Plot 2: barplot
ordre_def<-test_tmp %>% arrange(desc(mean_value))
test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def$feature)))

plot2<-ggplot(test_tmp,aes(y=feature,x=mean_value,fill=impact))+
  geom_bar(stat = "identity",col='black')+
  geom_text(aes(x=mean_val,y=feature,label=mean_SHAP),size=6,position=position_nudge(x = if_else(test_tmp$mean_val>=0,
                                                                                                 round(max(test_tmp$mean_val,na.rm=T)/10),
                                                                                                 -round(max(test_tmp$mean_val,na.rm=T)/10))))+
  scale_x_continuous("mean |SHAP value|",label=function(x) abs(x))+
  scale_y_discrete("")+
  theme_bw()+
  scale_fill_manual("Direction:",
                    na.value="white",
                    values=c("#A6DDCE","#F9CBC2"))+
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
       file=paste("SHAP_top 50_v2_",
                  Sys.Date(),".pdf",sep=""),
       dpi=600,width=60,height=30,units = "cm",limitsize=F)

#----------------------------------------------------------------
#### Saving results

# Best hyperparameters
fwrite(as_tibble(data.frame(Parameters=names(unlist(best_params)),Values=unlist(best_params))),
       paste("best_parameters_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)

# CV performance
fwrite(perf_xgboost,paste("perf_xgboost_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)

# C index from final model
fwrite(as_tibble(c_index_train_final),paste("c_index_train_final_mod_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE) # c index final model on training set
fwrite(as_tibble(c_index_test_final),paste("c_index_test_final_mod_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE) # c index final model on test set

# Feature importance
fwrite(importance_matrix,paste("importance_matrix_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)

# SHAP values
fwrite(Shap_val,paste("Shap_val_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
fwrite(test_tmp,paste("test_tmp_",Sys.Date(),".csv",sep=""), sep = ";", row.names=FALSE)
fwrite(test,paste("test_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)
fwrite(shap,paste("shap_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)
fwrite(shap_score_sub,paste("shap_score_sub_",Sys.Date(),".txt",sep=""), sep = ";", row.names=FALSE)

