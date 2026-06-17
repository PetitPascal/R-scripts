#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing packages if needed
pack_needed<-c("data.table","tidyverse","mllrnrs","broom","doParallel","foreach",
               "splitTools","conflicted","grid","gridExtra","RColorBrewer","mlbench",
               "mlexperiments","caret","MLmetrics","patchwork","yardstick","pROC","ordinal","performance",
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
seed<-123

## Setting cores
if(isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_")))){
  ncores<-2L
}else{
  ncores<-ifelse(test=parallel::detectCores() > 4,yes=4L,no=ifelse(test=parallel::detectCores() < 2L,yes=1L,no=parallel::detectCores()))
}

## Setting mlexperiments package options
options("mlexperiments.bayesian.max_init"=10L)
options("mlexperiments.optim.xgb.nrounds"=100L)
options("mlexperiments.optim.xgb.early_stopping_rounds"=10L)

#----------------------------------------------------------------
#### Creating functions ####

#- - - - - -
## Function for plotting SHAP values

plot.shap.summary<-function(data_long){
  x_bound<-max(abs(data_long$value))
  require('ggforce')
  plot1<-ggplot(data=data_long)+
    coord_flip() + 
    geom_sina(aes(x=variable, y=value, color=stdfvalue)) +
    geom_text(data=unique(data_long[, c("variable", "mean_value"), with=F]),
              aes(x=variable, y=-Inf, label=sprintf("%.3f", mean_value)),
              size=3, alpha=0.7,hjust=-0.2, fontface="bold") +
    scale_color_gradient(low="#FFCC33", high="#6600CC", breaks=c(0,1), labels=c("Low","High")) +
    theme_bw() + 
    theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), legend.position="bottom") + 
    geom_hline(yintercept=0) +
    scale_y_continuous(limits=c(-x_bound, x_bound)) +
    scale_x_discrete(limits=rev(levels(data_long$variable))) + 
    labs(y="SHAP value (impact on model output)", x="", color="Feature value") 
  return(plot1)
}

#- - - - - -
## Standardizing feature values into [0,1]

std1<-function(x){
  return((x - min(x, na.rm=T))/(max(x, na.rm=T) - min(x, na.rm=T)))
}

#- - - - - -
## Formatting summary statistics for display

test_format<-function(x){
  x<-as.numeric(x)
  sign_x<-if_else(x<0,"neg","pos")
  x<-abs(x)
  x_raw<-x
  if(is.na(x)|is.infinite(x)) x_raw<-0
  if(x_raw>=100){
    virg_pos<-str_locate(as.character(x_raw),"[.]")[1]
    if(!is.na(virg_pos)&as.numeric(substr(x_raw,virg_pos+1,virg_pos+1))>=5){
      x<-x+1
      x<-as.numeric(substr(x,1,virg_pos-1))
    }
  }
  if(is.na(x)|is.infinite(x)){
    x<-""
  }else{
    if(x<0.05|x>=10000){
      x<-format(signif(x,3), scientific=TRUE)
    }else{
      x_save<-x
      x<-signif(x,3)
      if(nchar(x)==6) x<-as.numeric(substr(x,1,5))
      if(nchar(x)==5){
        if(as.numeric(substr(x,5,5))>=5){
          x<-x+0.01
          x<-substr(x,1,4)
        }else{
          x<-substr(x,1,4)
        }
      }else{
        if(x>=1000) x<-as.character(signif(x_save,4)) else x<-as.character(x)
      }
    }
  }
  if(sign_x=="neg"&x_raw!=0) x<-paste("-",x,sep="")
  if(x=="0e+00") x<-"0"
  return(x)
}

#- - - - - -
## Custom ggplot theme

theme_Gaia<-function(){
  theme_bw() +
    theme(strip.text=element_text(size=12, colour="black", face="bold"),
          strip.background=element_rect(fill="#CAE1FF",colour="black"),
          axis.text=element_text(size=12, color="black"),
          axis.title=element_text(size=12, face="bold", color="black"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12, face="bold"),
          axis.line=element_line(color="black", linewidth=0.1))
}

#- - - - - -
## Functions for determining feature direction from SHAP values (i.e.: how a feature impacts the model)

# Approach 1: conventional (mean SHAP sign)

conventional_direction<-function(df, feature_col, shap_col){
  s<-df[[shap_col]]
  mean_s<-mean(s, na.rm=TRUE)
  if(abs(mean_s)<1e-10) return("neutral")
  ifelse(mean_s>0, "promoting", "mitigating")
}

# Approach 2: GAM derivative + Spearman fallback

gam_direction<-function(df, feature_col, shap_col){
  x<-df[[feature_col]]
  s<-df[[shap_col]]
  
  valid<-complete.cases(x, s)
  x<-x[valid]
  s<-s[valid]
  
  if(length(unique(x)) <= 1 || length(unique(s)) <= 1) return("undefined")
  
  # Binary or sparse
  if(length(unique(x)) <= 2 || quantile(x, 0.75, na.rm=TRUE) == 0){
    mean_diff<-mean(s[x > 0], na.rm=TRUE) - mean(s[x <= 0], na.rm=TRUE)
    return(ifelse(mean_diff > 0, "promoting",ifelse(mean_diff < 0, "mitigating", "neutral")))
  }
  
  # Continuous: GAM
  tryCatch({
    gam_m<-mgcv::gam(s ~ s(x, bs="cr", k=8), method="REML")
    derivs<-gratia::derivatives(gam_m, term="s(x)")
    mean_d<-mean(derivs$derivative, na.rm=TRUE)
    if(abs(mean_d) < 1e-10) return("neutral")
    return(ifelse(mean_d > 0, "promoting", "mitigating"))
  }, error=function(e){
    rho<-suppressWarnings(cor(x, s, method="spearman"))
    if(is.na(rho) || abs(rho) < 0.05) return("neutral")
    return(ifelse(rho > 0, "promoting", "mitigating"))
  })
}

# Approach 3: Pairwise bin concordance (Rcpp)

# Compiling C++ function
if(!exists("bin_pairwise_counts")){
  Rcpp::cppFunction('
Rcpp::List bin_pairwise_counts(NumericVector bx, NumericVector bs,
                               NumericVector bc){
  int B=bx.size();
  double pos=0.0, neg=0.0, total=0.0;
  for (int i=0; i < B; ++i){
    for (int j=i+1; j < B; ++j){
      double ci=bc[i], cj=bc[j];
      if(ci <= 0.0 || cj <= 0.0) continue;
      double pairs=ci * cj;
      if(bx[i] > bx[j]){
        total += pairs;
        if(bs[i] > bs[j]) pos += pairs;
        else if(bs[i] < bs[j]) neg += pairs;
      } else if(bx[j] > bx[i]){
        total += pairs;
        if(bs[j] > bs[i]) pos += pairs;
        else if(bs[j] < bs[i]) neg += pairs;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("pos")  =pos,
    Rcpp::Named("neg")  =neg,
    Rcpp::Named("total")=total);
}', depends="Rcpp")
}

pairwise_direction<-function(df, feature_col, shap_col,majority_threshold=0.55,n_quantile_bins=200){
  x<-df[[feature_col]]
  s<-df[[shap_col]]
  valid<-complete.cases(x, s)
  x<-x[valid]
  s<-s[valid]
  
  if(length(unique(x)) <= 1 || length(unique(s)) <= 1) return("undefined")
  
  # Binary / sparse
  if(length(unique(x)) <= 2 || quantile(x, 0.75, na.rm=TRUE) == 0){
    lv<-min(x)
    hv<-max(x)
    g0<-s[x==lv]
    g1<-s[x==hv]
    if(length(g0)==0 || length(g1)==0) return("neutral")
    g0s<-sort(g0)
    pos<-sum(findInterval(g1, g0s, left.open=TRUE))
    tot<-as.double(length(g0)) * as.double(length(g1))
    pp<-pos/tot
    if(pp>=majority_threshold) return("promoting")
    if(pp<=1-majority_threshold) return("mitigating")
    return("neutral")
  }
  
  # Continuous: binning
  probs<-seq(0, 1, length.out=n_quantile_bins+1)
  breaks<-unique(quantile(x, probs=probs, na.rm=TRUE, type=7))
  if(length(breaks) <= 2)
    breaks<-seq(min(x,na.rm=TRUE), max(x,na.rm=TRUE), length.out=3)
  bins<-cut(x, breaks=breaks, include.lowest=TRUE)
  
  bx<-as.numeric(tapply(x, bins, mean, na.rm=TRUE))
  bs<-as.numeric(tapply(s, bins, mean, na.rm=TRUE))
  bc<-as.numeric(tapply(s, bins, length))
  ok<-!is.na(bx) & !is.na(bs) & !is.na(bc)
  bx<-bx[ok]
  bs<-bs[ok]
  bc<-bc[ok]
  
  if(length(bx)<2) return("neutral")
  
  cnts<-bin_pairwise_counts(bx, bs, bc)
  pos_p<-as.numeric(cnts$pos)
  neg_p<-as.numeric(cnts$neg)
  tot_p<-as.numeric(cnts$total)
  if(tot_p <= 0) return("neutral")
  
  pp<-pos_p/tot_p
  np<-neg_p/tot_p
  if(pp>=majority_threshold) return("promoting")
  if(np>=majority_threshold) return("mitigating")
  return("neutral")
}

## Wrapper for computing SHAP directions from wide-format data
compute_shap_directions<-function(data_df, feature_cols,shap_prefix="shap_",
                                  methods=c("conventional","gam","pairwise"),
                                  threshold=0.55,n_bins=200){
  
  res<-lapply(feature_cols, function(f){
    sc<-paste0(shap_prefix, f)
    if(!sc %in% names(data_df)){
      message("SHAP column not found for: ", f); return(NULL)
    }
    
    row<-data.frame(feature=f,
                    mean_shap=round(mean(data_df[[sc]], na.rm=TRUE), 5),
                    median_shap=round(median(data_df[[sc]],na.rm=TRUE),5))
    
    if("conventional" %in% methods)
      row$conventional<-conventional_direction(data_df, f, sc)
    if("gam" %in% methods)
      row$gam_deriv<-gam_direction(data_df, f, sc)
    if("pairwise" %in% methods)
      row$pairwise_bins<-pairwise_direction(data_df, f, sc, majority_threshold=threshold,n_quantile_bins=n_bins)
    row
  })
  
  do.call(rbind, Filter(Negate(is.null), res))
}

## Wrapper for computing SHAP directions from long-format data
compute_shap_directions_long<-function(long_df,threshold,n_bins){
  long_df %>%
    group_split(feature) %>%
    map_dfr(function(df_feat){
      f<-as.character(df_feat$feature[1])
      tmp<-data.frame(x=df_feat$feature_value,  s=df_feat$shap_value)
      
      tibble(feature=f,
             n=nrow(df_feat),
             mean_shap=mean(df_feat$shap_value, na.rm=TRUE),
             median_shap=median(df_feat$shap_value, na.rm=TRUE),
             conventional=conventional_direction(tmp, "x", "s"),
             gam=gam_direction(tmp, "x", "s"),
             pairwise=pairwise_direction(tmp, "x", "s",majority_threshold=threshold,n_quantile_bins=n_bins))
    })
}

#- - - - - -
## not-in operator
`%ni%`<-Negate('%in%')

#----------------------------------------------------------------
#### Data import and preprocessing ####

# Simulating a dataset
n<-500
df<-data.frame(Age=rnorm(n, 50, 10),
  BMI=rnorm(n, 25, 5),
  Gender=factor(sample(c(0, 1), n, replace=TRUE)),
  Treatment=factor(sample(c("A", "B", "C"), n, replace=TRUE)))

# Ordinal outcome: 0,1,2
df$outcome<-factor(sample(0:2, n, replace=TRUE), ordered=TRUE)

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
dataset<-clean_data |> data.table::as.data.table()

## Creating vectors with the independent and dependent variable names
feature_cols<-colnames(dataset)[2:ncol(dataset)]
target_col<-"outcome"

## Identifying unordered factors
categorical_var<-var_type_tab %>% filter(feature %in% colnames(clean_data)) %>% 
  filter(Var_type %in% c("categorical")) %>% select(feature) %>% pull

## One-hot encoding of unordered factors
if(length(categorical_var) > 0) {
  dummies<-model.matrix(~ . -1, data=dataset[, ..categorical_var])
  dataset<-cbind(dataset[, !..categorical_var, with=FALSE], dummies)
}

## Updating feature columns after encoding
feature_cols<-setdiff(colnames(dataset), "outcome")

# Converting outcome to 0-based integers
dataset$outcome<-as.numeric(as.character(dataset$outcome))
num_class<-length(unique(dataset$outcome))

## Creating a stratified 70/30 train-test split
set.seed(seed)
data_split<-splitTools::partition(y=dataset[[target_col]],
  p=c(train=0.7, test=0.3),
  type="stratified",
  seed=seed)

## Preparing training and test datasets

# train dataset
X_train<-as.matrix(dataset[data_split$train, ..feature_cols])
X_train<-apply(X_train, 2, as.numeric)  # convert each column to numeric
y_train<-as.integer(dataset[data_split$train, get(target_col)])

# test dataset
X_test<-as.matrix(dataset[data_split$test, ..feature_cols])
X_test<-apply(X_test, 2, as.numeric)
y_test<-as.integer(dataset[data_split$test, get(target_col)])

#----------------------------------------------------------------
#### 10-fold nested cross-validation (CV) ####

#--------------------------
## Performing the outer CV on training set

# Creating the outer folds
outer_folds<-splitTools::create_folds(y_train, k=10, type="stratified", seed=seed)

# Creating empty lists to store results
outer_results<-list()
best_params_all<-list()

for(outer_idx in seq_along(outer_folds)){ # for each outer fold
  
  # selecting an outer fold
  val_idx<-outer_folds[[outer_idx]]
  train_idx_cv<-setdiff(seq_len(nrow(X_train)), val_idx)
  
  # creating training and test data sets for a given outer fold
  X_tr<-X_train[train_idx_cv, , drop=FALSE]
  y_tr<-y_train[train_idx_cv]
  X_val<-X_train[val_idx, , drop=FALSE]
  y_val<-y_train[val_idx]

  #--------------------------
  ## Performing the inner CV for tuning
  
  # Creating the inner folds
  inner_folds<-splitTools::create_folds(y_tr, k=10, type="stratified", seed=seed)
  
  # Creating random parameter space for hyperparameter tuning
  param_list<-expand.grid(subsample=seq(0.5, 1, 0.25),
    colsample_bytree=seq(0.5, 1, 0.25),
    min_child_weight=c(1, 5, 10),
    learning_rate=c(0.05, 0.1, 0.3),
    max_depth=c(3, 5, 7)) %>% dplyr::slice_sample(n=30, replace=TRUE) # Limiting space to 30 combinations for computational efficiency and environmental sustainability considerations
  
  # Initializing best parameters
  best_auc<--Inf
  best_params<-NULL
  
  for (i in 1:nrow(param_list)) {

    # Initializing cross validation
    xgb_cv<-mlexperiments::MLCrossValidation$new(learner=mllrnrs::LearnerXgboost$new(metric_optimization_higher_better=TRUE),
      fold_list=inner_folds,
      ncores=2,
      seed=0)

    # Defining the learner arguments
    xgb_cv$learner_args<-c( as.list(param_list[i, ]), list(objective="multi:softprob", eval_metric="mlogloss", nrounds=100L,num_class=num_class))
    xgb_cv$performance_metric<-mlexperiments::metric("mae")
    xgb_cv$set_data(x=X_tr, y=y_tr)
    
    # Executing CV tuning
    res_cv<-xgb_cv$execute()
    
    # Evaluating using MAE (ordinal-aware)
    fold_mae<-sapply(res_cv$results$folds, function(f){
      pred<-max.col(f[[5]]$predictions) - 1
      mean(abs(pred - f[[5]]$labels))
    })
    mean_mae<-mean(fold_mae)
    
    if(!is.na(mean_mae) && mean_mae > best_score){  # minimizing MAE
      best_score<-mean_mae
      best_params<-params
    }
  }
  
  # If no params are selected, fallback to default
  if (is.null(best_params)){
    best_params<-data.frame(subsample=1,
      colsample_bytree=1,
      min_child_weight=1,
      learning_rate=0.1,
      max_depth=3)
  }
  
  best_params_all[[outer_idx]]<-best_params
  
  ## Evaluating on outer fold
  
  # Building the model
  dtrain<-xgboost::xgb.DMatrix(data=X_tr, label=y_tr)
  dval<-xgboost::xgb.DMatrix(data=X_val, label=y_val)
  
  xgb_outer<-xgboost::xgb.train(params=c(as.list(best_params),list(objective="multi:softprob",num_class=num_class,eval_metric="mlogloss")),
    data=dtrain,
    nrounds=100,
    verbose=0)
  
  # Assessing the performance
  val_pred<-predict(xgb_outer, newdata=dval)
  val_pred<-matrix(val_pred, ncol=num_class, byrow=TRUE)
  pred_class<-max.col(val_pred) - 1
  MAE<-mean(abs(pred_class - y_val))
  
  ## Saving CV results
  outer_results[[outer_idx]]<-data.frame(Fold=outer_idx, MAE=MAE)
}

# Combining results from outer CV into a summary table
outer_summary<-dplyr::bind_rows(outer_results)
outer_summary

#--------------------------
## Retraining final model on full training set

# Picking the best params (e.g., from highest outer AUROC)
best_idx<-which.min(outer_summary$MAE)
final_params<-best_params_all[[best_idx]]

# Retraining final model
dtrain_full<-xgboost::xgb.DMatrix(data=X_train, label=y_train)
xgb_final<-xgboost::xgb.train(params=c(as.list(final_params),list(objective="multi:softprob",num_class=num_class,eval_metric="mlogloss")),
  data=dtrain_full,
  nrounds=100,
  verbose=1)

# Saving the final model
saveRDS(xgb_final, file="XGBoost_model.rds")
xgb.save(xgb_final,paste0("XGBoost_model_",Sys.Date(),".ubj"))

#--------------------------
## Final evaluation on holdout test dataset

# Predicted probabilities
pred_prob<-predict(xgb_final, newdata=X_test)
pred_prob<-matrix(pred_prob, ncol=num_class, byrow=TRUE)

# Predicted classes
pred_class<-max.col(pred_prob) - 1
pred_f<-factor(pred_class, levels=sort(unique(y_test)))
true_f<-factor(y_test, levels=sort(unique(y_test)))

# Calculating performance metrics
multi_acc<-yardstick::accuracy_vec(truth=true_f, estimate=pred_f)                 # Accuracy
multi_F1<-yardstick::f_meas_vec(truth=true_f, estimate=pred_f, estimator="macro") # Macro F1 score
multi_bal_acc<-yardstick::bal_accuracy_vec(truth=true_f, estimate=pred_f)         # Balanced accuracy
precision_macro<-yardstick::precision_vec(true_f, pred_f, estimator="macro")      # Macro precision
recall_macro<-yardstick::recall_vec(true_f, pred_f, estimator="macro")            # Macro recall
sens_macro<-yardstick::sens_vec(true_f, pred_f, estimator="macro")                # Macro sensitivity
spec_macro<-yardstick::spec_vec(true_f, pred_f, estimator="macro")                # Macro specificity
multi_kappa<-caret::confusionMatrix(pred_f, true_f)$overall["Kappa"]              # Cohen's kappa unweighted
QWK<-irr::kappa2(cbind(pred_class, y_test), weight="squared")$value               # Quadratic weighted kappa
multi_logloss<-MLmetrics::MultiLogLoss(y_pred=pred_prob, y_true=y_test)           # Multi-class Logloss
Multi_MAE<-yardstick::mae_vec(truth=y_test, estimate=pred_class)                  # MAE

colnames(pred_prob)<-levels(true_f)
multi_AUC<-as.vector(pROC::multiclass.roc(true_f, pred_prob)$auc)                 # Multi-class AUC

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
#### Association analysis across all classes
#----------------------------------------------------------------

#----------------------------------------------------------------
#### Univariate association analysis ####

var_test<-colnames(X_test)

Association_tab_res<-c()

for(n_var in 1:length(var_test)){
  
  test_tmpo<-as_tibble(X_test) %>% select(which(colnames(X_test)==var_test[n_var]))
  colnames(test_tmpo)<-"tmp"
  test_tmpo<-bind_cols(y_test,test_tmpo)
  colnames(test_tmpo)[1]<-"outcome"
  test_tmpo$outcome<-factor(test_tmpo$outcome,ordered=T,levels=sort(unique(test_tmpo$outcome)))
  
  clm_model<-ordinal::clm(data=test_tmpo, formula=outcome~.)
  pred_class<-predict(clm_model, type="class")$fit
  pred_class<-as.numeric(pred_class)
  true_class<-as.numeric(test_tmpo$outcome)
  
  mod_res<-broom::tidy(clm_model,exponentiate=T,conf.int=T) %>%
    filter(coef.type!="intercept") %>%
    mutate(logloss=-mean(rowSums(model.matrix(~ outcome - 1, data=test_tmpo) * log(pmax(pmin(predict(clm_model, type="prob")$fit, 1 - 1e-15), 1e-15)))),
           nom_test=as_tibble(data.frame(term=rownames(ordinal::nominal_test(clm_model)),
                                         ordinal::nominal_test(clm_model))) %>% filter(term=="tmp") %>% pull(Pr..Chi.),
           scale_test=as_tibble(data.frame(term=rownames(ordinal::scale_test(clm_model)),
                                           ordinal::scale_test(clm_model))) %>% filter(term=="tmp") %>% pull(Pr..Chi.),
           accuracy=mean(pred_class == true_class),
           MAE=mean(abs(pred_class - true_class)),
           kappa=irr::kappa2(data.frame(
             pred=as.numeric(pred_class),
             obs =as.numeric(test_tmpo$outcome)),
             weight="squared")$value) %>%
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
contr<-predict(xgb_final, as.matrix(X_test), predcontrib=TRUE)

# Reshaping predicted probabilities into n_samples x num_class
pred_prob_mat<-matrix(pred_prob, ncol=num_class, byrow=TRUE)

# Aggregating SHAP across classes (keeping SHAP distributions), weighted by their probability
if(is.list(contr)==T){
  shap<-list()
  for(i in 1:length(contr)){

    # Get SHAP for class i
    shap_i<-as_tibble(contr[[i]])
    
    # Multiply each row by predicted probability of that class
    shap_i_weighted<-shap_i * pred_prob_mat[, i]
    
    shap[[i]]<-shap_i_weighted
    
  }
  shap<-bind_rows(shap)
  shap_contrib<-as.data.table(shap)
}else{
  shap<-as_tibble(contr)
  shap_contrib<-as.data.table(contr)
}

colnames(shap_contrib)<-c("id","class","feature","shap_value")

#- - - - - - - - - -
## Removing BIAS term
shap_contrib<-shap_contrib[, !grepl("bias", names(shap_contrib), ignore.case=TRUE), with=FALSE]
shap_contrib<-shap_contrib[, !grepl("(Intercept)", names(shap_contrib), ignore.case=TRUE), with=FALSE]
shap_contrib<-shap_contrib %>% filter(feature %ni% c("(Intercept)","bias","BIAS","Bias"))

#- - - - - - - - - -
## Computing mean SHAP score
mean_shap_score<-shap_contrib %>% group_by(feature) %>%summarize(mean=mean(abs(shap_value),na.rm=T)) %>% pivot_wider(names_from=feature,values_from=mean)

#- - - - - - - - - -
## Computting summary statistics per feature
Shap_val<-as_tibble(shap_contrib) %>% 
  group_by(feature) %>% 
  mutate(IC_2.5=quantile(abs(shap_value),0.025),
         IC_97.5=quantile(abs(shap_value),0.975),
         mean_val=mean(abs(shap_value),na.rm=T)) %>%
  select(feature,mean_val,IC_2.5,IC_97.5) %>% 
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
tempopo<-X_test %>%
  as_tibble() %>%
  rowid_to_column("id") %>%
  pivot_longer(-id, names_to="feature", values_to="feature_value")

tempopo<-left_join(tempopo, shap_contrib, by=c("id", "feature"))

#- - - - - - - - - -
## SHAP summary and merging labels
tempopo<-left_join(tempopo,Shap_val,by=base::intersect(colnames(tempopo),colnames(Shap_val))) %>% distinct

test3<-Shap_val %>% arrange(desc(as.numeric(mean_val)))

#- - - - - - - - - -
## SHAP directionality
direction_impact<-compute_shap_directions_long(tempopo,threshold=0.55, n_bins=200) %>% select(-c(n,mean_shap,median_shap)) %>%
  filter(feature %ni% c("(Intercept)","BIAS","Bias")) %>%
  # Consensus: majority vote across three methods
  mutate(consensus=apply(cbind(conventional,gam, pairwise), 1,
                         function(x){
                           tbl<-sort(table(x), decreasing=TRUE)
                           if(tbl[1]>=2) names(tbl)[1] else "uncertain"
                         }))

test<-left_join(tempopo,direction_impact,by="feature")

test<-left_join(test,Shap_val,by=base::intersect(colnames(test),colnames(Shap_val)))

test_tmp<-test %>% 
  mutate(impact=consensus) %>% # selecting the consensus-based approach, adapt if needed 
  select(feature, mean_val,IC_2.5,IC_97.5,impact,mean_SHAP) %>% 
  distinct() %>%
  mutate(mean_val=as.numeric(mean_val),
         IC_2.5=as.numeric(IC_2.5),
         IC_97.5=as.numeric(IC_97.5)) %>%
  mutate(mean_val=if_else(impact=="mitigating",-mean_val,mean_val),
         IC_2.5=if_else(impact=="mitigating",-IC_2.5,IC_2.5),
         IC_97.5=if_else(impact=="mitigating",-IC_97.5,IC_97.5))

test_tmp<-test_tmp %>% arrange(desc(mean_val)) %>% filter(mean_val!=0) %>% filter(!is.na(mean_val))

#- - - - - - - - - -
## Plotting SHAP

# reordoring factors
ordre_def1<-test_tmp %>% filter(mean_val>=0) %>% arrange(desc(IC_2.5))
ordre_def2<-test_tmp %>% filter(mean_val<0) %>% arrange(desc(IC_97.5))
ordre_def<-bind_rows(ordre_def1,ordre_def2)

test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def$feature)))
test_tmp$impact<-factor(test_tmp$impact,
                        levels=c("mitigating","neutral","promoting","uncertain","undefined"),
                        labels=c("mitigating predictor","neutral","promoting predictor","uncertain","undefined"))

# Plot 1: error bar
plot1<-ggplot(test_tmp,aes(y=feature,x=mean_val,col=impact))+
  geom_point()+
  geom_errorbar(aes(xmin=IC_2.5, xmax=IC_97.5))+
  scale_x_continuous("mean |SHAP value|",label=function(x) abs(x))+
  geom_text(aes(x=IC_97.5,y=feature,label=abs(IC_97.5)),size=5,position=position_nudge(x=if_else(test_tmp$IC_97.5>=0,
                                                                                                 round(max(test_tmp$IC_97.5,na.rm=T)/10),
                                                                                                 -round(max(test_tmp$IC_97.5,na.rm=T)/10))))+
  scale_y_discrete("")+
  theme_bw()+
  scale_color_manual("Direction:",na.value="white",
                     values=rev(c("neutral"="grey80","promoting predictor"="#C35C33","mitigating predictor"="#40B696",
                                  "uncertain"="black","undefined"="grey50")))+
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="bottom",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

# Exporting plot 1
ggsave(plot1,file=paste("Across all classes_SHAP_errorbar_",Sys.Date(),".pdf",sep=""),dpi=600,width=60,height=30,units="cm",limitsize=F)

# Plot 2: barplot
ordre_def<-test_tmp %>% arrange(desc(mean_val))
test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def$feature)))

plot2<-ggplot(test_tmp,aes(y=feature,x=mean_val,fill=impact))+
  geom_bar(stat="identity",col='black')+
  geom_text(aes(x=mean_val,y=feature,label=mean_SHAP),size=6,position=position_nudge(x=if_else(test_tmp$mean_val>=0,
                                                                                               round(max(test_tmp$mean_val,na.rm=T)/10),
                                                                                               -round(max(test_tmp$mean_val,na.rm=T)/10))))+
  scale_x_continuous("mean |SHAP value|",label=function(x) abs(x))+
  scale_y_discrete("")+
  theme_bw()+
  scale_fill_manual("Direction:",na.value="white",
                    values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2","neutral"="white","uncertain"="grey","undefined"="grey50"))+
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text=element_text(size=16,color="black", face="bold"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="bottom",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

# Exporting plot 2
ggsave(plot2,file=paste("Across all classes_SHAP_barplot_",Sys.Date(),".pdf",sep=""),dpi=600,width=60,height=30,units="cm",limitsize=F)

# Plot 3: SHAP direction comparison
dir_long<-direction_impact %>%
  select(feature,conventional,gam,pairwise,consensus) %>%
  pivot_longer(-feature, names_to="method", values_to="direction") %>%
  mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","undefined","uncertain")),
         method=factor(method,
                       levels=c("conventional","gam","pairwise","consensus"),
                       labels=c("Conventional\n(mean sign)","GAM derivative","Pairwise bins","Consensus")))

plot3<-ggplot(dir_long,aes(x=method,y=feature,fill=direction))+
  geom_tile(color="black")+
  scale_fill_manual("Direction:",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","uncertain"="grey25",
                                          "mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
  scale_x_discrete("",expand=c(0,0))+
  scale_y_discrete("",expand=c(0,0))+
  theme_Gaia()+
  theme(legend.position="top")

# Exporting plot
ggsave(plot3,file=paste("Across all classes_SHAP direction comparison_",Sys.Date(),".pdf",sep=""),dpi=600,width=60,height=30,units="cm",limitsize=F)

## Plot 4: SHAP percentage contribution 

test_tmp<-test_tmp %>% mutate(SHAP_per=abs(mean_val)/sum(abs(mean_val))*100) # calculating the relative contribution to model output based on mean |SHAP|

ordre_def2<-test_tmp %>% arrange(desc(SHAP_per))
test_tmp$feature<-factor(test_tmp$feature, levels=rev(unique(ordre_def2$feature)))

plot4<-ggplot(test_tmp, aes(y=feature, x=SHAP_per, fill=impact))+
  geom_bar(stat="identity", col="black")+
  geom_text(aes(x=SHAP_per, y=feature, label=paste0(signif(abs(SHAP_per),3))), size=6,
            position=position_nudge(x=if_else(test_tmp$SHAP_per>=0,
                                              round(max(test_tmp$SHAP_per,na.rm=TRUE)/10),
                                              -round(max(test_tmp$SHAP_per,na.rm=TRUE)/10))))+
  scale_x_continuous("Relative contribution to model output (%)", label=function(x) paste0(abs(x)))+
  scale_y_discrete("Feature")+
  scale_fill_manual("Direction:", na.value="white",
                    values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2","neutral"="white"))+
  theme_Gaia()+
  theme(legend.position="bottom")

# Exporting plot 4
ggsave(plot4,file=paste("Across all classes_SHAP_contribution_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)

#- - - - - - - - - -
## SHAP and CLM comparison

tmp_CLM<-Association_tab_res %>% select(term,OR,direction,p.value,R2_Nagelkerke)
colnames(tmp_CLM)<-c("feature","OR [95% CI]","CLM - association direction","CLM - p-value","CLM - Pseudo-R2")

test_SHAP_tmp<-test %>%  
  select(feature,mean_SHAP,conventional,gam,pairwise,consensus,mean_val) %>% 
  distinct %>%
  arrange(desc(abs(as.numeric(mean_val)))) %>%
  rowwise %>%
  mutate_if(is.numeric,test_format) %>%
  ungroup %>%
  select(-c(mean_val))

colnames(test_SHAP_tmp)<-c("feature","mean_abs_SHAP [95% CI]",
                           "conventional SHAP direction","GAM-based SHAP direction","Pairwise-based SHAP direction","Consensus-based SHAP direction")

Compa_CLM_SHAP<-left_join(test_SHAP_tmp,tmp_CLM,by="feature")

#----------------------------------------------------------------
#### Saving SHAP results

# SHAP values
fwrite(Shap_val,paste("Across all classes_Shap_val_",Sys.Date(),".csv",sep=""), sep=";", row.names=FALSE)
fwrite(test_tmp,paste("Across all classes_test_tmp_",Sys.Date(),".csv",sep=""), sep=";", row.names=FALSE)
fwrite(test,paste("Across all classes_test_",Sys.Date(),".txt",sep=""), sep=";", row.names=FALSE)
fwrite(shap,paste("Across all classes_shap_",Sys.Date(),".txt",sep=""), sep=";", row.names=FALSE)
fwrite(shap_score_sub,paste("Across all classes_shap_score_sub_",Sys.Date(),".txt",sep=""), sep=";", row.names=FALSE)
fwrite(Association_tab_res,paste("Across all classes_Univariate CLM regression_",Sys.Date(),".csv",sep=""), sep=";", row.names=FALSE)
fwrite(Compa_CLM_SHAP,paste("Across all classes_Compa_CLM_SHAP_",Sys.Date(),".csv",sep=""), sep=";", row.names=FALSE)

#----------------------------------------------------------------
#### Association analysis for each class
#----------------------------------------------------------------

shap_contrib_save<-shap_contrib

for(n_clust in 1:num_class){
  
  #----------------------------------------------------------------
  #### Univariate association analysis ####  
  
  y_test_class<-if_else(y_test==sort(unique(y_test))[n_clust],1,0)
  
  var_test<-colnames(X_test)
  
  Association_tab_res<-c()
  
  for(n_var in 1:length(var_test)){
    
    test_tmpo<-as_tibble(X_test) %>% select(which(colnames(X_test)==var_test[n_var]))
    colnames(test_tmpo)<-"tmp"
    test_tmpo<-bind_cols(y_test_class,test_tmpo)
    colnames(test_tmpo)[1]<-"outcome"
    
    glm_model<-glm(outcome ~ tmp, data=test_tmpo,family="binomial")
    mod_res<-broom::tidy(glm_model,exponentiate=T,conf.int=T) %>%
      filter(term!="(Intercept)") %>%
      mutate(term=var_test[n_var],
             direction=case_when(conf.low>=1~"positive association",
                                 conf.high<1~"negative association",
                                 T~"uncertain association")) %>% 
      bind_cols(performance::performance(glm_model)) %>%
      rowwise %>% 
      mutate_if(is.numeric,test_format) %>%
      ungroup %>%
      mutate(OR=paste(estimate," [",conf.low,"; ",conf.high,"]",sep="")) %>%
      mutate(OR=if_else(OR==" [; ]","",OR))
    
    Association_tab_res<-bind_rows(Association_tab_res,mod_res)
    
  }
  
  #----------------------------------------------------------------
  #### SHAP analysis ####
  
  #- - - - - - - - - -
  ## Calculating SHAP values
  contr<-predict(xgb_final, as.matrix(X_test), predcontrib=TRUE)
  shap<-as_tibble(contr[, n_clust, ])
  shap_contrib<-as.data.table(contr[, n_clust, ])
  
  #- - - - - - - - - -
  ## Removing BIAS term
  shap_contrib<-shap_contrib[, !grepl("bias", names(shap_contrib), ignore.case=TRUE), with=FALSE]
  shap_contrib<-shap_contrib[, !grepl("(Intercept)", names(shap_contrib), ignore.case=TRUE), with=FALSE]
  
  #- - - - - - - - - -
  ## Computing mean SHAP score
  mean_shap_score<-colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing=T)]
  
  #- - - - - - - - - -
  ## Reshaping SHAP values to long format
  shap_score<-shap_contrib
  shap_score_sub<-as.data.table(shap_score)
  shap_score_sub<-shap_score_sub[, names(mean_shap_score), with=F]
  shap_score_long<-melt.data.table(shap_score_sub, measure.vars=colnames(shap_score_sub))
  
  #- - - - - - - - - -
  ## Matching feature values
  fv_sub<-as.data.table(X_test)[, names(mean_shap_score), with=F]
  fv_sub_long<-melt.data.table(fv_sub, measure.vars=colnames(fv_sub))
  fv_sub_long[, stdfvalue := std1(value), by="variable"]
  names(fv_sub_long)<-c("variable", "rfvalue", "stdfvalue" )
  
  #- - - - - - - - - -
  ## Merging  SHAP and feature values
  shap_long2<-cbind(shap_score_long, fv_sub_long[,c('rfvalue','stdfvalue')])
  shap_long2[, mean_value := mean(abs(value)), by=variable]
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
  shap_tempo<-shap %>% 
    rowid_to_column("id") %>%
    pivot_longer(-id, names_to="feature", values_to="shap_value")
  
  tempopo<-X_test %>%
    as_tibble() %>%
    rowid_to_column("id") %>%
    pivot_longer(-id, names_to="feature", values_to="feature_value")
  
  tempopo<-left_join(tempopo, shap_tempo, by=c("id", "feature"))
  
  #- - - - - - - - - -
  ## SHAP summary and merging labels
  tempopo<-left_join(tempopo,Shap_val,by=base::intersect(colnames(tempopo),colnames(Shap_val))) %>% distinct
  
  test3<-Shap_val %>% arrange(desc(as.numeric(mean_val)))
  
  #- - - - - - - - - -
  ## SHAP directionality

  direction_impact<-compute_shap_directions_long(tempopo,threshold=0.55, n_bins=200) %>% select(-c(n,mean_shap,median_shap)) %>%
    filter(feature %ni% c("(Intercept)","BIAS","Bias")) %>%
    # Consensus: majority vote across three methods
    mutate(consensus=apply(cbind(conventional,gam, pairwise), 1,
                           function(x){
                             tbl<-sort(table(x), decreasing=TRUE)
                             if(tbl[1]>=2) names(tbl)[1] else "uncertain"
                           }))
  
  test<-left_join(tempopo,direction_impact,by="feature") %>% 
    mutate(direction=consensus) # selecting the consensus-based approach, adapt if needed
  
  test<-left_join(test,Shap_val,by=base::intersect(colnames(test),colnames(Shap_val)))
  
  test_tmp<-test %>% mutate(impact=direction) %>%
    select(feature, mean_val,IC_2.5,IC_97.5,impact,mean_SHAP) %>% 
    distinct() %>%
    mutate(mean_val=as.numeric(mean_val),
           IC_2.5=as.numeric(IC_2.5),
           IC_97.5=as.numeric(IC_97.5)) %>%
    mutate(mean_val=if_else(impact=="mitigating",-mean_val,mean_val),
           IC_2.5=if_else(impact=="mitigating",-IC_2.5,IC_2.5),
           IC_97.5=if_else(impact=="mitigating",-IC_97.5,IC_97.5))
  
  test_tmp<-test_tmp %>% arrange(desc(mean_val)) %>% filter(mean_val!=0) %>% filter(!is.na(mean_val))
  
  test_tmp_save<-test_tmp
  test_tmp<-test_tmp %>% arrange(desc(abs(mean_val))) %>% slice(1:30)
  
  #- - - - - - - - - -
  ## Plotting SHAP
  
  # reordoring factors
  ordre_def1<-test_tmp %>% filter(mean_val>=0) %>% arrange(desc(IC_2.5))
  ordre_def2<-test_tmp %>% filter(mean_val<0) %>% arrange(desc(IC_97.5))
  ordre_def<-bind_rows(ordre_def1,ordre_def2)
  
  test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def$feature)))
  test_tmp$impact<-factor(test_tmp$impact,
                          levels=c("mitigating","neutral","promoting","uncertain","undefined"),
                          labels=c("mitigating predictor","neutral","promoting predictor","uncertain","undefined"))
  
  # Plot 1: error bar
  plot1<-ggplot(test_tmp,aes(y=feature,x=mean_val,col=impact))+
    geom_point()+
    geom_errorbar(aes(xmin=IC_2.5, xmax=IC_97.5))+
    scale_x_continuous("mean |SHAP value|",label=function(x) abs(x))+
    geom_text(aes(x=IC_97.5,y=feature,label=abs(IC_97.5)),size=5,position=position_nudge(x=if_else(test_tmp$IC_97.5>=0,
                                                                                                     round(max(test_tmp$IC_97.5,na.rm=T)/10),
                                                                                                     -round(max(test_tmp$IC_97.5,na.rm=T)/10))))+
    scale_y_discrete("")+
    theme_bw()+
  scale_color_manual("Direction:",na.value="white",
                     values=rev(c("neutral"="grey80","promoting predictor"="#C35C33","mitigating predictor"="#40B696",
                                  "uncertain"="black","undefined"="grey50")))+
    theme(strip.text.x=element_text(size=16, colour="black", angle=0),
          strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
          axis.text.y=element_text(size=16,color="black"),
          axis.text.x=element_text(size=16,color="black"),
          axis.title=element_text(size=16,face="bold",color="black"),
          legend.text=element_text(size=16, face="bold"),
          legend.title=element_text(size=16, face="bold"),
          legend.position="bottom",
          axis.line=element_line(color="black",size=0.1, linetype="solid"))
  
  # Exporting plot 1
  ggsave(plot1,file=paste("Class ",n_clust,"_SHAP_errorbar_",Sys.Date(),".pdf",sep=""),dpi=600,width=60,height=30,units="cm",limitsize=F)
  
  # Plot 2: barplot
  ordre_def<-test_tmp %>% arrange(desc(mean_val))
  test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def$feature)))
  
  plot2<-ggplot(test_tmp,aes(y=feature,x=mean_val,fill=impact))+
    geom_bar(stat="identity",col='black')+
    geom_text(aes(x=mean_val,y=feature,label=mean_SHAP),size=6,position=position_nudge(x=if_else(test_tmp$mean_val>=0,
                                                                                                   round(max(test_tmp$mean_val,na.rm=T)/10),
                                                                                                   -round(max(test_tmp$mean_val,na.rm=T)/10))))+
    scale_x_continuous("mean |SHAP value|",label=function(x) abs(x))+
    scale_y_discrete("")+
    theme_bw()+
    scale_fill_manual("Direction:",na.value="white",
                      values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2","neutral"="white","uncertain"="grey","undefined"="grey50"))+
    theme(strip.text.x=element_text(size=16, colour="black", angle=0),
          strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
          axis.text=element_text(size=16,color="black", face="bold"),
          axis.title=element_text(size=16,face="bold",color="black"),
          legend.text=element_text(size=16, face="bold"),
          legend.title=element_text(size=16, face="bold"),
          legend.position="bottom",
          axis.line=element_line(color="black",size=0.1, linetype="solid"))
  
  # Exporting plot 2
  ggsave(plot2,file=paste("Class ",n_clust,"_SHAP_barplot_",Sys.Date(),".pdf",sep=""),dpi=600,width=60,height=30,units="cm",limitsize=F)
  
  # Plot 3: direction comparison
  dir_long<-direction_impact %>%
    select(feature,conventional,gam,pairwise,consensus) %>%
    pivot_longer(-feature, names_to="method", values_to="direction") %>%
    mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","undefined","uncertain")),
           method=factor(method,
                         levels=c("conventional","gam","pairwise","consensus"),
                         labels=c("Conventional\n(mean sign)","GAM derivative","Pairwise bins","Consensus")))
  
  plot3<-ggplot(dir_long,aes(x=method,y=feature,fill=direction))+
    geom_tile(color="black")+
    scale_fill_manual("Direction:",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","uncertain"="grey25",
                                            "mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
    scale_x_discrete("",expand=c(0,0))+
    scale_y_discrete("",expand=c(0,0))+
    theme_Gaia()+
    theme(legend.position="top")
  
  # Exporting plot
  ggsave(plot3,file=paste("Class ",n_clust,"_SHAP direction comparison_",Sys.Date(),".pdf",sep=""),dpi=600,width=60,height=30,units="cm",limitsize=F)
  
  ## Plot 4: SHAP percentage contribution 
  
  test_tmp<-test_tmp %>% mutate(SHAP_per=abs(mean_val)/sum(abs(mean_val))*100) # calculating the relative contribution to model output based on mean |SHAP|
  
  ordre_def2<-test_tmp %>% arrange(desc(SHAP_per))
  test_tmp$feature<-factor(test_tmp$feature, levels=rev(unique(ordre_def2$feature)))
  
  plot4<-ggplot(test_tmp, aes(y=feature, x=SHAP_per, fill=impact))+
    geom_bar(stat="identity", col="black")+
    geom_text(aes(x=SHAP_per, y=feature, label=paste0(signif(abs(SHAP_per),3))), size=6,
              position=position_nudge(x=if_else(test_tmp$SHAP_per>=0,
                                                round(max(test_tmp$SHAP_per,na.rm=TRUE)/10),
                                                -round(max(test_tmp$SHAP_per,na.rm=TRUE)/10))))+
    scale_x_continuous("Relative contribution to model output (%)", label=function(x) paste0(abs(x)))+
    scale_y_discrete("Feature")+
    scale_fill_manual("Direction:", na.value="white",
                      values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2","neutral"="white"))+
    theme_Gaia()+
    theme(legend.position="bottom")
  
  # Exporting plot 4
  ggsave(plot4,file=paste("Class ",n_clust,"_SHAP_contribution_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)
  
  #- - - - - - - - - -
  ## SHAP and GLM comparison
  
  tmp_GLM<-Association_tab_res %>% select(term,OR,direction,p.value,Log_loss)
  colnames(tmp_GLM)<-c("feature","OR [95% CI]","GLM - association direction","GLM - p-value","GLM - log loss")
  
  test_SHAP_tmp<-test %>%  
    select(feature,mean_SHAP,conventional,gam,pairwise,consensus,mean_val) %>% 
    distinct %>%
    arrange(desc(abs(as.numeric(mean_val)))) %>%
    rowwise %>%
    mutate_if(is.numeric,test_format) %>%
    ungroup %>%
    select(-c(mean_val))
  
  colnames(test_SHAP_tmp)<-c("feature","mean_abs_SHAP [95% CI]",
                             "conventional SHAP direction","GAM-based SHAP direction","Pairwise-based SHAP direction","Consensus-based SHAP direction")
  
  Compa_GLM_SHAP<-left_join(test_SHAP_tmp,tmp_GLM,by="feature")
  
  #----------------------------------------------------------------
  #### Saving SHAP results
  
  # SHAP values
  fwrite(Shap_val,paste("Class ",n_clust,"_Shap_val_",Sys.Date(),".csv",sep=""), sep=";", row.names=FALSE)
  fwrite(test_tmp,paste("Class ",n_clust,"_test_tmp_",Sys.Date(),".csv",sep=""), sep=";", row.names=FALSE)
  fwrite(test,paste("Class ",n_clust,"_test_",Sys.Date(),".txt",sep=""), sep=";", row.names=FALSE)
  fwrite(shap,paste("Class ",n_clust,"_shap_",Sys.Date(),".txt",sep=""), sep=";", row.names=FALSE)
  fwrite(shap_score_sub,paste("Class ",n_clust,"_shap_score_sub_",Sys.Date(),".txt",sep=""), sep=";", row.names=FALSE)
  fwrite(Association_tab_res,paste("Class ",n_clust,"_Univariate CLM regression_",Sys.Date(),".csv",sep=""), sep=";", row.names=FALSE)
  fwrite(Compa_GLM_SHAP,paste("Cluster ",n_clust,"_Compa_GLM_SHAP_",Sys.Date(),".csv",sep=""), sep=";", row.names=FALSE)
  
}
