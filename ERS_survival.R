#-------------------------------------------------------------------------------
## Reproducible and generalisable script for calculating an exposomic risk score (ERS) for a survival outcome using XGBoost with nested cross-validation
#
# Method: Double Machine Learning (DML) - residualizing exposures on adjustment variables (cohort membership and/or covariates) to preserve the survival outcome structure.
#
# Steps:
# - Step 1 [optional, multi-cohort]: W_j = X_j - E[X_j | cohort] -> each exposure is cleaned of between-cohort differences
# - Step 2: W2_j = W_j - E[W_j | covariates] -> each exposure is additionally cleaned of covariate effects
# - Step 3: ERS = f(W2) -> ERS = predicted log partial hazard from doubly-adjusted exposures
#
# How to use the ERS in a final model:
# - coxph(Surv(time, event) ~ ers + age + sex + bmi, data=sim_data)
# - or stratify into tertiles for HR and KM curves.
#
# Settings to adjust before running:
# - use_cohort: TRUE/FALSE
# - cov_names: names of covariate adjustment variables
# - surv_time_col: name of the time-to-event column
# - surv_event_col: name of the event indicator column (1=event, 0=censored)
#-------------------------------------------------------------------------------

#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing and loading packages
pack_needed<-c("data.table","tidyverse","mllrnrs","mlsurvlrnrs","broom","doParallel","foreach","splitTools","conflicted",
               "grid","gridExtra","RColorBrewer","mlbench","mlexperiments","caret","MLmetrics","patchwork","performance",
               "xgboost","parallel","here","scales","dplyr","ggplot2","tidyr","tibble","mgcv","ggforce","gratia",
               "Rcpp","Metrics","MASS","survival","survminer","ggkm","R6","kdry","pROC","iml","ggsurvfit",
               "gtsummary","GGally","cmprsk","tidycmprsk","ggstats")

if(!("ggkm"%in%.packages(all.available=TRUE))){
  pak::pak("sachsmc/ggkm")
}

is_installed<-pack_needed %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed==FALSE)){
  install.packages(pack_needed[!is_installed], repos="http://cran.us.r-project.org")
}
invisible(lapply(pack_needed, library, character.only=TRUE))

## Preventing package conflicts
conflict_prefer("select","dplyr")
conflict_prefer("filter","dplyr")
conflict_prefer("slice","dplyr")
conflict_prefer("alpha","scales")
conflicts_prefer(caret::lift)

## Setting the working directory
here::here("XGBoost - ERS calculation - survival outcome")

## Setting seed
seed<-123

## Setting cores
if(isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_")))){
  ncores<-2L
}else{
  ncores<-ifelse(test=parallel::detectCores()>4,yes=4L,no=ifelse(test=parallel::detectCores()<2L,yes=1L,no=parallel::detectCores()))
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

#- - - - - -
## Computing concordance-based metrics between an exposure (or its residual) and the survival outcome
#
# Interpretation: a C-index > 0.5 means the exposure tends to rank higher-risk individuals correctly.

compute_survival_metrics<-function(x_vec, surv_time, surv_event){
  
  x_num<-as.numeric(x_vec)
  
  ## C-index (Harrell's concordance index)
  cindex<-tryCatch({
    srv<-Surv(time=surv_time, event=surv_event)
    round(survival::concordance(srv~x_num)$concordance, 3)
  }, error=function(e) NA)
  
  ## Log partial hazard association via univariate Cox (log hazard ratio + p-value)
  cox_res<-tryCatch({
    df_tmp<-data.frame(time=surv_time, event=surv_event, x=x_num)
    fit<-survival::coxph(Surv(time, event)~x, data=df_tmp)
    s<-summary(fit)
    list(logHR=round(coef(fit),3),
         HR=round(exp(coef(fit)),3),
         p=round(s$logtest["pvalue"],3))
  }, error=function(e) list(logHR=NA, HR=NA, p=NA))
  
  return(data.frame(C_index=cindex,
                    logHR=cox_res$logHR,
                    HR=cox_res$HR,
                    Cox_p=cox_res$p))
}

#- - - - - -
## Fitting XGBoost with nested CV for exposure adjustment (one exposure at a time)
#
# Rationale:
# - this function residualizes an exposure on adjustment variables (e.g., cohort and/or covariates), following the double machine learning framework.
# - the objective follows the nature of the exposure being modeled (the dependent variable of this function).
#
# Arguments:
# - X: matrix of adjustment variables (cohort dummies and/or covariates)
# - y: numeric vector of the exposure to adjust
# - is_binary_exposure: TRUE if y is binary (0/1), FALSE if continuous/ordinal
# - nb_outer_fold, nb_inner_fold: number of folds for nested CV
# - param_grid: optional custom hyperparameter grid
#
# Returns: OOF predicted values for the full dataset, on the original exposure scale. Residual = exposure - OOF_pred (computed outside this function)

fit_adjustment_oof<-function(X, y, is_binary_exposure=FALSE, nb_outer_fold=5, nb_inner_fold=5, param_grid=NULL){
  
  y<-as.numeric(y)
  
  ## If constant exposure, residual = 0 for all observations
  if(max(y, na.rm=TRUE)==min(y, na.rm=TRUE)) return(rep(y[1], length(y)))
  
  if(is_binary_exposure){
    
    ## Binary exposure
    y_model<-as.integer(y)
    objective_<-"binary:logistic"
    eval_metric_<-"logloss"
    
    ## scale_pos_weight to account for class imbalance within the exposure itself
    tbl<-table(y_model)
    spw<-if(length(tbl)==2) as.numeric(max(tbl)/min(tbl)) else 1
    fold_type<-"stratified"
    
  }else{
    
    ## Continuous or ordinal exposure
    y_model<-y
    objective_<-"reg:squarederror"
    eval_metric_<-"rmse"
    spw<-1
    fold_type<-"basic"
  }
  
  ## Default parameter grid if none provided
  if(is.null(param_grid)){
    param_grid<-expand.grid(subsample=seq(0.5,1,0.25),
                            colsample_bytree=seq(0.5,1,0.25),
                            min_child_weight=c(1,5,10),
                            learning_rate=c(0.05,0.1,0.3),
                            max_depth=c(3,5,7)) %>%
      dplyr::slice_sample(n=30, replace=TRUE) # Limiting space to 30 combinations for computational efficiency and environmental sustainability considerations
  }
  
  ## Creating outer folds (stratified for binary exposures, basic otherwise)
  outer_folds<-splitTools::create_folds(y_model, k=nb_outer_fold, type=fold_type, seed=seed)
  oof_preds<-rep(NA, length(y_model))
  best_params_all<-list()
  
  for(outer_idx in seq_along(outer_folds)){
    
    val_idx<-outer_folds[[outer_idx]]
    train_idx<-setdiff(seq_len(length(y_model)), val_idx)
    
    X_tr<-X[train_idx,,drop=FALSE]
    y_tr<-y_model[train_idx]
    X_val<-X[val_idx,,drop=FALSE]
    y_val<-y_model[val_idx]
    
    ## Computing scale_pos_weight per outer fold for binary exposures
    spw_fold<-if(is_binary_exposure){
      tbl_tr<-table(y_tr)
      if(length(tbl_tr)==2) as.numeric(max(tbl_tr)/min(tbl_tr)) else 1
    }else{ 1 }
    
    ## Inner CV
    best_perf<-Inf # both logloss and rmse are minimized
    best_params<-NULL
    
    for(i in seq_len(nrow(param_grid))){
      
      params_i<-list(objective=objective_,eval_metric=eval_metric_,subsample=param_grid$subsample[i],colsample_bytree=param_grid$colsample_bytree[i],
                     min_child_weight=param_grid$min_child_weight[i],eta=param_grid$learning_rate[i],max_depth=param_grid$max_depth[i],
                     scale_pos_weight=spw_fold)
      
      dtrain_inner<-xgboost::xgb.DMatrix(data=X_tr, label=y_tr)
      
      ## stratified=TRUE only meaningful (and used) for binary exposures
      cv_res<-tryCatch(suppressWarnings(xgboost::xgb.cv(params=params_i, data=dtrain_inner,nrounds=100, nfold=nb_inner_fold,
                                                        early_stopping_rounds=10, verbose=0,stratified=is_binary_exposure)),error=function(e) NULL)
      
      if(is.null(cv_res)) next
      
      ## Extracting best performance across rounds (lower = better for both logloss and rmse)
      metric_col<-paste0("test_",eval_metric_,"_mean")
      best_val<-tryCatch(min(cv_res$evaluation_log[[metric_col]], na.rm=TRUE),error=function(e) NA)
      
      if(!is.na(best_val) && best_val<best_perf){
        best_perf<-best_val
        best_params<-params_i
      }
    }
    
    ## Fallback to default params if all tuning attempts failed
    if(is.null(best_params)){
      best_params<-list(objective=objective_,eval_metric=eval_metric_,subsample=1, colsample_bytree=1, min_child_weight=1,
                        eta=0.1, max_depth=3, scale_pos_weight=spw_fold)
    }
    best_params_all[[outer_idx]]<-best_params
    
    ## Training outer model and predicting on validation fold
    dtrain_outer<-xgboost::xgb.DMatrix(data=X_tr, label=y_tr)
    mod<-xgboost::xgb.train(params=best_params, data=dtrain_outer, nrounds=100, verbose=0)
    
    oof_preds[val_idx]<-predict(mod, xgboost::xgb.DMatrix(X_val))
  }
  
  return(oof_preds)
}

#- - - - - -
## Custom survival XGBoost learner with sample weight support
#
# Method: Inherits from mlsurvlrnrs::LearnerSurvXgboostCox and adds per-observation weights to account for event imbalance (few events vs. many censored observations).

LearnerSurvXgboostCoxWeighted<-R6::R6Class(classname="LearnerSurvXgboostCoxWeighted",inherit=mlsurvlrnrs::LearnerSurvXgboostCox,
  public=list(sample_weight=NULL,
    set_sample_weight=function(w){ self$sample_weight<-w },
    train=function(x, y, ...){
      if(!is.null(self$sample_weight)){
        super$train(x=x, y=y, sample_weight=self$sample_weight, ...)
      }else{
        super$train(x=x, y=y, ...)
      }
    }))

#- - - - - -
## Fitting XGBoost ERS model on doubly-adjusted exposures with nested CV
#
# Method : Performance is measured by Harrell's C-index (concordance index), with higher C-index = better discrimination of event timing.
#
## Arguments:
# - dataset_dt: data.table with time and event columns + residualized exposures
# - surv_time_col: name of the time-to-event column
# - surv_event_col: name of the event indicator column (1=event, 0=censored)
# - feature_cols: vector of doubly-adjusted exposure column names
# - train_split, test_split: dataset split proportions
# - nb_inner_fold, nb_outer_fold: number of folds for nested CV
# - param_grid: optional custom hyperparameter grid
#
# Outputs (list):
# - oof_preds: OOF predicted log partial hazard scores for the full dataset
# - test_preds: predicted scores on holdout test set
# - metric_test: C-index and Cox performance on test set
# - final_model: model retrained on full training set with best hyperparameters
# - outer_summary: CV C-index summary
# - final_params: best hyperparameters selected

fit_ers_survival<-function(dataset_dt,
                           surv_time_col,
                           surv_event_col,
                           feature_cols,
                           train_split=0.7,
                           test_split=0.3,
                           nb_inner_fold=5,
                           nb_outer_fold=5,
                           param_grid=NULL){
  
  ## Ensuring splits are in correct format
  train_split<-as.numeric(train_split)
  test_split<-as.numeric(test_split)
  
  if((train_split+test_split>1)|is.na(test_split)|is.na(train_split)){
    train_split<-0.7
    test_split<-0.3
    }
  
  ## Ensuring fold numbers are in correct format
  nb_inner_fold<-as.numeric(nb_inner_fold)
  if(is.na(nb_inner_fold)|nb_inner_fold<2) nb_inner_fold<-5
  
  nb_outer_fold<-as.numeric(nb_outer_fold)
  if(is.na(nb_outer_fold)|nb_outer_fold<2) nb_outer_fold<-5
  
  ## Extracting survival outcome
  surv_time<-as.numeric(dataset_dt[[surv_time_col]])
  surv_event<-as.integer(dataset_dt[[surv_event_col]])
  
  ## Creating stratified train/test split using kmeans multi-strata on (time, event) to ensure both splits have similar event rates and time distributions
  surv_cols_dt<-dataset_dt[, .SD, .SDcols=c(surv_time_col, surv_event_col)]
  split_vector<-splitTools::multi_strata(df=surv_cols_dt, strategy="kmeans", k=4)
  data_split<-splitTools::partition(y=split_vector, p=c(train=train_split, test=test_split),type="stratified", seed=seed)
  
  ## Verifying no overlap between train and test
  stopifnot(length(intersect(data_split$train, data_split$test))==0)
  
  ## Creating training and test datasets
  X_train<-as.matrix(dataset_dt[data_split$train, ..feature_cols])
  t_train<-surv_time[data_split$train]
  e_train<-surv_event[data_split$train]
  y_train<-Surv(time=t_train, event=e_train)
  
  X_test<-as.matrix(dataset_dt[data_split$test, ..feature_cols])
  t_test<-surv_time[data_split$test]
  e_test<-surv_event[data_split$test]
  y_test<-Surv(time=t_test, event=e_test)
  
  ## Computing sample weights to account for event imbalance (observations with rare events receive higher weight)
  event_dt<-data.table(event=e_train)
  event_counts<-event_dt[, .N, by=event][, weight:=1/N]
  train_weights<-event_counts[event_dt, on="event", weight]
  
  ## Creating stratified CV folds on training set (stratified on survival strata)
  surv_cols_train<-dataset_dt[data_split$train, .SD, .SDcols=c(surv_time_col, surv_event_col)]
  split_vector_train<-splitTools::multi_strata(df=surv_cols_train, strategy="kmeans", k=4)
  outer_folds<-splitTools::create_folds(y=split_vector_train, k=nb_outer_fold, type="stratified",seed=seq(1,nb_outer_fold,1)*seed)
  
  ## Default parameter grid if none provided
  if(is.null(param_grid)){
    param_grid<-expand.grid(learning_rate=c(0.05,0.1),
                            max_depth=c(3,4,5),
                            min_child_weight=c(5,10,20),
                            subsample=c(0.7,0.9,1.0),
                            colsample_bytree=c(0.7,0.9,1.0),
                            nrounds=c(50,100,200)) %>%
      dplyr::slice_sample(n=30, replace=FALSE) # Limiting space to 30 combinations for computational efficiency and environmental sustainability considerations
  }
  
  ## Defining learner arguments
  learner_args<-list(objective="survival:cox", eval_metric="cox-nloglik", nthread=ncores)
  
  ## Outer CV
  outer_results<-list()
  best_params_all<-list()
  
  for(outer_idx in seq_along(outer_folds)){
    
    val_idx<-outer_folds[[outer_idx]]
    train_idx_cv<-setdiff(seq_len(nrow(X_train)), val_idx)
    
    X_tr<-X_train[train_idx_cv,,drop=FALSE]
    y_tr<-y_train[train_idx_cv]
    t_tr<-t_train[train_idx_cv]
    e_tr<-e_train[train_idx_cv]
    w_tr<-train_weights[train_idx_cv]
    
    X_val<-X_train[val_idx,,drop=FALSE]
    y_val<-y_train[val_idx]
    t_val<-t_train[val_idx]
    e_val<-e_train[val_idx]
    
    ## Skipping outer fold if it contains no events (C-index undefined)
    if(sum(e_val)==0){
      cat("Outer fold",outer_idx,": skipped (no events in validation set)\n")
      best_params_all[[outer_idx]]<-NULL
      next
    }
    
    ## Inner CV folds
    surv_cols_tr<-data.table(time=t_tr, event=e_tr)
    split_vec_tr<-splitTools::multi_strata(df=surv_cols_tr, strategy="kmeans", k=4)
    inner_folds<-splitTools::create_folds(y=split_vec_tr, k=nb_inner_fold, type="stratified",seed=seq(1,nb_inner_fold,1)*seed)
    
    ## Hyperparameter tuning using weighted survival Cox learner
    learner_inner<-LearnerSurvXgboostCoxWeighted$new(metric_optimization_higher_better=FALSE)
    learner_inner$set_sample_weight(w_tr)
    
    best_nloglik<-Inf
    best_params<-NULL
    
    for(i in seq_len(nrow(param_grid))){
      xgb_cv<-mlexperiments::MLCrossValidation$new(learner=learner_inner,fold_list=inner_folds,ncores=ncores,seed=seed)
      
      xgb_cv$learner_args<-c(as.list(param_grid[i,]), learner_args)
      xgb_cv$performance_metric<-c_index
      xgb_cv$set_data(x=X_tr, y=y_tr)
      
      res_cv<-tryCatch(xgb_cv$execute(), error=function(e) NULL)
      if(is.null(res_cv)) next
      
      ## Extracting mean cox-nloglik (lower = better) across inner folds
      mean_nloglik<-tryCatch(
        mean(sapply(res_cv$results$folds, function(f){
          val<-tryCatch(f$performance, error=function(e) NA)
          if(is.list(val)) val<-unlist(val)
          as.numeric(val[1])
        }), na.rm=TRUE),
        error=function(e) NA)
      
      if(!is.na(mean_nloglik)&&mean_nloglik<best_nloglik){
        best_nloglik<-mean_nloglik
        best_params<-as.list(param_grid[i,])
      }
    }
    
    ## Fallback to default params if tuning failed
    if(is.null(best_params)) best_params<-list(learning_rate=0.1,
                                               max_depth=3,
                                               min_child_weight=10,
                                               subsample=0.7,
                                               colsample_bytree=0.7,
                                               nrounds=100)
    best_params_all[[outer_idx]]<-best_params
    
    ## Training outer model and evaluating C-index on outer validation fold
    learner_outer<-LearnerSurvXgboostCoxWeighted$new(metric_optimization_higher_better=FALSE)
    learner_outer$set_sample_weight(w_tr)
    
    mod_outer<-learner_outer$fit(x=X_tr, y=y_tr, seed=seed, ncores=ncores,
                                 objective=learner_args$objective,
                                 eval_metric=learner_args$eval_metric,
                                 subsample=best_params$subsample,
                                 colsample_bytree=best_params$colsample_bytree,
                                 min_child_weight=best_params$min_child_weight,
                                 learning_rate=best_params$learning_rate,
                                 nrounds=best_params$nrounds,
                                 max_depth=best_params$max_depth)
    
    val_pred<-predict(mod_outer, xgboost::xgb.DMatrix(X_val))
    
    ## C-index on outer validation fold
    cindex_val<-tryCatch(round(survival::concordance(Surv(t_val,e_val)~val_pred)$concordance, 3),error=function(e) NA)
    
    outer_results[[outer_idx]]<-data.frame(Fold=outer_idx, C_index=cindex_val)
  }
  
  outer_summary<-dplyr::bind_rows(outer_results)
  cat("\nNested CV summary (C-index per outer fold):\n")
  print(outer_summary)
  cat("\nOuter mean C-index:", round(mean(outer_summary$C_index, na.rm=TRUE), 3),
      "SD:", round(sd(outer_summary$C_index, na.rm=TRUE), 3), "\n")
  
  ## Selecting best hyperparameters based on highest outer CV C-index
  best_idx<-which.max(outer_summary$C_index)
  final_params<-best_params_all[[best_idx]]
  
  ## Retraining final model on full training set
  learner_final<-LearnerSurvXgboostCoxWeighted$new(metric_optimization_higher_better=FALSE)
  learner_final$set_sample_weight(train_weights)
  
  final_model<-learner_final$fit(x=X_train, y=y_train, seed=seed, ncores=ncores,
                                 objective=learner_args$objective,
                                 eval_metric=learner_args$eval_metric,
                                 subsample=final_params$subsample,
                                 colsample_bytree=final_params$colsample_bytree,
                                 min_child_weight=final_params$min_child_weight,
                                 learning_rate=final_params$learning_rate,
                                 nrounds=final_params$nrounds,
                                 max_depth=final_params$max_depth)
  
  ## Computing OOF predicted log partial hazard scores on the full dataset
  # - train: predicted by per-fold models (each observation not seen during training)
  # - test: predicted by final model (no leakage)
  oof_preds_full<-rep(NA, nrow(dataset_dt))
  
  for(outer_idx in seq_along(outer_folds)){
    if(is.null(best_params_all[[outer_idx]])) next
    
    val_idx_global<-data_split$train[outer_folds[[outer_idx]]]
    train_idx_cv<-setdiff(seq_len(length(data_split$train)), outer_folds[[outer_idx]])
    X_tr_oof<-X_train[train_idx_cv,,drop=FALSE]
    y_tr_oof<-y_train[train_idx_cv]
    w_tr_oof<-train_weights[train_idx_cv]
    
    learner_oof<-LearnerSurvXgboostCoxWeighted$new(metric_optimization_higher_better=FALSE)
    learner_oof$set_sample_weight(w_tr_oof)
    bp<-best_params_all[[outer_idx]]
    
    mod_oof<-learner_oof$fit(x=X_tr_oof, y=y_tr_oof, seed=seed, ncores=ncores,
                             objective=learner_args$objective,
                             eval_metric=learner_args$eval_metric,
                             subsample=bp$subsample,
                             colsample_bytree=bp$colsample_bytree,
                             min_child_weight=bp$min_child_weight,
                             learning_rate=bp$learning_rate,
                             nrounds=bp$nrounds,
                             max_depth=bp$max_depth)
    
    oof_preds_full[val_idx_global]<-predict(mod_oof,xgboost::xgb.DMatrix(X_train[outer_folds[[outer_idx]],,drop=FALSE]))
  }
  
  ## For test observations: using final model
  test_preds<-predict(final_model, xgboost::xgb.DMatrix(X_test))
  oof_preds_full[data_split$test]<-test_preds
  
  ## Test set performance
  cindex_test<-round(survival::concordance(Surv(t_test,e_test)~test_preds)$concordance, 3)
  cox_test<-survival::coxph(Surv(t_test,e_test)~test_preds,data=data.frame(t_test=t_test, e_test=e_test, test_preds=test_preds))
  cox_test_s<-summary(cox_test)
  HR_test<-round(exp(coef(cox_test)), 3)
  HR_CI_test<-round(exp(confint(cox_test)), 3)
  HR_str<-paste(HR_test," (",HR_CI_test[1],"; ",HR_CI_test[2],")", sep="")
  p_test<-round(cox_test_s$logtest["pvalue"], 3)
  
  cat("\nTest set performance:\n")
  cat("C-index:", cindex_test, "\n")
  cat("HR (ERS):", HR_str, " LR p-value:", p_test, "\n")
  
  metric_test<-as_tibble(data.frame(C_index=cindex_test,
                                    HR_ERS=HR_test,
                                    HR_CI=HR_str,
                                    Cox_LR_pvalue=p_test))
  
  return(list(oof_preds=oof_preds_full,
              test_preds=test_preds,
              final_model=final_model,
              outer_summary=outer_summary,
              final_params=final_params,
              data_split=data_split,
              metric_test=metric_test,
              X_test=X_test,
              y_test=y_test,
              t_test=t_test,
              e_test=e_test))
}

#----------------------------------------------------------------
#### Settings ####

## Set to TRUE if data comes from multiple cohorts (activates Step 1)
use_cohort<-TRUE

## Names of survival outcome columns
surv_time_col<-"time" # name of the time-to-event column
surv_event_col<-"event" # name of the event indicator column (1=event, 0=censored)

## Names of covariates to adjust for in Step 2
cov_names<-c("age","sex","bmi")

#----------------------------------------------------------------
#### Data import and preprocessing (replace simulation with your own data) ####

## Simulating survival data for demonstration
set.seed(seed)
n<-800
K<-8

## Simulating continuous exposures
expo_df<-as.data.frame(exp(MASS::mvrnorm(n,rep(0,K),0.5^as.matrix(dist(1:K)))/3))

## Simulating a binary exposure
expo_df2<-data.frame(X9=sample(x=c(0,1),n,replace=TRUE))
expo_df<-bind_cols(expo_df,expo_df2)
colnames(expo_df)<-paste0("X",1:(K+1))
expo_names<-colnames(expo_df)

## Building exposure matrix: log1p for continuous, untransformed for binary
expo_mat<-as.matrix(cbind(log1p(expo_df[,paste0("X",1:K)]), expo_df[,"X9",drop=FALSE]))
colnames(expo_mat)<-expo_names

## Covariates and cohort
cov_df<-data.frame(age=rnorm(n,50,10), sex=rbinom(n,1,0.5), bmi=rnorm(n,26,4))
cohort<-sample(c("cohort_A","cohort_B","cohort_C"),n,replace=TRUE)

## Cohort effect on the log hazard scale
cohort_eff<-ifelse(cohort=="cohort_A",0,ifelse(cohort=="cohort_B",0.3,-0.3))

## True exposure effect (what ERS should capture after removing cohort + covariate effects)
h_z<-as.numeric(as.matrix(expo_mat)%*%c(0.5,0.3,-0.2,0.1,0.4,-0.1,0.05,0.05,1))

## Linear predictor on the log hazard scale, WITHOUT intercept
log_hazard_no_intercept<-0.02*cov_df$age + 0.3*cov_df$sex + 0.05*cov_df$bmi + cohort_eff + h_z

## Calibrating intercept so that the median event time (under Weibull) is centered around a target value (e.g., 10).
# - under Weibull(shape, scale), median = scale * (log(2))^(1/shape)
# - scale = exp(-log_hazard/shape) -> intercept shifts log_hazard so that median(scale) matches the target.
shape<-1.5
target_median_time<-10
target_scale<-target_median_time/(log(2)^(1/shape))
intercept_calibrated<- -shape*log(target_scale) - median(log_hazard_no_intercept)

log_hazard<-intercept_calibrated + log_hazard_no_intercept

## Simulating survival times from a Weibull model
scale<-exp(-log_hazard/shape)
surv_time_raw<-rweibull(n, shape=shape, scale=scale)

cat("Event time distribution (before censoring):\n")
cat("Median:", round(median(surv_time_raw),3),
    "Mean:", round(mean(surv_time_raw),3),
    "Max:", round(max(surv_time_raw),3), "\n")

## Simulating censoring times
target_censoring<-0.30

find_censoring_rate<-function(rate, surv_time_raw, seed_local){
  set.seed(seed_local)
  cens_time_tmp<-rexp(length(surv_time_raw), rate=rate)
  mean(surv_time_raw > cens_time_tmp)
}

rate_grid<-exp(seq(log(1/(target_median_time*50)), log(1/(target_median_time*0.01)), length.out=300))
cens_props<-sapply(rate_grid, find_censoring_rate, surv_time_raw=surv_time_raw, seed_local=seed)

best_rate<-rate_grid[which.min(abs(cens_props - target_censoring))]

set.seed(seed)
cens_time<-rexp(n, rate=best_rate)

## Observed time and event indicator
obs_time <-pmin(surv_time_raw, cens_time)
obs_event<-as.integer(surv_time_raw <= cens_time)

cat("Event rate:", round(mean(obs_event),3),
    "- Events:", sum(obs_event), "- Censored:", sum(obs_event==0), "\n")
cat("Median follow-up:", round(median(obs_time),3), "\n")

## Assembling full dataset
sim_data<-cbind(data.frame(time=obs_time, event=obs_event, cohort=cohort), cov_df, expo_mat)
dataset<-as.data.table(sim_data)

## Names of binary exposures in your dataset (all others treated as continuous/ordinal)
expo_binary<-c("X9") # adapt to your data

#----------------------------------------------------------------
#### Step 1 (optional): adjusting exposures for cohort membership (only if multi-cohort) ####

# Goal: to remove between-cohort differences from each exposure so that the ERS is not driven by which cohort a participant belongs to.
#       Residualization is applied to exposures to preserve the time-to-event structure throughout. DML framework: W_j = X_j - E[X_j | cohort].

# Method:
# - for each exposure X_j, fit XGBoost from cohort.
# - W_j = X_j - predicted(cohort): exposure cleaned of cohort effects.
# - OOF predictions used for all observations to avoid overfitting residuals.

# Output:
# - expo_mat_W: matrix of exposure residuals after cohort adjustment
# - Survival metrics (C-index, HR) before and after adjustment per exposure

if(use_cohort){
  
  ## One-hot encoding of cohort membership (no intercept to avoid collinearity)
  Z_mat<-model.matrix(~as.factor(sim_data$cohort)-1)
  colnames(Z_mat)<-paste0("cohort_",levels(as.factor(sim_data$cohort)))
  
  ## Adjusting each exposure for cohort using lapply
  cat("\nStep 1 - adjusting all exposures for cohort...\n")

  expo_hat_cohort_list<-lapply(expo_names, function(j){
    is_bin<-j %in% expo_binary
    cat("Adjusting:",j, if(is_bin) "[binary]" else "[continuous]","\n")
    fit_adjustment_oof(X=Z_mat, y=expo_mat[,j],is_binary_exposure=is_bin,nb_outer_fold=5, nb_inner_fold=5)
  })
  
  expo_hat_cohort<-do.call(cbind, expo_hat_cohort_list)
  colnames(expo_hat_cohort)<-expo_names
  
  ## W = X - E[X|cohort]: exposures cleaned of cohort effects
  expo_mat_W<-expo_mat-expo_hat_cohort
  
  ## Survival metrics (exposure ~ survival outcome) before and after cohort adjustment
  metrics_before_cohort<-do.call(rbind, lapply(expo_names, function(j)
    cbind(data.frame(exposure=j, adjustment="before"),compute_survival_metrics(expo_mat[,j], obs_time, obs_event))))
  
  metrics_after_cohort<-do.call(rbind, lapply(expo_names, function(j)
    cbind(data.frame(exposure=j, adjustment="after"),compute_survival_metrics(expo_mat_W[,j], obs_time, obs_event))))
  
  metrics_cohort<-bind_rows(metrics_before_cohort, metrics_after_cohort)
  
  cat("\nSurvival metrics (exposure ~ outcome) before and after cohort adjustment:\n")
  metrics_cohort %>% pivot_wider(names_from=adjustment,values_from=c(C_index,logHR,HR,Cox_p))
  
}else{
  expo_mat_W<-expo_mat
}

#----------------------------------------------------------------
#### Step 2: adjusting exposures for covariates ####

# Goal: to remove covariate effects from each exposure residual so that the ERS captures only the exposure-specific association with survival, independent of covariates (e.g., age, sex, BMI)
#       Second DML stage: W2_j = W_j - E[W_j | covariates].

cov_mat<-as.matrix(sim_data[,cov_names])

## Adjusting each exposure residual for covariates
cat("\nStep 2 - adjusting all exposures for covariates\n")

expo_hat_cov_list<-lapply(expo_names,function(j){
  is_bin<-j %in% expo_binary
  cat("Adjusting:",j, if(is_bin) "[binary]" else "[continuous]","\n")
  fit_adjustment_oof(X=cov_mat,y=expo_mat_W[,j],is_binary_exposure=is_bin,nb_outer_fold=5,nb_inner_fold=5)
})

expo_hat_cov<-do.call(cbind, expo_hat_cov_list)
colnames(expo_hat_cov)<-expo_names

## W2 = W - E[W|covariates]: exposures cleaned of both cohort and covariate effects
expo_mat_W2<-expo_mat_W-expo_hat_cov

## Survival metrics (exposure ~ survival outcome) before and after covariate adjustment
metrics_before_cov<-do.call(rbind, lapply(expo_names, function(j)
  cbind(data.frame(exposure=j, adjustment="before"),compute_survival_metrics(expo_mat_W[,j], obs_time, obs_event))))

metrics_after_cov<-do.call(rbind, lapply(expo_names, function(j)
  cbind(data.frame(exposure=j, adjustment="after"),compute_survival_metrics(expo_mat_W2[,j], obs_time, obs_event))))

metrics_cov<-bind_rows(metrics_before_cov, metrics_after_cov)

cat("\nSurvival metrics (exposure ~ outcome) before and after covariate adjustment:\n")
metrics_cov %>% pivot_wider(names_from=adjustment,values_from=c(C_index,logHR,HR,Cox_p))

#----------------------------------------------------------------
#### Step 3: Fitting ERS model on doubly-adjusted exposures (nested CV) ####

## Goal: to model the survival outcome using doubly-adjusted exposures W2 only.
#        The predicted log partial hazard f(W2) represents the exposure-driven component of the survival outcome (= ERS), cleanly separated from cohort and covariate effects

colnames(expo_mat_W2)<-expo_names
dataset_step3<-as.data.table(cbind(data.frame(time=obs_time, event=obs_event),as.data.frame(expo_mat_W2)))

step3_result<-fit_ers_survival(dataset_dt=dataset_step3,
                               surv_time_col=surv_time_col,
                               surv_event_col=surv_event_col,
                               feature_cols=expo_names,
                               train_split=0.7,
                               test_split=0.3,
                               nb_outer_fold=5,
                               nb_inner_fold=5)

## ERS = predicted log partial hazard from doubly-adjusted exposures on the full dataset

# Interpretation: higher ERS = higher predicted hazard -> higher risk
ers<-predict(step3_result$final_model, xgboost::xgb.DMatrix(expo_mat_W2))
sim_data$ers<-ers

## Evaluating ERS against the original survival outcome (full dataset)
ers_surv_metrics<-compute_survival_metrics(ers, obs_time, obs_event)
cindex_ers<-ers_surv_metrics$C_index

cat("\nERS performance against survival outcome (full dataset):\n")
cat("C-index:", cindex_ers, "\n")
cat("HR:", ers_surv_metrics$HR, " logHR:", ers_surv_metrics$logHR,
    "Cox p:", ers_surv_metrics$Cox_p, "\n")

## Saving ERS model
saveRDS(step3_result$final_model, paste0("ERS_model_survival_",Sys.Date(),".rds"))
xgb.save(step3_result$final_model,paste0("ERS_model_survival_",Sys.Date(),".ubj"))

#----------------------------------------------------------------
#### Step 4: Overall ERS direction (mixture effect) ####

## Goal: to determine whether the overall exposome (ERS) promotes or mitigates risk, using a same 3-method consensus framework
#        Here the outcome is the log hazard (ERS itself is on the log hazard scale), so we use the event indicator as the reference signal for directionality.

## Building a two-column data frame: ERS value vs. event indicator
ers_direction_df<-data.frame(ers=ers, y=as.numeric(obs_event))

## Approach 1: sign of mean association
ers_mean_sign<-conventional_direction(ers_direction_df, feature_col="ers", shap_col="y")

## Approach 2: GAM derivative + Spearman fallback
ers_gam<-gam_direction(ers_direction_df, feature_col="ers", shap_col="y")

## Approach 3: pairwise bin concordance
ers_pairwise<-pairwise_direction(ers_direction_df, feature_col="ers", shap_col="y",majority_threshold=0.55, n_quantile_bins=200)

## Consensus: majority vote across the three approaches
ers_consensus_vec<-c(ers_mean_sign, ers_gam, ers_pairwise)
ers_tbl<-sort(table(ers_consensus_vec), decreasing=TRUE)
ers_consensus<-if(ers_tbl[1]>=2) names(ers_tbl)[1] else "uncertain"

## Summary table: direction + survival metrics (no R² or binary metrics)
ers_direction_summary<-tibble(mean_sign_direction=ers_mean_sign,
                              gam_direction=ers_gam,
                              pairwise_direction=ers_pairwise,
                              consensus_direction=ers_consensus,
                              C_index=cindex_ers,
                              HR=ers_surv_metrics$HR,
                              logHR=ers_surv_metrics$logHR,
                              Cox_p=ers_surv_metrics$Cox_p)

cat("\nERS global direction and performance summary:\n")
ers_direction_summary

#----------------------------------------------------------------
#### Step 5: Using ERS in a final association model ####

# The ERS can be used as a predictor in standard survival models alongside covariates.
# Since W2 (doubly-adjusted exposures) were used to build the ERS, the ERS is by construction largely orthogonal to cohort and covariates.

## Cox model: ERS + covariates
final_cox_glm<-survival::coxph(Surv(time,event)~ers+age+sex+bmi, data=sim_data)
summary(final_cox_glm)
broom::tidy(final_cox_glm, exponentiate=TRUE, conf.int=TRUE)

## Stratification into tertiles
sim_data$ers_group<-ntile(sim_data$ers, 3)
sim_data$ers_group<-factor(sim_data$ers_group, labels=c("Low","Medium","High"))

ers_tertile_cox<-survival::coxph(Surv(time,event)~ers_group+age+sex+bmi, data=sim_data)
broom::tidy(ers_tertile_cox, exponentiate=TRUE, conf.int=TRUE)

#----------------------------------------------------------------
#### Step 6: Visualization ####

#- - - - -
## Figure 1: ERS distribution by event status
ggplot(sim_data, aes(x=ers, fill=factor(event)))+
  geom_histogram(position="identity", alpha=0.6, bins=40, color="black")+
  scale_fill_manual("Event", values=c("0"="#A6DDCE","1"="#F9CBC2"),
                    labels=c("Censored","Event"))+
  labs(title="Distribution of ERS by event status", x="ERS (log partial hazard scale)", y="Count")+
  theme_Gaia()

#- - - - -
## Figure 2: ERS by cohort
if(use_cohort){
  ggplot(sim_data, aes(x=cohort, y=ers, fill=cohort))+
    geom_violin(fill="transparent")+
    geom_boxplot(alpha=0.7, width=0.4)+
    scale_fill_brewer(palette="Set2")+
    labs(title="ERS by cohort",
         x="", y="ERS")+
    theme_Gaia()+
    theme(legend.position="none")
}

#- - - - -
## Figure 3: ERS by event status within each cohort
if(use_cohort){
  ggplot(sim_data, aes(x=factor(event), y=ers, fill=factor(event)))+
    geom_boxplot(alpha=0.7)+
    scale_fill_manual(values=c("0"="#A6DDCE","1"="#F9CBC2"))+
    facet_wrap(~cohort)+
    labs(title="ERS by event status within each cohort",
         x="Event (0=censored, 1=event)", y="ERS")+
    theme_Gaia()+
    theme(legend.position="none")
}else{
  ggplot(sim_data, aes(x=factor(event), y=ers, fill=factor(event)))+
    geom_boxplot(alpha=0.7)+
    scale_fill_manual(values=c("0"="#A6DDCE","1"="#F9CBC2"))+
    labs(title="ERS by event status",
         x="Event (0=censored, 1=event)", y="ERS")+
    theme_Gaia()+
    theme(legend.position="none")
}

#- - - - -
## Figure 4: Baseline survival (overall) for ERS

ers_cox_mod<-survival::coxph(Surv(time,event)~ers, data=sim_data)

# version 1
ggsurvplot(survfit(ers_cox_mod, data = sim_data),
           palette="#40B696",
           risk.table=TRUE,
           pval=TRUE,
           conf.int=TRUE,
           ggtheme=theme_Gaia(),
           xlab="Time",
           ylab="Overall survival probability")

# version 2
survfit2(ers_cox_mod) %>%
  ggsurvfit(color = "#40B696") +
  labs(x = "Time", y = "Overall survival probability") +
  add_confidence_interval(fill = "#40B696") +
  add_risktable()+
  theme_Gaia()           

# version 3
ggplot(sim_data, aes(time=time, status=event)) + geom_km() + geom_kmband() + theme_bw()+
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Survival probability")+
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="top",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

#- - - - -
## Figure 5: Cumulative events (overall) for ERS

# version 1
ggsurvplot(survfit(ers_cox_mod, data = sim_data),
           palette= "#40B696",
           risk.table= TRUE,
           pval= TRUE,
           conf.int= TRUE,
           ggtheme= theme_Gaia(),
           xlab="Time",
           ylab="Cumulative event probability",
           fun="event")

# version 2
survfit2(ers_cox_mod) %>%
  ggsurvfit(color="#40B696",type="risk") +
  labs(x = "Time", y = "Cumulative event probability") +
  add_confidence_interval(fill = "#40B696") +
  add_risktable()+
  theme_Gaia()

#- - - - -
## Figure 6: Cumulative hazards (overall) for ERS

# version 1
ggsurvplot(survfit(ers_cox_mod, data = sim_data),
           palette= "#40B696",
           risk.table= TRUE,
           pval= TRUE,
           conf.int= TRUE,
           ggtheme= theme_Gaia(),
           xlab="Time",
           ylab="Cumulative hazard",
           fun="cumhaz")

# version 2
survfit2(ers_cox_mod) %>%
  ggsurvfit(color="#40B696",type="cumhaz") +
  labs(x = "Time", y = "Cumulative hazard") +
  add_confidence_interval(fill = "#40B696") +
  add_risktable()+
  theme_Gaia()

#- - - - -
## Figure 7: Kaplan-Meier curves by ERS tertile group

# version 1
km_fit<-survfit(Surv(time,event)~ers_group, data=sim_data)

survminer::ggsurvplot(km_fit,
                      data=sim_data,
                      risk.table=TRUE,
                      pval=TRUE,
                      conf.int=TRUE,
                      palette=c("#40B696","#F1CB0E","red"),
                      legend.title="ERS Group",
                      legend.labs=c("Low","Medium","High"),
                      xlab="Time",
                      ylab="Survival probability",
                      title=paste0("KM curves by ERS tertile (C-index=",cindex_ers,")"))

# version 2
ggplot(sim_data, aes(time=time, status=event, fill=ers_group, color=ers_group)) + geom_km() + geom_kmband() + theme_bw()+
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Survival probability")+
  scale_fill_manual("",values=c(Low="#40B696", Medium="#F1CB0E",High="red"))+
  scale_color_manual("",values=c(Low="#40B696", Medium="#F1CB0E",High="red"))+
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="top",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

# version 3
ggplot(sim_data, aes(time=time, status=event, fill=ers_group, color=ers_group)) + geom_km() + theme_bw()+
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Survival probability")+
  scale_fill_manual("",values=c(Low="#40B696", Medium="#F1CB0E",High="red"))+
  scale_color_manual("",values=c(Low="#40B696", Medium="#F1CB0E",High="red"))+
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="top",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

#- - - - -
## Figure 8: Cumulative events by ERS tertile group

ers_group_cox_mod<-survival::coxph(Surv(time,event)~ers_group, data=sim_data)

# version 1
ggsurvplot(km_fit,
           ggtheme=theme_Gaia(),
           data=sim_data,
           risk.table=TRUE,
           pval=TRUE,
           conf.int=TRUE,
           palette=c("#40B696","#F1CB0E","red"),
           legend.title="ERS Group",
           legend.labs=c("Low","Medium","High"),
           xlab="Time",
           ylab="Cumulative event probability",
           fun="event")

# version 2
survfit2(Surv(time, event) ~ ers_group,data=sim_data) %>%
  ggsurvfit(aes(color=ers_group,fill=ers_group),type="risk") +
  labs(x = "Time", y = "Cumulative event probability") +
  add_confidence_interval() +
  add_risktable()+
  scale_fill_manual("",values=c(Low="#40B696", Medium="#F1CB0E", High="red"))+
  scale_color_manual("",values=c(Low="#40B696", Medium="#F1CB0E", High="red"))+
  theme_Gaia()

#- - - - -
## Figure 9: Cumulative hazard by ERS tertile group

# version 1
ggsurvplot(km_fit,
           ggtheme=theme_Gaia(),
           data=sim_data,
           risk.table=TRUE,
           pval=TRUE,
           conf.int=TRUE,
           palette=c("#40B696","#F1CB0E","red"),
           legend.title="ERS Group",
           legend.labs=c("Low","Medium","High"),
           xlab="Time",
           ylab="Cumulative hazard",
           fun="cumhaz")

# version 2
survfit2(Surv(time, event) ~ ers_group,data=sim_data) %>%
  ggsurvfit(aes(color=ers_group,fill=ers_group),type="cumhaz") +
  labs(x = "Time", y = "Cumulative hazard") +
  add_confidence_interval() +
  add_risktable()+
  scale_fill_manual("",values=c(Low="#40B696", Medium="#F1CB0E", High="red"))+
  scale_color_manual("",values=c(Low="#40B696", Medium="#F1CB0E", High="red"))+
  theme_Gaia()

#----------------------------------------------------------------
#### Step 7: Univariate Cox regression on doubly-adjusted exposures (test set) ####

# Interpretation: HR coefficients reflect exposure-specific associations with survival, cleaned of cohort and covariate effects via DML.

Cox_tab_res<-c()
X_test_ers<-step3_result$X_test
t_test_ers<-step3_result$t_test
e_test_ers<-step3_result$e_test

for(n_var in seq_along(expo_names)){
  
  test_tmpo<-as_tibble(X_test_ers) %>% select(all_of(expo_names[n_var]))
  colnames(test_tmpo)<-"tmp"
  test_tmpo<-bind_cols(data.frame(time=t_test_ers, event=e_test_ers), test_tmpo)
  
  cox_model<-tryCatch(survival::coxph(Surv(time,event)~tmp, data=test_tmpo), error=function(e) NULL)
  if(is.null(cox_model)) next
  
  mod_res<-broom::tidy(cox_model, exponentiate=TRUE, conf.int=TRUE) %>%
    mutate(term=expo_names[n_var],
           direction=case_when(conf.low>1~"positive association",
                               conf.high<1~"negative association",
                               TRUE~"uncertain association")) %>%
    bind_cols(broom::glance(cox_model)) %>%
    rowwise() %>% mutate_if(is.numeric, test_format) %>% ungroup() %>%
    mutate(HR=paste(estimate," [",conf.low,"; ",conf.high,"]",sep="")) %>%
    mutate(HR=if_else(HR==" [; ]","",HR))
  
  Cox_tab_res<-bind_rows(Cox_tab_res, mod_res)
}

#----------------------------------------------------------------
#### Step 8: SHAP analysis ####

## SHAP values are computed on the test set of Step 3. They show which doubly-adjusted exposures contribute most to the ERS, and in which direction:
# - promoting = increases predicted log partial hazard = 'risk-increasing' (shorter survival)
# - mitigating = decreases predicted log partial hazard = 'protective' (longer survival)

## Computing SHAP values
contr<-predict(step3_result$final_model, as.matrix(X_test_ers), predcontrib=TRUE)
shap<-as_tibble(contr)
shap_contrib<-as.data.table(contr)

## Removing BIAS term
shap_contrib<-shap_contrib[,!grepl("bias|Intercept|BIAS|Bias",names(shap_contrib),ignore.case=TRUE),with=FALSE]

## Computing mean absolute SHAP score per feature (= importance ranking)
mean_shap_score<-colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)),decreasing=TRUE)]

## Reshaping SHAP values to long format
shap_score_sub<-as.data.table(shap_contrib)[,names(mean_shap_score),with=FALSE]
shap_score_long<-melt.data.table(shap_score_sub, measure.vars=colnames(shap_score_sub))

## Matching standardized feature values for color scale
fv_sub<-as.data.table(X_test_ers)[,names(mean_shap_score),with=FALSE]
fv_sub_long<-melt.data.table(fv_sub, measure.vars=colnames(fv_sub))
fv_sub_long[,stdfvalue:=std1(value),by="variable"]
names(fv_sub_long)<-c("variable","rfvalue","stdfvalue")

## Merging SHAP and feature values
shap_long2<-cbind(shap_score_long, fv_sub_long[,c("rfvalue","stdfvalue")])
shap_long2[,mean_value:=mean(abs(value)),by=variable]
setkey(shap_long2,variable)

## Computing summary statistics per feature (mean |SHAP| with 95% CI)
Shap_val<-as_tibble(shap_long2) %>%
  group_by(variable) %>%
  mutate(IC_2.5=quantile(abs(value),0.025),
         IC_97.5=quantile(abs(value),0.975),
         mean_val=mean(abs(value),na.rm=TRUE)) %>%
  select(variable,mean_val,IC_2.5,IC_97.5) %>%
  distinct() %>%
  arrange(desc(mean_val)) %>%
  mutate(mean_SHAP=paste(signif(mean_val,3)," (",signif(IC_2.5,3),"; ",signif(IC_97.5,3),")",sep="")) %>%
  ungroup()
colnames(Shap_val)[1]<-"feature"

## Formatting SHAP values for display
Shap_val<-Shap_val %>%
  rowwise() %>%
  mutate(mean_val=test_format(mean_val), IC_2.5=test_format(IC_2.5), IC_97.5=test_format(IC_97.5)) %>%
  ungroup() %>%
  mutate(mean_SHAP=paste(mean_val," (",IC_2.5,"; ",IC_97.5,")",sep=""),
         mean_SHAP=if_else(mean_SHAP%in%c("0 (0; 0)","0e+00 (0e+00; 0e+00)"),"0",mean_SHAP))

## Pivoting to long format for direction analysis
shap_tempo<-shap %>%
  rowid_to_column("id") %>%
  pivot_longer(-id, names_to="feature", values_to="shap_value")

tempopo<-X_test_ers %>%
  as_tibble() %>%
  rowid_to_column("id") %>%
  pivot_longer(-id, names_to="feature", values_to="feature_value")

tempopo<-left_join(tempopo, shap_tempo, by=c("id","feature"))
tempopo<-left_join(tempopo, Shap_val, by=base::intersect(colnames(tempopo),colnames(Shap_val)))

## Computing SHAP directionality using 3-method consensus
direction_impact<-compute_shap_directions_long(tempopo, threshold=0.55, n_bins=200) %>%
  select(-c(n,mean_shap,median_shap)) %>%
  filter(feature %ni% c("(Intercept)","BIAS","Bias","bias")) %>%
  mutate(consensus=apply(cbind(conventional,gam,pairwise),1,function(x){
    tbl<-sort(table(x),decreasing=TRUE)
    if(tbl[1]>=2) names(tbl)[1] else "uncertain"
  }))

## Merging direction with SHAP long format
test<-left_join(tempopo, direction_impact, by="feature") %>% mutate(direction=consensus)
test<-left_join(test, Shap_val, by=base::intersect(colnames(test),colnames(Shap_val)))

## Preparing data for plots: applying sign to mean |SHAP| based on direction
test_tmp<-test %>%
  mutate(impact=direction) %>%
  select(feature,mean_val,IC_2.5,IC_97.5,impact,mean_SHAP) %>%
  distinct() %>%
  mutate(across(c(mean_val,IC_2.5,IC_97.5),as.numeric)) %>%
  mutate(mean_val=if_else(impact=="mitigating",-mean_val,mean_val),
         IC_2.5  =if_else(impact=="mitigating",-IC_2.5,  IC_2.5),
         IC_97.5 =if_else(impact=="mitigating",-IC_97.5, IC_97.5)) %>%
  arrange(desc(mean_val)) %>%
  filter(mean_val!=0,!is.na(mean_val))

## Ordering features for plots
ordre_def<-bind_rows(
  test_tmp %>% filter(mean_val>=0) %>% arrange(desc(IC_2.5)),
  test_tmp %>% filter(mean_val< 0) %>% arrange(desc(IC_97.5)))

test_tmp$feature<-factor(test_tmp$feature, levels=rev(unique(ordre_def$feature)))
test_tmp$impact<-factor(test_tmp$impact,
                        levels=c("mitigating","neutral","promoting","uncertain","undefined"),
                        labels=c("mitigating predictor","neutral","promoting predictor","uncertain","undefined"))

## Plot 1: error bar (mean |SHAP| with 95% CI, colored by direction)
plot1<-ggplot(test_tmp,aes(y=feature,x=mean_val,col=impact))+
  geom_point()+
  geom_errorbar(aes(xmin=IC_2.5,xmax=IC_97.5))+
  scale_x_continuous("mean |SHAP value|",label=function(x) abs(x))+
  geom_text(aes(x=IC_97.5,y=feature,label=abs(IC_97.5)),size=5,
            position=position_nudge(x=if_else(test_tmp$IC_97.5>=0,
                                              round(max(test_tmp$IC_97.5,na.rm=TRUE)/10),
                                              -round(max(test_tmp$IC_97.5,na.rm=TRUE)/10))))+
  scale_y_discrete("")+
  scale_color_manual("Direction:",na.value="white",
                     values=rev(c("neutral"="grey80","promoting predictor"="#C35C33","mitigating predictor"="#40B696",
                                  "uncertain"="black","undefined"="grey50")))+
  theme_Gaia()+
  theme(legend.position="bottom")

ggsave(plot1,file=paste0("ERS_SHAP_errorbar_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)

## Plot 2: barplot (mean |SHAP| with labels)
ordre_def2<-test_tmp %>% arrange(desc(mean_val))
test_tmp$feature<-factor(test_tmp$feature, levels=rev(unique(ordre_def2$feature)))

plot2<-ggplot(test_tmp,aes(y=feature,x=mean_val,fill=impact))+
  geom_bar(stat="identity",col="black")+
  geom_text(aes(x=mean_val,y=feature,label=mean_SHAP),size=6,
            position=position_nudge(x=if_else(test_tmp$mean_val>=0,
                                              round(max(test_tmp$mean_val,na.rm=TRUE)/10),
                                              -round(max(test_tmp$mean_val,na.rm=TRUE)/10))))+
  scale_x_continuous("mean |SHAP value|",label=function(x) abs(x))+
  scale_y_discrete("")+
  scale_fill_manual("Direction:",na.value="white",
                    values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2","neutral"="white","uncertain"="grey","undefined"="grey50"))+
  theme_Gaia()+
  theme(legend.position="bottom")

ggsave(plot2,file=paste0("ERS_SHAP_barplot_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)

## Plot 3: dispersion (beeswarm-like, colored by standardized feature value)
bornes_shap<-test %>% summarise(min=min(shap_value,na.rm=TRUE), max=max(shap_value,na.rm=TRUE))
test3<-test %>%
  filter(mean_SHAP!="0") %>%
  group_by(feature) %>%
  mutate(feature_value2=(feature_value-min(feature_value,na.rm=TRUE))/
           (max(feature_value,na.rm=TRUE)-min(feature_value,na.rm=TRUE))) %>%
  ungroup()
ordre_Def<-test3 %>% select(feature,mean_val) %>% distinct() %>%
  mutate(mean_val=as.numeric(mean_val)) %>% arrange(mean_val)
test3$feature<-factor(test3$feature, levels=ordre_Def$feature)
test2<-test3 %>% select(feature,mean_SHAP) %>% distinct()

plot3<-ggplot(test3,aes(x=feature,y=shap_value,colour=feature_value2))+
  geom_jitter(width=0.1)+
  scale_color_gradientn("Feature value",labels=c("Low","High"),breaks=c(0,1),limits=c(0,1),
                        colours=c("blue",scales::alpha("blue",0.75),scales::alpha("blue",0.5),
                                  scales::alpha("blue",0.35),scales::alpha("#FFA5001A",0.25),
                                  scales::alpha("red",0.35),scales::alpha("red",0.5),
                                  scales::alpha("red",0.75),"red"),
                        na.value="lightgrey",
                        guide=guide_colorbar(barheight=50,
                                             title.theme=element_text(angle=90,hjust=0.5,vjust=0,size=16,face="bold"),
                                             title.position="left"))+
  geom_abline(slope=0,intercept=0,colour="darkgrey")+
  scale_x_discrete("Feature")+
  scale_y_continuous("SHAP value",limits=c(bornes_shap$min,bornes_shap$max*1.5))+
  coord_flip()+
  geom_text(data=test2,aes(x=feature,y=bornes_shap$max*1.35,label=mean_SHAP),size=5,col="black")+
  theme_Gaia()

ggsave(plot3,file=paste0("ERS_SHAP_dispersion_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)


## Plot 4: SHAP percentage contribution within ERS

# Note: - percentage is computed as the share of each exposure in the total mean |SHAP|, i.e., its relative contribution to the ERS.
#       - this is not a percentage of outcome variance explained.

test_tmp<-test_tmp %>% mutate(SHAP_per=abs(mean_val)/sum(abs(mean_val))*100) # calculating the SHAP percentage contribution within ERS

ordre_def2<-test_tmp %>% arrange(desc(SHAP_per))
test_tmp$feature<-factor(test_tmp$feature, levels=rev(unique(ordre_def2$feature)))

plot4<-ggplot(test_tmp, aes(y=feature, x=SHAP_per, fill=impact))+
  geom_bar(stat="identity", col="black")+
  geom_text(aes(x=SHAP_per, y=feature, label=signif(abs(SHAP_per),3)), size=5,
            position=position_nudge(x=if_else(test_tmp$SHAP_per>=0,
                                              round(max(test_tmp$SHAP_per,na.rm=TRUE)/10),
                                              -round(max(test_tmp$SHAP_per,na.rm=TRUE)/10))))+
  scale_x_continuous("% contribution to ERS (based on mean |SHAP|)", label=function(x) paste0(abs(x)))+
  scale_y_discrete("")+
  scale_fill_manual("Direction:",na.value="white",
                    values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2","neutral"="white","uncertain"="grey","undefined"="grey50"))+
  theme_Gaia()+
  theme(legend.position="bottom")

ggsave(plot4,file=paste0("ERS_SHAP_contribution_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)

## Plot 5: direction comparison heatmap across 3 methods + consensus
dir_long<-direction_impact %>%
  select(feature,conventional,gam,pairwise,consensus) %>%
  pivot_longer(-feature, names_to="method", values_to="direction") %>%
  mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","undefined","uncertain")),
         method=factor(method,
                       levels=c("conventional","gam","pairwise","consensus"),
                       labels=c("Conventional\n(mean sign)","GAM derivative","Pairwise bins","Consensus")))

plot5<-ggplot(dir_long,aes(x=method,y=feature,fill=direction))+
  geom_tile(color="black")+
  scale_fill_manual("Direction:",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","uncertain"="grey25",
                                          "mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
  scale_x_discrete("",expand=c(0,0))+
  scale_y_discrete("",expand=c(0,0))+
  theme_Gaia()+
  theme(legend.position="top")

ggsave(plot5,file=paste0("ERS_SHAP_direction_comparison_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)

## Plot 6: SHAP dependence plots (one per exposure, if 15 or fewer exposures)
if(length(unique(tempopo$feature))<=15){
  
  plot_shap_dep<-function(data, feat, feature_value, shap_value){
    ggplot(data %>% filter(feature==feat), aes(x=feature_value, y=shap_value))+
      geom_point(alpha=0.3, color="#2166ac")+
      geom_smooth(method="loess", se=TRUE, color="#d6604d", linewidth=1.1)+
      geom_hline(yintercept=0, linetype="dashed", color="grey50")+
      labs(x=feat, y="SHAP value")+
      theme_Gaia()
  }
  
  dep_plots<-lapply(unique(tempopo$feature), function(f)
    plot_shap_dep(data=tempopo, feat=f, feature_value=feature_value, shap_value=shap_value))
  
  plot6<-patchwork::wrap_plots(dep_plots, ncol=3)+ patchwork::plot_annotation(title="SHAP dependence plots")
  
  ggsave(plot6,file=paste0("ERS_SHAP_dependence_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)
}

#----------------------------------------------------------------
#### SHAP and Cox regression comparison table ####

tmp_Cox<-Cox_tab_res %>% select(term,HR,direction,concordance,p.value.log)
colnames(tmp_Cox)<-c("feature","HR [95% CI]","Cox direction","Cox C-index","Cox LR p-value")

test_SHAP_tmp<-test %>%
  select(feature,mean_SHAP,conventional,gam,pairwise,consensus,mean_val) %>%
  distinct() %>%
  arrange(desc(abs(as.numeric(mean_val)))) %>%
  rowwise() %>% mutate_if(is.numeric,test_format) %>% ungroup() %>%
  select(-mean_val)

colnames(test_SHAP_tmp)<-c("feature","mean |SHAP| [95% CI]",
                           "Conventional direction","GAM-based direction","Pairwise-based direction","Consensus direction")

Compa_Cox_SHAP<-left_join(test_SHAP_tmp, tmp_Cox, by="feature")

#----------------------------------------------------------------
#### Saving results ####

fwrite(step3_result$outer_summary, paste0("ERS_CV_summary_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(step3_result$metric_test, paste0("ERS_test_performance_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(Shap_val, paste0("ERS_Shap_val_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(test_tmp, paste0("ERS_test_tmp_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(test, paste0("ERS_test_",Sys.Date(),".txt"), sep=";",row.names=FALSE)
fwrite(shap, paste0("ERS_shap_",Sys.Date(),".txt"), sep=";",row.names=FALSE)
fwrite(ers_direction_summary, paste0("ERS_overall_direction_summary_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(Cox_tab_res, paste0("ERS_Univariate_Cox_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(Compa_Cox_SHAP, paste0("ERS_Compa_Cox_SHAP_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(as_tibble(step3_result$final_params), paste0("ERS_best_params_",Sys.Date(),".csv"), sep=";",row.names=FALSE)

## Saving exposure-level adjustment metrics
if(use_cohort){
  fwrite(metrics_cohort %>% pivot_wider(names_from=adjustment,values_from=c(C_index,logHR,HR,Cox_p)),
         paste0("ERS_exposure_metrics_cohort_adjustment_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
}
fwrite(metrics_cov %>% pivot_wider(names_from=adjustment,values_from=c(C_index,logHR,HR,Cox_p)),
       paste0("ERS_exposure_metrics_covariate_adjustment_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
