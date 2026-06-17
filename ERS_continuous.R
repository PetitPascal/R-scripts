#-------------------------------------------------------------------------------
## Reproducible and generalisable script for calculating an exposomic risk score (ERS) for a continuous outcome using XGBoost with nested cross-validation
#
# Method: residualizing outcome on adjustment variables (cohort membership and/or covariates)
#
# Steps:
# - Step 1 [optional, multi-cohort]: outcome is cleaned of between-cohort differences
# - Step 2: outcome is additionally cleaned of covariate effects
# - Step 3: modeling the doubly-residualized outcome using exposure variables only -> the predicted values f(X) represent the exposure-driven component of Y (i.e., ERS), after removing cohort and covariate effects
#
#
# How to use the ERS in a final model:
# - glm(y ~ ers + age + sex + bmi, data=sim_data, family='gaussian')
# - or stratify into tertiles
#
# Settings to adjust before running:
# - use_cohort: TRUE/FALSE
# - cov_names: names of covariate adjustment variables
#-------------------------------------------------------------------------------

#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing and loading packages
pack_needed<-c("data.table","tidyverse","mllrnrs","broom","doParallel","foreach","splitTools","conflicted","grid","gridExtra","RColorBrewer","mlbench",
                 "mlexperiments","caret","MLmetrics","patchwork","performance","xgboost","parallel","here","scales","dplyr","ggplot2","tidyr",
                 "tibble","mgcv","ggforce","gratia","Rcpp","Metrics","MASS","shapr","iml")

is_installed<-pack_needed %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed==FALSE)){
  install.packages(pack_needed[!is_installed],repos="http://cran.us.r-project.org")
}
invisible(lapply(pack_needed, library, character.only=TRUE))

## Preventing package conflicts
conflict_prefer("select","dplyr")
conflict_prefer("filter","dplyr")
conflict_prefer("slice","dplyr")
conflict_prefer("alpha","scales")

## Setting the working directory
here::here("XGBoost - ERS calculation - linear outcome")

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

#- - - - - -
## SHAP interaction summary
#
# Goal: to summarize SHAP interaction values
#
# Outputs: dataset with columns:
# - feature: feature name
# - interaction_strength: mean absolute off-diagonal SHAP interaction -> high values indicate this feature interacts strongly with others
# - main_effect_mean: mean absolute diagonal (main effect) SHAP -> useful to compare interaction vs. main effect magnitude
# - features with high interaction_strength may have non-monotone effects due to interactions
# - if main_effect_mean >> interaction_strength, the feature is mostly additive

shap_interaction_summary<-function(xgb_fit, X_mat, feature_cols){
  inter<-predict(xgb_fit, X_mat, predinteraction = TRUE)
  p<-length(feature_cols)
  
  arr<-if(length(dim(inter)) == 3){
    inter
  }else{
    array(inter, dim = c(nrow(X_mat), p + 1, p + 1))
  }
  
  dplyr::bind_rows(lapply(seq_len(p), function(j){
    off<-arr[, j, 1:p, drop = FALSE]
    if(length(dim(off)) == 3) off<-off[, 1, , drop = FALSE]
    main<-arr[, j, j]
    offdiag<-rowSums(abs(off), na.rm = TRUE) - abs(main)
    
    data.frame(feature = feature_cols[j],
               interaction_strength = mean(offdiag, na.rm = TRUE),
               main_effect_mean = mean(main, na.rm = TRUE))
  }))
}

#- - - - - -
## Fitting XGBoost with nested CV and returning OOF predictions + final model

# Arguments:
#   - dataset_dt: data.table with outcome in col 1 and features in remaining cols
#   - target_col: name of the outcome column
#   - feature_cols: vector of feature column names
#   - objective: XGBoost objective (e.g., "reg:squarederror")
#   - eval_metric: XGBoost eval metric (e.g., "rmse")
#   - train_split: training set split (e.g., 0.7)
#   - test_split: holdout test set split (e.g., 0.3)
#   - higher_better: TRUE if higher metric = better (e.g. AUROC), FALSE if lower (e.g., RMSE)
#   - select_best: function to pick best outer fold index (e.g., which.min or which.max)
#   - nb_inner_fold: number of inner folds for hyperparameter tuning (e.g., 5)
#   - nb_outer_fold: number of outer folds for nested CV (e.g., 5)
#   - perf_extractor: function to extract scalar performance from one outer fold result
#
# Returns a list with:
#   - oof_preds: out-of-fold (OOF) predictions on the full dataset (needed for residuals)
#   - test_preds: predictions on holdout test set
#   - metric_test: performance table on holdout test set
#   - final_model: model retrained on full training set with best hyperparameters
#   - outer_summary: cross-validation (CV) performance table
#   - final_params: best hyperparameters selected

fit_ers_step<-function(dataset_dt,
                        target_col,
                        feature_cols,
                        objective,
                        eval_metric,
                        higher_better,
                        select_best,
                        perf_extractor,
                        extra_learner_args=list(),
                       train_split,
                       test_split,
                       nb_inner_fold,
                       nb_outer_fold,
                        param_grid=NULL){
  
  ## ensuring training and test splits are in correct format
  train_split<-as.numeric(train_split)
  test_split<-as.numeric(test_split)
  
  if((train_split+test_split>1)|is.na(test_split)|is.na(train_split)){
    train_split<-0.7
    test_split<-0.3
  }
  
  ## Ensuring the number of outer and inner folds, for nested CV and hyperparameter tuning, are in correct format
  nb_inner_fold<-as.numeric(nb_inner_fold)
  if(is.na(nb_inner_fold)|nb_inner_fold<2){
    nb_inner_fold<-5
  }
  nb_outer_fold<-as.numeric(nb_outer_fold)
  if(is.na(nb_outer_fold)|nb_outer_fold<2){
    nb_outer_fold<-5
  }
  
  ## Creating a 70/30 stratified split
  data_split<-splitTools::partition(y=dataset_dt[[target_col]],
    p=c(train=train_split, test=test_split),
    type="stratified",
    seed=seed)
  
  ## Creating training and test data sets
  X_train<-as.matrix(dataset_dt[data_split$train, ..feature_cols])
  y_train<-dataset_dt[[target_col]][data_split$train]
  X_test<-as.matrix(dataset_dt[data_split$test,  ..feature_cols])
  y_test<-dataset_dt[[target_col]][data_split$test]
  
  ## Default parameter grid if none provided
  if(is.null(param_grid)){
    param_grid<-expand.grid(subsample=seq(0.5,1,0.25),
                    colsample_bytree=seq(0.5,1,0.25),
      min_child_weight=c(1,5,10),
      learning_rate=c(0.05,0.1,0.3),
      max_depth=c(3,5,7)) %>%
      dplyr::slice_sample(n=30, replace=TRUE) # Limiting space to 30 combinations for computational efficiency and environmental sustainability considerations
  }
  
  ## Creating outer folds for nested CV
  outer_folds<-splitTools::create_folds(y_train, k=nb_outer_fold, type="stratified", seed=seed)
  outer_results<-list()
  best_params_all<-list()
  
  for(outer_idx in seq_along(outer_folds)){
    
    val_idx<-outer_folds[[outer_idx]]
    train_idx_cv<-setdiff(seq_len(nrow(X_train)), val_idx)
    
    X_tr<-X_train[train_idx_cv,,drop=FALSE]
    y_tr<-y_train[train_idx_cv]
    X_val<-X_train[val_idx,,drop=FALSE]
    y_val<-y_train[val_idx]
    
    ## Inner folds for hyperparameter tuning
    inner_folds<-splitTools::create_folds(y_tr, k=nb_inner_fold, type="stratified", seed=seed)
    
    best_perf<-ifelse(higher_better, -Inf, Inf)
    best_params<-NULL
    
    for(i in seq_len(nrow(param_grid))){
      xgb_cv<-mlexperiments::MLCrossValidation$new(learner=mllrnrs::LearnerXgboost$new(metric_optimization_higher_better=higher_better),
        fold_list=inner_folds,
        ncores=2,
        seed=123)
      
      xgb_cv$learner_args<-c(as.list(param_grid[i,]),list(objective=objective, eval_metric=eval_metric, nrounds=100L),
        extra_learner_args)
      
      ## Setting performance metric and data
      xgb_cv$performance_metric<-mlexperiments::metric("rmse")
      xgb_cv$set_data(x=X_tr, y=y_tr)
      
      res_cv<-xgb_cv$execute()
      
      ## Extracting mean performance across inner folds
      mean_perf<-tryCatch(mean(sapply(res_cv$results$folds, perf_extractor), na.rm=TRUE),
        error=function(e) NA)
      
      if(!is.na(mean_perf)){
        better<-if(higher_better) mean_perf>best_perf else mean_perf<best_perf
        if(better){
          best_perf<-mean_perf
        best_params<-param_grid[i,]
        }
      }
    }
    
    ## Fallback to default params if tuning failed
    if(is.null(best_params)){
      best_params<-data.frame(subsample=1,
                              colsample_bytree=1,
                              min_child_weight=1,
                              learning_rate=0.1,
                              max_depth=3)
    }
    best_params_all[[outer_idx]]<-best_params
    
    ## Training outer model and evaluating on outer validation fold
    dtrain<-xgboost::xgb.DMatrix(data=X_tr, label=y_tr)
    xgb_outer<-xgboost::xgb.train(params=c(as.list(best_params),list(objective=objective, eval_metric=eval_metric),
               extra_learner_args),
      data=dtrain,
      nrounds=100,
      verbose=0)
    
    val_pred<-predict(xgb_outer, xgboost::xgb.DMatrix(X_val))
    
    ## Computing performance metrics on outer fold
    R2<-1-sum((y_val-val_pred)^2)/sum((y_val-mean(y_val))^2)
    RMSE<-Metrics::rmse(y_val, val_pred)
    MAE<-Metrics::mae(y_val, val_pred)
    r<-cor(y_val, val_pred)
    CCC<-(2*cov(y_val,val_pred,use="complete.obs"))/(var(y_val,na.rm=TRUE)+var(val_pred,na.rm=TRUE)+
         (mean(y_val,na.rm=TRUE)-mean(val_pred,na.rm=TRUE))^2)
    
    outer_results[[outer_idx]]<-data.frame(Fold=outer_idx, R2=R2, RMSE=RMSE, MAE=MAE, r=r, CCC=CCC)
  }
  
  outer_summary<-dplyr::bind_rows(outer_results)
  cat("\nNested CV summary:\n")
  print(outer_summary)
  
  ## Selecting best hyperparameters based on outer CV performance
  best_idx<-select_best(outer_summary$RMSE)
  final_params<-best_params_all[[best_idx]]
  
  ## Retraining final model on full training set
  dtrain_full<-xgboost::xgb.DMatrix(data=X_train, label=y_train)
  final_model<-xgboost::xgb.train(params=c(as.list(final_params),list(objective=objective, eval_metric=eval_metric),
             extra_learner_args),
    data=dtrain_full,
    nrounds=100,
    verbose=0)
  
  ## Computing OOF predictions on the full dataset (train + test)
  # - OOF are needed so that residuals in the next step are not overfitted.
  # - Each observation is predicted by a model that was not trained on it.
  # - Strategy: use the final model on test, and per-fold models on train OOF
 
  X_full<-as.matrix(dataset_dt[,..feature_cols])
  oof_preds_full<-rep(NA, nrow(dataset_dt))
  
  ## For training observations: use OOF from outer CV.We recompute per-fold predictions to cover all training observations
  for(outer_idx in seq_along(outer_folds)){
    val_idx_global<-data_split$train[outer_folds[[outer_idx]]]
    train_idx_cv<-setdiff(seq_len(length(data_split$train)),outer_folds[[outer_idx]])
    X_tr_oof<-X_train[train_idx_cv,,drop=FALSE]
    y_tr_oof<-y_train[train_idx_cv]
    dtrain_oof<-xgboost::xgb.DMatrix(data=X_tr_oof, label=y_tr_oof)
    
    mod_oof<-xgboost::xgb.train(params=c(as.list(best_params_all[[outer_idx]]),list(objective=objective, eval_metric=eval_metric),
               extra_learner_args),
      data=dtrain_oof,
      nrounds=100,
      verbose=0)
    
    oof_preds_full[val_idx_global]<-predict(mod_oof,xgboost::xgb.DMatrix(X_train[outer_folds[[outer_idx]],,drop=FALSE]))
  }
  
  ## For test observations: use final model (no leakage since it was trained on train only)
  oof_preds_full[data_split$test]<-predict(final_model,xgboost::xgb.DMatrix(X_test))
  
  ## Test set performance
  test_preds<-predict(final_model, xgboost::xgb.DMatrix(X_test))
  R2_test<-1-sum((y_test-test_preds)^2)/sum((y_test-mean(y_test))^2)
  RMSE_test<-Metrics::rmse(y_test, test_preds)
  MAE_test<-Metrics::mae(y_test, test_preds)
  r_test<-cor(y_test, test_preds)
  CCC_test<-(2*cov(y_test,test_preds,use="complete.obs"))/(var(y_test,na.rm=TRUE)+var(test_preds,na.rm=TRUE)+
       (mean(y_test,na.rm=TRUE)-mean(test_preds,na.rm=TRUE))^2)
  
  cat("\nTest set performance:\n")
  cat("R²:", round(R2_test,3), " RMSE:", round(RMSE_test,3),
      " MAE:", round(MAE_test,3), " r:", round(r_test,3), " CCC:", round(CCC_test,3), "\n")
  
  metric_test<-as_tibble(data.frame(R2_test=R2_test,
                       RMSE_test=RMSE_test,
                       MAE_test=MAE_test,
                       r_test=r_test,
                       CCC_test=CCC_test))
  
  return(list(oof_preds=oof_preds_full,  # OOF predictions for all observations -> used to compute residuals
    test_preds=test_preds,
    final_model=final_model,
    outer_summary=outer_summary,
    final_params=final_params,
    data_split=data_split,
    metric_test=metric_test,
    X_test=X_test,
    y_test=y_test))
}

#- - - - - -
## Safely extracting scalar performance from one inner fold
perf_extractor_metric<-function(f){
  val<-tryCatch(f$performance, error=function(e) NA)
  if(is.list(val)) val<-unlist(val)
  if(is.null(val)||length(val)==0) return(NA)
  return(as.numeric(val[1]))
}

#----------------------------------------------------------------
#### Setting ####

## Set to TRUE if data comes from multiple cohorts (activates Step 1)
use_cohort<-TRUE

## Names of covariates to adjust for in Step 3
cov_names<-c("age","sex","bmi")

#----------------------------------------------------------------
#### Data import and preprocessing (replace simulation with your own data) ####

## Simulating data for demonstration
set.seed(seed)
n<-800
K<-8

## Building exposure matrix

# Simulating continuous exposures
expo_df<-as.data.frame(exp(MASS::mvrnorm(n, rep(0, K), 0.5^as.matrix(dist(1:K)))/3))

# Simulating binary exposure
expo_df2<-data.frame(X9=sample(x=c(0,1), n, replace = TRUE))
expo_df<-bind_cols(expo_df,expo_df2)
colnames(expo_df)<-paste0("X", 1:(K+1))
expo_names<-colnames(expo_df)

# Building exposure matrix
expo_mat<-as.matrix(cbind(log1p(expo_df[,paste0("X", 1:K)]), # log-transform continuous variables only 
                          expo_df[,"X9"])) # binary/ordinal untransformed

colnames(expo_mat)<-expo_names

## Adding cohort                           
cov_df<-data.frame(age=rnorm(n, 50, 10), sex=rbinom(n, 1, 0.5), bmi=rnorm(n, 26, 4))
cohort<-sample(c("cohort_A","cohort_B","cohort_C"), n, replace=TRUE)

## Cohort effect: cohort_A = reference, cohort_B adds 1.5, cohort_C subtracts 1
cohort_eff<-ifelse(cohort=="cohort_A", 0, ifelse(cohort=="cohort_B", 1.5, -1))

## True exposure effect (what ERS should capture after removing cohort + covariate effects)
h_z<-as.numeric(as.matrix(expo_mat)%*%c(0.5,0.3,-0.2,0.1,0.4,-0.1,0.05,0.05,1))

## Continuous outcome
y<-2 + 0.03*cov_df$age + 0.5*cov_df$sex + 0.1*cov_df$bmi + cohort_eff + h_z + rnorm(n, 0, 1)

## Assembling dataset
sim_data<-cbind(data.frame(y=y, cohort=cohort), cov_df, expo_mat)

#----------------------------------------------------------------
#### Step 1 (optional): residualizing outcome (y) on cohort membership (z) (only if multi-cohort) ####

## Goal: to remove between-cohort differences from Y so that the ERS is not driven by which cohort a participant belongs to.

## Method:
# - fitting XGBoost predicting Y from one-hot encoded cohort indicators.
# - The residuals U = Y - predicted(cohort) are an outcome "cleaned" of cohort effects.
# - OOF predictions are used for all observations (train + test) to avoid overfitting the residuals. 
# - Each observation's residual is computed using a model that was not trained on that observation.

## Outputs:
# - U=residuals after removing cohort effect
# - R²=proportion of Y variance explained by cohort membership (e.g., R²=0.15 means 15% of outcome variance was cohort-driven)

if(use_cohort){
  
  ## One-hot encoding of cohort membership (no intercept to avoid collinearity)
  Z_mat<-model.matrix(~as.factor(sim_data$cohort)-1)
  colnames(Z_mat)<-paste0("cohort_",levels(as.factor(sim_data$cohort)))
  
  dataset_step1<-as.data.table(cbind(data.frame(Y=y), as.data.frame(Z_mat)))
  cohort_cols<-colnames(Z_mat)
  
  step1_result<-fit_ers_step(dataset_dt=dataset_step1,
    target_col="Y",
    feature_cols=cohort_cols,
    objective="reg:squarederror",
    eval_metric="rmse",
    higher_better=FALSE,
    select_best=which.min,
    train_split=0.7,
    test_split=0.3,
    nb_outer_fold=5,
    nb_inner_fold=5,
    perf_extractor=perf_extractor_metric)
  
  ## U = Y - g(cohort): outcome cleaned of cohort effects
  U<-y - step1_result$oof_preds
  
}else{ ## If no cohort adjustment needed, U = Y unchanged
  U<-y
}

#----------------------------------------------------------------
#### Step 2: residualizing U on covariates ####

## Goal: to remove the effects of age, sex, BMI (and any other covariate) from U, so that the ERS captures only the exposure-specific signal.

## Method: 
# - fitting XGBoost predicting U from covariates.
# - the residuals V = U - predicted(covariates).
# - OOF predictions used for the same reason as Step 1.

## Outputs:
# - V = residuals after removing covariate effects from U
# - R² = proportion of U variance explained by covariates (e.g., R²=0.20 means 20% of U was covariate-driven)

dataset_step2<-as.data.table(cbind(data.frame(U=U), sim_data[,cov_names]))

step2_result<-fit_ers_step(dataset_dt=dataset_step2,
  target_col="U",
  feature_cols=cov_names,
  objective="reg:squarederror",
  eval_metric="rmse",
  higher_better=FALSE,
  train_split=0.7,
  test_split=0.3,
  nb_outer_fold=5,
  nb_inner_fold=5,
  select_best=which.min,
  perf_extractor=perf_extractor_metric)

## V = U - f(covariates): outcome cleaned of both cohort and covariate effects
V<-U - step2_result$oof_preds

#----------------------------------------------------------------
#### Step 3: Fitting ERS model on exposures (nested CV) ####

## Goal: to model the doubly-residualized outcome V using exposure variables only.
#        The predicted values f(X) represent the exposure-driven component of Y (i.e., ERS), after removing cohort and covariate effects.

## Method: nested CV

## Outputs:
# - ers: f(X) = predicted V from exposures = Exposomic Risk Score (ERS)
# - R²: proportion of V variance explained by exposures. This is the "pure" exposure effect size, cleanly separated from cohort and covariate effects.
# - correlation(ers, Y) = overall association between ers and raw outcome
# - correlation(ers, U) = association between ers and outcome after cohort effect removal
# - correlation(ers, V) = association between ers and outcome after cohort and covariate effect removal

expo_mat<-as.matrix(sim_data[,expo_names])

dataset_step3<-as.data.table(cbind(data.frame(V=V), as.data.frame(expo_mat)))

step3_result<-fit_ers_step(dataset_dt=dataset_step3,
  target_col="V",
  feature_cols=expo_names,
  objective="reg:squarederror",
  eval_metric="rmse",
  higher_better=FALSE,
  select_best=which.min,
  train_split=0.7,
  test_split=0.3,
  nb_outer_fold=5,
  nb_inner_fold=5,
  perf_extractor=perf_extractor_metric)

# Using final model predictions on full dataset to compute ERS
ers<-predict(step3_result$final_model, xgboost::xgb.DMatrix(expo_mat))
sim_data$ers<-ers

cat("\nERS summary:\n")
cat("Correlation(ERS, Y):", round(cor(ers,y),3)," R²(ERS~Y):", round(cor(ers,y)^2,3), "\n")
if(use_cohort){
  cat("Correlation(ERS, U):", round(cor(ers,U),3)," R²(ERS~U):", round(cor(ers,U)^2,3), "\n")
}
cat("Correlation(ERS, V):", round(cor(ers,V),3)," R²(ERS~V):", round(cor(ers,V)^2,3),
    "(= pure exposure effect after cleaning)\n")

## Saving ERS model
saveRDS(step3_result$final_model, paste0("ERS_model_continuous_",Sys.Date(),".rds"))
xgb.save(step3_result$final_model,paste0("ERS_model_continuous_",Sys.Date(),".ubj"))

#----------------------------------------------------------------
#### Step 4: Overall ERS direction (mixture effect) ####

# Goal: to determine whether the overall ERS promotes or mitigates the outcome using a 3-method consensus framework

## Building a two-column data frame: ERS value vs. residualized outcome V (doubly-residualized outcome)
ers_direction_df<-data.frame(ers=ers, V=V, U=U, y=y)

## Approach 1: sign of mean association
ers_mean_sign<-conventional_direction(ers_direction_df, feature_col="ers", shap_col="V")

## Approach 2: GAM derivative + Spearman fallback
ers_gam<-gam_direction(ers_direction_df, feature_col="ers", shap_col="V")

## Approach 3: pairwise bin concordance
ers_pairwise<-pairwise_direction(ers_direction_df, feature_col="ers", shap_col="V",majority_threshold=0.55, n_quantile_bins=200)

## Consensus: majority vote across the three approaches
ers_consensus_vec<-c(ers_mean_sign, ers_gam, ers_pairwise)
ers_tbl<-sort(table(ers_consensus_vec), decreasing=TRUE)
ers_consensus<-if(ers_tbl[1]>=2) names(ers_tbl)[1] else "uncertain"

## Summary table
ers_direction_summary<-tibble(mean_sign_direction=ers_mean_sign,
                              gam_direction=ers_gam,
                              pairwise_direction=ers_pairwise,
                              consensus_direction=ers_consensus,
                              cor_ERS_Y=round(cor(ers, y), 3),
                              cor_ERS_U=round(cor(ers, U), 3),
                              cor_ERS_V=round(cor(ers, V), 3),
                              R2_ERS_Y=round(cor(ers, y)^2, 3),
                              R2_ERS_U=round(cor(ers, U)^2, 3),
                              R2_ERS_V=round(cor(ers, V)^2, 3))


cat("\nERS global direction:\n")
ers_direction_summary

#----------------------------------------------------------------
#### Step 5: Using ERS in a final association model ####

## The ERS can be used as a predictor in standard models alongside covariates.
## Since V (doubly-residualized outcome) was used to build the ERS, the ERS is by construction largely orthogonal to cohort and covariates. 
## The coefficient of ERS reflects the association between the exposure-driven component of Y and the raw outcome Y, on the residual scale of V.

## Linear regression: ERS + covariates -> raw outcome Y
final_model_lm<-lm(y~ers+age+sex+bmi, data=sim_data)
summary(final_model_lm)

## Coefficient table: beta = mean difference in Y per unit increase in ERS
broom::tidy(final_model_lm, conf.int=TRUE)

## Model diagnostics (linearity, homoscedasticity, normality of residuals, influential points)
performance::check_model(final_model_lm)
check_autocorrelation(final_model_lm) # independence of residuals
check_normality(final_model_lm) # normality of residuals
check_heteroscedasticity(final_model_lm) # homoscedasticity of residuals
check_outliers(final_model_lm) # outliers

## Saving model summary
fwrite(broom::tidy(final_model_lm, conf.int=TRUE),paste0("ERS_final_model_coefficients_",Sys.Date(),".csv"),sep=";", row.names=FALSE)

## Stratification into tertiles, with beta = mean difference in Y between medium/high vs. low ERS tertile
sim_data$ers_group<-ntile(sim_data$ers, 3)
sim_data$ers_group<-factor(sim_data$ers_group, labels=c("Low","Medium","High"))

ers_tertile_model<-lm(y~ers_group+age+sex+bmi, data=sim_data)
summary(ers_tertile_model)
broom::tidy(ers_tertile_model, conf.int=TRUE)

fwrite(broom::tidy(ers_tertile_model, conf.int=TRUE),paste0("ERS_tertile_model_coefficients_",Sys.Date(),".csv"),sep=";", row.names=FALSE)

## Figure: mean Y by ERS tertile group
sim_data %>%
  group_by(ers_group) %>%
  summarise(mean_y=mean(y, na.rm=TRUE),
            se_y=sd(y, na.rm=TRUE)/sqrt(n()),
            ci_lo=mean_y - 1.96*se_y,
            ci_hi=mean_y + 1.96*se_y) %>%
  ggplot(aes(x=ers_group, y=mean_y, ymin=ci_lo, ymax=ci_hi, color=ers_group))+
  geom_point(size=4)+
  geom_errorbar(width=0.2, linewidth=1)+
  scale_color_manual(values=c(Low="#40B696", Medium="#F1CB0E", High="#C35C33"))+
  labs(title="Mean outcome by ERS tertile",
       x="ERS group", y="Mean outcome (Y)")+
  theme_Gaia()+
  theme(legend.position="none")

#----------------------------------------------------------------
#### Step 6: Visualization ####

# Figure 1: ERS distribution
ggplot(sim_data, aes(x=ers)) +
  geom_histogram(fill="#A6DDCE", color="black", bins=40) +
  labs(title="Distribution of ERS", x="ERS", y="Count") +
  theme_Gaia()

# Figure 2: ERS vs. outcome
ggplot(sim_data, aes(x=ers, y=y)) +
  geom_point(alpha=0.3, color="#2166ac") +
  geom_smooth(method="gam", formula=y ~ s(x, bs="cs"), color="#d6604d", linewidth=1.5) +
  labs(title=paste0("R2=", round(cor(ers,y)^2,3)),
       x="ERS", y="Outcome") +
  theme_Gaia()

## Figure 3: ERS vs. residualized outcome V
ggplot(ers_direction_df, aes(x=ers, y=U)) +
  geom_point(alpha=0.3, color="#2166ac") +
  geom_smooth(method="gam", formula=y ~ s(x, bs="cs"), color="#d6604d", linewidth=1.5) +
  labs(title=paste0("R2=",round(cor(ers,U)^2,3)),
       x="ERS", y="Residualized outcome U") +
  theme_Gaia()

## Figure 4: ERS vs. residualized outcome V
ggplot(ers_direction_df, aes(x=ers, y=V)) +
  geom_point(alpha=0.3, color="#2166ac") +
  geom_smooth(method="gam", formula=y ~ s(x, bs="cs"), color="#d6604d", linewidth=1.5) +
  labs(title=paste0("R2=", round(cor(ers,V)^2,3)),
       x="ERS", y="Residualized outcome V") +
  theme_Gaia()

# Figure 5: ERS by cohort (if applicable)
if(use_cohort){
  ggplot(sim_data, aes(x=cohort, y=ers, fill=cohort)) +
    geom_violin(fill="transparent")+
    geom_boxplot(alpha=0.7,width=0.4) +
    scale_fill_brewer(palette="Set2") +
    labs(title="ERS by cohort", x="", y="ERS") +
    theme_Gaia() + 
    theme(legend.position="none")
}

#----------------------------------------------------------------
#### Step 7: Univariate GLM on exposures (test set) — for comparison with SHAP ####

Association_tab_res<-c()
X_test_ers<-step3_result$X_test

for(n_var in seq_along(expo_names)){
  
  test_tmpo<-as_tibble(X_test_ers) %>% select(all_of(expo_names[n_var]))
  colnames(test_tmpo)<-"tmp"
  test_tmpo<-bind_cols(data.frame(outcome=step3_result$y_test), test_tmpo)
  glm_model<-glm(outcome~tmp, data=test_tmpo, family="gaussian")
  mod_res<-broom::tidy(glm_model, exponentiate=FALSE, conf.int=TRUE) %>%
    filter(term!="(Intercept)") %>%
    mutate(term=expo_names[n_var],
           direction=case_when(conf.low>=0~"positive association",
                               conf.high<0~"negative association",
                               TRUE~"uncertain association")) %>%
    bind_cols(performance::performance(glm_model)) %>%
    rowwise() %>% mutate_if(is.numeric,test_format) %>% ungroup() %>%
    mutate(beta=paste(estimate," [",conf.low,"; ",conf.high,"]",sep="")) %>%
    mutate(beta=if_else(beta==" [; ]","",beta))
  
  Association_tab_res<-bind_rows(Association_tab_res,mod_res)
}

#----------------------------------------------------------------
#### Step 8: SHAP analysis ####

# - SHAP values are computed on the test set of Step 3.
# - SHAP values show which exposures contribute most to the ERS, and in which direction (promoting = higher ERS, mitigating = lower ERS).

## Computing SHAP values
contr<-predict(step3_result$final_model, as.matrix(X_test_ers), predcontrib=TRUE)
shap<-as_tibble(contr)
shap_contrib<-as.data.table(contr)

## Removing BIAS term
shap_contrib<-shap_contrib[,!grepl("bias|Intercept|BIAS|Bias",names(shap_contrib),ignore.case=TRUE),with=FALSE]

## Computing mean absolute SHAP score per feature (= importance ranking)
mean_shap_score<-colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing=TRUE)]

## Reshaping SHAP values to long format
shap_score_sub<-as.data.table(shap_contrib)[, names(mean_shap_score), with=FALSE]
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

tempopo<-left_join(tempopo,shap_tempo,by=c("id","feature"))
tempopo<-left_join(tempopo,Shap_val,by=base::intersect(colnames(tempopo),colnames(Shap_val)))

## Computing SHAP directionality using 3-method consensus
direction_impact<-compute_shap_directions_long(tempopo,threshold=0.55,n_bins=200) %>%
  select(-c(n,mean_shap,median_shap)) %>%
  filter(feature %ni% c("(Intercept)","BIAS","Bias","bias")) %>%
  mutate(consensus=apply(cbind(conventional,gam,pairwise),1,function(x){
    tbl<-sort(table(x),decreasing=TRUE)
    if(tbl[1]>=2) names(tbl)[1] else "uncertain"
  }))

## Merging direction with SHAP long format
test<-left_join(tempopo,direction_impact,by="feature") %>% mutate(direction=consensus) # selecting the consensus-based approach, adapt if needed
test<-left_join(test,Shap_val,by=base::intersect(colnames(test),colnames(Shap_val)))

## Preparing data for plots: applying sign to mean |SHAP| based on direction
test_tmp<-test %>%
  mutate(impact=direction) %>%
  select(feature,mean_val,IC_2.5,IC_97.5,impact,mean_SHAP) %>%
  distinct() %>%
  mutate(across(c(mean_val,IC_2.5,IC_97.5),as.numeric)) %>%
  mutate(mean_val=if_else(impact=="mitigating",-mean_val,mean_val),
         IC_2.5=if_else(impact=="mitigating",-IC_2.5, IC_2.5),
         IC_97.5=if_else(impact=="mitigating",-IC_97.5, IC_97.5)) %>%
  arrange(desc(mean_val)) %>%
  filter(mean_val!=0, !is.na(mean_val))

## Ordering features for plots
ordre_def<-bind_rows(test_tmp %>% filter(mean_val>=0) %>% arrange(desc(IC_2.5)),
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

## Plot 3: SHAP percentage contribution within ERS

# Note: - percentage is computed as the share of each exposure in the total mean |SHAP|, i.e., its relative contribution to the ERS.
#       - this is not a percentage of outcome variance explained.

test_tmp<-test_tmp %>% mutate(SHAP_per=abs(mean_val)/sum(abs(mean_val))*100) # calculating the SHAP percentage contribution within ERS

ordre_def2<-test_tmp %>% arrange(desc(SHAP_per))
test_tmp$feature<-factor(test_tmp$feature, levels=rev(unique(ordre_def2$feature)))

plot3<-ggplot(test_tmp, aes(y=feature, x=SHAP_per, fill=impact))+
  geom_bar(stat="identity", col="black")+
  geom_text(aes(x=SHAP_per, y=feature, label=signif(abs(SHAP_per),3)), size=5,
            position=position_nudge(x=if_else(test_tmp$SHAP_per>=0,
                                              round(max(test_tmp$SHAP_per,na.rm=TRUE)/10),
                                              -round(max(test_tmp$SHAP_per,na.rm=TRUE)/10))))+
  scale_x_continuous("% contribution to ERS (based on mean |SHAP|)", label=function(x) paste0(abs(x),"%"))+
  scale_y_discrete("")+
  scale_fill_manual("Direction:",na.value="white",
                    values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2","neutral"="white","uncertain"="grey","undefined"="grey50"))+
  theme_Gaia()+
  theme(legend.position="bottom")

ggsave(plot3,file=paste0("ERS_SHAP_contribution_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)

## Plot 4: dispersion (beeswarm-like, colored by standardized feature value)

bornes_shap<-test %>% summarise(min=min(shap_value,na.rm=TRUE), max=max(shap_value,na.rm=TRUE))
test3<-test %>%
  filter(mean_SHAP!="0") %>%
  group_by(feature) %>%
  mutate(feature_value2=(feature_value-min(feature_value,na.rm=TRUE))/(max(feature_value,na.rm=TRUE)-min(feature_value,na.rm=TRUE))) %>%
  ungroup()
ordre_Def<-test3 %>% select(feature,mean_val) %>% distinct() %>% mutate(mean_val=as.numeric(mean_val)) %>% arrange(mean_val)
test3$feature<-factor(test3$feature, levels=ordre_Def$feature)
test2<-test3 %>% select(feature,mean_SHAP) %>% distinct()

plot4<-ggplot(test3,aes(x=feature,y=shap_value,colour=feature_value2))+
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

ggsave(plot4,file=paste0("ERS_SHAP_dispersion_",Sys.Date(),".pdf"), dpi=600,width=60,height=30,units="cm",limitsize=FALSE)

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

## Plot 6: SHAP dependence plots

if(length(unique(tempopo$feature))<=15){
  
  plot_shap_dep<-function(data, feat,feature_value, shap_value){
    
    ggplot(data %>% filter(feature==feat), aes(x=feature_value, y=shap_value)) +
      geom_point(alpha = 0.3, color = "#2166ac") +
      geom_smooth(method = "loess", se = TRUE,color = "#d6604d", linewidth = 1.1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(x = feat, y = "SHAP value") +
      theme_Gaia()
    
  }
  
  dep_plots<-lapply(unique(tempopo$feature), function(f)
    plot_shap_dep(data=tempopo, feat=f,feature_value=feature_value,shap_value=shap_value))
  
  plot6<-patchwork::wrap_plots(dep_plots, ncol = 3) + patchwork::plot_annotation(title = "SHAP dependence plots")
  
  ggsave(plot6,file=paste0("ERS_SHAP_dependence_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)
}

## Computing interaction strengths
interaction_df<-shap_interaction_summary(step3_result$final_model, as.matrix(X_test_ers), expo_names)
interaction_df

#----------------------------------------------------------------
#### SHAP and GLM comparison table ####

tmp_GLM<-Association_tab_res %>% select(term,beta,direction,R2,p.value)
colnames(tmp_GLM)<-c("feature","Beta [95% CI]","GLM direction","GLM R2","GLM p-value")

test_SHAP_tmp<-test %>%
  select(feature,mean_SHAP,conventional,gam,pairwise,consensus,mean_val) %>%
  distinct() %>%
  arrange(desc(abs(as.numeric(mean_val)))) %>%
  rowwise() %>% mutate_if(is.numeric,test_format) %>% ungroup() %>%
  select(-mean_val)

colnames(test_SHAP_tmp)<-c("feature","mean |SHAP| [95% CI]","Conventional direction","GAM-based direction","Pairwise-based direction","Consensus direction")

Compa_GLM_SHAP<-left_join(test_SHAP_tmp,tmp_GLM,by="feature")

#----------------------------------------------------------------
#### Saving results ####

fwrite(step3_result$outer_summary, paste0("ERS_CV_summary_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(Shap_val,paste0("ERS_Shap_val_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(test_tmp,paste0("ERS_test_tmp_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(test,paste0("ERS_test_",Sys.Date(),".txt"), sep=";",row.names=FALSE)
fwrite(ers_direction_summary,paste0("ERS_verall direction summary_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(shap,paste0("ERS_shap_",Sys.Date(),".txt"), sep=";",row.names=FALSE)
fwrite(Association_tab_res, paste0("ERS_Univariate_GLM_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(Compa_GLM_SHAP, paste0("ERS_Compa_GLM_SHAP_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(interaction_df, paste0("ERS_SHAP interaction summary_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(as_tibble(step3_result$final_params), paste0("ERS_best_params_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
