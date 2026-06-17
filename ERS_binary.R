#-------------------------------------------------------------------------------
## Reproducible and generalisable script for calculating an exposomic risk score (ERS) for a binary outcome using XGBoost with nested cross-validation
#
# Method: Double Machine Learning (DML) - residualizing exposures on adjustment variables (cohort membership and/or covariates)
#
# Steps:
# - Step 1 [optional, multi-cohort]: W_j = X_j - E[X_j | cohort] -> each exposure is cleaned of between-cohort differences
# - Step 2: W2_j = W_j - E[W_j | covariates] -> each exposure is additionally cleaned of covariate effects
# - Step 3: ERS = f(W2) -> ERS = predicted log-odds from doubly-adjusted exposures
#
#
# How to use the ERS in a final model:
# - glm(y ~ ers + age + sex + bmi, data=sim_data, family=binomial())
# - or stratify into tertiles for OR and AUROC
#
# Settings to adjust before running:
# - use_cohort: TRUE/FALSE
# - cov_names: names of covariate adjustment variables
# - ers_metric: metric for ERS model tuning (step 3), with one of: "AUROC", "F1", "LogLoss"
#-------------------------------------------------------------------------------

#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing and loading packages
pack_needed<-c("data.table","tidyverse","mllrnrs","broom","doParallel","foreach","splitTools","conflicted","grid","gridExtra","RColorBrewer","mlbench",
               "mlexperiments","caret","MLmetrics","patchwork","performance","xgboost","parallel","here","scales","dplyr","ggplot2","tidyr",
               "tibble","mgcv","ggforce","gratia","Rcpp","Metrics","MASS","pROC","iml")

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
conflict_prefer("auc","pROC")
conflict_prefer("roc","pROC")
conflicts_prefer(caret::lift)

## Setting the working directory
here::here("XGBoost - ERS calculation - binary outcome")

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
## Finding optimal classification threshold by maximizing the chosen metric on a dataset
find_optimal_threshold<-function(probs, obs, metric){
  thresholds<-seq(0.01, 0.99, by=0.01)
  perf<-sapply(thresholds, function(thr){
    pred_bin<-as.integer(probs >= thr)
    if(length(unique(pred_bin))<2) return(NA)
    tryCatch({
      if(metric=="AUROC"){
        ## AUROC is threshold-independent; use 0.5 as default
        return(NA)
      }else if(metric=="F1"){
        MLmetrics::F1_Score(pred_bin, obs)
      }else if(metric=="LogLoss"){
        ## LogLoss is threshold-independent; use 0.5 as default
        return(NA)
      }
    }, error=function(e) NA)
  })
  if(all(is.na(perf))) return(0.5)
  ## For threshold-independent metrics (AUROC, LogLoss): use prevalence as threshold
  if(metric %in% c("AUROC","LogLoss")) return(mean(obs, na.rm=TRUE))
  thresholds[which.max(perf)]
}

#- - - - - -
## Computing all binary classification metrics between an exposure (or its residual) and Y

compute_binary_metrics<-function(x_vec, y_bin){
  x_num<-as.numeric(x_vec)
  y_int<-as.integer(y_bin)
  
  ## AUROC
  AUROC<-round(tryCatch(as.numeric(pROC::auc(pROC::roc(y_int,x_num,quiet=TRUE))),error=function(e) NA),3)
  
  ## PRAUC
  PRAUC<-round(tryCatch(MLmetrics::PRAUC(x_num,y_int),error=function(e) NA),3)
  
  ## LogLoss: rank-based probability proxy to handle raw continuous values outside [0,1]
  prob_proxy<-rank(x_num,na.last="keep")/(sum(!is.na(x_num))+1)
  LogLoss<-round(tryCatch(MLmetrics::LogLoss(prob_proxy,y_int),error=function(e) NA),3)
  
  ## Binarizing at median for threshold-dependent metrics
  thresh<-median(x_num,na.rm=TRUE)
  x_bin<-as.integer(x_num>thresh)
  F1<-round(tryCatch(MLmetrics::F1_Score(x_bin,y_int),error=function(e) NA),3)
  Sensitivity<-round(tryCatch(MLmetrics::Sensitivity(x_bin,y_int),error=function(e) NA),3)
  Specificity<-round(tryCatch(MLmetrics::Specificity(x_bin,y_int),error=function(e) NA),3)
  Balanced_accuracy<-round((Sensitivity+Specificity)/2,3)
  
  return(data.frame(AUROC=AUROC,PRAUC=PRAUC,LogLoss=LogLoss,F1=F1,Sensitivity=Sensitivity,Specificity=Specificity,Balanced_accuracy=Balanced_accuracy))
}

#- - - - - -
## Fitting XGBoost with nested CV for exposure adjustment (one exposure at a time)
#
# Rationale:
# - this function residualizes an exposure (not the outcome Y) on adjustment variables (e.g., cohort and/or covariates), following the double machine learning framework.
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
## Fitting XGBoost ERS model on doubly-adjusted exposures with nested CV
#
# scale_pos_weight is computed per outer fold to account for class imbalance.
#
# Arguments:
# - dataset_dt: dataset with binary outcome in col 1, residualized exposures in remaining cols
# - target_col: name of the binary outcome column (0/1)
# - feature_cols: vector of residualized exposure column names
# - ers_metric: metric for hyperparameter tuning and best fold selection (one of: "AUROC" (higher=better), "F1" (higher=better), "LogLoss" (lower=better))
# - train_split, test_split: dataset split proportions (training and holdout test)
# - nb_inner_fold, nb_outer_fold: number of folds for nested CV
# - param_grid: optional custom hyperparameter grid

fit_ers_binary<-function(dataset_dt,
                         target_col,
                         feature_cols,
                         ers_metric="AUROC",
                         train_split=0.7,
                         test_split=0.3,
                         nb_inner_fold=5,
                         nb_outer_fold=5,
                         param_grid=NULL){
  
  ## Validating ers_metric
  valid_metrics<-c("AUROC","F1","LogLoss")
  if(!ers_metric %in% valid_metrics){
    stop("ers_metric not recognised. Please choose one of:\n",
         "'AUROC': area under the ROC curve (recommended for balanced or moderately imbalanced datasets)\n",
         "'F1': harmonic mean of precision and recall (recommended when class imbalance is a concern)\n",
         "'LogLoss': log-loss / cross-entropy (recommended when predicted probability calibration matters)",
         call.=FALSE)
  }
  
  ## Mapping ers_metric to the corresponding native XGBoost eval_metric: AUROC -> "auc", LogLoss -> "logloss", F1 -> "aucpr" (best native XGBoost proxy for F1)
  xgb_eval_metric<-switch(ers_metric,
                          "AUROC"="auc",
                          "LogLoss"="logloss",
                          "F1"="aucpr")
  
  ## Whether higher or lower is better for the chosen eval_metric
  higher_better_inner<-switch(ers_metric,
                              "AUROC"=TRUE,
                              "LogLoss"=FALSE,
                              "F1"=TRUE)
  
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
  
  ## Creating a stratified train/test split (stratified on Y to preserve class balance)
  data_split<-splitTools::partition(y=dataset_dt[[target_col]],
                                    p=c(train=train_split, test=test_split),
                                    type="stratified", seed=seed)
  
  ## Creating training and test datasets
  X_train<-as.matrix(dataset_dt[data_split$train, ..feature_cols])
  y_train<-as.integer(dataset_dt[[target_col]][data_split$train])
  X_test<-as.matrix(dataset_dt[data_split$test,  ..feature_cols])
  y_test<-as.integer(dataset_dt[[target_col]][data_split$test])
  
  ## Computing global scale_pos_weight (ratio of majority to minority class)
  spw<-max(table(y_train))/min(table(y_train))
  
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
    
    X_tr<-X_train[train_idx_cv,, drop=FALSE]
    y_tr<-y_train[train_idx_cv]
    X_val<-X_train[val_idx,, drop=FALSE]
    y_val<-y_train[val_idx]
    
    ## Skipping outer fold if it contains only one class (cannot evaluate)
    if(length(unique(y_val))<2){
      cat("Outer fold",outer_idx,": skipped (single class in validation set)\n")
      best_params_all[[outer_idx]]<-NULL
      next
    }
    
    ## Computing scale_pos_weight per outer fold
    spw_fold<-max(table(y_tr))/min(table(y_tr))
    
    ## Initializing best performance tracker for inner CV
    best_inner_perf<-if(higher_better_inner) -Inf else Inf
    best_params<-NULL
    
    for(i in seq_len(nrow(param_grid))){
      
      ## Building parameter list with the eval_metric matching the chosen ers_metric
      params_i<-list(objective="binary:logistic",
        eval_metric=xgb_eval_metric,
        subsample=param_grid$subsample[i],
        colsample_bytree=param_grid$colsample_bytree[i],
        min_child_weight=param_grid$min_child_weight[i],
        learning_rate=param_grid$learning_rate[i],
        max_depth=param_grid$max_depth[i],
        scale_pos_weight=spw_fold)
      
      dtrain_inner<-xgboost::xgb.DMatrix(data=X_tr, label=y_tr)
      
      ## Inner CV via xgb.cv with stratified=TRUE to guarantee both classes in every inner fold
      cv_res<-tryCatch(suppressWarnings(xgboost::xgb.cv(params=params_i,data=dtrain_inner,nrounds=100,nfold=nb_inner_fold,early_stopping_rounds=10,verbose=0,stratified=TRUE)),
        error=function(e) NULL)
      
      if(is.null(cv_res)) next
      
      ## Extracting best performance across rounds for the chosen metric
      metric_col<-paste0("test_",xgb_eval_metric,"_mean")
      mean_perf<-tryCatch({
        vals<-cv_res$evaluation_log[[metric_col]]
        if(higher_better_inner) max(vals, na.rm=TRUE) else min(vals, na.rm=TRUE)
      }, error=function(e) NA)
      
      if(!is.na(mean_perf)){
        better<-if(higher_better_inner) mean_perf>best_inner_perf else mean_perf<best_inner_perf
        if(better){
          best_inner_perf<-mean_perf
          best_params<-params_i
          }
      }
    }
    
    ## Fallback to default params if all tuning attempts failed
    if(is.null(best_params)){
      best_params<-list(objective="binary:logistic",
        eval_metric=xgb_eval_metric,
        subsample=1,
        colsample_bytree=1,
        min_child_weight=1,
        learning_rate=0.1,
        max_depth=3,
        scale_pos_weight=spw_fold)
    }
    best_params_all[[outer_idx]]<-best_params
    
    ## Training outer model and evaluating on outer validation fold
    dtrain_outer<-xgboost::xgb.DMatrix(data=X_tr, label=y_tr)
    xgb_outer<-xgboost::xgb.train(params=best_params,data=dtrain_outer,nrounds=100,verbose=0)
    
    val_pred<-predict(xgb_outer, xgboost::xgb.DMatrix(X_val))
    
    ## Finding optimal classification threshold on the outer training fold predictions
    tr_pred<-predict(xgb_outer, xgboost::xgb.DMatrix(X_tr))
    opt_threshold<-find_optimal_threshold(tr_pred, y_tr, ers_metric)
    
    val_pred_bin<-as.integer(val_pred >= opt_threshold)
    
    ## Computing all performance metrics on outer fold
    AUROC<-tryCatch(as.numeric(pROC::auc(pROC::roc(y_val, val_pred, quiet=TRUE))), error=function(e) NA)
    F1<-tryCatch(MLmetrics::F1_Score(val_pred_bin, y_val),error=function(e) NA)
    logloss<-tryCatch(MLmetrics::LogLoss(val_pred, y_val), error=function(e) NA)
    PRAUC<-tryCatch(MLmetrics::PRAUC(val_pred, y_val), error=function(e) NA)
    Sensitivity<-tryCatch(MLmetrics::Sensitivity(val_pred_bin, y_val), error=function(e) NA)
    Specificity<-tryCatch(MLmetrics::Specificity(val_pred_bin, y_val), error=function(e) NA)
    Balanced_accuracy<-tryCatch((Sensitivity+Specificity)/2, error=function(e) NA)
    
    outer_results[[outer_idx]]<-data.frame(Fold=outer_idx, AUROC=AUROC, PRAUC=PRAUC,F1=F1, LogLoss=logloss,Sensitivity=Sensitivity, Specificity=Specificity,
      Balanced_accuracy=Balanced_accuracy)
  }
  
  outer_summary<-dplyr::bind_rows(outer_results)
  cat("\nNested CV summary (all metrics):\n")
  print(outer_summary)
  cat("\nMetric used for hyperparameter selection:",ers_metric,"(XGBoost eval_metric:",xgb_eval_metric,")\n")
  
  ## Selecting best outer fold based on the chosen metric
  best_idx<-switch(ers_metric,
                   "AUROC"=which.max(outer_summary$AUROC),
                   "F1"=which.max(outer_summary$F1),
                   "LogLoss"=which.min(outer_summary$LogLoss))
  final_params<-best_params_all[[best_idx]]
  
  ## Retraining final model on full training set with best hyperparameters
  final_params$scale_pos_weight<-spw
  dtrain_full<-xgboost::xgb.DMatrix(data=X_train, label=y_train)
  final_model<-xgboost::xgb.train(params=final_params,data=dtrain_full,nrounds=100,verbose=0)
  
  ## Finding optimal threshold on full training set predictions
  train_preds_for_thr<-predict(final_model, xgboost::xgb.DMatrix(X_train))
  final_threshold<-find_optimal_threshold(train_preds_for_thr, y_train, ers_metric)

  ## Computing OOF predictions on the full dataset (train + test)
  # - Train: predicted by per-fold models (each observation predicted by model not trained on it)
  # - Test: predicted by final model (trained on train only, no leakage)
  oof_preds_full<-rep(NA, nrow(dataset_dt))
  
  for(outer_idx in seq_along(outer_folds)){
    if(is.null(best_params_all[[outer_idx]])) next
    val_idx_global<-data_split$train[outer_folds[[outer_idx]]]
    train_idx_cv<-setdiff(seq_len(length(data_split$train)), outer_folds[[outer_idx]])
    X_tr_oof<-X_train[train_idx_cv,, drop=FALSE]
    y_tr_oof<-y_train[train_idx_cv]
    spw_oof <-max(table(y_tr_oof))/min(table(y_tr_oof))
    
    ## Updating scale_pos_weight for this OOF fold while keeping all other best params
    params_oof<-best_params_all[[outer_idx]]
    params_oof$scale_pos_weight<-spw_oof
    
    dtrain_oof<-xgboost::xgb.DMatrix(data=X_tr_oof, label=y_tr_oof)
    mod_oof<-xgboost::xgb.train(params=params_oof, data=dtrain_oof, nrounds=100, verbose=0)
    oof_preds_full[val_idx_global]<-predict(mod_oof,xgboost::xgb.DMatrix(X_train[outer_folds[[outer_idx]],, drop=FALSE]))
  }
  
  ## For test observations: use final model (no leakage since it was trained on train only)
  oof_preds_full[data_split$test]<-predict(final_model, xgboost::xgb.DMatrix(X_test))
  
  ## Test set performance (all metrics computed regardless of ers_metric)
  test_preds<-predict(final_model, xgboost::xgb.DMatrix(X_test))
  test_preds_bin<-as.integer(test_preds >= final_threshold)
  
  test_preds_fac<-factor(test_preds_bin,levels=c(0,1))
  y_test_fac<-factor(y_test,levels=c(0,1))
  
  ## Computing performance metrics
  roc_curve_test<-pROC::roc(y_test, test_preds, quiet=TRUE)
  AUROC_test<-as.numeric(pROC::auc(roc_curve_test))
  AUC_CI<-signif(pROC::ci.auc(roc_curve_test), 3)
  AUC_CI_str<-paste(AUC_CI[2]," (",AUC_CI[1],"; ",AUC_CI[3],")", sep="")
  F1_test<-tryCatch(MLmetrics::F1_Score(test_preds_bin, y_test), error=function(e) NA)
  logloss_test<-tryCatch(MLmetrics::LogLoss(test_preds, y_test), error=function(e) NA)
  PRAUC_test<-tryCatch(MLmetrics::PRAUC(test_preds, y_test),error=function(e) NA)
  conf_mat<-caret::confusionMatrix(factor(test_preds_bin), factor(y_test), mode="everything", positive="1")
  Sensitivity_test<-conf_mat$byClass[["Sensitivity"]]
  Specificity_test<-conf_mat$byClass[["Specificity"]]
  Balanced_accuracy_test<-(Sensitivity_test+Specificity_test)/2
  PPV_test<-conf_mat$byClass[["Pos Pred Value"]]
  NPV_test<-conf_mat$byClass[["Neg Pred Value"]]
  
  cat("\nTest set performance (all metrics):\n")
  cat("AUROC:",round(AUROC_test,3),"[",AUC_CI_str,"]\n")
  cat("PRAUC:",round(PRAUC_test,3)," F1:",round(F1_test,3)," LogLoss:",round(logloss_test,3),"\n")
  cat("Sensitivity:",round(Sensitivity_test,3)," Specificity:",round(Specificity_test,3)," Balanced accuracy:",round(Balanced_accuracy_test,3),"\n")
  cat("PPV:",round(PPV_test,3)," NPV:",round(NPV_test,3),"\n")
  cat("Metric used for model selection:",ers_metric,"\n")
  
  metric_test<-as_tibble(data.frame(metric_used_for_selection=ers_metric,AUROC=AUROC_test, AUROC_CI=AUC_CI_str,PRAUC=PRAUC_test, F1=F1_test,
                                    LogLoss=logloss_test,Sensitivity=Sensitivity_test, Specificity=Specificity_test,Balanced_accuracy=Balanced_accuracy_test,
                                    PPV=PPV_test, NPV=NPV_test))
  
  return(list(oof_preds=oof_preds_full,
              test_preds=test_preds,
              final_model=final_model,
              outer_summary=outer_summary,
              final_params=final_params,
              data_split=data_split,
              metric_test=metric_test,
              X_test=X_test,
              y_test=y_test))
}

#----------------------------------------------------------------
#### Settings ####

## Set to TRUE if data comes from multiple cohorts (activates Step 1)
use_cohort<-TRUE

## Names of covariates to adjust for in Step 2
cov_names<-c("age","sex","bmi")

## Metric used for ERS model (Step 3) hyperparameter tuning and best fold selection. Choose one of:
# - "AUROC": area under the ROC curve -> recommended for balanced or moderately imbalanced datasets
# - "F1": harmonic mean of precision/recall -> recommended when class imbalance is a concern
# - "LogLoss": log-loss / cross-entropy -> recommended when predicted probability calibration matters

ers_metric<-"AUROC"

#----------------------------------------------------------------
#### Data import and preprocessing (replace simulation with your own data) ####

## Simulating data
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
expo_mat<-as.matrix(cbind(log1p(expo_df[,paste0("X",1:K)]), # log-transform continuous only
  expo_df[,"X9",drop=FALSE])) # binary untransformed
colnames(expo_mat)<-expo_names

## Covariates and cohort
cov_df<-data.frame(age=rnorm(n,50,10),sex=rbinom(n,1,0.5),bmi=rnorm(n,26,4))
cohort<-sample(c("cohort_A","cohort_B","cohort_C"),n,replace=TRUE)

## Cohort effect on the log-odds scale
cohort_eff<-ifelse(cohort=="cohort_A",0,ifelse(cohort=="cohort_B",1.5,-1))

## True exposure effect (i.e., what ERS should capture after removing cohort + covariate effects)
h_z<-as.numeric(as.matrix(expo_mat)%*%c(0.5,0.3,-0.2,0.1,0.4,-0.1,0.05,0.05,1))

## Computing the linear predictor without intercept
log_odds_no_intercept<-0.03*cov_df$age + 0.5*cov_df$sex + 0.1*cov_df$bmi + cohort_eff + h_z

## Calibrating intercept so that prevalence is around 50%
intercept_calibrated<- -median(log_odds_no_intercept)

## Binary outcome via logistic link
y<-rbinom(n, 1, plogis(intercept_calibrated + log_odds_no_intercept))

cat("Observed prevalence:", round(mean(y),3), "- Cases:", sum(y==1), "- Controls:", sum(y==0), "\n")

## Assembling full dataset
sim_data<-cbind(data.frame(y=y,cohort=cohort),cov_df,expo_mat)

## Names of binary exposures in your dataset (all others treated as continuous/ordinal)
expo_binary<-c("X9") # adapt to your data

#----------------------------------------------------------------
#### Step 1 (optional): adjusting exposures for cohort membership (only if multi-cohort) ####

# Goal: to remove between-cohort differences from each exposure so that the ERS is not driven by which cohort a participant belongs to.
# - Residualization is applied to exposures (not Y) to preserve the binary nature of Y.
# - DML framework: W_j = X_j - E[X_j | cohort].
# - All exposures use binary:logistic after [0,1] rescaling (see fit_adjustment_oof).

# Method:
# - for each exposure X_j, fit XGBoost (binary:logistic on rescaled X_j) predicting X_j from cohort.
# - W_j = X_j - predicted(cohort): exposure cleaned of cohort effects.
# - OOF predictions used for all observations to avoid overfitting residuals.

# Output:
# - expo_mat_W: matrix of exposure residuals after cohort adjustment
# - All binary metrics (exposure ~ Y) before and after adjustment

if(use_cohort){
  
  ## One-hot encoding of cohort membership (no intercept to avoid collinearity)
  Z_mat<-model.matrix(~as.factor(sim_data$cohort)-1)
  colnames(Z_mat)<-paste0("cohort_",levels(as.factor(sim_data$cohort)))
  
  ## Residualizing each exposure on cohort. W_j = X_j - E[X_j | cohort]
  cat("\nStep 1 - adjusting all exposures for cohort\n")

  expo_hat_cohort_list<-lapply(expo_names, function(j){
    is_bin<-j %in% expo_binary
    cat("Adjusting:",j, if(is_bin) "[binary]" else "[continuous]","\n")
    fit_adjustment_oof(X=Z_mat, y=expo_mat[,j],is_binary_exposure=is_bin,nb_outer_fold=5, nb_inner_fold=5)
  })
  
  expo_hat_cohort<-do.call(cbind,expo_hat_cohort_list)
  colnames(expo_hat_cohort)<-expo_names
  
  ## W = X - E[X|cohort]: exposures cleaned of cohort effects
  expo_mat_W<-expo_mat-expo_hat_cohort
  
  ## All binary metrics (exposure ~ Y) before and after cohort adjustment
  metrics_before_cohort<-do.call(rbind,lapply(expo_names,function(j)
    cbind(data.frame(exposure=j,adjustment="before"),compute_binary_metrics(expo_mat[,j],y))))
  
  metrics_after_cohort<-do.call(rbind,lapply(expo_names,function(j)
    cbind(data.frame(exposure=j,adjustment="after"),compute_binary_metrics(expo_mat_W[,j],y))))
  
  metrics_cohort<-bind_rows(metrics_before_cohort,metrics_after_cohort)
  
  cat("\nAll binary metrics (exposure ~ Y) before and after cohort adjustment:\n")
  metrics_cohort %>% pivot_wider(names_from=adjustment, values_from=c(AUROC,PRAUC,F1,LogLoss,Sensitivity,Specificity,Balanced_accuracy))
  
}else{
  expo_mat_W<-expo_mat
}

#----------------------------------------------------------------
#### Step 2: adjusting exposures for covariates ####

# Second DML stage: W2_j = W_j - E[W_j | covariates].
# Goal: to remove the effects of age, sex, BMI or any other covariates from each exposure so that the ERS captures only the exposure-specific signal independent of classical covariates.

cov_mat<-as.matrix(sim_data[,cov_names])

## Adjusting each exposure residual for covariates. W2_j = W_j - E[W_j | covariates]
cat("\nStep 2 - adjusting all exposures for covariates\n")

expo_hat_cov_list<-lapply(expo_names,function(j){
  is_bin<-j %in% expo_binary
  cat("Adjusting:",j, if(is_bin) "[binary]" else "[continuous]","\n")
  fit_adjustment_oof(X=cov_mat,y=expo_mat_W[,j],is_binary_exposure=is_bin,nb_outer_fold=5,nb_inner_fold=5)
})

expo_hat_cov<-do.call(cbind,expo_hat_cov_list)
colnames(expo_hat_cov)<-expo_names

## W2 = W - E[W|covariates]: exposures cleaned of both cohort and covariate effects
expo_mat_W2<-expo_mat_W-expo_hat_cov

## All binary metrics (exposure ~ Y) before and after covariate adjustment
metrics_before_cov<-do.call(rbind,lapply(expo_names,function(j)
  cbind(data.frame(exposure=j,adjustment="before"),compute_binary_metrics(expo_mat_W[,j],y))))

metrics_after_cov<-do.call(rbind,lapply(expo_names,function(j)
  cbind(data.frame(exposure=j,adjustment="after"),compute_binary_metrics(expo_mat_W2[,j],y))))

metrics_cov<-bind_rows(metrics_before_cov,metrics_after_cov)

cat("\nAll binary metrics (exposure ~ Y) before and after covariate adjustment:\n")
metrics_cov %>% pivot_wider(names_from=adjustment, values_from=c(AUROC,PRAUC,F1,LogLoss,Sensitivity,Specificity,Balanced_accuracy))

#----------------------------------------------------------------
#### Step 3: Fitting ERS model on doubly-adjusted exposures (nested CV) ####

# Goal: to model the binary outcome Y using doubly-adjusted exposures W2 only. Because exposures have been cleaned of cohort and covariate effects (Steps 1-2),
#       the predicted log-odds f(W2) represent the exposure-driven component of Y (= ERS), cleanly separated from cohort and covariate effects.

colnames(expo_mat_W2)<-expo_names
dataset_step3<-as.data.table(cbind(data.frame(y=y),as.data.frame(expo_mat_W2)))

step3_result<-fit_ers_binary(dataset_dt=dataset_step3,
                             target_col="y",
                             feature_cols=expo_names,
                             ers_metric=ers_metric,
                             train_split=0.7,
                             test_split=0.3,
                             nb_outer_fold=5,
                             nb_inner_fold=5)

## ERS = predicted log-odds from doubly-adjusted exposures on the full dataset
ers<-predict(step3_result$final_model,xgboost::xgb.DMatrix(expo_mat_W2))
sim_data$ers<-ers

## Evaluating ERS against the original binary outcome
ers_metrics_full<-compute_binary_metrics(ers,y)
roc_ers<-pROC::roc(y,ers,quiet=TRUE)
AUROC_ers<-round(as.numeric(pROC::auc(roc_ers)),3)
AUC_CI_ers<-signif(pROC::ci.auc(roc_ers),3)

cat("\nERS performance against binary outcome Y (full dataset):\n")
cat("AUROC:",AUROC_ers,"[",AUC_CI_ers[1],";",AUC_CI_ers[3],"]\n")
cat("PRAUC:",ers_metrics_full$PRAUC," F1:",ers_metrics_full$F1," LogLoss:",ers_metrics_full$LogLoss,"\n")
cat("Sensitivity:",ers_metrics_full$Sensitivity," Specificity:",ers_metrics_full$Specificity,
    "Balanced accuracy:",ers_metrics_full$Balanced_accuracy,"\n")

## Saving ERS model
saveRDS(step3_result$final_model,paste0("ERS_model_binary_",Sys.Date(),".rds"))
xgb.save(step3_result$final_model,paste0("ERS_model_binary_",Sys.Date(),".ubj"))

#----------------------------------------------------------------
#### Step 4: Overall ERS direction (mixture effect) ####

# Goal: to determine whether the overall exposome (ERS) promotes or mitigates the binary outcome, using a 3-method consensus framework

## Building a two-column data frame: ERS value vs. binary outcome Y
ers_direction_df<-data.frame(ers=ers,y=as.numeric(y))

## Approach 1: sign of mean association
ers_mean_sign<-conventional_direction(ers_direction_df,feature_col="ers",shap_col="y")

## Approach 2: GAM derivative + Spearman fallback
ers_gam<-gam_direction(ers_direction_df,feature_col="ers",shap_col="y")

## Approach 3: pairwise bin concordance
ers_pairwise<-pairwise_direction(ers_direction_df,feature_col="ers",shap_col="y",majority_threshold=0.55,n_quantile_bins=200)

## Consensus: majority vote across the three approaches
ers_consensus_vec<-c(ers_mean_sign,ers_gam,ers_pairwise)
ers_tbl<-sort(table(ers_consensus_vec),decreasing=TRUE)
ers_consensus<-if(ers_tbl[1]>=2) names(ers_tbl)[1] else "uncertain"

## Summary table: direction + all binary metrics
ers_direction_summary<-tibble(mean_sign_direction=ers_mean_sign,
  gam_direction=ers_gam,
  pairwise_direction=ers_pairwise,
  consensus_direction=ers_consensus,
  AUROC=AUROC_ers,
  AUROC_CI=paste(AUC_CI_ers[1],";",AUC_CI_ers[3]),
  PRAUC=ers_metrics_full$PRAUC,
  F1=ers_metrics_full$F1,
  LogLoss=ers_metrics_full$LogLoss,
  Sensitivity=ers_metrics_full$Sensitivity,
  Specificity=ers_metrics_full$Specificity,
  Balanced_accuracy=ers_metrics_full$Balanced_accuracy)

cat("\nERS global direction and performance summary:\n")
ers_direction_summary

#----------------------------------------------------------------
#### Step 5: Using ERS in a final association model ####

# The ERS can be used as a predictor in standard models alongside covariates.
# Since W2 (doubly-adjusted exposures) were used to build the ERS, the ERS is by construction largely orthogonal to cohort and covariates.

## Logistic regression: ERS + covariates
final_model_glm<-glm(y~ers+age+sex+bmi,data=sim_data,family=binomial())
summary(final_model_glm)
broom::tidy(final_model_glm,exponentiate=TRUE,conf.int=TRUE)

## Stratification into tertiles
sim_data$ers_group<-ntile(sim_data$ers,3)
sim_data$ers_group<-factor(sim_data$ers_group,labels=c("Low","Medium","High"))

ers_tertile_model<-glm(y~ers_group+age+sex+bmi,data=sim_data,family=binomial())
broom::tidy(ers_tertile_model,exponentiate=TRUE,conf.int=TRUE)

#----------------------------------------------------------------
#### Step 6: Visualization ####

## Figure 1: ERS distribution by outcome group
ggplot(sim_data,aes(x=ers,fill=factor(y)))+
  geom_histogram(position="identity",alpha=0.6,bins=40,color="black")+
  scale_fill_manual("Outcome",values=c("0"="#A6DDCE","1"="#F9CBC2"),labels=c("Control","Case"))+
  labs(title="Distribution of ERS by outcome",x="ERS (log-odds scale)",y="Count")+
  theme_Gaia()

## Figure 2: ROC curve of ERS vs. binary outcome
roc_df<-data.frame(specificity=rev(roc_ers$specificities),sensitivity=rev(roc_ers$sensitivities))
ggplot(roc_df,aes(x=1-specificity,y=sensitivity))+
  geom_line(color="#d6604d",linewidth=1.5)+
  geom_abline(linetype="dashed",color="grey50")+
  labs(title=paste0("ROC curve - ERS vs. outcome  (AUROC=",AUROC_ers,")"),x="1 - Specificity",y="Sensitivity")+
  theme_Gaia()

## Figure 3: ERS by cohort
if(use_cohort){
  ggplot(sim_data,aes(x=cohort,y=ers,fill=cohort))+
    geom_violin(fill="transparent")+
    geom_boxplot(alpha=0.7,width=0.4)+
    scale_fill_brewer(palette="Set2")+
    labs(x="",y="ERS")+
    theme_Gaia()+
    theme(legend.position="none")
}

## Figure 4: ERS by outcome group within each cohort
if(use_cohort){
  ggplot(sim_data,aes(x=factor(y),y=ers,fill=factor(y)))+
    geom_boxplot(alpha=0.7)+
    scale_fill_manual(values=c("0"="#A6DDCE","1"="#F9CBC2"))+
    facet_wrap(~cohort)+
    labs(title="ERS by outcome group within each cohort",x="Outcome (0=control, 1=case)",y="ERS")+
    theme_Gaia()+
    theme(legend.position="none")
}else{
  ggplot(sim_data, aes(x=factor(y), y=ers, fill=factor(y)))+
    geom_boxplot(alpha=0.7)+
    scale_fill_manual(values=c("0"="#A6DDCE","1"="#F9CBC2"))+
    labs(title="ERS by outcome group",
         x="Outcome (0=control, 1=case)", y="ERS")+
    theme_Gaia()+
    theme(legend.position="none")
}

#----------------------------------------------------------------
#### Step 7: Univariate logistic regression on doubly-adjusted exposures (test set) ####

Association_tab_res<-c()
X_test_ers<-step3_result$X_test
y_test_ers<-step3_result$y_test

for(n_var in seq_along(expo_names)){
  
  test_tmpo<-as_tibble(X_test_ers) %>% select(all_of(expo_names[n_var]))
  colnames(test_tmpo)<-"tmp"
  test_tmpo<-bind_cols(data.frame(outcome=y_test_ers),test_tmpo)
  
  glm_model<-glm(outcome~tmp,data=test_tmpo,family="binomial")
  mod_res<-broom::tidy(glm_model,exponentiate=TRUE,conf.int=TRUE) %>%
    filter(term!="(Intercept)") %>%
    mutate(term=expo_names[n_var],
           direction=case_when(conf.low>=1~"positive association",
                               conf.high<1~"negative association",
                               TRUE~"uncertain association")) %>%
    bind_cols(performance::performance(glm_model)) %>%
    rowwise() %>% mutate_if(is.numeric,test_format) %>% ungroup() %>%
    mutate(OR=paste(estimate," [",conf.low,"; ",conf.high,"]",sep="")) %>%
    mutate(OR=if_else(OR==" [; ]","",OR))
  
  Association_tab_res<-bind_rows(Association_tab_res,mod_res)
}

#----------------------------------------------------------------
#### Step 8: SHAP analysis ####

## SHAP values are computed on the test set of Step 3. They show which doubly-adjusted exposures contribute most to the ERS, and in which direction:
# - promoting  = increases predicted log-odds = risk-increasing effect on Y
# - mitigating = decreases predicted log-odds = protective effect on Y

## Computing SHAP values
contr<-predict(step3_result$final_model,as.matrix(X_test_ers),predcontrib=TRUE)
shap<-as_tibble(contr)
shap_contrib<-as.data.table(contr)

## Removing BIAS term
shap_contrib<-shap_contrib[,!grepl("bias|Intercept|BIAS|Bias",names(shap_contrib),ignore.case=TRUE),with=FALSE]

## Computing mean absolute SHAP score per feature (= importance ranking)
mean_shap_score<-colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)),decreasing=TRUE)]

## Reshaping SHAP values to long format
shap_score_sub<-as.data.table(shap_contrib)[,names(mean_shap_score),with=FALSE]
shap_score_long<-melt.data.table(shap_score_sub,measure.vars=colnames(shap_score_sub))

## Matching standardized feature values for color scale
fv_sub<-as.data.table(X_test_ers)[,names(mean_shap_score),with=FALSE]
fv_sub_long<-melt.data.table(fv_sub,measure.vars=colnames(fv_sub))
fv_sub_long[,stdfvalue:=std1(value),by="variable"]
names(fv_sub_long)<-c("variable","rfvalue","stdfvalue")

## Merging SHAP and feature values
shap_long2<-cbind(shap_score_long,fv_sub_long[,c("rfvalue","stdfvalue")])
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
  mutate(mean_val=test_format(mean_val),IC_2.5=test_format(IC_2.5),IC_97.5=test_format(IC_97.5)) %>%
  ungroup() %>%
  mutate(mean_SHAP=paste(mean_val," (",IC_2.5,"; ",IC_97.5,")",sep=""),
         mean_SHAP=if_else(mean_SHAP%in%c("0 (0; 0)","0e+00 (0e+00; 0e+00)"),"0",mean_SHAP))

## Pivoting to long format for direction analysis
shap_tempo<-shap %>%
  rowid_to_column("id") %>%
  pivot_longer(-id,names_to="feature",values_to="shap_value")

tempopo<-X_test_ers %>%
  as_tibble() %>%
  rowid_to_column("id") %>%
  pivot_longer(-id,names_to="feature",values_to="feature_value")

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
test<-left_join(tempopo,direction_impact,by="feature") %>%  mutate(direction=consensus)
test<-left_join(test,Shap_val,by=base::intersect(colnames(test),colnames(Shap_val)))

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
ordre_def<-bind_rows(test_tmp %>% filter(mean_val>=0) %>% arrange(desc(IC_2.5)),
  test_tmp %>% filter(mean_val< 0) %>% arrange(desc(IC_97.5)))

test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def$feature)))
test_tmp$impact<-factor(test_tmp$impact,
                        levels=c("mitigating","neutral","promoting","uncertain","undefined"),
                        labels=c("mitigating predictor","neutral","promoting predictor","uncertain","undefined"))

#- - - -
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

#- - - -
## Plot 2: barplot (mean |SHAP| with labels)

ordre_def2<-test_tmp %>% arrange(desc(mean_val))
test_tmp$feature<-factor(test_tmp$feature,levels=rev(unique(ordre_def2$feature)))

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

#- - - -
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
  scale_x_continuous("% contribution to ERS (based on mean |SHAP|)", label=function(x) paste0(abs(x)))+
  scale_y_discrete("")+
  scale_fill_manual("Direction:",na.value="white",
                    values=c("mitigating predictor"="#A6DDCE","promoting predictor"="#F9CBC2","neutral"="white","uncertain"="grey","undefined"="grey50"))+
  theme_Gaia()+
  theme(legend.position="bottom")

ggsave(plot3,file=paste0("ERS_SHAP_contribution_",Sys.Date(),".pdf"),dpi=600,width=60,height=30,units="cm",limitsize=FALSE)

#- - - -
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

#- - - -
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

#- - - -
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
#### SHAP and logistic regression comparison table ####

tmp_GLM<-Association_tab_res %>% select(term,OR,direction,p.value,Log_loss)
colnames(tmp_GLM)<-c("feature","OR [95% CI]","GLM direction","GLM p-value","GLM log-loss")

test_SHAP_tmp<-test %>%
  select(feature,mean_SHAP,conventional,gam,pairwise,consensus,mean_val) %>%
  distinct() %>%
  arrange(desc(abs(as.numeric(mean_val)))) %>%
  rowwise() %>% mutate_if(is.numeric,test_format) %>% ungroup() %>%
  select(-mean_val)

colnames(test_SHAP_tmp)<-c("feature","mean |SHAP| [95% CI]",
                           "Conventional direction","GAM-based direction",
                           "Pairwise-based direction","Consensus direction")

Compa_GLM_SHAP<-left_join(test_SHAP_tmp,tmp_GLM,by="feature")

#----------------------------------------------------------------
#### Saving results ####

fwrite(step3_result$outer_summary,paste0("ERS_CV_summary_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
fwrite(step3_result$metric_test,paste0("ERS_test_performance_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
fwrite(Shap_val,paste0("ERS_Shap_val_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
fwrite(test_tmp,paste0("ERS_test_tmp_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
fwrite(test,paste0("ERS_test_",Sys.Date(),".txt"),sep=";",row.names=FALSE)
fwrite(shap,paste0("ERS_shap_",Sys.Date(),".txt"),sep=";",row.names=FALSE)
fwrite(ers_direction_summary,paste0("ERS_overall_direction_summary_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
fwrite(Association_tab_res,paste0("ERS_Univariate_logistic_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
fwrite(interaction_df, paste0("ERS_SHAP interaction summary_",Sys.Date(),".csv"), sep=";",row.names=FALSE)
fwrite(Compa_GLM_SHAP,paste0("ERS_Compa_logistic_SHAP_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
fwrite(as_tibble(step3_result$final_params),paste0("ERS_best_params_",Sys.Date(),".csv"),sep=";",row.names=FALSE)

## Saving exposure-level adjustment metrics (before/after cohort and covariate adjustment)
if(use_cohort){
  fwrite(metrics_cohort %>% pivot_wider(names_from=adjustment,values_from=c(AUROC,PRAUC,F1,LogLoss,Sensitivity,Specificity,Balanced_accuracy)),
         paste0("ERS_exposure_metrics_cohort_adjustment_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
}
fwrite(metrics_cov %>% pivot_wider(names_from=adjustment,values_from=c(AUROC,PRAUC,F1,LogLoss,Sensitivity,Specificity,Balanced_accuracy)),
       paste0("ERS_exposure_metrics_covariate_adjustment_",Sys.Date(),".csv"),sep=";",row.names=FALSE)
