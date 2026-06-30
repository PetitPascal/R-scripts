#-------------------------------------------------------------------------------
## Reproducible and generalisable global SHAP direction script
#
# Covers:
#   - Approach 1: mean SHAP sign
#   - Approach 2: pairwise correlation  
#   - Approach 3: overall trend of the sign of SHAP conditional on feature based on parametric or non-parametric linear regression
#   - Approach 4: overall trend of the sign of SHAP conditional on feature based on the derivative sign of a smooth regression, evaluated on a clean synthetic grid
#   - Approach 5: Theil-Sen slope estimator
#   - Approach 6: number of zero-crossings
#   - Approach 7: pairwise bin concordance (probabilistic dominance) based on Kendall's tau (Rcpp)
#   - Approach 8: stochastic dominance
#   - Approach 9: signed-AUC dominance index
#   - Approach 10: consensus of all approaches (unweighted majority vote)
#   - Approach 11: consensus of all approaches (weighted majority vote)
#   - Approach 12: model-agnostic confidence measure based on bootstrap
#   - Approach 13: model-refit stability diagnostic

# Examples for SHAP and ALE values:
#   - Comparison across all approaches
#   - Visualization: importance plot, and dependence plots
#   - Sensitivity analyses: threshold and bin number sensitivity for approach 7; majority threshold, quantile threshold, and minimum observations per group sensitivity for approach 8; signed-AUC dominance threshold sensitivity (approach 9)
#   - SHAP interaction
#-------------------------------------------------------------------------------


#--------------------------------------------
#### Step 1: Setup #### 

rm(list=ls())
set.seed(123)

## Installing and loading packages

required_pkgs<-c("tidyverse", "mgcv", "gratia", "patchwork", "xgboost", "shapr", "iml", "Rcpp", "MASS", "stats", "splines", "np", "ale", "shapviz")

is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(required_pkgs[!is_installed],repos="http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only=TRUE))

## Preventing package conflicts
conflict_prefer("select", "dplyr") 
conflict_prefer("filter", "dplyr") 
conflict_prefer("slice", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("rename", "dplyr")

#--------------------------------------------
#### Step 2: Creating useful helpers ####

#- - - - - -
## Computing numerical derivative from predictions evaluated on a clean, evenly-spaced grid
#
# Parameters:
# - x: vector of feature values
# - yhat: vector of attribution/explanation/effect (e.g., SHAP, ALE) values

#
# Goal: to estimate the average slope of a fitted curve, guarding against the Inf values that arise when x contains duplicate values (diff(x)==0)
#       this is why this is always called on a synthetic grid (see make_grid) rather than on raw, possibly-duplicated x.
#
# Parameters:
# - x: vector of feature values (expected to come from make_grid, strictly increasing)
# - yhat: vector of fitted attribution/explanation/effect values (e.g., SHAP, ALE) evaluated at each x
#
# Returns: mean finite slope (numeric), or NA if no finite slopes exist

mean_derivative<-function(x, yhat){
  o<-order(x)
  x<-x[o]
  yhat<-yhat[o]
  d<-diff(yhat)/diff(x)
  d<-d[is.finite(d)]
  if(length(d)==0) return(NA_real_)
  mean(d, na.rm=TRUE)
}

#- - - - - -
## Building a clean, strictly increasing synthetic grid over the observed range of x
#
# Goal: used by every derivative/area-based approach to avoid duplicate-x issues and to give every smoother a common, 
#       comparable evaluation grid regardless of the original (possibly irregular or duplicated) distribution of observed feature values
#
# Parameters:
# - x: vector of feature values
# - n_grid: number of grid points (default: 200)
#
# Returns: a numeric vector of n_grid evenly-spaced values spanning range(x), or NULL if x is constant (range of 0)

make_grid<-function(x, n_grid=200){
  rng<-range(x, na.rm=TRUE)
  if(diff(rng)==0) return(NULL)
  seq(rng[1], rng[2], length.out=n_grid)
}

#- - - - - -
## Predict function for ALE computation with XGBoost model
#
# Goal: to wraps xgboost's predict() in the signature expected by ale::ALE(), converting the newdata data frame into the xgb.DMatrix format the model requires
#
# Parameters:
# - object: fitted xgboost model (passed by ale::ALE() internally)
# - newdata: data frame of feature values at which to predict (passed by ale::ALE())
# - ...: unused, present only to match the generic predict() signature ale::ALE() expects
#
# Returns: numeric vector of model predictions
#
# Note: relies on feature_cols existing in the calling environment (set in Step 5)

predict_xgb_ale<-function(object, newdata, ...){
  X_new<-as.matrix(newdata[, feature_cols, drop = FALSE])
  colnames(X_new)<-feature_cols
  predict(object, xgboost::xgb.DMatrix(X_new))
}

#- - - - - -
## Converting numeric score to direction label
#
# Parameters:
# - z: numeric score (e.g., mean slope, correlation)
# - eps: numerical tolerance threshold

direction_from_value<-function(z, eps=1e-10){
  if(is.na(z) || abs(z) < eps) return("neutral")
  if(z > 0) return("promoting")
  "mitigating"
}

#- - - - -
## Defining safe null operator

`%||%`<-function(a, b) if (!is.null(a)) a else b

#- - - - - -
## Defining not-in operator

`%ni%`<-Negate('%in%')

#- - - - - -
## Creating dependence plots
#
# Parameters:
# - x_col: vector of feature values
# - y_col: vector of SHAP values
# - df: dataframe in a long format, containing a column with feature names, a column with feature values and a column with attribution/explanation/effect values (e.g., SHAP, ALE) 
# - y_name: attribution/explanation/effect metric name
# - ncol: number of columns to use for display

plot_dep_all<-function(df, feature_col, x_col, y_col, y_name=NULL, ncol=3){
  
  plots<-lapply(split(df, df[[feature_col]]), function(d) {
    ggplot2::ggplot(d, ggplot2::aes(x=.data[[x_col]], y=.data[[y_col]])) +
      ggplot2::geom_point(alpha=0.3, color="#2166ac") +
      ggplot2::geom_smooth(method="loess", se=TRUE, color="#d6604d", linewidth=1.1) +
      ggplot2::geom_hline(yintercept=0, linetype="dashed", color="grey50") +
      ggplot2::labs(x=unique(d[[feature_col]]),y=if(is.null(y_name)) y_col else y_name) +
      theme_Gaia()
  })
  
  patchwork::wrap_plots(plots, ncol=ncol)
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
## SHAP interaction summary
#
# Goal: to summarize SHAP interaction values
#
# Parameters:
# - xgb_fit: XGBoost model
# - X_mat: matrix of feature values
# - feature_cols: vector of feature names
#
# Returns:
# - feature: feature name
# - interaction_strength: mean absolute off-diagonal SHAP interaction -> high values indicate this feature interacts strongly with others
# - main_effect_mean: mean absolute diagonal (main effect) SHAP -> to compare interaction vs. main effect magnitude

shap_interaction_summary<-function(xgb_fit, X_mat, feature_cols){
  inter<-predict(xgb_fit, X_mat, predinteraction = TRUE)
  p<-length(feature_cols)
  
  arr<-if(length(dim(inter))==3){
    inter
  }else{
    array(inter, dim = c(nrow(X_mat), p + 1, p + 1))
  }
  
  dplyr::bind_rows(lapply(seq_len(p), function(j){
    off<-arr[, j, 1:p, drop = FALSE]
    if(length(dim(off))==3) off<-off[, 1, , drop = FALSE]
    main<-arr[, j, j]
    offdiag<-rowSums(abs(off), na.rm = TRUE) - abs(main)
    
    data.frame(feature=feature_cols[j],
               interaction_strength=mean(offdiag, na.rm = TRUE),
               main_effect_mean=mean(main, na.rm = TRUE))
  }))
}

#- - - - - -
## Checking cardinality
#
# Goal: to decide whether a feature is "low-cardinality" (binary/ordinal-sparse), in which case fitting splines/LOESS is meaningless and a direct two-group contrast is used
#
# Parameters:
# - x: vector of feature values
# - threshold: cutoff value (default: 5), below which a smooth curve is considered to have no interpretable shape to estimate

is_low_cardinality<-function(x, threshold=5){
  length(unique(x))<=threshold
}

#- - - - - -
## Adaptive basis dimension for any mgcv smooth term
#
# Goal: to avoid errors when using mgcv fuction by capping k at n_unique-1 and additionally capping it at the user-requested default
#
# Parameters:
# - x: vector of feature values
# - default_k: dimension for the smooth term (default: 10)

safe_k<-function(x, default_k=10){
  nun<-length(unique(x))
  max(3, min(default_k, nun-1))
}

#- - - - - -
## Direct two-group contrast for low-cardinality features (binary or few-level ordinal)
#
# Goal: to compare mean attribution/explanation/effect (e.g., SHAP, ALE) value at the highest level vs. lowest level
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - alpha: alpha risk threshold (default: 0.05)

two_group_contrast<-function(x, y, alpha=0.05){
  ok<-complete.cases(x,y); 
  x<-x[ok]
  y<-y[ok]
  
  lv<-min(x, na.rm=TRUE)
  hv<-max(x, na.rm=TRUE)
  g0<-y[x==lv]
  g1<-y[x==hv]
  
  if(length(g0)<2 || length(g1)<2) return(list(dir="undefined", p=NA))
  if(var(g0, na.rm=TRUE)==0 && var(g1, na.rm=TRUE)==0){
    if(mean(g1)==mean(g0)) return(list(dir="neutral", p=NA))
    return(list(dir=if_else(mean(g1)>=mean(g0), "promoting", "mitigating"), p=NA))
  }
  tt<-tryCatch(suppressWarnings(t.test(g1, g0)), error=function(e) NULL)
  if(is.null(tt) || is.na(tt$p.value) || tt$p.value > alpha) return(list(dir="neutral", p=tt$p.value %||% NA))
  list(dir=if_else(mean(g1,na.rm=TRUE) >= mean(g0,na.rm=TRUE), "promoting", "mitigating"), p=tt$p.value)
}

#------------------------------------------------------------------------------
#### Step 3: Creating helpers to determine the feature's direction based on SHAP values ####

#- - - - - -
## Approach 1: mean SHAP sign
#
# Parameters:
# - shap_vec: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - eps: numerical tolerance (default: 1e-10)
#
# Returns: "promoting", "mitigating", or "neutral"
#
# Caveat: this collapses to "neutral" whenever positive and negative contributions  cancel out on average

mean_sign_direction<- function(shap_vec, eps=1e-10){
  direction_from_value(mean(shap_vec, na.rm=TRUE), eps=eps)
}

#- - - - - -
## Approach 2: pairwise correlation
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - alpha: alpha risk threshold (default: 0.05)
# - cor.method: "pearson", "spearman" (default), or "kendall"
#
# Returns: "promoting", "mitigating", "neutral", or "undefined"

cor_dir<-function(x, y, alpha=0.05, cor.method="spearman"){
  ok<-complete.cases(x, y)
  x<-x[ok]
  y<-y[ok]
  if(length(unique(x)) <= 2 || length(unique(y)) <= 2) return("undefined")
  
  ct<-suppressWarnings(cor.test(x, y, method=cor.method))
  
  if(is.na(ct$p.value) || ct$p.value > alpha){
    return("neutral")
  }
  if_else(ct$estimate >= 0, "promoting", "mitigating")
}

#- - - - - -
## Approach 3: overall trend of the sign of SHAP conditional on feature based on parametric or non-parametric linear regression
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - alpha: alpha risk threshold (default: 0.05)
# - glm.type: "LR", "RLR", "spline", "GAM_l", "GAM_s" (default), "GAM_p", "GAM_tp", "GAM_cr", "GAM_ds", "GAM_bs", "GAM_ad", "GAM_gp", "GAM_re", "GAM_sz", "GAM_fs"
#
# Returns: "promoting", "mitigating", "neutral", or "undefined"

GLM_trend<-function(x, y, alpha=0.05, glm.type="GAM_s"){
  ok<-complete.cases(x, y)
  x<-x[ok]
  y<-y[ok]
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return("undefined")
  
  if(is_low_cardinality(x)){
    res<-two_group_contrast(x, y, alpha)
    return(res$dir)
  }
  
  grid<-make_grid(x)
  if(is.null(grid)) return("undefined")
  beta_glm<-NA_real_
  p_glm<-NA_real_
  k_use<-safe_k(x)
  
  tryCatch({
    if(glm.type=="LR"){
      
      fit<-suppressWarnings(glm(y ~ x, family="gaussian"))
      beta_glm<-summary(fit)$coeff[2,1]
      p_glm<-summary(fit)$coeff[2,4]
      
    }else if(glm.type=="RLR"){
      
      fit<-suppressWarnings(MASS::rlm(y ~ x, method="MM"))
      sm<-summary(fit)
      beta_glm<-sm$coefficients[2,1]
      se_glm<-sm$coefficients[2,2]
      p_glm<-2*pnorm(-abs(beta_glm/se_glm))
      
    }else if(glm.type=="spline"){
      
      fit<-suppressWarnings(lm(y ~ splines::ns(x, df=min(3, length(unique(x))-1))))
      beta_glm<-summary(fit)$coeff[2,1]
      p_glm<-summary(fit)$coeff[2,4]
      
    }else if(glm.type=="GAM_l"){
      
      fit<-suppressWarnings(gam(y ~ x, method="REML"))
      beta_glm<-summary(fit)$p.table[2,1]
      p_glm<-summary(fit)$p.table[2,4]
      
    }else{
      
      bs_map<-c(GAM_s="tp", GAM_p="ps", GAM_tp="tp", GAM_cr="cr", GAM_ds="ds", GAM_bs="bs", GAM_re="re", GAM_gp="gp", GAM_ad="ad", GAM_sz="sz", GAM_fs="fs")
      bs_choice<-bs_map[[glm.type]] %||% "tp"
      fit<-suppressWarnings(gam(y ~ s(x, bs=bs_choice, k=k_use), method="REML"))
      p_glm<-summary(fit)$s.table[4]
      pred_on_grid<-as.numeric(predict(fit, newdata=data.frame(x=grid)))
      beta_glm<-mean_derivative(grid, pred_on_grid)
    }
  }, error=function(e) NULL)
  
  if(is.na(p_glm) || is.na(beta_glm) || p_glm > alpha) return("neutral")
  if_else(beta_glm >= 0, "promoting", "mitigating")
}

#- - - - - -
## Approach 4: overall trend of the sign of SHAP conditional on feature based on the derivative sign of a smooth regression, evaluated on a clean synthetic grid
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - method.type: "kernel", "LOESS" (default), "supsmu" (Friedman's SuperSmoother), or "GAM"
#
# Returns: "promoting", "mitigating", "neutral", "uncertain", or "undefined"

reg_deriv<-function(x, y, method.type="LOESS"){
  ok<-complete.cases(x, y)
  x<-x[ok]
  y<-y[ok]
  
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return("undefined")
  
  if(is_low_cardinality(x)){
    res<-two_group_contrast(x, y)
    return(res$dir)
  }
  
  grid<-make_grid(x)
  if(is.null(grid)) return("undefined")
  sign_deriv<-NA_real_
  k_use<-safe_k(x)
  
  sign_deriv<-tryCatch({
    if(method.type=="kernel"){
      
      fit<-suppressWarnings(npreg(y ~ x))
      mean_derivative(grid, predict(fit, newdata=data.frame(x=grid)))
      
    }else if(method.type=="LOESS"){
      
      fit<-suppressWarnings(loess(y ~ x, span=0.75, degree=2))
      mean_derivative(grid, predict(fit, newdata=data.frame(x=grid)))
      
    }else if(method.type=="supsmu"){
      
      o<-order(x)
      sm<-suppressWarnings(supsmu(x[o], y[o]))
      mean_derivative(sm$x, sm$y)
      
    }else if(method.type=="GAM"){
      
      fit<-suppressWarnings(gam(y ~ s(x, k=k_use), method="REML"))
      mean_derivative(grid, predict(fit, newdata=data.frame(x=grid)))
      
    }else NA_real_
  }, error=function(e) NA_real_)
  
  if(is.na(sign_deriv)) return("uncertain")
  if(sign_deriv==0) return("neutral")
  if_else(sign_deriv > 0, "promoting", "mitigating")
}

#- - - - - -
## Approach 5: Theil-Sen robust slope 
#
# Goal: to estimate the median pairwise slope using a robust non-parametric method insensitive to outliers
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - alpha: alpha risk threshold for the associated Kendall's tau test (default: 0.05)
#
# Returns: "promoting", "mitigating", "neutral", "uncertain", or "undefined"

theil_sen_direction<-function(x, y, alpha=0.05){
  ok<-complete.cases(x, y)
  x<-x[ok]
  y<-y[ok]
  
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return("undefined")
  if(is_low_cardinality(x)) return(two_group_contrast(x, y, alpha)$dir)
  
  n<-length(x)
  max_pairs<-200000
  
  if(choose(n,2) > max_pairs){
    idx<-sample(seq_len(n), size=min(n, 2000))
    x_s<-x[idx]
    y_s<-y[idx]
  }else{
    x_s<-x
    y_s<-y
  }
  
  pairs<-expand.grid(i=seq_along(x_s), j=seq_along(x_s))
  pairs<-pairs[pairs$i<pairs$j,]
  dx<-x_s[pairs$j]-x_s[pairs$i]
  dy<-y_s[pairs$j]-y_s[pairs$i]
  valid<-dx!=0
  slopes<-dy[valid]/dx[valid]
  if(length(slopes)==0) return("uncertain")
  med_slope<-median(slopes, na.rm=TRUE)
  
  kt<-suppressWarnings(cor.test(x, y, method="kendall"))
  if(is.na(kt$p.value) || kt$p.value > alpha) return("neutral")
  if_else(med_slope >= 0, "promoting", "mitigating")
}

#- - - - - -
## Approach 6: number of zero-crossings (area-weighted)
#
# Goal: to determine whether a feature's smoothed attribution/explanation/effect metric-vs.-feature curve has a dominant sign once integrated
#
# Method: the curve is fitted with a chosen smoother evaluated on a synthetic grid, split into segments at each sign change, and each segment is weighted by its trapezoidal area before comparing total positive vs. total negative area
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - glm.type: regression type used to obtain the smooth fitted curve, including "GAM_s" (default), "LR", "RLR", "spline", "GAM_l", "GAM_p", "GAM_tp", "GAM_cr", "GAM_ds",
#             "GAM_bs", "GAM_ad", "GAM_gp", "GAM_re", "GAM_sz", "GAM_fs"
# - majority_threshold: proportion of total integrated area required for one sign to be declared dominant (default: 0.55). Below this threshold in both directions, the curve is reported as "non-monotonic"
#
# Returns: "promoting", "mitigating", "non-monotonic", "neutral" (zero total area), or "undefined" (insufficient unique values, or curve fitting failed)
#
# Note: for binary/low-cardinality x, this bypasses curve-fitting entirely and falls back to a direct two-group contrast, since a spline/GAM has no meaningful "shape" to estimate from so few distinct x values.

n_crossing<-function(x, y, glm.type="GAM_s",majority_threshold=0.55){
  ok<-complete.cases(x, y)
  x<-x[ok]
  y<-y[ok]
  
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return("undefined")
  if(is_low_cardinality(x)) return(two_group_contrast(x, y)$dir)
  
  grid<-make_grid(x)
  if(is.null(grid)) return("undefined")
  k_use<-safe_k(x)
  
  fit<-tryCatch({
    if(glm.type=="LR") suppressWarnings(lm(y ~ x))
    else if(glm.type=="RLR") suppressWarnings(MASS::rlm(y ~ x, method="MM"))
    else if(glm.type=="spline") suppressWarnings(lm(y ~ splines::ns(x, df=min(3, length(unique(x))-1))))
    else if(glm.type=="GAM_l") suppressWarnings(gam(y ~ x, method="REML"))
    else{
      bs_map<-c(GAM_s="tp", GAM_p="ps", GAM_tp="tp", GAM_cr="cr", GAM_ds="ds", GAM_bs="bs", GAM_re="re", GAM_gp="gp", GAM_ad="ad", GAM_sz="sz", GAM_fs="fs")
      suppressWarnings(gam(y ~ s(x, bs=bs_map[[glm.type]] %||% "tp", k=k_use), method="REML"))
    }
  }, error=function(e) NULL)
  
  if(is.null(fit)) return("uncertain")
  pred_on_grid<-as.numeric(predict(fit, newdata=data.frame(x=grid)))
  if(anyNA(pred_on_grid)) return("uncertain")
  
  zero_indices<-which(diff(sign(pred_on_grid)) != 0)
  if(length(zero_indices)==0){
    return(if_else(mean(pred_on_grid, na.rm=TRUE) >= 0, "promoting", "mitigating"))
  }
  
  seg_bounds<-unique(c(1, zero_indices, length(grid)))
  seg_areas<-numeric(length(seg_bounds)-1)
  for(k in seq_len(length(seg_bounds)-1)){
    rng_idx<-seg_bounds[k]:seg_bounds[k+1]
    xg<-grid[rng_idx]
    yg<-pred_on_grid[rng_idx]
    dxv<-diff(xg)
    ymid<-(yg[-1]+yg[-length(yg)])/2
    seg_areas[k]<-sum(dxv*ymid)
  }
  pos_area<-sum(seg_areas[seg_areas>0])
  neg_area<-sum(-seg_areas[seg_areas<0])
  total<-pos_area+neg_area
  
  if(total==0) return("neutral")
  prop_pos<-pos_area/total
  if(prop_pos>=majority_threshold) return("promoting")
  if(prop_pos<=(1-majority_threshold)) return("mitigating")
  "non-monotonic"
}

#- - - - -
# Approach 7: pairwise bin concordance (probabilistic dominance) based on Kendall's tau (Rcpp)
#
# Goal: to estimate, in a fully non-parametric way, the probability that a higher feature value is associated with a higher attribution/explanation/effect value, by comparing all pairs of quantile bins of x and counting how often the bin with the larger x also has the larger mean y.
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - majority_threshold: proportion of concordant pairs (weighted by bin counts) required for a confident directional call (default: 0.55)
# - n_quantile_bins: number of quantile bins used to discretize x before pairwise comparison (default: 200)
#
# Returns: "promoting", "mitigating", "neutral" (no sign dominates), "uncertain" (insufficient comparable pairs), or "undefined" (x or y has only one unique value)
#
# Note: for binary/sparse x, this bypasses quantile binning and instead directly compares the two existing groups via a Mann-Whitney-style rank comparison.

# Compile C++ function
if (!exists("bin_pairwise_counts")) {
  Rcpp::cppFunction('
Rcpp::List bin_pairwise_counts(NumericVector bx, NumericVector bs, NumericVector bc) {
  int B=bx.size();
  double pos=0.0, neg=0.0, total=0.0;
  for (int i=0; i < B; ++i) {
    for (int j=i+1; j < B; ++j) {
      double ci=bc[i], cj=bc[j];
      if (ci <= 0.0 || cj <= 0.0) continue;
      double pairs=ci*cj;
      if (bx[i] > bx[j]) {
        total += pairs;
        if (bs[i] > bs[j]) pos += pairs;
        else if (bs[i] < bs[j]) neg += pairs;
      } else if (bx[j] > bx[i]) {
        total += pairs;
        if (bs[j] > bs[i]) pos += pairs;
        else if (bs[j] < bs[i]) neg += pairs;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("pos")=pos,
    Rcpp::Named("neg")=neg,
    Rcpp::Named("total")=total);
}', depends="Rcpp")
}

pairwise_direction<-function(x, y, majority_threshold=0.55, n_quantile_bins=200){
  ok<-complete.cases(x, y)
  x<-x[ok]
  y<-y[ok]
  
  if(length(unique(x)) <= 1 || length(unique(y)) <= 1) return("undefined")
  
  # Binary / sparse
  if(length(unique(x)) <= 2 || quantile(x, 0.75, na.rm=TRUE) == 0){
    lv<-min(x)
    hv<-max(x)
    g0<-y[x==lv]
    g1<-y[x==hv]
    
    if(length(g0)==0 || length(g1)==0) return("uncertain")
    g0s<-sort(g0)
    pos<-sum(findInterval(g1, g0s, left.open=TRUE))
    tot<-as.double(length(g0))*as.double(length(g1))
    pp<-pos/tot
    if(pp >= majority_threshold) return("promoting")
    if(pp <= 1-majority_threshold) return("mitigating")
    return("neutral")
  }
  
  # Continuous: binning
  probs<-seq(0, 1, length.out=n_quantile_bins+1)
  breaks<-unique(quantile(x, probs=probs, na.rm=TRUE, type=7))
  if(length(breaks) <= 2)
    breaks<-seq(min(x,na.rm=TRUE), max(x,na.rm=TRUE), length.out=3)
  
  bins<-cut(x, breaks=breaks, include.lowest=TRUE)
  bx<-as.numeric(tapply(x, bins, mean, na.rm=TRUE))
  bs<-as.numeric(tapply(y, bins, mean, na.rm=TRUE))
  bc<-as.numeric(tapply(y, bins, length))
  ok2<-!is.na(bx) & !is.na(bs) & !is.na(bc)
  
  bx<-bx[ok2]
  bs<-bs[ok2]
  bc<-bc[ok2]
  
  if(length(bx) < 2) return("neutral")
  
  cnts<-bin_pairwise_counts(bx, bs, bc)
  if(cnts$total<=0) return("uncertain")
  pos_p<-as.numeric(cnts$pos)
  neg_p<-as.numeric(cnts$neg)
  tot_p<-as.numeric(cnts$total)
  
  pp<-pos_p/tot_p
  np<-neg_p/tot_p
  if(pp>=majority_threshold) return("promoting")
  if(np>=majority_threshold) return("mitigating")
  return("neutral")
}

#- - - - - -
## Approach 8: stochastic dominance
#
# Goal: to test, directly on the raw (unbinned) attribution/explanation/effect metric distributions, whether observations with feature values in the top quantile have systematically higher (or lower) attribution/explanation/effect values than observations in the bottom quantile (i.e., whether the upper-x distribution stochastically dominates the lower-x distribution)
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - q: quantile threshold defining the low/high split (default: 0.20, i.e. bottom vs. top quintile)
# - min_n: minimum number of observations required in each extreme group for the comparison to be considered reliable (default: 30)
# - majority_threshold: proportion of pairwise comparisons required for a confident directional call (default: 0.55)
#
# Returns: "promoting", "mitigating", "neutral" (no sign dominates), "uncertain" (groups too small), or "undefined" (sample too small overall to form two groups of size >= min_n)

stochastic_dominance_direction<-function(x, y, q=0.2, min_n=30, majority_threshold=0.55){
  ok<-complete.cases(x, y)
  x<-x[ok]
  y<-y[ok]
  
  if(length(x) < 2 * min_n){
    return("undefined")
  }
  
  lo_thr<-quantile(x, q, na.rm=TRUE)
  hi_thr<-quantile(x, 1 - q, na.rm=TRUE)
  low<-y[x <= lo_thr]
  high<-y[x >= hi_thr]
  
  if(length(low) < min_n || length(high) < min_n){
    return("uncertain")
  }
  
  p_sup<-mean(outer(high, low, `>`))
  p_inf<-mean(outer(high, low, `<`))
  
  if(p_sup>=majority_threshold) return("promoting")
  if(p_inf>=majority_threshold) return("mitigating")
  return("neutral")
}

#- - - - - -
## Approach 9: signed-AUC dominance index
#
# Goal: to determine the feature's net direction by comparing the total integrated magnitude of positive vs. negative attribution/explanation/effect values across the feature's range
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - smoother: smoothing method used to obtain the fitted curve over a clean evaluation grid before integrating: "GAM" (default), "LOESS", "supsmu" (Friedman's SuperSmoother), or "kernel".
# - dominance_threshold: proportion of total integrated area required for one sign to be declared dominant (default: 0.55)
# - n_grid: number of points in the synthetic evaluation grid used for integration (default: 300)
#
# Returns: "promoting", "mitigating", "non-monotonic" (neither sign reaches dominance_threshold), "neutral" (zero total integrated area), "uncertain" (the chosen smoother failed to fit), or "undefined" (insufficient unique values in x or y)
#
# Note: 
# - for binary/low-cardinality x, this bypasses curve-fitting entirely and falls back to a direct two-group contrast, for the same reason as Approach 6
# - different smoothers make different bias/variance trade-offs, so comparing results across smoothers is a useful robustness check for borderline features

signed_auc_direction<-function(x, y, smoother="GAM", dominance_threshold=0.55, n_grid=300){
  ok<-complete.cases(x, y)
  x<-x[ok]
  y<-y[ok]
  
  if(length(unique(x)) < 2 || length(unique(y)) < 2) return("undefined")
  
  if(is_low_cardinality(x)){
    res<-two_group_contrast(x, y)
    return(res$dir)
  }
  
  grid<-make_grid(x, n_grid=n_grid)
  if(is.null(grid)) return("undefined")
  k_use<-safe_k(x)
  
  pred_on_grid<-tryCatch({
    if(smoother=="GAM"){
      
      fit<-suppressWarnings(gam(y ~ s(x, k=k_use), method="REML"))
      as.numeric(predict(fit, newdata=data.frame(x=grid)))
      
    }else if(smoother=="LOESS"){
      
      fit<-suppressWarnings(loess(y ~ x, span=0.75, degree=2))
      as.numeric(predict(fit, newdata=data.frame(x=grid)))
      
    }else if(smoother=="kernel"){
      
      fit<-suppressWarnings(npreg(y ~ x))
      as.numeric(predict(fit, newdata=data.frame(x=grid)))
      
    }else if(smoother=="supsmu"){
      
      o<-order(x)
      sm<-suppressWarnings(supsmu(x[o], y[o]))
      approx(sm$x, sm$y, xout=grid)$y
      
    }else stop("unknown smoother")
  }, error=function(e) rep(NA_real_, length(grid)))
  
  if(all(is.na(pred_on_grid))) return("uncertain")
  
  dxv<-diff(grid)
  ymid<-(pred_on_grid[-1]+pred_on_grid[-length(pred_on_grid)])/2
  seg_area<-dxv*ymid
  
  pos_area<-sum(seg_area[seg_area>0], na.rm=TRUE)
  neg_area<-sum(-seg_area[seg_area<0], na.rm=TRUE)
  total<-pos_area+neg_area
  
  if(total==0 || is.na(total)) return("neutral")
  prop_pos<-pos_area/total
  if(prop_pos>= dominance_threshold) return("promoting")
  if(prop_pos<= 1-dominance_threshold) return("mitigating")
  "non-monotonic"
}

#- - - - - -
## Pattern descriptor
#
# Rationale: a single direction label (promoting/mitigating/non-monotonic) cannot distinguish, for example, a simple monotonic increase from a relationship that 
# rises, plateaus, and then mildly declines, both might receive the same "promoting" label, but they describe very different clinical pictures and should be reported together, not as a substitute for one another.
#
# Goal: to separately characterize the shape of the attribution/explanation/effect metric-vs.-feature relationship, independently of its net direction.
#
# Parameters:
# - x: vector of feature values
# - y: vector of attribution/explanation/effect (e.g., SHAP, ALE) values
# - smoother: smoothing method used to fit the curve before classifying its shape: "GAM" (default), "LOESS", or "supsmu"
# - n_grid: grid resolution for evaluating the fitted curve (default: 300)
# - flat_tol: proportion of the attribution/explanation/effect metric's overall standard deviation below which the fitted curve's range is considered too small to represent a real effect (default: 0.05)
#
# Returns one of:
# - "monotonic-increasing", "monotonic-decreasing" (single inflection-free trend)
# - "U-shaped", "inverted-U-shaped" (one inflection near the middle of the range)
# - "threshold/plateau" (one inflection near either edge of the range)
# - "complex-non-monotonic" (two or more inflections)
# - "flat/no-effect" (fitted curve's range is negligible)
# - "binary" / "ordinal-low-cardinality" (shape is not meaningfully defined by a continuous curve)
# - "undefined" (curve fitting failed or insufficient data)

pattern_descriptor<-function(x, y, smoother="GAM", n_grid=300, flat_tol=0.05){
  ok<-complete.cases(x, y)
  x<-x[ok]
  y<-y[ok]
  
  if(length(unique(x)) < 2) return("undefined")
  
  if(is_low_cardinality(x)){
    if(length(unique(x))==2) return("binary")
    return("ordinal-low-cardinality")
  }
  
  grid<-make_grid(x, n_grid=n_grid)
  if(is.null(grid)) return("undefined")
  k_use<-safe_k(x)
  
  pred<-tryCatch({
    if(smoother=="GAM"){
      
      fit<-suppressWarnings(gam(y ~ s(x, k=k_use), method="REML"))
      as.numeric(predict(fit, newdata=data.frame(x=grid)))
      
    }else if(smoother=="LOESS"){
      
      fit<-suppressWarnings(loess(y ~ x, span=0.75, degree=2))
      as.numeric(predict(fit, newdata=data.frame(x=grid)))
      
    }else if(smoother=="supsmu"){
      
      o<-order(x)
      sm<-suppressWarnings(supsmu(x[o], y[o]))
      approx(sm$x, sm$y, xout=grid)$y
      
    }else NULL
  }, error=function(e) NULL)
  
  if(is.null(pred) || all(is.na(pred))) return("undefined")
  
  rng<-diff(range(pred, na.rm=TRUE))
  if(is.na(rng) || rng < flat_tol*sd(y, na.rm=TRUE)) return("flat/no-effect")
  
  d<-diff(pred)/diff(grid)
  d<-d[is.finite(d)]
  sign_d<-sign(d)
  sign_d<-sign_d[sign_d!=0]
  
  if(length(sign_d)<2) return("undefined")
  n_slope_changes<-sum(diff(sign_d)!=0)
  
  if(n_slope_changes==0){
    return(if_else(mean(d, na.rm=TRUE)>=0, "monotonic-increasing", "monotonic-decreasing"))
  }
  if(n_slope_changes==1){
    turn_idx<-which(diff(sign_d)!=0)[1]
    turn_pos<-turn_idx/length(sign_d)
    if(turn_pos > 0.2 && turn_pos < 0.8){
      return(if_else(d[1]<0, "U-shaped", "inverted-U-shaped"))
    } else return("threshold/plateau")
  }
  "complex-non-monotonic"
}

#--------------------------------------------
#### Step 4 - Creating the overall wrappers ####

#- - - - - -
## Global direction
#
# Goals: 
# - to compute, for every feature, all individual direction-estimating approaches above as well as a pattern descriptor
# - to combine direction and pattern descriptor into two complementary consensus:
#      ¤ a shape-aware majority vote (Approach 10)
#      ¤ a weighted vote (Approach 11) computed independently as a sensitivity check on Approach 10
#
# Parameters:
# - feature: vector of feature names, repeated once per observation (long format)
# - feature_val: vector of feature values, same length/order as `feature`
# - shap_val: vector of attribution/explanation/effect values (e.g., SHAP, ALE), same length/order as `feature`
# - eps: numerical tolerance for Approach 1 (default: 1e-10)
# - cor.method: correlation method for Approach 2: "pearson", "spearman" (default), "kendall"
# - alpha: alpha risk threshold used by Approaches 2, 3, 5 (default: 0.05)
# - glm.type: regression/smoother type used by Approaches 3 and 6 (default: "GAM_s")
# - method.type: smoother used by Approach 4: "kernel", "LOESS" (default), "supsmu", "GAM"
# - smoother: smoother used by Approach 9 and the pattern descriptor: "GAM" (default), "LOESS", "supsmu", "kernel"
# - majority_threshold: proportion required for a confident call in Approaches 7 and 8 (default: 0.55)
# - n_quantile_bins: number of bins used by Approach 7 (default: 200)
# - q: quantile threshold for the low/high split in Approach 8 (default: 0.20)
# - min_n: minimum group size required by Approach 8 (default: 30)
# - dominance_threshold: proportion of integrated area required for a confident call in Approaches 6 and 9 (default: 0.55)
# - n_grid: resolution of the synthetic evaluation grid used by Approach 9 and the pattern descriptor (default: 300)
# - flat_tol: flatness tolerance used by the pattern descriptor (default: 0.05)
#
# Returns: a tibble with 
# - one row per feature
# - one column per individual approach
# - a pattern descriptor column
# - two consensus columns (consensus, consensus_weighted)

global_direction<-function(feature=feature, feature_val=feature_val, shap_val=shap_val,
                           eps=1e-10, cor.method="spearman", alpha=0.05, glm.type="GAM_s",method.type="LOESS", smoother="GAM",
                           majority_threshold=0.55, n_quantile_bins=200, q=0.2, min_n=30, dominance_threshold=0.55, n_grid=300, flat_tol=0.05){
  
  data_df<-as_tibble(data.frame(feature_name=feature, feature_val=feature_val, shap_val=shap_val))
  
  res_1<-suppressWarnings(dplyr::group_by(data_df, feature_name) |>
                            dplyr::group_modify(~{
                              tibble::tibble(mean_sign_dir=mean_sign_direction(.x$shap_val, eps),
                                             correlation_dir=cor_dir(.x$feature_val, .x$shap_val, alpha, cor.method),
                                             GLM_trend_dir=GLM_trend(.x$feature_val, .x$shap_val, alpha, glm.type),
                                             deriv_dir=reg_deriv(.x$feature_val, .x$shap_val, method.type),
                                             theil_sen_dir=theil_sen_direction(.x$feature_val, .x$shap_val, alpha),
                                             n_crossing_dir=n_crossing(.x$feature_val, .x$shap_val, glm.type, majority_threshold),
                                             pairwise_dir=pairwise_direction(.x$feature_val, .x$shap_val, majority_threshold, n_quantile_bins),
                                             SD_dir=stochastic_dominance_direction(.x$feature_val, .x$shap_val, q, min_n, majority_threshold),
                                             AUC_dir=signed_auc_direction(.x$feature_val, .x$shap_val, smoother, dominance_threshold, n_grid),
                                             pattern=pattern_descriptor(.x$feature_val, .x$shap_val, smoother, n_grid, flat_tol))
                            })) %>% rename(feature=feature_name) %>% dplyr::ungroup()
  
  vote_cols<-c("mean_sign_dir","correlation_dir","GLM_trend_dir","deriv_dir","theil_sen_dir", "n_crossing_dir","pairwise_dir","SD_dir","AUC_dir")
  
  vote_mat<-as.matrix(res_1[, vote_cols])
  
  #- - - -
  ## Approach 10: unweighted majority vote
  
  res_1$consensus<-mapply(function(row_idx){
    x<-vote_mat[row_idx, ]
    directional<-x[x %in% c("promoting","mitigating")]
    
    if(length(directional) >= ceiling(length(vote_cols)/2)){
      tbl<-sort(table(directional), decreasing=TRUE)
      return(names(tbl)[1])
    }
    
    shape_flags<-c(x[["AUC_dir"]], x[["n_crossing_dir"]])
    pattern_flag<-res_1$pattern[row_idx] %in% c("U-shaped","inverted-U-shaped","complex-non-monotonic","threshold/plateau")
    
    if(sum(shape_flags=="non-monotonic", na.rm=TRUE) >= 2 || pattern_flag){
      return("non-monotonic")
    }
    
    tbl<-sort(table(x), decreasing=TRUE)
    if(tbl[1] >= ceiling(length(vote_cols)/2)) return(names(tbl)[1])
    "uncertain"
  }, seq_len(nrow(vote_mat)))

  #- - - -
  ## Approach 11: weighted consensus (sensitivity analysis of approach 10)
  
  approach_weights<-c(mean_sign_dir=0.5, # most fragile approach under non-monotonicity -> weight 0.5
                      deriv_dir=0.5, # most fragile approach under non-monotonicity -> weight 0.5
                      correlation_dir=0.75, # assumes monotonic association, blind to non-monotonic patterns -> weight 0.75
                      GLM_trend_dir=1.0, # rely on a single smoother specification -> weight 1.0
                      n_crossing_dir=1.0, # rely on a single smoother specification -> weight 1.0
                      AUC_dir=1.25, # robust to outliers/non-monotonicity -> weight 1.25
                      theil_sen_dir=1.25, # robust to outliers/non-monotonicity -> weight 1.25
                      pairwise_dir=1.5, # most assumption-light -> highest weight
                      SD_dir=1.5) # most assumption-light -> highest weight
  
  res_1$consensus_weighted<-apply(vote_mat, 1, function(x){
    
    w<-approach_weights[vote_cols]
    pos_w<-sum(w[x=="promoting"], na.rm=TRUE)
    neg_w<-sum(w[x=="mitigating"], na.rm=TRUE)
    tot_w<-sum(w, na.rm=TRUE)
    if(pos_w==0 && neg_w==0) return("uncertain")
    if(pos_w/tot_w >= 0.5*sum(w)/tot_w && pos_w > neg_w && pos_w/(pos_w+neg_w) >= 0.6) return("promoting")
    if(neg_w/(pos_w+neg_w) >= 0.6) return("mitigating")
    "uncertain"
    
  })
  
  res_1
}

#- - - - - -
## Approach 12: Bootstrap stability wrapper (model-agnostic confidence measure)
#
# Rationale: a single deterministic direction label hides whether the underlying estimate is stable. This is especially relevant for non-monotonic, low-signal, or low-cardinality
#            features, or when correlated predictors distort SHAP credit allocation
#
# Goal: to resample observations with replacement, recompute the full global_direction() pipeline on each replicate, and report the proportion of replicates whose consensus label matches
#       the consensus label obtained on the original (non-resampled) sample
#
# Interpretation: 
# - boot_agreement close to 1 means the direction label is stable under resampling
# - lower values mean the direction label is sensitive to which observations are drawn and should be reported with that caveat
#
# Caveat specific to ALE: ALE values are already aggregated onto a small number of grid bins (as few as 2 for a binary feature). Resampling at the bin level, as done here, does
#                         not carry the same statistical meaning as resampling raw observations. For a rigorous ALE stability estimate, prefer resampling the raw data and 
#                         recomputing ALE from scratch on each replicate; bin-level resampling here is a quick diagnostic, not a final result.
#
# Parameters:
# - feature: vector of feature names, repeated once per observation/bin (long format)
# - feature_val: vector of feature values, same length/order as `feature`
# - shap_val: vector of attribution/explanation/effect values (e.g., SHAP, ALE), same length/order as `feature`
# - n_boot: number of bootstrap replicates (default: 100; reduce for quick checks)
# - ...: additional parameters forwarded unchanged to global_direction() on every replicate
#
# Returns: a tibble with 
# - one row per feature
# - the original consensus direction label
# - the bootstrap agreement proportion (boot_agreement, between 0 and 1)

bootstrap_direction_stability<-function(feature=feature, shap_val=shap_val, feature_val=feature_val, n_boot=100, ...){
  
  full_df<-data.frame(feature=feature, shap_val=shap_val, feature_val=feature_val)
  
  base_res<-global_direction(full_df$feature, full_df$shap_val, full_df$feature_val, ...)
  base_labels<-base_res %>% select(feature, consensus)
  
  boot_labels<-dplyr::bind_rows(lapply(seq_len(n_boot), function(b){
    
    boot_idx<-sample(seq_len(nrow(full_df)), size=nrow(full_df), replace=TRUE)
    boot_df<-full_df[boot_idx, ]
    boot_res<-global_direction(boot_df$feature, boot_df$shap_val, boot_df$feature_val, ...)
    boot_res %>% select(feature, consensus) %>% mutate(boot_id=b)
    
  }))
  
  boot_labels %>%
    left_join(base_labels, by="feature", suffix=c("_boot","_base")) %>%
    group_by(feature) %>%
    summarise(base_consensus=unique(consensus_base),
              boot_agreement=mean(consensus_boot==consensus_base, na.rm=TRUE),
              .groups="drop")
}

#- - - - - -
## Approach 13: Model-refit stability wrapper
#
# Goal: to distinguish a genuine effect from a spurious one picked up by overfitting in a single model fit
#
# Method: unlike bootstrap_direction_stability, which resamples already-computed SHAP rows from one fixed model, this refits the underlying model from 
#         scratch on each bootstrap resample of the raw training data, recomputes SHAP, and checks whether the feature's direction is consistent across 
#         independently-fitted models
#
# Parameters:
# - raw_df: original training data frame (must include the outcome column)
# - feature_cols: vector of feature column names used in the model
# - y_col: name of the outcome column
# - n_boot: number of model refits (default: 20)
# - xgb_params: list of XGBoost hyperparameters matching the original fit
# - nrounds: number of boosting rounds matching the original fit
#
# Returns: a tibble with the proportion of refits in which the direction was "promoting", "mitigating", or non-directional for each variable
#
# Interpretation: low max-proportion indicates the original fit's direction for that feature is not trustworthy

model_refit_stability<-function(raw_df, feature_cols, y_col="y", n_boot=20,
                                xgb_params=list(objective="reg:squarederror", max_depth=4,learning_rate=0.05, subsample=0.8),nrounds=200, ...){
  
  results<-dplyr::bind_rows(lapply(seq_len(n_boot), function(b){
    boot_idx<-sample(seq_len(nrow(raw_df)), nrow(raw_df), replace=TRUE)
    boot_data<-raw_df[boot_idx, ]
    
    X_b<-as.matrix(boot_data[, feature_cols, drop=FALSE])
    colnames(X_b)<-feature_cols
    
    dtrain_b<-xgboost::xgb.DMatrix(X_b, label=boot_data[[y_col]])
    
    # model fitting
    fit_b<-xgboost::xgb.train(params=xgb_params, data=dtrain_b, nrounds=nrounds, verbose=0)
    
    # computing SHAP
    shap_b<-predict(fit_b, dtrain_b, predcontrib=TRUE)
    shap_df_b<-as.data.frame(shap_b[, feature_cols])
    shap_long_b<-shap_df_b %>% mutate(id=row_number()) %>% pivot_longer(cols=-id, names_to="feature", values_to="shap_val")
    feat_long_b<-as.data.frame(X_b) %>% mutate(id=row_number()) %>% pivot_longer(cols=-id, names_to="feature", values_to="feature_val")
    merged_b<-left_join(feat_long_b, shap_long_b, by=c("id","feature"))
    
    dir_b<-global_direction(feature=merged_b$feature, feature_val=merged_b$feature_val, shap_val=merged_b$shap_val, ...)
    dir_b %>% select(feature, consensus) %>% mutate(boot_id=b)
  }))
  
  results %>% group_by(feature, consensus) %>% summarise(n=n(), .groups="drop") %>% group_by(feature) %>% mutate(prop=n/sum(n)) %>% ungroup()
}

#--------------------------------------------
#### Step 5: Simulating a dataset with known feature behaviors ####

set.seed(123)
n<-800
K<-8

# Feature types with known true directions
x_pos<-rnorm(n, 5, 2) # monotone positive
x_neg<-rnorm(n, 5, 2) # monotone negative
x_ushape<-rnorm(n, 0, 2) # U-shaped
x_thresh<-rnorm(n, 0, 1) # threshold effect
x_bin<-rbinom(n, 1, 0.4) # binary
x_noise<-rnorm(n, 0, 1) # no effect
non_monotone<-as.data.frame(exp(MASS::mvrnorm(n,rep(0,8),0.5^as.matrix(dist(1:8)))/3))[,c(3,6,7,8)] # other non-monotone distributions
x_ord<-sample(1:5, n, replace=TRUE, prob=c(0.3,0.25,0.2,0.15,0.1)) # ordinal feature (e.g., disease stage / severity scale 1-5), monotone increasing effect
x_exp<-rexp(n, rate=0.5) # exponential-distributed feature (e.g., time-to-event-like exposure, always positive, right-skewed)
x_gamma<-rgamma(n, shape=2, rate=0.5) # gamma-distributed feature (e.g., biomarker concentration, right-skewed, bounded at 0)
x_negbin<-rnbinom(n, size=2, mu=3) # negative-binomial-distributed feature (e.g., count of comorbidities/prior hospitalisations, over-dispersed count data)

# Simulating an "age-like" feature
x_age<-runif(n, 50, 100)
age_effect<-function(a){
  dplyr::case_when(a<70~1.8*(a-50)/20, # steep rise: from 50 to 70
                   a>=70&a<=80~1.8, # plateau: from 70 to 80
                   a>80~1.8-0.6*(a-80)/20) # decline: from 80 to 100
}

# Feature names
feature_cols<-c("x_pos","x_neg","x_ushape","x_thresh","x_bin","x_noise","x_NM1","x_NM2","x_NM3","x_NM4","x_age","x_ord","x_exp","x_gamma","x_negbin")

# Adding the feature simulated into a dataset
sim_df<-data.frame(x_pos, x_neg, x_ushape, x_thresh, x_bin, x_noise, non_monotone, x_age, x_ord, x_exp, x_gamma, x_negbin)
colnames(sim_df)<-feature_cols

# Creating the true outcome
y<-1.5*x_pos -
  1.2*x_neg +
  0.8*x_ushape^2 +
  1.5*(x_thresh > 0.8) +
  1.0*x_bin +
  rnorm(n, 0, 1.5) +
  1.0*non_monotone[,1] +
  1.2*non_monotone[,2] -
  0.8*non_monotone[,3] +
  1.9*non_monotone[,1] +
  age_effect(x_age) +
  0.6*x_ord +
  0.3*log1p(x_exp) +
  0.4*x_gamma +
  0.25*x_negbin

# Adding y to the simulated dataset
sim_df$y<-y

#--------------------------------------------
#### Step 5: Fitting XGBoost and computing SHAP and ALE values ####

# Creating the feature matrix
X_mat<-as.matrix(sim_df[, feature_cols, drop=FALSE])
colnames(X_mat)<-feature_cols

# Creating the xgb.DMatrix
dtrain<-xgboost::xgb.DMatrix(X_mat, label=y)

# Fitting XGBoost
xgb_fit<-xgboost::xgb.train(params=list(objective="reg:squarederror",max_depth=4,learning_rate=0.05,subsample=0.8),data=dtrain,nrounds=200,verbose=0)

# Computing SHAP values
shap_mat<-predict(xgb_fit, dtrain, predcontrib=TRUE)
shap_df<-as.data.frame(shap_mat[, feature_cols]) # drop BIAS/intercept column

shap_df2<-shap_df %>% dplyr::mutate(id=row_number()) %>% pivot_longer(cols=-id)
colnames(shap_df2)<-c('id','feature_names','shap_val')
sim_df_long<-sim_df %>% mutate(id=row_number()) %>% dplyr::select(-y) %>% pivot_longer(cols=-c(id)) %>% distinct
colnames(sim_df_long)<-c('id','feature_names','feature_val')

shap_val<-left_join(sim_df_long,shap_df2,by=c('id','feature_names'))
shap_val<-shap_val %>% select(-c(id))
colnames(shap_val)[1]<-"feature"

# computing ALE values
ale_obj<-ale::ALE(xgb_fit, y_col="y", x_cols=feature_cols, data=sim_df, pred_fun=predict_xgb_ale)

ale_list<-ale_obj@effect$y$ale$d1

ale_val<-as_tibble(dplyr::bind_rows(lapply(names(ale_list), function(f){
  dat<-ale_list[[f]]
  grid_col<-setdiff(names(dat), c(".n", ".y", ".y_lo", ".y_mean", ".y_median", ".y_hi"))[1]
  data.frame(feature=f,
             feature_val=as.numeric(as.character(dat[[grid_col]])),
             ale_val=as.numeric(dat$.y),
             stringsAsFactors=FALSE)
})))

#--------------------------------------------
#### Step 6: Visualization ####

sv<-shapviz(xgb_fit, X_pred=X_mat, X=sim_df[, feature_cols, drop = FALSE], interactions=TRUE)

#- - - - -
### SHAP local explanation

sv_force(sv,fill_colors = c("#F9CBC2", "#40B696"),row_id=1) # for observation 1
sv_waterfall(sv,fill_colors = c("#F9CBC2", "#40B696"),row_id=1) # for observation 1

#- - - - -
### SHAP importance plot

sv_importance(sv,fill="#40B696")+theme_Gaia()

#- - - - -
### Dependence plots

## ALE
ale_val$feature<-factor(ale_val$feature,
                        levels=feature_cols,
                        labels=c("monotonic positive","monotonic negative","U-shaped","threshold effect","binary","noise","non-monotonic 1",
                                  "non-monotonic 2","non-monotonic 3","non-monotonic 4","non-monotonic 5","ordinal","exponential","gamma",
                                  "negative binomial"))

plot_dep_all(df=ale_val, feature_col="feature", x_col="feature_val", y_col="ale_val", y_name="ALE values", ncol=3)

## SHAP
shap_val$feature<-factor(shap_val$feature,
                         levels=feature_cols,
                         labels=c("monotonic positive","monotonic negative","U-shaped","threshold effect","binary","noise","non-monotonic 1",
                                  "non-monotonic 2","non-monotonic 3","non-monotonic 4","non-monotonic 5","ordinal","exponential","gamma",
                                  "negative binomial"))

plot_dep_all(df=shap_val, feature_col="feature", x_col="feature_val", y_col="shap_val", y_name="SHAP values", ncol=3)

#--------------------------------------------
#### Step 7: Determining the direction ####

#- - - - -
## SHAP

# Computing direction
direction_shap<-global_direction(feature=shap_val$feature, shap_val=shap_val$shap_val, feature_val=shap_val$feature_val,
                                 eps=1e-10, cor.method="spearman", alpha=0.05, glm.type="GAM_s", method.type="LOESS", smoother="GAM",
                                 majority_threshold=0.55, n_quantile_bins=200, q=0.2, min_n=30, dominance_threshold=0.55)

# Adding true expected directions for comparison
truth<-data.frame(feature=unique(direction_shap$feature), true_direction=c("promoting","mitigating","non-monotonic","promoting","promoting","uncertain",
                                                                           "promoting","promoting","mitigating","mitigating","promoting","promoting","promoting",
                                                                           "promoting","promoting"))
truth<-left_join(direction_shap, truth, by="feature")
truth %>% select(feature, consensus, consensus_weighted, true_direction, pattern)

# Direction heatmap

dir_long<-truth %>%
  dplyr::select(feature,mean_sign_dir,correlation_dir,GLM_trend_dir,deriv_dir,theil_sen_dir,n_crossing_dir,pairwise_dir,SD_dir,AUC_dir,consensus,
                consensus_weighted,true_direction) %>%
  tidyr::pivot_longer(-feature, names_to="method", values_to="direction") %>%
  mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","non-monotonic","undefined","uncertain")),
         method=factor(method,
                       levels=c("true_direction","mean_sign_dir","correlation_dir","GLM_trend_dir","deriv_dir",
                                "theil_sen_dir","n_crossing_dir","pairwise_dir","SD_dir","AUC_dir","consensus","consensus_weighted"),
                       labels=c("True direction",
                                "Approach 1\n(Mean sign)",
                                "Approach 2\n(Pairwise correlation)",
                                "Approach 3\n(Overall trend)",
                                "Approach 4\n(Derivative)",
                                "Approach 5\n(Theil-Sen slope)",
                                "Approach 6\n(Zero crossings)",
                                "Approach 7\n(Pairwise bin concordance)",
                                "Approach 8\n(Stochastic dominance)",
                                "Approach 9\n(Signed-AUC)",
                                "Approach 10\n(Unweighted consensus)",
                                "Approach 11\n(Weighted consensus)")))

ggplot(dir_long, aes(x=method, y=feature, fill=direction)) +
  geom_tile(color="black") +
  scale_fill_manual("Direction",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","non-monotonic"="#FDDDA0",
                                          "uncertain"="grey25","mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
  scale_x_discrete("",expand=c(0,0))+
  scale_y_discrete("Feature",expand=c(0,0))+
  labs(fill="Direction") +
  theme_Gaia()

# Bootstrap stability

set.seed(123)
stability_direction_shap<-bootstrap_direction_stability(feature=shap_val$feature, shap_val=shap_val$shap_val, feature_val=shap_val$feature_val,
                                                        n_boot=50, eps=1e-10, cor.method="spearman", alpha=0.05, glm.type="GAM_s",
                                                        method.type="LOESS", smoother="GAM", majority_threshold=0.55, n_quantile_bins=200,
                                                        q=0.2, min_n=30, dominance_threshold=0.55)
stability_direction_shap

# model-refit stability

set.seed(123)
mode_refit_stab<-model_refit_stability(sim_df, feature_cols, y_col="y", n_boot=20)
mode_refit_stab

#- - - - -
## ALE

# Computing direction
direction_ale<-global_direction(feature=ale_val$feature, shap_val=ale_val$ale_val, feature_val=ale_val$feature_val,
                                eps=1e-10, cor.method="spearman", alpha=0.05, glm.type="GAM_s", method.type="LOESS", smoother="GAM",
                                majority_threshold=0.55, n_quantile_bins=200, q=0.2, min_n=30, dominance_threshold=0.55)

# Adding true expected directions for comparison
truth<-data.frame(feature=unique(direction_ale$feature), true_direction=c("promoting","mitigating","non-monotonic","promoting","promoting","uncertain",
                                                                          "promoting","promoting","mitigating","mitigating","promoting","promoting","promoting",
                                                                          "promoting","promoting"))
truth<-left_join(direction_ale, truth, by="feature")
truth %>% select(feature, consensus, consensus_weighted, true_direction, pattern)

# Direction heatmap

dir_long<-truth %>%
  dplyr::select(feature,mean_sign_dir,correlation_dir,GLM_trend_dir,deriv_dir,theil_sen_dir,n_crossing_dir,pairwise_dir,SD_dir,AUC_dir,consensus,
                consensus_weighted,true_direction) %>%
  tidyr::pivot_longer(-feature, names_to="method", values_to="direction") %>%
  mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","non-monotonic","undefined","uncertain")),
         method=factor(method,
                       levels=c("true_direction","mean_sign_dir","correlation_dir","GLM_trend_dir","deriv_dir",
                                "theil_sen_dir","n_crossing_dir","pairwise_dir","SD_dir","AUC_dir","consensus","consensus_weighted"),
                       labels=c("True direction",
                                "Approach 1\n(Mean sign)",
                                "Approach 2\n(Pairwise correlation)",
                                "Approach 3\n(Overall trend)",
                                "Approach 4\n(Derivative)",
                                "Approach 5\n(Theil-Sen slope)",
                                "Approach 6\n(Zero crossings)",
                                "Approach 7\n(Pairwise bin concordance)",
                                "Approach 8\n(Stochastic dominance)",
                                "Approach 9\n(Signed-AUC)",
                                "Approach 10\n(Unweighted consensus)",
                                "Approach 11\n(Weighted consensus)")))

ggplot(dir_long, aes(x=method, y=feature, fill=direction)) +
  geom_tile(color="black") +
  scale_fill_manual("Direction",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","non-monotonic"="#FDDDA0",
                                          "uncertain"="grey25","mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
  scale_x_discrete("",expand=c(0,0))+
  scale_y_discrete("Feature",expand=c(0,0))+
  labs(fill="Direction") +
  theme_Gaia()

# Bootstrap stability

set.seed(123)
stability_direction_ale<-bootstrap_direction_stability(feature=ale_val$feature, shap_val=ale_val$ale_val, feature_val=ale_val$feature_val,
                                                        n_boot=50, eps=1e-10, cor.method="spearman", alpha=0.05, glm.type="GAM_s",
                                                        method.type="LOESS", smoother="GAM", majority_threshold=0.55, n_quantile_bins=200,
                                                        q=0.2, min_n=30, dominance_threshold=0.55)
stability_direction_ale

#--------------------------------------------
#### Step 8: Sensitivity analysis

#- - - -
### Approach 7 - pairwise bin concordance (probabilistic dominance) based on Kendall's tau (Rcpp)

## Threshold sensitivity

# ALE
thresh_df<-dplyr::bind_rows(lapply(seq(0.50, 0.75, by=0.05), function(thr){
  dplyr::bind_rows(lapply(split(ale_val, ale_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               pairwise_dir=pairwise_direction(d$feature_val, d$ale_val, majority_threshold=thr, n_quantile_bins=200))
  }))
}))

thresh_df %>% pivot_wider(names_from=threshold,values_from=pairwise_dir)

ALE_thresh<-thresh_df %>% mutate(metric="ALE")
colnames(ALE_thresh)[3]<-"direction"

# SHAP
thresh_df<-dplyr::bind_rows(lapply(seq(0.50, 0.75, by=0.05), function(thr){
  dplyr::bind_rows(lapply(split(shap_val, shap_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               pairwise_dir=pairwise_direction(d$feature_val, d$shap_val, majority_threshold=thr, n_quantile_bins=200))
  }))
}))

thresh_df %>% pivot_wider(names_from=threshold,values_from=pairwise_dir)

SHAP_thresh<-thresh_df %>% mutate(metric="SHAP")
colnames(SHAP_thresh)[3]<-"direction"

SHAP_thresh<-bind_rows(SHAP_thresh,ALE_thresh) %>%
             mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","non-monotonic","undefined","uncertain")),
                    metric=factor(metric,levels=c("SHAP","ALE")),
                    threshold=factor(threshold,levels=seq(0.50, 0.75, by=0.05)))

ggplot(SHAP_thresh, aes(x=threshold, y=feature, fill=direction)) +
  geom_tile(color="black") +
  scale_fill_manual("Direction",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","non-monotonic"="#FDDDA0",
                                         "uncertain"="grey25","mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
  scale_x_discrete("Majority threshold",expand=c(0,0))+
  scale_y_discrete("Feature",expand=c(0,0))+
  labs(fill="Direction") +
  theme_Gaia()+
  facet_grid(.~metric)

## Bin number sensitivity

# ALE
bin_df<-dplyr::bind_rows(lapply(seq(50,500,50), function(thr){
  dplyr::bind_rows(lapply(split(ale_val, ale_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               pairwise_dir=pairwise_direction(d$feature_val, d$ale_val, majority_threshold=0.55, n_quantile_bins=thr))
  }))
}))

bin_df %>% pivot_wider(names_from=threshold,values_from=pairwise_dir)

ALE_bin<-bin_df %>% mutate(metric="ALE")
colnames(ALE_bin)[3]<-"direction"

# SHAP
bin_df<-dplyr::bind_rows(lapply(seq(50,500,50), function(thr){
  dplyr::bind_rows(lapply(split(shap_val, shap_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               pairwise_dir=pairwise_direction(d$feature_val, d$shap_val, majority_threshold=0.55, n_quantile_bins=thr))
  }))
}))

bin_df %>% pivot_wider(names_from=threshold,values_from=pairwise_dir)

SHAP_bin<-bin_df %>% mutate(metric="SHAP")
colnames(SHAP_bin)[3]<-"direction"

SHAP_bin<-bind_rows(SHAP_bin,ALE_bin) %>%
  mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","non-monotonic","undefined","uncertain")),
         metric=factor(metric,levels=c("SHAP","ALE")),
         threshold=factor(threshold,levels=seq(50,500,50)))

ggplot(SHAP_bin, aes(x=threshold, y=feature, fill=direction)) +
  geom_tile(color="black") +
  scale_fill_manual("Direction",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","non-monotonic"="#FDDDA0",
                                         "uncertain"="grey25","mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
  scale_x_discrete("Bin number",expand=c(0,0))+
  scale_y_discrete("Feature",expand=c(0,0))+
  labs(fill="Direction") +
  theme_Gaia()+
  facet_grid(.~metric)


#- - - -
### Approach 8: Stochastic dominance

## Majority threshold sensitivity

# ALE
thresh_df_sd<-dplyr::bind_rows(lapply(seq(0.50, 0.75, by=0.05), function(thr){
  dplyr::bind_rows(lapply(split(ale_val, ale_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               direction=stochastic_dominance_direction(x=d$feature_val, y=d$ale_val, majority_threshold=thr, q=0.2, min_n=30))
  }))
}))

thresh_df_sd %>% pivot_wider(names_from=threshold,values_from=direction)

ALE_thresh<-thresh_df_sd %>% mutate(metric="ALE")
colnames(ALE_thresh)[3]<-"direction"

# SHAP
thresh_df_sd<-dplyr::bind_rows(lapply(seq(0.50, 0.75, by=0.05), function(thr){
  dplyr::bind_rows(lapply(split(shap_val, shap_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               direction=stochastic_dominance_direction(x=d$feature_val, y=d$shap_val, majority_threshold=thr, q=0.2, min_n=30))
  }))
}))

thresh_df_sd %>% pivot_wider(names_from=threshold,values_from=direction)

SHAP_thresh<-thresh_df_sd %>% mutate(metric="SHAP")
colnames(SHAP_thresh)[3]<-"direction"

SHAP_thresh<-bind_rows(SHAP_thresh,ALE_thresh) %>%
  mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","non-monotonic","undefined","uncertain")),
         metric=factor(metric,levels=c("SHAP","ALE")),
         threshold=factor(threshold,levels=seq(0.50, 0.75, by=0.05)))

ggplot(SHAP_thresh, aes(x=threshold, y=feature, fill=direction)) +
  geom_tile(color="black") +
  scale_fill_manual("Direction",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","non-monotonic"="#FDDDA0",
                                         "uncertain"="grey25","mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
  scale_x_discrete("Majority threshold",expand=c(0,0))+
  scale_y_discrete("Feature",expand=c(0,0))+
  labs(fill="Direction") +
  theme_Gaia()+
  facet_grid(.~metric)

## Quantile threshold sensitivity

# ALE
q_thresh_df<-dplyr::bind_rows(lapply(seq(0.10, 0.80, by=0.05), function(thr){
  dplyr::bind_rows(lapply(split(ale_val, ale_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               direction=stochastic_dominance_direction(x=d$feature_val, y=d$ale_val, majority_threshold=0.55, q=thr, min_n=30))
  }))
}))

q_thresh_df %>% pivot_wider(names_from=threshold,values_from=direction)

ALE_thresh<-q_thresh_df %>% mutate(metric="ALE")
colnames(ALE_thresh)[3]<-"direction"

# SHAP
q_thresh_df<-dplyr::bind_rows(lapply(seq(0.10, 0.80, by=0.05), function(thr){
  dplyr::bind_rows(lapply(split(shap_val, shap_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               direction=stochastic_dominance_direction(x=d$feature_val, y=d$shap_val, majority_threshold=0.55, q=thr, min_n=30))
  }))
}))

q_thresh_df %>% pivot_wider(names_from=threshold,values_from=direction)

SHAP_thresh<-q_thresh_df %>% mutate(metric="SHAP")
colnames(SHAP_thresh)[3]<-"direction"

SHAP_thresh<-bind_rows(SHAP_thresh,ALE_thresh) %>%
  mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","non-monotonic","undefined","uncertain")),
         metric=factor(metric,levels=c("SHAP","ALE")),
         threshold=factor(threshold,levels=seq(0.10, 0.80, by=0.05)))

ggplot(SHAP_thresh, aes(x=threshold, y=feature, fill=direction)) +
  geom_tile(color="black") +
  scale_fill_manual("Direction",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","non-monotonic"="#FDDDA0",
                                         "uncertain"="grey25","mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
  scale_x_discrete("Quantile threshold",expand=c(0,0))+
  scale_y_discrete("Feature",expand=c(0,0))+
  labs(fill="Direction") +
  theme_Gaia()+
  facet_grid(.~metric)

## Minimum observations per group sensitivity

# ALE
n_thresh_df<-dplyr::bind_rows(lapply(seq(5, 50, by=5), function(thr){
  dplyr::bind_rows(lapply(split(ale_val, ale_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               direction=stochastic_dominance_direction(x=d$feature_val, y=d$ale_val, majority_threshold=0.55, q=0.2, min_n=thr))
  }))
}))

n_thresh_df %>% pivot_wider(names_from=threshold,values_from=direction)

ALE_thresh<-n_thresh_df %>% mutate(metric="ALE")
colnames(ALE_thresh)[3]<-"direction"

# SHAP
n_thresh_df<-dplyr::bind_rows(lapply(seq(5, 50, by=5), function(thr){
  dplyr::bind_rows(lapply(split(shap_val, shap_val$feature), function(d){
    data.frame(threshold=thr,
               feature=unique(d$feature),
               direction=stochastic_dominance_direction(x=d$feature_val, y=d$shap_val, majority_threshold=0.55, q=0.2, min_n=thr))
  }))
}))

n_thresh_df %>% pivot_wider(names_from=threshold,values_from=direction)

SHAP_thresh<-n_thresh_df %>% mutate(metric="SHAP")
colnames(SHAP_thresh)[3]<-"direction"

SHAP_thresh<-bind_rows(SHAP_thresh,ALE_thresh) %>%
  mutate(direction=factor(direction,levels=c("promoting","neutral","mitigating","non-monotonic","undefined","uncertain")),
         metric=factor(metric,levels=c("SHAP","ALE")),
         threshold=factor(threshold,levels=seq(5, 50, by=5)))

ggplot(SHAP_thresh, aes(x=threshold, y=feature, fill=direction)) +
  geom_tile(color="black") +
  scale_fill_manual("Direction",values=c("promoting"="#F9CBC2","neutral"="#f7f7f7","non-monotonic"="#FDDDA0",
                                         "uncertain"="grey25","mitigating"="#A6DDCE","undefined"="grey75"),na.value="grey80")+
  scale_x_discrete("Minimum observations per group",expand=c(0,0))+
  scale_y_discrete("Feature",expand=c(0,0))+
  labs(fill="Direction") +
  theme_Gaia()+
  facet_grid(.~metric)

#--------------------------------------------
#### SHAP interaction

# feature-level average of absolute SHAP interaction magnitude
shap_interaction_summary(xgb_fit, X_mat, feature_cols)

# global interaction strength for the model
predictor<-Predictor$new(model=xgb_fit,
                         data=data.frame(X_mat),
                         y=y,
                         predict.function=function(model, newdata){
                           predict(model, xgboost::xgb.DMatrix(as.matrix(newdata)))
                         })

interact<-Interaction$new(predictor)
plot(interact)+theme_Gaia()

# pairwise SHAP interaction strengths
sv_interaction(sv, kind="no")

# SHAP dependence plot
sv_dependence(sv, v=feature_cols, interactions = TRUE)

# SHAP interaction plot
sv_interaction(sv, kind="bar",fill="#40B696")+theme_Gaia()
sv_interaction(sv, kind="beeswarm")+theme_Gaia()
