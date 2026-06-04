#-------------------------------------------------------------------------------
## Reproducible & generalisable SHAP direction script
# Covers:
#   - Conventional approach (mean SHAP sign)
#   - Approach 1: GAM derivative (+ Spearman fallback)
#   - Approach 2: pairwise bin concordance
#   - Approach 3: consensus of all approaches (majority voting)
#   - Simulated dataset with monotonic, U-shaped, binary, threshold, and noise features
#   - Comparison across all three approaches
#   - Visualization: SHAP dependence plots, direction heatmap
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs<-c("dplyr", "ggplot2", "tidyr", "tibble","mgcv", "gratia", "patchwork","xgboost", "shapr","iml", "Rcpp")

is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(required_pkgs[!is_installed],repos = "http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

#--------------------------------------------
## Step 2: Simulating a dataset with known feature behaviors
n<-800

# Feature types with known true directions
x_pos<-rnorm(n, 5, 2) # monotone positive
x_neg<-rnorm(n, 5, 2) # monotone negative
x_ushape<-rnorm(n, 0, 2) # U-shaped (non-monotone)
x_thresh<-rnorm(n, 0, 1) # threshold effect
x_bin<-rbinom(n, 1, 0.4) # binary
x_noise<-rnorm(n, 0, 1) # no effect

sim_df<-data.frame(x_pos, x_neg, x_ushape, x_thresh, x_bin, x_noise)

# True outcome
y<-1.5 * x_pos -1.2 * x_neg +
  0.8 * x_ushape^2 + # U-shape: ambiguous mean sign
  1.5 * (x_thresh > 0.8) + # threshold
  1.0 * x_bin + rnorm(n, 0, 1.5)

sim_df$y<-y

feature_cols<-c("x_pos","x_neg","x_ushape","x_thresh","x_bin","x_noise")

#--------------------------------------------
## Step 3: Fitting XGBoost and computing SHAP values

X_mat<-as.matrix(sim_df[, feature_cols])
dtrain<-xgboost::xgb.DMatrix(X_mat, label = y)

xgb_fit<-xgboost::xgb.train(
  params=list(objective = "reg:squarederror",
    max_depth=4,
    learning_rate=0.05,
    subsample=0.8),
  data=dtrain,
  nrounds=200,
  verbose=0)

# SHAP values using xgboost's built-in SHAP
shap_mat<-predict(xgb_fit, dtrain, predcontrib = TRUE)
shap_df<-as.data.frame(shap_mat[, feature_cols])  # drop BIAS column
colnames(shap_df)<-paste0("shap_", feature_cols)

# Combining
full_df<-cbind(sim_df[, feature_cols], shap_df)

cat("Mean SHAP per feature:\n")
round(colMeans(shap_df), 4)

#--------------------------------------------
## Step 4: SHAP dependence plots

plot_shap_dep<-function(df, feature, shap_col, title = NULL) {
  ggplot(df, aes(x = .data[[feature]], y = .data[[shap_col]])) +
    geom_point(alpha = 0.3, color = "#2166ac") +
    geom_smooth(method = "loess", se = TRUE,color = "#d6604d", linewidth = 1.1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = title %||% paste("SHAP dependence:", feature),
         x = feature, y = "SHAP value") +
    theme_bw(base_size = 13)
}

dep_plots<-lapply(feature_cols, function(f)
  plot_shap_dep(full_df, f, paste0("shap_", f)))

patchwork::wrap_plots(dep_plots, ncol = 3) +
patchwork::plot_annotation(title = "SHAP dependence plots")

#--------------------------------------------
## Determining the direction

# Approach 1 - conventional (mean SHAP sign)

conventional_direction<-function(df, feature_col, shap_col){
  s<-df[[shap_col]]
  mean_s<-mean(s, na.rm = TRUE)
  if(abs(mean_s)< 1e-10) return("neutral")
  ifelse(mean_s>0, "promoting", "mitigating")
}

# Approach 2 - GAM derivative (+ Spearman fallback)

gam_direction<-function(df, feature_col, shap_col){
  x<-df[[feature_col]]
  s<-df[[shap_col]]
  
  valid<-complete.cases(x, s)
  x<-x[valid]
  s<-s[valid]
  
  if (length(unique(x)) <= 1 || length(unique(s)) <= 1) return("undefined")
  
  # Binary or sparse
  if (length(unique(x)) <= 2 || quantile(x, 0.75, na.rm=TRUE) == 0) {
    mean_diff <- mean(s[x > 0], na.rm=TRUE) - mean(s[x <= 0], na.rm=TRUE)
    return(ifelse(mean_diff > 0, "promoting",
                  ifelse(mean_diff < 0, "mitigating", "neutral")))
  }
  
  # Continuous: GAM
  tryCatch({
    gam_m  <- mgcv::gam(s ~ s(x, bs="cr", k=8), method="REML")
    derivs <- gratia::derivatives(gam_m, term="s(x)")
    mean_d <- mean(derivs$derivative, na.rm=TRUE)
    if (abs(mean_d) < 1e-10) return("neutral")
    return(ifelse(mean_d > 0, "promoting", "mitigating"))
  }, error = function(e) {
    rho <- suppressWarnings(cor(x, s, method="spearman"))
    if (is.na(rho) || abs(rho) < 0.05) return("neutral")
    return(ifelse(rho > 0, "promoting", "mitigating"))
  })
}

# Approach 3- Pairwise bin concordance (Rcpp)

# Compile C++ function
if (!exists("bin_pairwise_counts")) {
  Rcpp::cppFunction('
Rcpp::List bin_pairwise_counts(NumericVector bx, NumericVector bs,
                               NumericVector bc) {
  int B = bx.size();
  double pos = 0.0, neg = 0.0, total = 0.0;
  for (int i = 0; i < B; ++i) {
    for (int j = i+1; j < B; ++j) {
      double ci = bc[i], cj = bc[j];
      if (ci <= 0.0 || cj <= 0.0) continue;
      double pairs = ci * cj;
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
    Rcpp::Named("pos")   = pos,
    Rcpp::Named("neg")   = neg,
    Rcpp::Named("total") = total);
}', depends = "Rcpp")
}

pairwise_direction<-function(df, feature_col, shap_col,majority_threshold = 0.55,n_quantile_bins=200) {
  x<-df[[feature_col]]
  s<-df[[shap_col]]
  valid<-complete.cases(x, s)
  x<-x[valid]
  s<-s[valid]
  
  if(length(unique(x)) <= 1 || length(unique(s)) <= 1) return("undefined")
  
  # Binary / sparse
  if(length(unique(x)) <= 2 ||
      quantile(x, 0.75, na.rm=TRUE) == 0) {
    lv <- min(x); hv <- max(x)
    g0 <- s[x == lv]; g1 <- s[x == hv]
    if (length(g0)==0 || length(g1)==0) return("neutral")
    g0s<-sort(g0)
    pos<-sum(findInterval(g1, g0s, left.open=TRUE))
    tot<-as.double(length(g0)) * as.double(length(g1))
    pp<-pos/tot
    if(pp >= majority_threshold) return("promoting")
    if(pp <= 1-majority_threshold) return("mitigating")
    return("neutral")
  }
  
  # Continuous: binning
  probs<-seq(0, 1, length.out = n_quantile_bins+1)
  breaks<-unique(quantile(x, probs=probs, na.rm=TRUE, type=7))
  if (length(breaks) <= 2)
    breaks<-seq(min(x,na.rm=TRUE), max(x,na.rm=TRUE), length.out=3)
  bins<-cut(x, breaks=breaks, include.lowest=TRUE)
  
  bx<-as.numeric(tapply(x, bins, mean, na.rm=TRUE))
  bs<-as.numeric(tapply(s, bins, mean, na.rm=TRUE))
  bc<-as.numeric(tapply(s, bins, length))
  ok<-!is.na(bx) & !is.na(bs) & !is.na(bc)
  bx<-bx[ok]; bs <- bs[ok]; bc <- bc[ok]
  
  if(length(bx) < 2) return("neutral")
  
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

#--------------------------------------------
## Applying all three approaches to all features

results<-lapply(feature_cols, function(f){
  sc<-paste0("shap_", f)
  conv<-conventional_direction(full_df, f, sc)
  gam<-gam_direction(full_df, f, sc)
  pair<-pairwise_direction(full_df, f, sc)
  
  # Numeric SHAP summaries for context
  mean_shap<-round(mean(full_df[[sc]], na.rm=TRUE), 4)
  median_shap<-round(median(full_df[[sc]], na.rm=TRUE), 4)
  spearman_r<-round(cor(full_df[[f]], full_df[[sc]], method="spearman", use="complete.obs"), 4)
  
  data.frame(feature=f,
             mean_shap=mean_shap,
             median_shap=median_shap,
             spearman_r=spearman_r,
             conventional=conv,
             gam_deriv=gam,
             pairwise_bins=pair)
})

direction_df<-do.call(rbind, results) %>% 
  # Consensus: majority vote across three methods
  mutate(consensus=apply(cbind(conventional,gam_deriv, pairwise_bins), 1,
                         function(x) {
                           tbl<-sort(table(x), decreasing = TRUE)
                           if(tbl[1]>=2) names(tbl)[1] else "uncertain"
                         }))
direction_df

# True expected directions for comparison
truth<-data.frame(feature=feature_cols,
  true_direction=c("promoting","mitigating","neutral","promoting","promoting","neutral"),
  note=c("monotone +","monotone -","U-shaped: ambiguous mean","threshold","binary +","noise"))
left_join(direction_df,truth,by="feature")

#--------------------------------------------
## Visualisation 1 — Direction heatmap

dir_long<-direction_df %>%
  dplyr::select(feature, conventional, gam_deriv, pairwise_bins,consensus) %>%
  tidyr::pivot_longer(-feature, names_to="method", values_to="direction") %>%
  mutate(direction=factor(direction,levels = c("promoting","neutral","mitigating","undefined")),
    method=factor(method,
                       levels = c("conventional","gam_deriv","pairwise_bins","consensus"),
                       labels = c("Conventional\n(mean sign)",
                                  "Approach 1\n(GAM derivative)",
                                  "Approach 2\n(pairwise bins)",
                                  "Approach 3\n(consensus)")))

ggplot(dir_long, aes(x=method, y=feature, fill=direction)) +
  geom_tile(color="white", linewidth=0.8) +
  scale_fill_manual(values = c("promoting"="#A6DDCE",
               "neutral"="#f7f7f7",
               "mitigating"="#F9CBC2",
               "undefined"="grey80"),
    na.value="grey80") +
  labs(title="", x = "", y = "Feature", fill = "Direction") +
  theme_bw(base_size=14) +
  theme(axis.text=element_text(color="black"))

# #--------------------------------------------
# ## Visualisation 2 — Detailed dependence with GAM curve
# 
# plot_shap_gam<-function(df, feature, shap_col){
#   x<-df[[feature]]
#   s<-df[[shap_col]]
#   gd<-tryCatch({
#     gm<-mgcv::gam(s ~ s(x, bs="cr", k=8), method="REML")
#     pf<-data.frame(x=seq(min(x),max(x),length.out=200))
#     pp<-predict(gm, newdata=pf, se.fit=TRUE)
#     data.frame(x=pf$x, fit=pp$fit,
#                lo=pp$fit-1.96*pp$se.fit,
#                hi=pp$fit+1.96*pp$se.fit)
#   }, error=function(e) NULL)
#   
#   p<-ggplot(data.frame(x=x,s=s), aes(x=x,y=s)) +
#     geom_point(alpha=0.25, color="#2166ac", size=1) +
#     geom_hline(yintercept=0, linetype="dashed", color="grey50")
#   
#   if(!is.null(gd))
#     p<-p +
#     geom_ribbon(data=gd, aes(x=x, ymin=lo, ymax=hi), inherit.aes=FALSE,
#                 fill="#d6604d", alpha=0.2) +
#     geom_line(data=gd, aes(x=x, y=fit), inherit.aes=FALSE,
#               color="#d6604d", linewidth=1.2)
#   
#   p + labs(title=paste("SHAP dependence + GAM:", feature),
#            x=feature, y="SHAP") +
#     theme_bw(base_size=13)
# }
# 
# gam_plots<-lapply(feature_cols, function(f)
#   plot_shap_gam(full_df, f, paste0("shap_",f)))
# 
# patchwork::wrap_plots(gam_plots, ncol=3) + patchwork::plot_annotation(title="SHAP dependence plots with GAM smooth (Approach 1)")

#--------------------------------------------
## Sensitivity analysis: threshold sensitivity for approach 2

thresh_vals<-seq(0.50, 0.75, by=0.05)
thresh_res<-lapply(thresh_vals, function(thr) {
  dirs<-sapply(feature_cols, function(f)
    pairwise_direction(full_df, f, paste0("shap_",f),majority_threshold=thr))
  data.frame(threshold=thr, t(dirs))
})
thresh_df<-do.call(rbind, thresh_res)
thresh_df

#--------------------------------------------
## Sensitivity analysis: bin number sensitivity for approach 2

bin_vals<-c(50, 100, 200, 500)
bin_res<-lapply(bin_vals, function(nb){
  dirs<-sapply(feature_cols, function(f)
    pairwise_direction(full_df, f, paste0("shap_",f),n_quantile_bins=nb))
  data.frame(n_bins=nb, t(dirs))
})
bin_df<-do.call(rbind, bin_res)
bin_df

#--------------------------------------------
## Reusable wrapper: computing directions for all features

compute_shap_directions<-function(data_df, feature_cols,shap_prefix="shap_",
                                    methods=c("conventional","gam","pairwise"),
                                    threshold=0.55,n_bins=200) {
  
  res<-lapply(feature_cols, function(f){
    sc<-paste0(shap_prefix, f)
    if(!sc %in% names(data_df)){
      message("SHAP column not found for: ", f); return(NULL)
    }
    
    row<-data.frame(feature=f,
                      mean_shap=round(mean(data_df[[sc]], na.rm=TRUE), 5),
                      median_shap=round(median(data_df[[sc]],na.rm=TRUE),5))
    
    if ("conventional" %in% methods)
      row$conventional<-conventional_direction(data_df, f, sc)
    if ("gam" %in% methods)
      row$gam_deriv<-gam_direction(data_df, f, sc)
    if ("pairwise" %in% methods)
      row$pairwise_bins<-pairwise_direction(data_df, f, sc, majority_threshold=threshold,n_quantile_bins=n_bins)
    row
  })
  
  do.call(rbind, Filter(Negate(is.null), res)) %>%
    # Consensus: majority vote across three methods
    mutate(consensus=apply(cbind(conventional, gam_deriv, pairwise_bins), 1,
                           function(x) {
                             tbl<-sort(table(x), decreasing = TRUE)
                             if(tbl[1]>=2) names(tbl)[1] else "uncertain"
                           }))
}

# Applying the pipeline
direction_results<-compute_shap_directions(data_df=full_df,
  feature_cols=feature_cols,
  methods=c("conventional","gam","pairwise"),
  threshold=0.55,
  n_bins=200)

direction_results
