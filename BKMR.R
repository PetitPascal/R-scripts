#-------------------------------------------------------------------------------
## Reproducible and generalizable BKMR analysis script
# - Simulates a correlated exposure mixture, covariates, and a binary outcome
# - Runs BKMR (kmbayes) with variable selection, knots, and posterior summaries
# - Produces diagnostics and plots
# - Support continuous and binary outcomes
# - Use it as a template for your own dataset by replacing the simulated data
#-------------------------------------------------------------------------------

## -------------------------
## Step 1 - Setup
## -------------------------

## Packages
required_pkgs<-c("bkmr", "fields", "openxlsx", "ggplot2", "corrplot", "MASS", "dplyr", "tidyr",
                 "scales", "gridExtra","gratia","mgcv")
is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(required_pkgs[!is_installed],repos="http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only=TRUE))

## Seed
set.seed(123)  # reproducibility seed (change as needed)

## Helpers

# Small helper to standardize a matrix (center and scale)
standardize_mat<-function(m){
  m<-as.matrix(m)
  apply(m, 2, function(x) if(sd(x, na.rm=TRUE) > 0) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) else x)
}

# Effect direction function based on the mean sign of est
bkmr_dir_mean_sign<-function(z_vec, est_vec) {
  m<- mean(est_vec, na.rm=TRUE)
  if(abs(m)<1e-10) return("uncertain")
  ifelse(m>0, "risk factor", "protective factor")
}

# Effect direction function based on GAM derivative on (z, est) curve
bkmr_dir_gam<-function(z_vec, est_vec){
  df_tmp<-data.frame(z=z_vec, est=est_vec)
  df_tmp<-df_tmp[complete.cases(df_tmp), ]
  if(length(unique(df_tmp$z))<5) return("undefined")
  
  tryCatch({
    gam_m<-mgcv::gam(est ~ s(z, bs="cr", k=8), data=df_tmp, method="REML")
    derivs<-gratia::derivatives(gam_m, term="s(z)")
    mean_d<-mean(derivs$derivative, na.rm=TRUE)
    if (abs(mean_d)<1e-10) return("uncertain")
    return(ifelse(mean_d>0, "risk factor", "protective factor"))
  }, error=function(e){
    rho<-suppressWarnings(cor(df_tmp$z, df_tmp$est, method="spearman"))
    if(is.na(rho) || abs(rho) < 0.05) return("uncertain")
    return(ifelse(rho>0, "risk factor", "protective factor"))
  })
}

# Effect direction function based on pairwise bin concordance

# Compiling C++ function
if (!exists("bin_pairwise_counts")) {
  Rcpp::cppFunction('
Rcpp::List bin_pairwise_counts(NumericVector bx, NumericVector bs,
                               NumericVector bc) {
  int B=bx.size();
  double pos=0.0, neg=0.0, total=0.0;
  for (int i=0; i < B; ++i) {
    for (int j=i+1; j < B; ++j) {
      double ci=bc[i], cj=bc[j];
      if (ci <= 0.0 || cj <= 0.0) continue;
      double pairs=ci * cj;
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
    Rcpp::Named("pos")  =pos,
    Rcpp::Named("neg")  =neg,
    Rcpp::Named("total")=total);
}', depends="Rcpp")
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
    if (length(g0)==0 || length(g1)==0) return("uncertain")
    g0s<-sort(g0)
    pos<-sum(findInterval(g1, g0s, left.open=TRUE))
    tot<-as.double(length(g0)) * as.double(length(g1))
    pp<-pos/tot
    if(pp >= majority_threshold) return("risk factor")
    if(pp <= 1-majority_threshold) return("protective factor")
    return("uncertain")
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
  
  if(length(bx) < 2) return("uncertain")
  
  cnts<-bin_pairwise_counts(bx, bs, bc)
  pos_p<-as.numeric(cnts$pos)
  neg_p<-as.numeric(cnts$neg)
  tot_p<-as.numeric(cnts$total)
  if(tot_p <= 0) return("uncertain")
  
  pp<-pos_p/tot_p
  np<-neg_p/tot_p
  if(pp>=majority_threshold) return("risk factor")
  if(np>=majority_threshold) return("protective factor")
  return("uncertain")
}

bkmr_dir_pairwise<-function(z_vec, est_vec, majority_threshold=0.55,
                            n_quantile_bins=50) {
  
  df_tmp<-data.frame(feature = z_vec, shap = est_vec)
  
  pairwise_direction(df_tmp, "feature", "shap",
                     majority_threshold=majority_threshold,
                     n_quantile_bins=n_quantile_bins)
}

## -------------------------
## Step 2 - Simulating a dataset
#
# - replace this block with code to load your own data
# - or set outcome_type to "gaussian" (continuous) or "binomial" (binary)
#
## -------------------------

outcome_type<-"binomial" # can be set to "gaussian" or "binomial"

n<-500
p_exposures<-10

# Creating correlated exposures using multivariate normal, then exponentiating to resemble concentrations
Sigma<-0.6 ^ as.matrix(dist(1:p_exposures)) # exponential decay correlation
raw_expos<-MASS::mvrnorm(n=n, mu=rep(0, p_exposures), Sigma=Sigma)

# Making exposures positive and on a scale like concentrations
expo_df<-as.data.frame(exp(raw_expos/3)) # scale down variation

# Naming exposures with a prefix 'X_'
colnames(expo_df)<-paste0("X_", seq_len(p_exposures))

# Covariates: continuous age, BMI, and a binary sex variable
cov_df<-data.frame(age_years=round(rnorm(n, mean=50, sd=12)),
  bmi=rnorm(n, mean=27, sd=5),
  sex_male=rbinom(n, 1, p=0.48))

# Creating an outcome (binary) where a non-linear combination of exposures influences probability:
# true underlying function for simulation: h(Z)=sum(w_j * log(X_j+1)) + interaction between X1 and X2
weights<-c(0.4, 0.25, rep(0.05, p_exposures-2))
hZ<-as.numeric(as.matrix(log(expo_df + 1)) %*% weights)
hZ<-hZ + 0.6 * (log(expo_df$X_1 + 1) * log(expo_df$X_2 + 1))  # some interaction

# Linear predictor includes covariate effects
linpred<- -2 + 0.02 * cov_df$age_years + 0.06 * cov_df$bmi + 0.4 * cov_df$sex_male + hZ

# Creating the outcome

if(outcome_type=="binomial"){
  set.seed(123)
  outcome_y<-rbinom(n, 1, plogis(linpred))
}else{
  set.seed(123)
  outcome_y<-linpred + rnorm(n, 0, 1)
}

# Combining into one data.frame
sim_data<-bind_cols(data.frame(id=seq_len(n), outcome=outcome_y), cov_df, expo_df)

## -------------------------
## Step 3 - EDA: correlation matrix of exposures + pairs plot (small)
## -------------------------

expo_mat<-as.matrix(sim_data %>% select(starts_with("X_")))
corrplot::corrplot.mixed(cor(expo_mat, method="spearman"), 
                         number.cex=0.7, tl.cex=0.8, 
                         main="Spearman correlations: exposures")

## optional pairs: (comment out if too large)
# pairs(log(expo_mat + 1), main="Pairwise log-exposure relationships")

## -------------------------
## Step 4 - Preparing matrices for BKMR (standardized exposures)
## -------------------------

expo_std<-standardize_mat(log(expo_mat + 1))  # log + standardization
covs_matrix<-as.matrix(sim_data %>% select(age_years, bmi, sex_male))
outcome_vector<-sim_data$outcome
expo_names<-colnames(expo_std) # names for exposures for plotting labels

## -------------------------
## Step 5 - Determining knots (cover.design) and groups (optional)
#  
#  knots: choose number based on n and complexity (rule of thumb: 25-100)
## -------------------------

# Choosing 50 knots by default (change as needed)
num_knots<-50
knots_coords<-fields::cover.design(expo_std, nd=num_knots)$design

## Demonstrating a possible group structure (example: group some exposures)
# you can create a vector length p_exposures specifying group ids or leave NULL
group_vector<-c(1,1,2,3,3,4,4,5,6,7)  # arbitrary example; change if you have domain knowledge

## -------------------------
## Step 6 - Fitting BKMR model
## -------------------------
## NOTES:
# - family="binomial" for binary outcome
# - iter: set to 5000 for demonstration; consider 10000-30000 for final inference
# - varsel=TRUE enables variable selection (PIPs)
# - est.h=TRUE will estimate h(z) for predictor-response plots
# - control.params can be tuned; r.jump2 controls MCMC mixing for h

set.seed(123)
bkmr_fit<-kmbayes(y=outcome_vector,
  Z=expo_std,
  X=covs_matrix,
  iter=5000,
  family=outcome_type,
  varsel=TRUE,
  knots=knots_coords,
  groups=group_vector,
  est.h=TRUE,
  verbose=TRUE,
  control.params=list(r.jump2=0.5))

## -------------------------
## Step 7 - Posterior inclusion probabilities (PIPs) and variable selection
## -------------------------

## Extracting PIPs
pip_df<-ExtractPIPs(bkmr_fit)
pip_df$PIP<-pip_df$groupPIP*pip_df$condPIP # computing overall variable-level PIP
pip_df

## Estimating direction of each exposure effect using univariate predictor-response
pred_univar<-PredictorResponseUnivar(fit=bkmr_fit)

direction_df<-pred_univar %>%
  group_by(variable) %>%
  summarise(dir_mean_sign=bkmr_dir_mean_sign(z, est),
    dir_gam=bkmr_dir_gam(z, est),
    dir_pairwise=bkmr_dir_pairwise(z, est),
    .groups = "drop") %>%
  # Consensus: majority vote across three methods
  mutate(direction_consensus=apply(cbind(dir_mean_sign, 
                                         dir_gam, dir_pairwise), 1,
      function(x) {
        tbl<-sort(table(x), decreasing = TRUE)
        if(tbl[1]>=2) names(tbl)[1] else "uncertain"
      }))

# Sensitivity check - direction from SingVarRiskSummaries (mixture context)
# Uses the sign of the estimate when other exposures are at their median (q.fixed=0.5)
single_risks_dir<-single_risks %>%
  filter(q.fixed==0.5) %>%
  mutate(dir_mixture=case_when(est>0~"risk factor",
    est<0~"protective factor",
    TRUE~"uncertain")) %>%
  select(variable, dir_mixture)

direction_df<-left_join(direction_df,single_risks_dir,by="variable")

## comparing direction between approaches
dir_long<-direction_df %>%
  dplyr::select(variable, dir_mean_sign, dir_gam, dir_pairwise,
                direction_consensus,dir_mixture) %>%
  tidyr::pivot_longer(-variable, names_to="method", values_to="direction") %>%
  mutate(direction=factor(direction,levels = c("risk factor","uncertain","protective factor","undefined")),
         method=factor(method,
                       levels = c("dir_mean_sign","dir_gam","dir_pairwise",
                                  "direction_consensus","dir_mixture"),
                       labels = c("Mean sign",
                                  "GAM derivative",
                                  "Pairwise bins",
                                  "Consensus",
                                  "Mixture context")))

ggplot(dir_long, aes(x=method, y=variable, fill=direction)) +
  geom_tile(color="white", linewidth=0.8) +
  scale_fill_manual(values = c("protective factor"="#A6DDCE",
                               "uncertain"="#f7f7f7",
                               "risk factor"="#F9CBC2",
                               "undefined"="grey80"),
                    na.value="grey80") +
  labs(title="", x = "", y = "Variable", fill = "Direction") +
  theme_minimal()+
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="top",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

## Visualizing PIPs with direction

# Preparing data for plotting
pip_plot_df<-pip_df %>%
  left_join(direction_df, by="variable")

# Plotting
ggplot(pip_plot_df, aes(x=reorder(variable, PIP), y=PIP, 
                        fill=direction_consensus)) + # other direction can be used
  geom_col(col="black") + 
  geom_text(aes(x=reorder(variable, PIP), y=PIP+0.05,label=as.character(signif(PIP,3))),size=5)+
  coord_flip() + 
  ylab("Posterior Inclusion Probability") +
  xlab("")+
  scale_fill_manual('',values=c("risk factor"="#E87D7D",
    "protective factor"="#A6DDCE",
    "uncertain"="darkgrey",
    "undefined"="lightgrey"))+
  theme_bw()+
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="top",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

## -------------------------
## Step 8 - Univariate predictor-response functions: effect of each exposure holding others at median
## -------------------------

ggplot(pred_univar, aes(z, est, ymin=est - 1.96 * se, ymax=est + 1.96 * se)) +
  geom_ribbon(alpha=0.18) +
  geom_line() +
  geom_hline(yintercept=0,linetype="dashed")+
  facet_wrap(~ variable, scales="free_x") +
  labs(x="Exposure (standardized)", y="Estimated h(z)", title="Univariate predictor-response estimates") +
  theme_bw()+
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="bottom",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

## -------------------------
## Step 9 - Overall mixture risk summary
## -------------------------

qs_seq<-seq(0.25, 0.75, by=0.05)
overall_risks<-OverallRiskSummaries(fit=bkmr_fit, qs=qs_seq, q.fixed=0.5)

## Testing whether the 95% CI at q=0.75 vs. q=0.25 exclude 0?

# Extracting h samples at q=0.75 and q=0.25 (joint shift)
sel<-seq(2501, 5000, by=5) # post burn-in thinned samples

h_q25<-ComputePostmeanHnew(fit=bkmr_fit,
                           Znew=matrix(apply(expo_std, 2, quantile, 0.25), nrow=1),
                           sel=sel)

h_q75<-ComputePostmeanHnew(fit=bkmr_fit,
  Znew=matrix(apply(expo_std, 2, quantile, 0.75), nrow=1),
  sel=sel)

h_samples_q25<-SamplePred(fit=bkmr_fit,
                          Znew=matrix(apply(expo_std, 2, quantile, 0.25), nrow = 1),
                          Xnew=matrix(colMeans(covs_matrix), nrow = 1),
                          sel= sel)

h_samples_q75<-SamplePred(fit=bkmr_fit,
  Znew=matrix(apply(expo_std, 2, quantile, 0.75), nrow = 1),
  Xnew=matrix(colMeans(covs_matrix), nrow = 1),
  sel=sel)

# contrast = h(Q75) - h(Q25) across posterior samples
mixture_contrast<-as.numeric(h_samples_q75)-as.numeric(h_samples_q25)
contrast_mean<-mean(mixture_contrast)
contrast_ci<-quantile(mixture_contrast, c(0.025, 0.975))
contrast_mean
contrast_ci

# Posterior probability that mixture effect > 0 (risk factor direction)
pp_positive<-mean(mixture_contrast>0)
pp_negative<-mean(mixture_contrast<0)
pp_positive
pp_negative

# Plotting overall mixture risk
ggplot(overall_risks, aes(quantile, est, ymin=est-1.96*sd, ymax=est+1.96*sd)) +
  geom_pointrange() +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  labs(title="Overall mixture risk: change from median (q.fixed=0.5)",
       x="Quantile of joint change", y="Estimated change in h(z)") +
  theme_minimal()+
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="bottom",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

ggplot(overall_risks, aes(quantile, est,
                          ymin=est-1.96*sd,
                          ymax=est+1.96*sd)) +
  geom_ribbon(alpha=0.2) +
  geom_line(linewidth=1) +
  geom_point(size=2) +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  annotate("text", x=min(qs_seq) + 0.02, y=max(overall_risks$est),
           label=paste0("PP(effect>0)=", round(pp_positive, 3)),
           hjust=0, size=4, color="darkred") +
  labs(title="Overall mixture risk (Q75 vs Q25 of joint shift)",
       x="Quantile of joint change", y="Estimated change in h(z)") +
  theme_bw(base_size=14)

## -------------------------
## Step 10 - Single-exposure risk summaries (pairwise comparisons)
## -------------------------

single_risks<-SingVarRiskSummaries(fit=bkmr_fit,
  y=outcome_vector,
  Z=expo_std,
  X=covs_matrix,
  qs.diff=c(0.25, 0.75),        # compare these quantiles
  q.fixed=c(0.25, 0.50, 0.75))  # settings to fix other exposures

# Ploting single-exposure risks
ggplot(single_risks, aes(x=variable, y=est,
                         ymin=est-1.96*sd,
                         ymax=est+1.96*sd,
                         color=factor(q.fixed))) +
  geom_pointrange(position=position_dodge(width=0.7)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_brewer(palette="Set2") +
  labs(title="Single-exposure risk summaries (others fixed)",
       color="q.fixed", x="", y="Estimated change in h(z)") +
  theme_bw(base_size=14) +
  theme(strip.text.x=element_text(size=16, colour="black", angle=0),
        strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold"),
        legend.position="bottom",
        axis.line=element_line(color="black",size=0.1, linetype="solid"))

# ## -------------------------
# ## Step 11 - Interaction plot: bivariate predictor-response (i.e., effect of X_i at different levels of X_j)
# ## -------------------------
# 
# # Choosing pairs to investigate (either from domain knowledge or highest PIPs)
# top_vars<-pip_plot_df %>% arrange(desc(PIP)) %>% pull(variable)
# pairs_to_plot<-list(c(top_vars[1], top_vars[2]),
#                     c(top_vars[1], top_vars[3]))  # add more pairs as needed
# 
# for(pair in pairs_to_plot){
#   var1<-pair[1]
#   var2<-pair[2]
#   
#   pred_bivar<-PredictorResponseBivar(fit=bkmr_fit,min.plot.dist=0.5,sel=sel)
#   
#   # Filtering for the pair of interest
#   bivar_sub<-pred_bivar %>% filter(variable1==var1, variable2==var2)
#   
#   if(nrow(bivar_sub)==0) next
#   
#   # Slicing at low/medium/high levels of var2
#   z2_quantiles<-quantile(bivar_sub$z2, probs=c(0.1, 0.5, 0.9))
#   bivar_sub$z2_group<-cut(bivar_sub$z2,
#                             breaks=c(-Inf, z2_quantiles, Inf),
#                             labels=c("Q10", "Q25-Q50", "Q50-Q90", "Q90"),
#                             include.lowest=TRUE)
#   
#   p<-ggplot(bivar_sub, aes(x=z1, y=est,
#                              ymin=est-1.96*se,
#                              ymax=est+1.96*se,
#                              color=z2_group, fill=z2_group)) +
#     geom_ribbon(alpha=0.15, color=NA) +
#     geom_line(linewidth=1) +
#     scale_color_brewer(palette="RdYlBu", direction=-1) +
#     scale_fill_brewer(palette="RdYlBu", direction=-1) +
#     labs(title=paste("Interaction:", var1, "×", var2),
#       subtitle=paste("Effect of", var1, "at different levels of", var2),
#       x=paste(var1, "(standardized)"),
#       y="Estimated h(z)",
#       color=paste(var2, "level"),
#       fill=paste(var2, "level")) +
#     theme_bw(base_size=14)+
#     theme(strip.text.x=element_text(size=16, colour="black", angle=0),
#           strip.background=element_rect(fill="#A6DDCE", colour="black", size=1),
#           axis.text.y=element_text(size=16,color="black"),
#           axis.text.x=element_text(size=16,color="black"),
#           axis.title=element_text(size=16,face="bold",color="black"),
#           legend.text=element_text(size=16, face="bold"),
#           legend.title=element_text(size=16, face="bold"),
#           legend.position="bottom",
#           axis.line=element_line(color="black",size=0.1, linetype="solid"))
#   
#   print(p)
# }

## -------------------------
## Step 12 - Saving the model
## -------------------------

# # Basic trace plot for kernel parameters if available
# if(!is.null(bkmr_fit$chain$h)){
#   r_df<-as.data.frame(bkmr_fit$chain$r)
#   colnames(r_df)<-expo_names
#   r_long<-r_df %>%
#     mutate(iter=row_number()) %>%
#     pivot_longe(-iter, names_to="exposure", values_to="r")
#   
#   ggplot(r_long, aes(x=iter, y=r)) +
#     geom_line(alpha=0.6) +
#     facet_wrap(~ exposure, scales="free_y") +
#     labs(title="Trace plots: kernel bandwidth parameters (r)", x="Iteration", y="r") +
#     theme_bw(base_size=11)
# }

# Saving the fitted object for later use
saveRDS(bkmr_fit, file="bkmr_fit_example.rds")
message("Saved fitted BKMR object to: bkmr_fit_example.rds")

## -------------------------
## Notes & recommendations
## -------------------------
cat("Notes:
- This script simulates data for demonstration. To use your own data, replace the simulation block
  and set 'expo_std', 'covs_matrix', and 'outcome_vector' accordingly.
- For final inference increase 'iter' (e.g., 10000-30000) and check mixing (trace plots) and effective
  sample sizes. MCMC heavy computations can take many minutes to hours depending on data and iter.
- Consider sensitivity analyses: different knot numbers, different priors, or excluding influential observations.
- If you have domain knowledge for grouping exposures (e.g., chemical classes), use the 'groups' argument.")
