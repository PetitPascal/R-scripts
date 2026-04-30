# ==============================================================================
# Main script information
# ==============================================================================

# Description: A generic, reproducible Cox PH regression script covering:
#              - Data import and preprocessing
#              - Model fitting and characterization
#              - Full assumption testing:
#                  (1) Independence assumption (cluster-robust SE)
#                  (2) Proportional hazards assumption (Schoenfeld residuals)
#                  (3) Linearity assumption (Martingale residuals / restricted cubic splines)
#                  (4) No multicollinearity (VIF)
#                  (5) No effect modification (interaction tests, e.g., by sex)
#              - PH violation remediation:
#                  Approach 1: covariate × time interaction terms
#                  Approach 2: time-dependent covariates (counting-process format, with optional timeSplitter from Greg package)
#              - Competing risks analysis (Fine-Gray subdistribution hazard model)
#              - Interaction / effect modification test (e.g., sex)
#              - Visualization (forest plots, KM curves, cumulative hazard/events)
#
# NOTE: All user-defined parameters are clearly flagged with <<USER>> comments.
#       To adapt this script to your own data, search for <<USER>> and edit only those sections.

# ==============================================================================
#### Configurations ####
# ==============================================================================

## Disabling memory torture
gctorture(FALSE)

## Installing packages if needed
pack_needed<-c("data.table","tidyverse","survival","survminer","conflicted","kdry","broom","here","ggsurvfit","ggkm",
               "gtsummary","GGally","cmprsk","tidycmprsk","car","rms","Hmisc","Greg","patchwork","ggstats")

for (i in 1:length(pack_needed)){
  if(pack_needed[i]%in%.packages(all.available=TRUE)){
  }else{
    install.packages(pack_needed[i])
  }
}

## Loading packages
library(survival)
library(survminer)
library(cmprsk)
library(tidycmprsk)
library(car)
library(rms)
library(Hmisc)
library(Greg)
library(data.table)
library(tidyverse)
library(broom)
library(conflicted)
library(ggsurvfit)
library(kdry)
library(here)
library(ggkm)
library(gtsummary)
library(GGally)
library(patchwork)

## Preventing package conflicts
conflict_prefer("select", "dplyr") 
conflict_prefer("filter", "dplyr") 
conflict_prefer("slice", "dplyr")
conflict_prefer("alpha", "scales")
conflict_prefer("summarize", "dplyr")

## Setting the working directory
here::here("Cox proportional hazards regression analysis")

## Output folder - <<USER>> Change the output folder path if needed
output_dir<-"Results"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
safe_path<-function(filename) file.path(output_dir, filename)

# ==============================================================================
#### Creating functions ####
# ==============================================================================

## Function for formatting the display of summary statistics
test_format<-function(x){
  x<-as.numeric(x)
  sign_x<-if_else(x < 0, "neg", "pos")
  x<-abs(x)
  x_raw<-x
  
  if(is.na(x)||is.infinite(x)) x_raw<-0
  
  if(x_raw>=100){
    virg_pos<-str_locate(as.character(x_raw), "[.]")[1]
    if (!is.na(virg_pos) && as.numeric(substr(x_raw, virg_pos + 1, virg_pos + 1)) >= 5){
      x<-x+1
      x<-as.numeric(substr(x, 1, virg_pos - 1))
    }
  }
  
  if(is.na(x) || is.infinite(x)){
    x<-""
  }else if (x < 0.05 || x >= 10000){
    x<-format(signif(x, 3), scientific = TRUE)
    if(nchar(x) == 8 && substr(x, 4, 4) != "e"){
      if(as.numeric(substr(x, 4, 4)) >= 5 && x < 0.01){
        x1<-as.numeric(substr(x, 1, 4)) + 0.1
        x2<-substr(x, 5, 8)
        x<-paste(substr(x1, 1, 4), x2, sep = "")
      }else{
        x<-paste(substr(x, 1, 4), substr(x, 5, 8), sep = "")
      }
    }
  }else{
    x_save<-x
    x<-signif(x, 3)
    if(nchar(x)==6) x<-as.numeric(substr(x, 1, 5))
    if(nchar(x)==5){
      x<-if(as.numeric(substr(x, 5, 5)) >= 5) substr(x + 0.01, 1, 4) else substr(x, 1, 4)
    }else{
      x<-if(x>= 1000) as.character(signif(x_save,4)) else as.character(x)
    }
  }

  if(sign_x == "neg" && x_raw != 0) x<-paste("-", x, sep = "")
  if(x == "0e+00") x<-"0"
  return(x)
}

## Function returning a tidy tibble from car::vif() regardless of whether GVIF or VIF is returned
tidy_vif<-function(model,VIF_set){
  v<-car::vif(model)
  if(is.matrix(v)){
    # GVIF form (factors with >2 levels)
    as_tibble(v, rownames = "term") %>%
      rename(GVIF = 1, Df = 2, `GVIF^(1/(2*Df))` = 3) %>%
      mutate(VIF_flag = `GVIF^(1/(2*Df))` > sqrt(VIF_set))
  }else{
    tibble(term = names(v), VIF = v) %>%
      mutate(VIF_flag = VIF > VIF_set)
  }
}

## Function for generating a Schoenfeld-based PH summary
tidy_ph_test<-function(ph_test_obj) {
  as_tibble(data.frame(term = rownames(ph_test_obj$table),
                       ph_test_obj$table)) %>%
    rename_with(~ c("term", "chisq", "df", "p_value"), everything()) %>%
    mutate(PH_violated = p_value < 0.05)
}

## Function for customizing plots
theme_cox<-function(){
  theme_bw() +
    theme(strip.text=element_text(size=12, colour="black", face="bold"),
          strip.background=element_rect(fill="#CAE1FF",colour="black"),
          axis.text=element_text(size=12, color="black"),
          axis.title=element_text(size=12, face="bold", color="black"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12, face="bold"),
          axis.line=element_line(color="black", linewidth=0.1))
}

## Function for safely creating filename
to_filename<-function(label){
  label %>%
    str_replace_all("[^[:alnum:] ]", "") %>% # remove special characters
    str_replace_all("\\s+", "_") %>%         # spaces to underscores
    str_trunc(60, ellipsis = "")             # cap length for OS safety
}

# ==============================================================================
#### Data import and preprocessing ####
# ==============================================================================

#----------------------------------------------------------------
## <<USER>> Replace the block below with your own data import
##          Ensure your final data frame is named `clean_data` and contains:
##            - time_var: the follow-up time variable (numeric)
##            - event_var: the event indicator (0 = censored, 1 = event)
##            - covariates: all predictor variables
##            - (optional): competing event indicator (0/1/2, etc.)
##            - (optional): a cluster/frailty variable if independence is violated
#----------------------------------------------------------------

## Importing the clean dataset - example: survival::colon dataset
clean_data<-as_tibble(na.omit(survival::colon) %>% filter(etype == 2) %>% select(-c(etype,study)) %>%
                      mutate(rx=case_when(rx=="Obs"~0,rx=="Lev"~1,T~2)))

## Renaming outcome variables - <<USER>> Adjust variable names to match your dataset
clean_data<-clean_data %>%
            rename(Disease_status=status,
                   time_to_diagnosis=time)

## Analysis setting - <<USER>> Set all key parameters here

# Outcome
time_var<-"time_to_diagnosis"   # name of the time variable
event_var<-"Disease_status"      # name of the event indicator (1 = event)

# Covariates to include in the main Cox model
covariates<-c("sex", "age", "rx", "nodes")

# Effect modification / interaction variable (set to NULL to skip)
interaction_var<-"sex"

# Reference level label for the interaction variable (for plots)
interaction_var_labels<-c("0"="Female", "1"="Male")

# Competing event indicator variable name (set to NULL if no competing risks), it should be a factor: 0=censored, 1=main event, 2=competing event
competing_risk_var<-NULL

# Cluster variable for shared frailty / non-independence (set to NULL if inapplicable)
cluster_var<-NULL

# Observation ID (set to NULL if not existing, or name of the corresponding column/variable otherwise)
ID_obs<-"id"

# Time point for survival probability estimation (in the same unit as time_var)
X_time<-c(150,365.25)  # e.g., 1 year when time is in days

# Alpha-risk level for PH violation detection
alpha_ph<-0.05

# VIF cutoff value to identify problematic collinearity
VIF_set<-10

# Whether to use Greg::timeSplitter
use_time_splitter<-TRUE

# Number of knots for restricted cubic splines (linearity check)
rcs_knots<-4

# ==============================================================================
#### Descritive statistics
# ==============================================================================

desc_table<-clean_data %>%
            select(all_of(c(time_var, event_var, covariates))) %>%
            tbl_summary(by= all_of(event_var),
                        missing="no",
                        statistic=list(all_continuous()  ~ "{mean} ({sd}) / median {median} [{p25}; {p75}]",
                                       all_categorical() ~ "{n} ({p}%)")) %>%
            add_overall() %>%
            add_p() %>%
            bold_labels()

fwrite(as_tibble(as.data.frame(desc_table)),
       safe_path(paste0("Descriptive_statistics_", Sys.Date(), ".csv")),
       sep = ";")

# ==============================================================================
#### Main Cox proportional hazards regression model
# ==============================================================================

## Building the survival formula
surv_formula_str<-paste0("Surv(", time_var, ", ", event_var, ") ~ ",paste(covariates, collapse = " + "))
surv_formula<-as.formula(surv_formula_str)

## Fitting the model
if(!is.null(cluster_var)){ # if a cluster variable is provided, use cluster() for robust standard errors
  surv_formula_cluster<-update(surv_formula,as.formula(paste(". ~ . + cluster(", cluster_var, ")")))
  cox_mod<-coxph(surv_formula_cluster, data = clean_data, ties = "efron")
}else{
  cox_mod<-coxph(surv_formula, data = clean_data, ties = "efron")
}

summary(cox_mod)

## Tidy model characteristics
cox_mod_charact<-broom::tidy(cox_mod, exponentiate = TRUE, conf.int = TRUE) %>%
                  mutate(direction = case_when(conf.low  > 1 ~ "positive association",
                                               conf.high < 1 ~ "negative association",
                                               T~ "uncertain association")) %>%
                  bind_cols(broom::glance(cox_mod)) %>%
                  rowwise() %>%
                  mutate(across(where(is.numeric), test_format)) %>%
                  ungroup() %>%
                  mutate(HR = paste0(estimate, " [", conf.low, "; ", conf.high, "]")) %>%
                  mutate(HR = if_else(HR == " [; ]", "", HR))

## Survival probability at X_time
X_time_surv<-summary(survfit(cox_mod), times=X_time)
X_time_surv_df<-tibble(surv=X_time_surv$surv,
                       lower=X_time_surv$lower,
                       upper=X_time_surv$upper) %>%
                rowwise() %>%
                mutate(across(everything(), test_format)) %>%
                ungroup() %>%
                mutate(X_time_result = paste0(surv, " (95% CI: ", lower, "; ", upper, ")"),
                       times=X_time) %>%
                select(times,X_time_result)

## Estimating the median survival time
survfit(cox_mod)

## Likelihood-ratio test (LRT) for each term
lrt_results <- drop1(cox_mod, test = "Chisq")
lrt_results

## Saving results
fwrite(cox_mod_charact,safe_path(paste0("Model_characteristics_", Sys.Date(), ".csv")),sep = ";")
fwrite(X_time_surv_df,safe_path(paste0("Survival_at_X_time_", Sys.Date(), ".csv")),sep = ";")
fwrite(as_tibble(data.frame(lrt_results), rownames = "term"),safe_path(paste0("LR test_", Sys.Date(), ".csv")),sep = ";")

# ==============================================================================
#### Residual diagnostics
# ==============================================================================

## Deviance residuals (should be ~N(0,1) under correct model)
dev_resi<-ggcoxdiagnostics(cox_mod, type = "deviance",   ggtheme = theme_bw())
dev_resi

ggsave(dev_resi, file = safe_path(paste0("Deviance residuals_", Sys.Date(), ".png")),
       dpi = 600, width = 20, height = 10, units = "cm", limitsize = FALSE)

## Martingale residuals (non-linearity check)
marti_desi<-ggcoxdiagnostics(cox_mod, type = "martingale", ggtheme = theme_bw())
marti_desi

ggsave(marti_desi, file = safe_path(paste0("Martingale residuals_", Sys.Date(), ".png")),
       dpi = 600, width = 20, height = 10, units = "cm", limitsize = FALSE)

## Score residuals
score_resi<-ggcoxdiagnostics(cox_mod, type = "score",    ggtheme = theme_bw())
score_resi

ggsave(score_resi, file = safe_path(paste0("Score residuals_", Sys.Date(), ".png")),
       dpi = 600, width = 20, height = 10, units = "cm", limitsize = FALSE)

## dfbeta (influence on coefficient vector)
dfbeta_resi<-ggcoxdiagnostics(cox_mod, type = "dfbeta",  ggtheme = theme_bw())
dfbeta_resi

ggsave(dfbeta_resi, file = safe_path(paste0("dfbeta residuals_", Sys.Date(), ".png")),
       dpi = 600, width = 20, height = 10, units = "cm", limitsize = FALSE)

## dfbetas (dfbeta scaled)
dfbetas_resi<-ggcoxdiagnostics(cox_mod, type = "dfbetas", ggtheme = theme_bw())
dfbetas_resi

ggsave(dfbetas_resi, file = safe_path(paste0("dfbeta scaled residuals_", Sys.Date(), ".png")),
       dpi = 600, width = 20, height = 10, units = "cm", limitsize = FALSE)

## Schoenfeld residuals
sch_resi<-ggcoxdiagnostics(cox_mod, type = "schoenfeld", ggtheme = theme_bw())
sch_resi

ggsave(sch_resi, file = safe_path(paste0("Schoenfeld residuals_", Sys.Date(), ".png")),
       dpi = 600, width = 20, height = 10, units = "cm", limitsize = FALSE)

## Scaled Schoenfeld residuals
schs_resi<-ggcoxdiagnostics(cox_mod, type = "scaledsch", ggtheme = theme_bw())
schs_resi

ggsave(schs_resi, file = safe_path(paste0("Scaled Schoenfeld residuals_", Sys.Date(), ".png")),
       dpi = 600, width = 20, height = 10, units = "cm", limitsize = FALSE)

# ==============================================================================
#### Assumption testing
# ==============================================================================

#- - - - - - - - - - - - - - - -
## Independence assumption

# The Cox model assumes that subjects are independent. Violations arise from:
#   - clustered data (patients within hospitals, families, etc.)
#   - repeated / matched observations
# Strategy:
#   (a) If cluster_var is specified then cluster-robust SEs already applied above
#   (b) Optionally fit a shared frailty model as a sensitivity check
# A frailty model explicitly models unobserved heterogeneity within clusters

if(!is.null(cluster_var)){ # with shared frailty model (sensitivity check)
  
  frailty_formula <- update(surv_formula,as.formula(paste(". ~ . + frailty(", cluster_var, ", distribution = 'gamma')")))
  cox_frailty <- coxph(frailty_formula, data = clean_data, ties = "efron")
  summary(cox_frailty)
  cox_frailty$history[[1]]$theta # if frailty variance (theta), unobserved clustering is present
  
}

#- - - - - - - - - - - - - - - -
## Proportional hazards (PH) assumption

# Approach 1 : Scaled Schoenfeld residuals (cox.zph) – global and per-covariate
# Approach 2 : log(-log(S(t))) vs. log(t) plot for categorical covariates
# p-value < alpha_ph indicates a PH violation

PH_test<-cox.zph(cox_mod)
ph_summary<-tidy_ph_test(PH_test)
ph_summary

# Identifying violated covariates
ph_violated_vars<-ph_summary %>%
                  filter(PH_violated & term != "GLOBAL") %>%
                  pull(term)

# Plotting scaled Schoenfeld residuals
png(file = safe_path(paste0("PH test_scaled Schoenfeld residuals_", Sys.Date(), ".png")), res = 600, width = 800, height = 700)
ggcoxzph(PH_test, font.main = 11)
dev.off()
ph_plots<-ggcoxzph(PH_test, font.main = 11)
ph_plots

# log(-log(S(t))) vs. log(t) for categorical covariates
cat_covariates<-covariates[sapply(clean_data[covariates], function(x) is.factor(x) | length(unique(x)) <= 5)]
if(length(cat_covariates) > 0){
  for(cv in cat_covariates){
    km_fit <- do.call(survfit,list(formula=as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", cv)),data=clean_data))
    p_loglog <- ggsurvplot(km_fit,
                           fun="cloglog",
                           ggtheme=theme_bw(),
                           title=paste0("log(-log(S(t))) vs. log(t): ", cv),
                           xlab="log(Time)",
                           ylab="log(-log(S(t)))")
    print(p_loglog)
  }
}

# Saving results
fwrite(ph_summary,safe_path(paste0("PH_test_results_", Sys.Date(), ".csv")),sep = ";")

#- - - - - - - - - - - - - - - -
## Linearity assumption (for continuous covariates only)

# Approach 1 : Martingale residuals from a null model plotted against each continuous covariate (non-linear relationship=curve)
# Approach 2 : Restricted cubic splines (RCS) via rms::cph() and rms::anova(). A non-linear (p<alpha_ph) term suggests linearity is violated

# Creating the null model
null_mod  <- coxph(as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ 1")),data = clean_data)

# Calculating the martingale residuals
mart_resid <- residuals(null_mod, type = "martingale")

# Selecting continuous covariates
cont_covariates <- covariates[-which(covariates %in% cat_covariates)]

# Linearity check
linearity_plots <- list()
for(cv in cont_covariates){
  df_plot<-tibble(x=clean_data[[cv]], resid=mart_resid)
  p<-ggplot(df_plot, aes(x=x,y=resid)) +
    geom_point(alpha=0.3, size=0.8, color="#40B696") +
    geom_smooth(method="loess",se=T,color= "red", linewidth = 1) +
    geom_hline(yintercept=0, linetype="dashed", color="black") +
    labs(title= paste0("Linearity check: ", cv),
         subtitle="Martingale residuals vs. covariate (loess smooth)",
         x=cv,y="Martingale residual") +
    theme_bw()
  linearity_plots[[cv]] <- p
  print(p)
}

# Testing non-linearity via restricted cubic splines (RCS)
dd<-datadist(clean_data) # distribution summary
options(datadist = "dd")

# Building RCS formula
rcs_terms<-sapply(cont_covariates, function(cv) paste0("rcs(", cv, ", ", rcs_knots, ")"))
non_cont<-setdiff(covariates, cont_covariates)
all_terms<-c(rcs_terms, non_cont)
rcs_formula<-as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", paste(all_terms, collapse = " + ")))

cox_rcs <- rms::cph(rcs_formula, data = clean_data, x = TRUE, y = TRUE, surv = TRUE) # RCS model ANOVA (non-linear terms)
anova_rcs <- anova(cox_rcs)
anova_rcs

#- - - - - - - - - - - - - - - -
## Multicollinearity

# Fitting a temporary linear model for VIF calculation
lm_temp  <-lm(as.formula(paste0(time_var, " ~ ", paste(covariates, collapse = " + "))),data = clean_data)

# VIF calculation
vif_result<-tidy_vif(lm_temp,VIF_set)
vif_result

# Saving results
fwrite(vif_result, safe_path(paste0("Multicollinearity check_", Sys.Date(), ".csv")), sep = ";")

# Correlation matrix plot for continuous covariates (if <=10 for readability)
if(length(cont_covariates)>1&length(cont_covariates)<=10) {
  cor_mat<-GGally::ggpairs(clean_data %>% select(all_of(cont_covariates)))+theme_cox()
  
  ggsave(cor_mat, file = safe_path(paste0("Correlation matrix_continuous var_", Sys.Date(), ".png")),
         dpi = 600, width = 30, height = 15, units = "cm", limitsize = FALSE)
}

#- - - - - - - - - - - - - - - -
## Interaction test (effect modification)

# Testing whether the effect of covariates differs across strata of the interaction_var (default: sex)
# Approach: Add covariate × interaction_var terms and compare via likelihood-ratio test 

if(!is.null(interaction_var) && interaction_var %in% covariates){ # if interaction_var is not NULL
  
  other_covs<-setdiff(covariates, interaction_var)
  
  # Testing each covariate × interaction_var individually
  interaction_results<-map_dfr(other_covs, function(cv){
    f_no_int<-surv_formula
    f_int<-update(surv_formula,as.formula(paste(". ~ . +", cv, "*", interaction_var)))
    m_no_int<-coxph(f_no_int, data = clean_data)
    m_int<-coxph(f_int,    data = clean_data)
    lrt<-anova(m_no_int, m_int)
    p_value <- lrt$`Pr(>|Chi|)`[2]
    tibble(covariate=cv,
           interaction=paste0(cv, " × ", interaction_var),
           LRT_chisq=lrt$Chisq[2],
           df=lrt$Df[2],
           p_value=p_value,
           modification=p_value < 0.05)})
  interaction_results
  
  # Global testing: all interactions simultaneously
  all_int_terms<-paste(paste0(other_covs, " * ", interaction_var), collapse = " + ")
  f_full_int<-update(surv_formula,as.formula(paste(". ~ . +", all_int_terms)))
  cox_full_int<-coxph(f_full_int, data=clean_data)
  global_lrt<-anova(cox_mod, cox_full_int)
  global_lrt
  
  # Stratified models for visual inspection
  strata_levels<-unique(clean_data[[interaction_var]])
  strat_models<-map(strata_levels, function(lvl){
    sub_data<-clean_data %>% filter(.data[[interaction_var]] == lvl)
    other_formula<-as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ",paste(other_covs, collapse = " + ")))
    coxph(other_formula, data = sub_data)
  })
  names(strat_models)<-paste0(interaction_var, "=", strata_levels)
  
  walk2(strat_models, names(strat_models), function(m, nm){
    cat("\nStratum:", nm, "\n")
    print(broom::tidy(m, exponentiate = TRUE, conf.int = TRUE))
  })
  
  # Saving interaction test results
  fwrite(interaction_results,safe_path(paste0("Interaction_test_", interaction_var, "_", Sys.Date(), ".csv")),sep = ";")
}

# ==============================================================================
#### In case of PH violation
# ==============================================================================

if(length(ph_violated_vars)>1){
  
  #- - - - - - - - - - - - - - - -
  ## Approach 1: covariate × log(time) interaction terms - <<USER>> covariate x time is also possible
  
  # Building time-transformation using tt() terms for violated covariates
  tt_terms<-paste0("tt(", ph_violated_vars, ")")
  non_tt_covs<-setdiff(covariates, ph_violated_vars)
  all_tt_terms<-c(non_tt_covs, tt_terms)
  formula_tt<-as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ",paste(all_tt_terms, collapse = " + ")))
  
  # Fitting the model
  cox_tt<- coxph(formula_tt,data = clean_data,tt= function(x, t, ...) x * log(t + 1))  # +1 avoids log(0)
  summary(cox_tt)
  
  # Comparing with original model via AIC/BIC
  AIC(cox_mod);AIC(cox_tt)
  BIC(cox_mod);BIC(cox_tt)
  
  # Tidy model characteristics
  cox_tt_charact<-broom::tidy(cox_tt, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(direction = case_when(conf.low  > 1 ~ "positive association",
                                 conf.high < 1 ~ "negative association",
                                 T~ "uncertain association")) %>%
    bind_cols(broom::glance(cox_tt)) %>%
    rowwise() %>%
    mutate(across(where(is.numeric), test_format)) %>%
    ungroup() %>%
    mutate(HR = paste0(estimate, " [", conf.low, "; ", conf.high, "]")) %>%
    mutate(HR = if_else(HR == " [; ]", "", HR))
  
  # Saving results
  fwrite(cox_tt_charact,safe_path(paste0("Time-dependent covariate model characteristics_", Sys.Date(), ".csv")),sep = ";")
  
  #- - - - - - - - - - - - - - - -
  ## Approach 2: time-dependent covariates (counting-process format)
  
  # Reformating the data into (tstart, tstop, event) counting-process format
  # PH-violating covariates are interacted with a piecewise-constant function of time (binary indicator: t > time_split_point), making them truly time-dependent
  #
  # Two sub-approaches:
  #   A: manual time splitting at the median follow-up time (or user-defined)
  #   B: using Greg::timeSplitter (flexible interval splits)
  
  # <<USER>> Set time split points (e.g., quartiles of follow-up)
  time_split_points<-quantile(clean_data[[time_var]], probs = c(0.25, 0.5, 0.75))
  
  # A - Manual time splitting
  cp_data<-survSplit(as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ .")),
                     data=clean_data,
                     cut=time_split_points,
                     episode="time_interval")
  
  # Creating time-dependent interaction terms for PH-violating covariates
  median_split<-median(clean_data[[time_var]])
  for(v in ph_violated_vars){
    cp_data[[paste0(v, "_td")]]<-cp_data[[v]] * as.integer(cp_data[["tstart"]] >= median_split)
  }
  
  # Building the counting-process formula
  td_int_terms<-paste0(ph_violated_vars, "_td")
  cp_formula_str<-paste0("Surv(tstart, ", time_var, ", ", event_var, ") ~ ",paste(c(covariates, td_int_terms), collapse = " + "))
  cox_td_A<-coxph(as.formula(cp_formula_str), data = cp_data, ties = "efron")
  summary(cox_td_A)
  
  # Tidy model characteristics
  cox_td_A_charact<-broom::tidy(cox_td_A, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(direction = case_when(conf.low  > 1 ~ "positive association",
                                 conf.high < 1 ~ "negative association",
                                 T~ "uncertain association")) %>%
    bind_cols(broom::glance(cox_td_A)) %>%
    rowwise() %>%
    mutate(across(where(is.numeric), test_format)) %>%
    ungroup() %>%
    mutate(HR = paste0(estimate, " [", conf.low, "; ", conf.high, "]")) %>%
    mutate(HR = if_else(HR == " [; ]", "", HR))
  
  # Saving the results results
  fwrite(cox_td_A_charact,safe_path(paste0("Time-dependent covariate_counting_manual_", Sys.Date(), ".csv")),sep = ";")
  
  # B - using Greg::timeSplitter
  if(use_time_splitter){
    ## <<USER>> Set split_by to your time unit (e.g., 365.25 = years from days)
    ## <<USER>> Adjust t.start and t.event to the correct column names after data re-coding if they differ from tstart / time_to_diagnosis
    
    # timeSplitter requires an ID variable – adding one if absent
    if(!ID_obs %in% colnames(clean_data)){
      clean_data_ts<-clean_data %>% mutate(id = row_number())
    }else{
      clean_data_ts<-clean_data
      colnames(clean_data_ts)[which(colnames(clean_data_ts)==ID_obs)]<-"id"
    }
    
    # Splitting time axis into intervals of equal length - <<USER>> Adjust split_by and max.follow to match your time scale
    split_by<-median(clean_data_ts[[time_var]])/4  # quarterly splits
    max_follow<-max(clean_data_ts[[time_var]])
    
    ts_data<-Greg::timeSplitter(data=clean_data_ts,
                                by=split_by,
                                event_var=event_var,
                                time_var=time_var,
                                event_start_status=0)
    
    # Creating time-dependent interaction for PH-violating covariates
    ts_median_split<-median(clean_data[[time_var]])
    for(v in ph_violated_vars){
      ts_data[[paste0(v, "_td")]]<-ts_data[[v]]*as.integer(ts_data[["Start_time"]]>=ts_median_split)
    }
    
    # Building and fitting the model
    td_terms_ts<-paste0(ph_violated_vars, "_td")
    ts_formula<-as.formula(paste0("Surv(Start_time, Stop_time, ", event_var, ") ~ ",paste(c(covariates, td_terms_ts), collapse = " + ")))
    cox_td_B <- coxph(ts_formula, data = ts_data, cluster = id, ties = "efron")
    summary(cox_td_B)
    
    # Tidy model characteristics
    cox_td_B_charact<-broom::tidy(cox_td_B, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(direction = case_when(conf.low  > 1 ~ "positive association",
                                   conf.high < 1 ~ "negative association",
                                   T~ "uncertain association")) %>%
      bind_cols(broom::glance(cox_td_B)) %>%
      rowwise() %>%
      mutate(across(where(is.numeric), test_format)) %>%
      ungroup() %>%
      mutate(HR = paste0(estimate, " [", conf.low, "; ", conf.high, "]")) %>%
      mutate(HR = if_else(HR == " [; ]", "", HR))
    
    # Saving the results results
    fwrite(cox_td_B_charact,safe_path(paste0("Time-dependent covariate_counting_Greg_", Sys.Date(), ".csv")),sep = ";")
    
  }
  
  # Comparing models with AIC
  AIC(cox_mod)
  AIC(cox_tt)
  AIC(cox_td_A)
  if (use_time_splitter && exists("cox_td_B")) {
    print(AIC(cox_td_B))
  }
  
}

# ==============================================================================
#### Competing risk analysis
# ==============================================================================

# Approach: competing_risk_var is specified, fitting a Fine-Gray subdistribution hazard
# Requirements: competing_risk_var must be a factor with levels: 0 = censored, 1 = event of interest, 2 = competing event

if(!is.null(competing_risk_var) && competing_risk_var %in% colnames(clean_data)){
  
  # Ensuring the competing risk variable is a factor with correct levels
  clean_data[[competing_risk_var]]<-factor(clean_data[[competing_risk_var]],levels=c(0,1,2),
                                           labels = c("censored","event","competing"))
  
  # Creating the model formula
  cr_formula<-as.formula(paste0("Surv(", time_var, ", ", competing_risk_var, ") ~ ",paste(covariates, collapse = " + ")))
  
  # Fitting a Fine-Gray subdistribution hazard model
  crr_mod<-tidycmprsk::crr(cr_formula, data = clean_data)
  
  # Tidy model characteristics
  crr_mod_charact<-broom::tidy(crr_mod, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(direction = case_when(conf.low  > 1 ~ "positive association",
                                 conf.high < 1 ~ "negative association",
                                 T~ "uncertain association")) %>%
    bind_cols(broom::glance(crr_mod)) %>%
    rowwise() %>%
    mutate(across(where(is.numeric), test_format)) %>%
    ungroup() %>%
    mutate(HR = paste0(estimate, " [", conf.low, "; ", conf.high, "]")) %>%
    mutate(HR = if_else(HR == " [; ]", "", HR))
  
  # Fitting a cause-specific Cox model for the event of interest (censoring competing event)
  cs_data_event<-clean_data %>% mutate(cs_event = as.integer(.data[[competing_risk_var]] == "event"))
  cs_formula_event<-as.formula(paste0("Surv(", time_var, ", cs_event) ~ ",paste(covariates, collapse = " + ")))
  cox_cs_event<-coxph(cs_formula_event, data = cs_data_event)
  summary(cox_cs_event)
  
  # Tidy model characteristics
  cox_cs_event_charact<-broom::tidy(cox_cs_event, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(direction = case_when(conf.low  > 1 ~ "positive association",
                                 conf.high < 1 ~ "negative association",
                                 T~ "uncertain association")) %>%
    bind_cols(broom::glance(cox_cs_event)) %>%
    rowwise() %>%
    mutate(across(where(is.numeric), test_format)) %>%
    ungroup() %>%
    mutate(HR = paste0(estimate, " [", conf.low, "; ", conf.high, "]")) %>%
    mutate(HR = if_else(HR == " [; ]", "", HR))
  
  # Fitting a cause-specific Cox model for the competing event
  cs_data_comp<-clean_data %>% mutate(cs_comp = as.integer(.data[[competing_risk_var]] == "competing"))
  cs_formula_comp<-as.formula(paste0("Surv(", time_var, ", cs_comp) ~ ",paste(covariates, collapse = " + ")))
  cox_cs_comp<-coxph(cs_formula_comp, data = cs_data_comp)
  summary(cox_cs_comp)
  
  # Tidy model characteristics
  cox_cs_comp_charact<-broom::tidy(cox_cs_comp, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(direction = case_when(conf.low  > 1 ~ "positive association",
                                 conf.high < 1 ~ "negative association",
                                 T~ "uncertain association")) %>%
    bind_cols(broom::glance(cox_cs_comp)) %>%
    rowwise() %>%
    mutate(across(where(is.numeric), test_format)) %>%
    ungroup() %>%
    mutate(HR = paste0(estimate, " [", conf.low, "; ", conf.high, "]")) %>%
    mutate(HR = if_else(HR == " [; ]", "", HR))
  
  # Cumulative Incidence Function plot (CIF)
  cif_formula<-as.formula(paste0("Surv(", time_var, ", ", competing_risk_var, ") ~ 1"))
  cif_fit<-tidycmprsk::cuminc(cif_formula, data = clean_data)
  p_cif<-ggcuminc(cif_fit,color = "#40B696") +
    labs(x="Time",y="Cumulative Incidence") +
    add_confidence_interval(fill="#40B696") +
    add_risktable() +
    theme_bw()
  p_cif
  
  ggsave(p_cif, file = safe_path(paste0("Cumulative Incidence Function_",competing_risk_var,"_", Sys.Date(), ".png")),
         dpi = 600, width = 30, height = 15, units = "cm", limitsize = FALSE)
  
  # CIF stratified by interaction_var (if available)
  if (!is.null(interaction_var)) {
    cif_strat_formula<-as.formula(paste0("Surv(", time_var, ", ", competing_risk_var, ") ~ ", interaction_var))
    cif_strat<-tidycmprsk::cuminc(cif_strat_formula, data = clean_data)
    p_cif_strat<-ggcuminc(cif_strat) +
      labs(x="Time",y="Cumulative Incidence") +
      add_confidence_interval() +
      add_risktable() +
      theme_bw()+
      scale_color_manual("Sex",values=c('0'="orange",'1'="#40B696"),labels=c('0'="female",'1'="male"))+
      scale_fill_manual("Sex",values=c('0'="orange",'1'="#40B696"),labels=c('0'="female",'1'="male"))
    p_cif_strat
    
    ggsave(p_cif_strat, file = safe_path(paste0("CIF_",competing_risk_var,"stratified by_",interaction_var,"_", Sys.Date(), ".png")),
           dpi = 600, width = 30, height = 15, units = "cm", limitsize = FALSE)
    
  }
  
  # Saving results
  fwrite(cox_cs_event_charact,safe_path(paste0("Competing_Risk_Cause_Specific_event_", Sys.Date(), ".csv")),sep = ";")
  fwrite(cox_cs_comp_charact,safe_path(paste0("Competing_Risk_Cause_Specific_competing_", Sys.Date(), ".csv")),sep = ";")
  
}

# ==============================================================================
#### Visualisation
# ==============================================================================

## Collecting all models in a named list
model_list<-list("Cox model"=if(exists("cox_mod"))cox_mod,
                 "Time dependent Cox model with interaction term"=if(exists("cox_tt"))cox_tt,
                 "Time dependent Cox model with manual counting process"=if(exists("cox_td_A"))cox_td_A,
                 "Time dependent Cox model Greg counting process"=if(exists("cox_td_B"))cox_td_B,
                 "Competing risk Fine-Gray subdistribution hazard model"= if (exists("crr_mod"))crr_mod,
                 "Competing risk cause-specific Cox model for the event of interest"=if(exists("cox_cs_event"))cox_cs_event,
                 "Competing risk cause-specific Cox model for the competing event"=if(exists("cox_cs_comp"))cox_cs_comp) %>%
  purrr::compact()   # drops any NULL entries (models that were never created)


## Generating and saving all plots
imap(model_list, function(mod, label){
  
  #- - - - - - - - - - - - - - - -
  ## Forest plot
  
  if (inherits(mod, "coxph")) {
    
    # version 1
    p<-ggforest(mod, main = paste("Hazard Ratios –", label))
    
    ggsave(plot=p, filename=safe_path(paste0("Forest_plot_v1_", to_filename(label), "_", Sys.Date(), ".png")),
           dpi=600, width= 35, height= 20, units= "cm", limitsize=F)
  }
  
  # version 2
  p2<-broom::tidy(mod,exponentiate=T,conf.int=T) %>%
    mutate(direction=case_when(conf.low>1~"positive",
                               conf.high<1~"negative",
                               T~"uncertain"),
           HR=cox_mod_charact$HR) %>%
    ggplot() + geom_point(aes(x=term,y=estimate,color=factor(direction)),size=2) +
    geom_errorbar(aes(x=term,ymin=conf.low,ymax=conf.high,color=factor(direction)))+
    scale_x_discrete("")+
    scale_y_log10(breaks=c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,10,50,100),
                  labels = c("0.01","0.05","0.1","0.2","0.3","0.4","0.5","1",
                             "2","3","4","5","10","50","100"),
                  limits=c(0.01,100))+
    geom_hline(yintercept = 1,linetype=2, col="darkblue")+
    scale_color_manual("Association:",
                       values=c('positive'="red",'negative'="#40B696",'uncertain'="black")) +
    geom_text(aes(y=50,x=term,label=HR),color="black",size=4)+
    theme_cox()+
    coord_flip()
  
  ggsave(plot=p2, filename=safe_path(paste0("Forest_plot_v2_", to_filename(label), "_", Sys.Date(), ".png")),
         dpi=600, width= 35, height= 20, units= "cm", limitsize=F)
  
  # version 3
  
  if (!grepl("interaction term", label)) {
  p3<-ggstats::ggcoef_table(mod, exponentiate=T, conf.int=T)
  
  ggsave(plot=p3, filename=safe_path(paste0("Forest_plot_v3_", to_filename(label), "_", Sys.Date(), ".png")),
         dpi=600, width= 35, height= 20, units= "cm", limitsize=F)
  }
  
  #- - - - - - - - - - - - - - - -
  ## Baseline survival (overall)
  
  if (!grepl("interaction term", label) && !grepl("Fine-Gray", label)) {
    
  # Version 1
  png(filename = safe_path(paste0("Overall_baseline_survival_v1_", to_filename(label), "_", Sys.Date(), ".png")),
      res = 600, width = 35, height = 20, units = "cm")
  ggsurvplot(survfit(mod, data = clean_data),
                   palette="#40B696",
                   risk.table=TRUE,
                   pval=TRUE,
                   conf.int=TRUE,
                   ggtheme=theme_cox(),
                   xlab="Time",
                   ylab="Overall survival probability")
  dev.off()
  
  # Version 2
  p2<-survfit2(mod) %>%
    ggsurvfit(color = "#40B696") +
    labs(x = "Time", y = "Overall survival probability") +
    add_confidence_interval(fill = "#40B696") +
    add_risktable()+
    theme_cox()
  
  ggsave(plot=p2, filename=safe_path(paste0("Overall baseline survival_v2_", to_filename(label), "_", Sys.Date(), ".png")),
         dpi=600, width= 35, height= 20, units= "cm", limitsize=F)
  
  }
  
  # Version 3
  if(grepl("Cox model", label)){
    p3<-ggplot(clean_data, aes(time = .data[[time_var]], status = .data[[event_var]])) +
      geom_km(fill = "#40B696", color = "#40B696") +
      geom_kmband(fill = "#40B696", color = "#40B696") +
      scale_x_continuous("Time") +
      scale_y_continuous("Overall survival probability") +
      theme_cox()
    
    ggsave(plot=p3, filename=safe_path(paste0("Overall baseline survival_v3_", to_filename(label), "_", Sys.Date(), ".png")),
           dpi=600, width= 35, height= 20, units= "cm", limitsize=F)
  }
  
  #- - - - - - - - - - - - - - - -
  ## Cumulative events and cumulative hazards (overall)
  
  if (!grepl("interaction term", label) && !grepl("Fine-Gray", label)) {
    
    walk(list(
      list(fun = "event",  ylabel = "Cumulative event probability",title = "Cumulative events"),
      list(fun = "cumhaz", ylabel = "Cumulative hazard",title = "Cumulative hazard")), function(cfg){
        
        # version 1
        png(filename = safe_path(paste0(cfg$title,"_v1_", to_filename(label), "_", Sys.Date(), ".png")),
            res = 600, width = 35, height = 20, units = "cm")
        
        p1<-ggsurvplot(survfit(mod, data = clean_data),
                       palette= "#40B696",
                       risk.table= TRUE,
                       pval= TRUE,
                       conf.int= TRUE,
                       ggtheme= theme_cox(),
                       xlab="Time",
                       ylab=cfg$ylabel,
                       fun=cfg$fun)
        
        print(p1)
        dev.off()
        
        # version 2
        p2<-survfit2(mod) %>%
          ggsurvfit(color = "#40B696",
                    type  = if (cfg$fun == "event") "risk" else "cumhaz") +
          labs(x = "Time", y = cfg$ylabel) +
          add_confidence_interval(fill = "#40B696") +
          add_risktable()+
          theme_cox()
        
        ggsave(plot=p2, filename=safe_path(paste0(cfg$title,"_v2_", to_filename(label), "_", Sys.Date(), ".png")),
               dpi=600, width= 35, height= 20, units= "cm", limitsize=F)
      })
  }
 
})

#- - - - - - - - - - - - - - - -
## Stratified plots (by interaction_var, e.g., sex)

if(!is.null(interaction_var)){
  strat_km_formula<-as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", interaction_var))
  
  # Log-rank test
  logrank_test<-survdiff(strat_km_formula, data = clean_data)
  fwrite(as_tibble(logrank_test$pvalue),safe_path(paste0("p-value log-rank test_", Sys.Date(), ".csv")),sep = ";")
  
  # Pre-building the survfit object
  strat_km_fit <- do.call(survfit,list(formula = strat_km_formula, data = clean_data))
  
  # Plotting and saving
  walk(list(
    list(fun=NULL, ylabel="Survival probability", title="Survival by "),
    list(fun="event", ylabel="Cumulative event probability", title = "Cumulative events by "),
    list(fun="cumhaz", ylabel="Cumulative hazard", title="Cumulative hazard by ")), function(cfg){
      
      #version 1
      png(filename = safe_path(paste0(cfg$title,interaction_var,"_v1_", Sys.Date(), ".png")),
          res = 600, width = 35, height = 20, units = "cm")
      
      p1<-ggsurvplot(strat_km_fit,
                     ggtheme=theme_cox(),
                     data=clean_data,
                     risk.table=TRUE,
                     pval=TRUE,
                     conf.int=TRUE,
                     palette=c("orange", "#40B696"),
                     legend.title=interaction_var,
                     legend.labs=names(interaction_var_labels),
                     xlab="Time",
                     ylab=cfg$ylabel,
                     fun=cfg$fun)
      
      print(p1)
      dev.off()
      
      # version 2
      p2<-survfit2(strat_km_formula, data = clean_data) %>%
        ggsurvfit(type = if (is.null(cfg$fun)) "surv" else if (cfg$fun == "event") "risk" else "cumhaz") +
        labs(x = "Time", y = cfg$ylabel) +
        add_confidence_interval() +
        add_risktable() +
        scale_color_manual(interaction_var,
                           values = c("orange", "#40B696"),
                           labels = interaction_var_labels) +
        scale_fill_manual(interaction_var,
                          values = c("orange", "#40B696"),
                          labels = interaction_var_labels)+
        theme_cox()
      
      ggsave(plot=p2, filename=safe_path(paste0(cfg$title,interaction_var,"_v2_", Sys.Date(), ".png")),
             dpi=600, width= 35, height= 20, units= "cm", limitsize=F)
      
      # version 3
      p3<-ggplot(clean_data,
                 aes(time   = .data[[time_var]],
                     status = .data[[event_var]],
                     fill   = factor(.data[[interaction_var]]),
                     color  = factor(.data[[interaction_var]]))) +
        geom_km() +
        geom_kmband() +
        scale_fill_manual(interaction_var,
                          values = c("orange", "#40B696"),
                          labels = interaction_var_labels) +
        scale_color_manual(interaction_var,
                           values = c("orange", "#40B696"),
                           labels = interaction_var_labels) +
        scale_x_continuous("Time") +
        scale_y_continuous("Survival probability") +
        theme_cox()
      
      ggsave(plot=p3, filename=safe_path(paste0(cfg$title,interaction_var,"_v3_", Sys.Date(), ".png")),
             dpi=600, width= 35, height= 20, units= "cm", limitsize=F)
      
    })
}
