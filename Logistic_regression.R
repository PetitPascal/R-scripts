#-------------------------------------------------------------------------------
## Reproducible & Generalisable Logistic Regression Script
# Covers:
#   - Simple binary logistic regression
#   - Multiple binary logistic regression
#   - Assumption testing (overdispersion, separation, multicollinearity, influential observations, Hosmer-Lemeshow, ROC/AUC)
#   - Ordinal logistic regression (proportional odds model)
#   - Odds ratio plots, predicted probability plots, model comparison
#   - Reusable pipeline functions
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs <- c("MASS", "dplyr", "tidyr", "ggplot2", "broom","gtsummary", "car", "lmtest", "performance",
                   "ResourceSelection", "pROC", "ordinal","ggeffects", "nnet", "broom.helpers",
                   "forestmodel", "ggstats", "questionr","conflicted", "tibble","ggstats")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# for avoiding package conflicts
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("tidy","broom")
conflicted::conflicts_prefer(stats::predict)

#--------------------------------------------
## Step 2: Simulating dataset

n<-600
age<-round(rnorm(n, 50, 12))
sex<-factor(rbinom(n, 1, 0.5), labels = c("Female", "Male"))
bmi<-rnorm(n, 26, 4)
smoking<-factor(rbinom(n, 1, 0.3), labels = c("No", "Yes"))
exposure<-rnorm(n, 5, 2) # continuous exposure (e.g. pollutant)
job<-factor(sample(c("Office","Manual","Agriculture"), n,replace = TRUE, prob = c(0.4, 0.4, 0.2)))

# Linear predictor for binary outcome
lp_bin<- -3 + 0.04*age + 0.5*(sex=="Male") +0.08*bmi + 0.7*(smoking=="Yes") + 0.35*exposure

# Binary outcomes
disease<-rbinom(n, 1, plogis(lp_bin)) # primary binary outcome
cancer<-rbinom(n, 1, plogis(lp_bin - 0.5)) # secondary binary outcome

# Ordinal outcome (disease severity: None < Low < Moderate < High)
sev_lp<- lp_bin
disease_sev<-cut(sev_lp,
                   breaks = quantile(sev_lp, c(0, 0.35, 0.60, 0.82, 1)),
                   labels = c("None","Low","Moderate","High"),
                   include.lowest = TRUE)
disease_sev<-factor(disease_sev,
                      levels = c("None","Low","Moderate","High"),
                      ordered = TRUE)

sim_data <- tibble(id=seq_len(n),
  disease=disease,
  cancer=cancer,
  disease_sev=disease_sev,
  age=age,
  sex=sex,
  bmi=bmi,
  smoking=smoking,
  exposure=exposure,
  job=job)


#--------------------------------------------
## Step 3: EDA
sim_data %>%
  gtsummary::tbl_summary(by= disease,include = c(age, sex, bmi, smoking, exposure, job)) %>%
  gtsummary::add_p() %>%
  gtsummary::bold_labels()

#--------------------------------------------
## Simple binary logistic regression

# Fitting the model
mod_simple<-glm(disease ~ exposure,data   = sim_data,family = binomial(link = "logit"))
summary(mod_simple)

# Odds ratio
or_simple<-broom::tidy(mod_simple, conf.int = TRUE, exponentiate = TRUE)
or_simple

# Assumption check - minimum events per variable (rule of thumb: >= 10 per predictor)
table(sim_data$disease)

# Assumption check - Overdispersion (residual deviance/df; should be ~1)
disp_simple<-mod_simple$deviance / mod_simple$df.residual
disp_simple

# Assumption check - Linearity of log-odds (Box-Tidwell for continuous predictor)
sim_data$exp_log<-sim_data$exposure * log(sim_data$exposure + abs(min(sim_data$exposure)) + 1)
bt_test<-glm(disease ~ exposure + exp_log, data = sim_data, family = binomial())
bt_p <-broom::tidy(bt_test) %>% filter(term == "exp_log") %>% pull(p.value)
round(bt_p, 4) # Box-Tidwell linearity test (p < 0.05 = non-linearity)
sim_data$exp_log <- NULL   # remove temporary variable

# Visualization

ggplot(sim_data, aes(x = exposure, y = disease)) +
  geom_jitter(height = 0.03, alpha = 0.3) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE, color = "#2166ac") +
  labs(title = "Simple logistic regression: disease ~ exposure",
       x     = "Exposure", y     = "P(disease = 1)") +
  theme_bw(base_size = 14)

# Predicted probabilities
plot(ggeffects::ggeffect(mod_simple, terms = "exposure")) +
  labs(title = "Marginal predicted probabilities (simple logistic)")

# Formatted table
mod_simple %>%
  gtsummary::tbl_regression(intercept = TRUE, exponentiate = TRUE) %>%
  gtsummary::bold_p()

#--------------------------------------------
## Multiple binary logistic regression

#- - - - 
## Model fitting
mod1<-glm(disease ~ exposure + age, data = sim_data, family = binomial(link = "logit"))

mod2<-glm(disease ~ exposure + age + sex + bmi + smoking + job,
            data = sim_data, family = binomial(link = "logit"))

summary(mod1)
summary(mod2)

#- - - - 
## Assumption checks
check_logistic_assumptions <- function(model, label = "") {
  df_mod<-model$data
  y_resp<-model$y
  n_obs<-length(y_resp)
  p_vars<-length(coef(model)) - 1
  
  # Events per variable (EPV >= 10 recommended)
  n_events<-sum(y_resp)
  epv<-n_events / p_vars
  cat("EPV:", round(epv, 1),
      if (epv < 10) "*** LOW EPV — results may be unreliable ***" else "OK", "\n")
  
  # Overdispersion
  disp<-model$deviance / model$df.residual
  cat("Overdispersion (dev/df):", round(disp, 3),
      if (disp > 1.5) "*** OVERDISPERSION ***" else "OK", "\n")
  
  # Multicollinearity (VIF)
  vif_res <- tryCatch(car::vif(model), error = function(e) NULL)
  if (!is.null(vif_res)) {
    cat("Max VIF:", round(max(vif_res), 3),
        if (max(vif_res) > 5) "*** HIGH VIF ***" else "OK", "\n")
  }
  
  # Hosmer-Lemeshow goodness-of-fit
  hl<-tryCatch(
    ResourceSelection::hoslem.test(y_resp, fitted(model), g = 10),
    error = function(e) NULL)
  if (!is.null(hl))
    cat("Hosmer-Lemeshow: chi2 =", round(hl$statistic, 3),
        "| p =", round(hl$p.value, 4),
        "(p > 0.05 = adequate fit)\n")
  
  # Influential observations (Cook's distance)
  cooks<-cooks.distance(model)
  n_infl<-sum(cooks > 4 / n_obs)
  cat("Influential obs (Cook > 4/n):", n_infl, "\n")
  
  # Complete separation check (very high coefficients = warning sign)
  max_coef <- max(abs(coef(model)[-1]), na.rm = TRUE)
  cat("Max |coefficient|:", round(max_coef, 2),
      if (max_coef > 10) "*** POSSIBLE SEPARATION ***" else "OK", "\n")
  
  # AUC / ROC
  roc_obj <- pROC::roc(y_resp, fitted(model), quiet = TRUE)
  cat("AUC:", round(pROC::auc(roc_obj), 4), "\n")
  
  invisible(list(epv = epv, disp = disp, vif = vif_res,
                 hl = hl, auc = pROC::auc(roc_obj)))
}

checks_mod1<-check_logistic_assumptions(mod1, "Model 1")
checks_mod2<-check_logistic_assumptions(mod2, "Model 2")

# ROC curves
roc1<-pROC::roc(sim_data$disease, fitted(mod1), quiet = TRUE)
roc2<-pROC::roc(sim_data$disease, fitted(mod2), quiet = TRUE)

pROC::ggroc(list("Model 1" = roc1, "Model 2" = roc2),
            linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 1,
              linetype = "dashed", color = "grey50") +
  labs(title = "",
       subtitle = paste0("AUC Mod1 = ", round(pROC::auc(roc1), 3),
                         " | AUC Mod2 = ", round(pROC::auc(roc2), 3)),
       color = "Model") +
  theme_bw(base_size = 14)

#- - - - 
## Odds ratios table and plot
extract_or <- function(model) {
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(sig= case_when(
        conf.low  > 1 ~ "OR > 1",
        conf.high < 1 ~ "OR < 1",
        TRUE          ~ "NS (1 in CI)"),
      sig    = factor(sig, levels = c("OR < 1","NS (1 in CI)","OR > 1")))
}

plot_or<- function(or_df, title = "Odds Ratios (95% CI)") {
  ggplot(or_df, aes(x = reorder(term, estimate),
                    y = estimate,
                    ymin = conf.low, ymax = conf.high,
                    color = sig)) +
    geom_pointrange(size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "darkblue") +
    coord_flip() +
    scale_color_manual(
      values = c("OR < 1"       = "darkgreen",
                 "NS (1 in CI)" = "black",
                 "OR > 1"       = "red"),
      name = "") +
    labs(title = title, x = "", y = "Odds Ratio (95% CI)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
}

# OR table
or_mod1 <- extract_or(mod1)
or_mod2 <- extract_or(mod2)
or_mod1
or_mod2

# Forest plot
ggstats::ggcoef_model(mod1, exponentiate = TRUE)
ggstats::ggcoef_model(mod2, exponentiate = TRUE)

ggstats::ggcoef_table(mod1, exponentiate = TRUE)
ggstats::ggcoef_table(mod2, exponentiate = TRUE)

forestmodel::forest_model(mod1)
forestmodel::forest_model(mod2)

# Predicted probabilities - marginal effect of exposure at representative values
plot(ggeffects::ggeffect(mod2, terms = "exposure")) +
  labs(title = "Predicted P(disease) by exposure (Model 2)")

plot(ggeffects::ggeffect(mod2, terms = c("exposure", "sex"))) +
  labs(title = "Predicted P(disease): exposure × sex")

#- - - - 
## Model comparison
mod_null <- glm(disease ~ 1, data = sim_data, family = binomial())

# Likelihood ratio tests
anova(mod_null, mod1, mod2, test = "Chisq")
lmtest::lrtest(mod_null, mod1, mod2)

# AIC / BIC
AIC(mod_null, mod1, mod2) %>% arrange(AIC)
BIC(mod_null, mod1, mod2) %>% arrange(BIC)

# McFadden pseudo-R2
pseudo_r2 <- function(mod, null_mod)
  1 - as.numeric(logLik(mod)) / as.numeric(logLik(null_mod))
cat("McFadden R2 — Mod1:", round(pseudo_r2(mod1, mod_null), 4),
    "| Mod2:", round(pseudo_r2(mod2, mod_null), 4), "\n")

# Summary table
gtsummary::tbl_merge(
  tbls = list(
    gtsummary::tbl_regression(mod1,exponentiate = TRUE),
    gtsummary::tbl_regression(mod2,exponentiate = TRUE)),
  tab_spanner = c("Model 1", "Model 2"))

ggcoef_compare(list("Model 1" = mod1, "Model 2" = mod2),
               type = "faceted",
               intercept=F,
               exponentiate = TRUE)

#- - - - 
## Overdispersion correction (quasibinomial) (use if overdispersion > 1.5)

mod2_quasi<-glm(disease ~ exposure + age + sex + bmi + smoking + job,
                  data= sim_data,
                  family= quasibinomial(link = "logit"))

# Type II Wald test using F-statistic
car::Anova(mod2_quasi, test = "F")

#- - - - 
## Influential observation diagnostics
infl_df <- data.frame(obs=seq_len(nrow(sim_data)),
  cooks=cooks.distance(mod2),
  leverage=hatvalues(mod2),
  std_res=rstandard(mod2))

ggplot(infl_df, aes(x = leverage, y = std_res, size = cooks)) +
  geom_point(alpha = 0.4, color = "#2166ac") +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") +
  labs(title    = "Influence plot: leverage vs standardised residuals",
       subtitle = "Point size = Cook's distance",
       x = "Leverage", y = "Standardised residual", size = "Cook's D") +
  theme_bw(base_size = 13)

# Flag highly influential observations
high_infl <- infl_df %>% filter(cooks > 4 / nrow(sim_data))
cat("\nHighly influential observations (Cook > 4/n):", nrow(high_infl), "\n")
if (nrow(high_infl) > 0) print(head(high_infl))

#--------------------------------------------
## Ordinal logistic regression

## Model fitting
mod_ord <- ordinal::clm(disease_sev ~ exposure + age + bmi + smoking,
  data = sim_data,
  link = "logit")   # proportional odds

summary(mod_ord)

## Odds ratios
or_ord <- broom::tidy(mod_ord, conf.int = TRUE, exponentiate = TRUE) %>%
          filter(coef.type == "location")   # exclude threshold parameters
or_ord

## Assumption check — proportional odds

# Nominal test: tests if odds ratios differ across response levels
nom_test <- ordinal::nominal_test(mod_ord)
nom_test # p < 0.05 for a term = PO assumption violated for that term

# Scale test: tests proportionality of scale
scale_test <- tryCatch(ordinal::scale_test(mod_ord),
                       error = function(e) NULL)
if (!is.null(scale_test)) {
  cat("\n--- Scale test ---\n")
  print(scale_test)
}

# Visual check: log-odds of cumulative probabilities should be parallel
cum_probs <- sim_data %>%
  mutate(age_group = cut(age, breaks = 4)) %>%
  group_by(age_group) %>%
  summarise(p_none=mean(disease_sev == "None"),
    p_low=mean(disease_sev %in% c("None","Low")),
    p_mod=mean(disease_sev %in% c("None","Low","Moderate")),
    .groups= "drop") %>%
  mutate(logit_none=log(p_none/(1 - p_none   + 1e-6)),
    logit_low=log(p_low/(1 - p_low    + 1e-6)),
    logit_mod=log(p_mod/(1 - p_mod    + 1e-6)))

cum_probs_long<-cum_probs %>%
  tidyr::pivot_longer(cols = starts_with("logit_"),
                      names_to = "threshold", values_to = "logit")

ggplot(cum_probs_long, aes(x = age_group, y = logit,
                           group = threshold, color = threshold)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
  labs(title = "Proportional odds check: cumulative logits by age group",
       subtitle = "Parallel lines = PO assumption satisfied",
       x = "Age group", y = "Cumulative log-odds", color = "Threshold") +
  theme_bw(base_size = 13)

## Type II Wald tests for ordinal model
car::Anova(mod_ord)

## Formatted table
mod_ord %>%
  gtsummary::tbl_regression(exponentiate = TRUE) %>%
  gtsummary::bold_p()

## Coefficient plot
ggstats::ggcoef_model(mod_ord, exponentiate = TRUE) +
  labs(title = "Ordinal logistic regression: OR plot")

## Predicted probabilities of each category across exposure levels
exp_seq<-seq(min(sim_data$exposure), max(sim_data$exposure), length.out = 100)
pred_frame<-data.frame(
  exposure=exp_seq,
  age=median(sim_data$age),
  sex=factor("Female", levels = levels(sim_data$sex)),
  bmi=median(sim_data$bmi),
  smoking=factor("No", levels = levels(sim_data$smoking)))

pred_ord<-predict(mod_ord, newdata = pred_frame, type = "prob")
pred_ord_df<-cbind(pred_frame, as.data.frame(pred_ord$fit)) %>%
  tidyr::pivot_longer(cols = c("None","Low","Moderate","High"),
                      names_to = "category", values_to = "prob") %>%
  mutate(category = factor(category,levels = c("None","Low","Moderate","High")))

ggplot(pred_ord_df, aes(x = exposure, y = prob, color = category)) +
  geom_line(linewidth = 1.2) +
  labs(title="Ordinal logistic: predicted probabilities by exposure",
       x="Exposure",
       y="Predicted probability",
       color="Severity") +
  scale_color_brewer(palette = "RdYlGn", direction = -1) +
  theme_bw(base_size = 14)

## Model comparison
mod_ord_null<-ordinal::clm(disease_sev ~ 1, data = sim_data, link = "logit")
lrt_ord<-anova(mod_ord_null, mod_ord)
cat("\nOrdinal LRT vs. null\n")
lrt_ord
cat("AIC — Null:", round(AIC(mod_ord_null), 2),
    "| Full:", round(AIC(mod_ord), 2), "\n")

#--------------------------------------------
## Reusable pipeline functions

## Binary logistic pipeline
run_logistic_pipeline <- function(data, outcome_col, predictor_cols,
                                  family=binomial(link = "logit"),
                                  ref_levels=NULL,
                                  run_roc=TRUE,
                                  run_hl=TRUE,
                                  seed=123) {
  set.seed(seed)
  df <- data
  
  # Apply reference levels if provided
  if (!is.null(ref_levels)) {
    for (v in names(ref_levels))
      df[[v]] <- relevel(factor(df[[v]]), ref = ref_levels[[v]])
  }
  
  form<-as.formula(paste(outcome_col, "~",
                           paste(predictor_cols, collapse = " + ")))
  fit <-glm(form, data = df, family = family)
  
  cat("\nLogistic regression summary\n")
  print(summary(fit))
  
  or_df <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  cat("\nOdds Ratios\n"); print(or_df)
  
  # EPV
  n_ev <- sum(fit$y)
  cat("\nEPV:", round(n_ev / length(predictor_cols), 1), "\n")
  
  # Overdispersion
  disp <- fit$deviance / fit$df.residual
  cat("Overdispersion:", round(disp, 3),
      if (disp > 1.5) "— consider quasibinomial" else "— OK", "\n")
  
  # VIF
  vif_r <- tryCatch(car::vif(fit), error = function(e) NULL)
  if (!is.null(vif_r))
    cat("Max VIF:", round(max(vif_r), 3), "\n")
  
  # HL test
  if (run_hl) {
    hl <- tryCatch(ResourceSelection::hoslem.test(fit$y, fitted(fit)),
                   error = function(e) NULL)
    if (!is.null(hl))
      cat("Hosmer-Lemeshow p:", round(hl$p.value, 4), "\n")
  }
  
  # AUC
  if (run_roc) {
    auc_v <- pROC::auc(pROC::roc(fit$y, fitted(fit), quiet = TRUE))
    cat("AUC:", round(auc_v, 4), "\n")
  }
  
  # Plot
  print(plot_or(extract_or(fit), paste("OR:", outcome_col)))
  
  return(fit)
}

# Example
bin_res <- run_logistic_pipeline(data= sim_data,
  outcome_col= "disease",
  predictor_cols= c("exposure","age","bmi","smoking"))

## Ordinal logistic pipeline
run_ordinal_pipeline <- function(data, outcome_col, predictor_cols,link = "logit", seed = 123) {
  set.seed(seed)
  df<-data
  df[[outcome_col]]<-factor(df[[outcome_col]], ordered = TRUE)
  
  form <- as.formula(paste(outcome_col, "~", paste(predictor_cols, collapse = " + ")))
  fit<- ordinal::clm(form, data = df, link = link)
  
  cat("\nOrdinal CLM summary\n")
  print(summary(fit))
  
  # PO test
  cat("\nProportional odds test (nominal_test)\n")
  print(ordinal::nominal_test(fit))
  
  # OR table
  or_df <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(coef.type == "location")
  cat("\nOdds Ratios\n"); print(or_df)
  
  # AIC
  cat("AIC:", round(AIC(fit), 2), "\n")
  
  return(fit)
}

ord_res<-run_ordinal_pipeline(data= sim_data,
  outcome_col= "disease_sev",
  predictor_cols= c("exposure","age","bmi","smoking"))
