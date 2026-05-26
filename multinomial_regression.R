#-------------------------------------------------------------------------------
## Reproducible & Generalisable Multinomial Regression Script
# - Simulates data with a 4-category outcome
# - Fits multinomial logistic regression via nnet::multinom
# - Tests assumptions: IIA (Hausman-McFadden), multicollinearity, model fit
# - Computes relative risk ratios, marginal effects, visualises results
# - Includes run_multinom_pipeline() wrapper
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs<-c("nnet", "MASS", "dplyr", "ggplot2", "broom","tidyr", "car", "mlogit", "marginaleffects","gtsummary")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

#--------------------------------------------
## Step 2: Simulating dat
n<-800

age<-rnorm(n, 50, 12)
sex<-rbinom(n, 1, 0.5)
bmi<-rnorm(n, 26, 4)
score<-rnorm(n, 0, 1)
smoke<-rbinom(n, 1, 0.3)

# True log-odds for 4 outcome categories (ref = category 1)
# Using category-specific linear predictors
lp2<- -1.0 + 0.02*age + 0.4*sex - 0.05*bmi + 0.6*score
lp3<-0.5 - 0.01*age + 0.2*sex + 0.08*bmi - 0.3*score + 0.8*smoke
lp4<- -2.0 + 0.04*age - 0.3*sex + 0.10*bmi + 0.5*score + 1.2*smoke

# Softmax probabilities
denom<-1 + exp(lp2) + exp(lp3) + exp(lp4)
p1<-1/denom
p2<-exp(lp2)/denom
p3<-exp(lp3)/denom
p4<-exp(lp4)/denom

probs<-cbind(p1, p2, p3, p4)
y_cat<-apply(probs, 1, function(p) sample(1:4, 1, prob = p))

sim_data<-data.frame(
  outcome=factor(y_cat, levels = 1:4,
                   labels = c("cat1","cat2","cat3","cat4")),
  age=age, sex=factor(sex, labels = c("F","M")),
  bmi=bmi, score=score, smoke=factor(smoke, labels = c("No","Yes")))

cat("Outcome distribution:\n")
prop.table(table(sim_data$outcome))

#--------------------------------------------
## Step 3: Assumption checks

# 3a. Outcome category frequencies (each should be >= ~5% for stable estimation)
freq_tab <- prop.table(table(sim_data$outcome))
round(freq_tab, 3)
if (any(freq_tab < 0.05))
  warning("At least one category has < 5% of observations — estimates may be unstable.")

# 3b. Multicollinearity among predictors
vif_lm <- car::vif(lm(as.numeric(outcome) ~ age + sex + bmi + score + smoke,data = sim_data))
round(vif_lm, 3)

# 3c. Cell sizes: cross-tab predictors vs outcome (for categorical predictors)
table(sim_data$sex, sim_data$outcome)

#--------------------------------------------
## Step 4: Fit multinomial logistic regression
# Reference category = "cat1" (first level)
sim_data$outcome <- relevel(sim_data$outcome, ref = "cat1")

multinom_fit <- nnet::multinom(outcome ~ age + sex + bmi + score + smoke,
  data=sim_data,
  trace=FALSE,
  maxit=300)

summary(multinom_fit)

#--------------------------------------------
## Step 5: Wald tests and p-values
z_scores<-summary(multinom_fit)$coefficients/summary(multinom_fit)$standard.errors
p_values<-2 * (1 - pnorm(abs(z_scores)))
round(p_values, 4)

#--------------------------------------------
## Step 6: Relative Risk Ratios (RRR = exp(coef))
rrr_df <- broom::tidy(multinom_fit, exponentiate = TRUE, conf.int = TRUE)
rrr_df

# Forest-style plot of RRRs
rrr_df %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high,
             color = y.level)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_x_log10() +
  labs(title = "Multinomial regression: Relative Risk Ratios (vs. cat1)",
       x = "RRR (log scale)", y = "Predictor", color = "Category") +
  theme_bw(base_size = 13)

#--------------------------------------------
## Step 7: Model fitting

# Null model
null_fit<-nnet::multinom(outcome ~ 1, data = sim_data, trace = FALSE)

# Likelihood ratio test vs null
lr_stat<-2 * (logLik(multinom_fit) - logLik(null_fit))
lr_df<-attr(logLik(multinom_fit), "df") - attr(logLik(null_fit), "df")
lr_pval<-pchisq(as.numeric(lr_stat), df = lr_df, lower.tail = FALSE)

cat("\nLikelihood ratio test vs null model\n")
cat("Chi2 =", round(lr_stat, 3), "| df =", lr_df,
    "| p =", round(lr_pval, 6), "\n")

# McFadden pseudo-R2
mcfadden_r2 <- 1 - (as.numeric(logLik(multinom_fit)) /
                      as.numeric(logLik(null_fit)))
round(mcfadden_r2, 4)

# AIC / BIC
cat("AIC:", round(AIC(multinom_fit), 2), "\n")
cat("BIC:", round(BIC(multinom_fit), 2), "\n")

#--------------------------------------------
## Step 8: Test of Independence of Irrelevant Alternatives (IIA)

# Hausman-McFadden test: re-estimate dropping one category at a time
# If IIA holds, coefficients should not change significantly

iia_test<-function(full_model, data, drop_cat) {
  sub_data<-data %>% filter(outcome != drop_cat) %>%
    mutate(outcome = droplevels(outcome))
  sub_fit<-nnet::multinom(outcome ~ age + sex + bmi + score + smoke, data = sub_data, trace = FALSE, maxit = 300)
  
  # Comparing overlapping coefficients
  coef_full<-coef(full_model)
  coef_sub<-coef(sub_fit)
  
  # Keeping common rows
  common_cats<-intersect(rownames(coef_full), rownames(coef_sub))
  if (length(common_cats) == 0) return(NA)
  
  diff_coef<-coef_full[common_cats, ] - coef_sub[common_cats, ]
  cat("Drop category:", drop_cat,
      "| Max coef diff:", round(max(abs(diff_coef), na.rm=TRUE), 4), "\n")
}

for (cat_drop in c("cat2","cat3","cat4")) {
  tryCatch(iia_test(multinom_fit, sim_data, cat_drop),
           error = function(e) cat("IIA test failed for:", cat_drop, "\n"))
}

#--------------------------------------------
## Step 9: Marginal effects

# Average marginal effects
if (requireNamespace("marginaleffects", quietly = TRUE)) {
  me <- marginaleffects::avg_slopes(multinom_fit,
                                    variables = c("age","bmi","score"))
  cat("\nAverage marginal effects\n")
  print(me)
}

#--------------------------------------------
## Step 10: Predicted probabilities plot

# Predicting over range of 'score' holding others at median/mode
score_seq<-seq(min(sim_data$score), max(sim_data$score), length.out = 100)
pred_frame<-data.frame(age=median(sim_data$age),
  sex=factor("F", levels = levels(sim_data$sex)),
  bmi=median(sim_data$bmi),
  score=score_seq,
  smoke=factor("No", levels = levels(sim_data$smoke)))
pred_probs<-predict(multinom_fit, newdata = pred_frame, type = "probs")
pred_df<-cbind(pred_frame, as.data.frame(pred_probs)) %>%
  tidyr::pivot_longer(cols = starts_with("cat"),
                      names_to = "category", values_to = "prob")

ggplot(pred_df, aes(x = score, y = prob, color = category)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Multinomial: predicted probabilities vs. score",
       x = "Score", y = "Predicted probability", color = "Category") +
  theme_bw(base_size = 14)

#--------------------------------------------
## Step 11: summary table
tbl_multinom <- sim_data %>%
  gtsummary::tbl_uvregression(
    method= nnet::multinom,
    y =outcome,
    exponentiate=TRUE,
    hide_n= TRUE)
print(tbl_multinom)

#--------------------------------------------
## Step 12: Reusable pipeline
run_multinom_pipeline <- function(data, outcome_col, predictor_cols,
                                  ref_cat = NULL, seed = 2025) {
  set.seed(seed)
  
  df<-data
  df[[outcome_col]]<-factor(df[[outcome_col]])
  if (!is.null(ref_cat))
    df[[outcome_col]]<-relevel(df[[outcome_col]], ref = ref_cat)
  
  form<-as.formula(paste(outcome_col, "~",
                           paste(predictor_cols, collapse = " + ")))
  
  fit<-nnet::multinom(form, data = df, trace = FALSE, maxit = 300)
  null<-nnet::multinom(as.formula(paste(outcome_col, "~ 1")),
                         data = df, trace = FALSE)
  
  # Frequency check
  freq<-prop.table(table(df[[outcome_col]]))
  if(any(freq < 0.05)) message("WARNING: category with < 5% observations.")
  
  # VIF
  lm_proxy<-lm(as.numeric(as.factor(df[[outcome_col]])) ~.,
                 data = df[, predictor_cols])
  vif_res<-tryCatch(car::vif(lm_proxy), error = function(e) NULL)
  
  # Results
  rrr<-broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  lr_stat<-2 * (logLik(fit) - logLik(null))
  pseudo_r2<-1 - as.numeric(logLik(fit)) / as.numeric(logLik(null))
  
  cat("\nMcFadden R2:", round(pseudo_r2, 4),
      " | AIC:", round(AIC(fit), 2),
      " | BIC:", round(BIC(fit), 2), "\n")
  
  return(list(fit = fit, rrr = rrr, vif = vif_res,pseudo_r2 = pseudo_r2))
}

res <- run_multinom_pipeline(data=sim_data,
  outcome_col="outcome",
  predictor_cols=c("age","sex","bmi","score","smoke"),
  ref_cat="cat1")
