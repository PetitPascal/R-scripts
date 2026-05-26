#-------------------------------------------------------------------------------
## Reproducible & Generalisable Linear Mixed-Effect Model (LMM) Script
# Covers: LMM, robust LMM, one-way mixed ANOVA, two-way mixed ANOVA
# - Uses simulated RCT-style longitudinal dataset
# - Tests assumptions: normality, homoscedasticity, random effect structure
# - Computes emmeans, pairwise contrasts, effect sizes
# - Includes run_lmm_pipeline() wrapper
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup
rm(list = ls())
set.seed(123)

required_pkgs <- c("lme4", "lmerTest", "emmeans", "broom.mixed","performance", "robustlmm", "dplyr", "ggplot2",
                   "tidyr", "car", "effectsize", "gtsummary","conflicted")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# for avoiding package conflict
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(lme4::lmer)

#--------------------------------------------
## Step 2: Simulating longitudinal repeated data
n_id<-180
n_time<-4
arms<-c("A","B")

ID<-seq_len(n_id)
arm<-sample(arms, n_id, replace = TRUE)
dat<-expand.grid(ID = ID, time = 0:(n_time-1))
dat<-merge(dat, data.frame(ID=ID, arm=arm), by="ID")

# Random intercept per subject
u_i<-rnorm(n_id, 0, 2)[match(dat$ID, ID)]

# Continuous covariate
age_id<-rnorm(n_id, 50, 10)[match(dat$ID, ID)]
dat$age<-age_id

# Primary continuous outcome
dat$y<-with(dat,10 +
                ifelse(arm=="B", 1.5, 0) +
                0.4 * time +
                ifelse(arm=="B", 0.3, 0) * time +
                0.05 * age +
                u_i +
                rnorm(nrow(dat), 0, 1.5))

# Factor variables
dat$ID<-factor(dat$ID)
dat$arm<-factor(dat$arm, levels = arms)
dat$time<-factor(dat$time, levels = 0:3,
                   labels = paste0("T", 0:3))

df <- dplyr::as_tibble(dat)

#--------------------------------------------
## Step 3: EDA

df %>%
  group_by(arm, time) %>%
  summarise(mean = mean(y), se = sd(y)/sqrt(n()), .groups="drop") %>%
  ggplot(aes(x=time, y=mean, color=arm, group=arm)) +
  geom_line(linewidth=1.1) + geom_point(size=3) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.15) +
  labs(title="Observed means ± SE by arm and time",
       x="Time", y="Mean outcome") +
  theme_bw(base_size=14)

#--------------------------------------------
## Step 4 - Model fitting

#- - - -
## Model 1 - One-Way Mixed ANOVA (between-subject factor: arm/group; within-subject factor: time)

lmm_1way<-lmerTest::lmer(y ~ arm + (1 | ID), data= df, REML= FALSE)
summary(lmm_1way)

# ANOVA table (Type III)
anova(lmm_1way, type="III")

# Pairwise contrasts between groups using marginal means (emmeans)
emm_arm <- emmeans::emmeans(lmm_1way, ~ arm)
pairs(emm_arm, adjust="BH")

# Effect size (Cohen's d)
effectsize::cohens_d(y ~ arm, data=df)

#- - - -
## Model 2 - Two-Way Mixed ANOVA (arm/group × time interaction))

lmm_2way<-lmerTest::lmer(y ~ arm * time + (1 | ID),data  = df,REML  = FALSE)

# ANOVA table (Type III)
anova(lmm_2way, type="III")

# Interaction p-value interpretation
anova_tab<-as.data.frame(anova(lmm_2way, type="III"))
int_p<-anova_tab["arm:time", "Pr(>F)"]
int_p

#- - - -
## Model 3 - Full LMM with covariates

lmm_full<-lmerTest::lmer(y ~ arm * time + age + (1 | ID),data  = df,REML  = FALSE)
summary(lmm_full)
broom.mixed::tidy(lmm_full, effects="fixed", conf.int=TRUE)

## Assumption checks

# normality of residuals
shapiro.test(residuals(lmm_full))

# normality of random effects
re_vals<-unlist(lme4::ranef(lmm_full)$ID)
shapiro.test(re_vals)

# homoscedasticity
performance::check_heteroscedasticity(lmm_full)

# independence
performance::check_autocorrelation(lmm_full)

# visual diagnostics
par(mfrow=c(1,2))
qqnorm(residuals(lmm_full), main="QQ-plot: residuals")
qqline(residuals(lmm_full), col="red")
qqnorm(re_vals, main="QQ-plot: random effects")
qqline(re_vals, col="red")
par(mfrow=c(1,1))

# random effects structure comparison (random slope vs. intercept only)
lmm_rs <- tryCatch(
  lmerTest::lmer(y ~ arm*time + age + (time | ID), data=df, REML=FALSE),
  error = function(e) {
    message("Random slope model failed: ", e$message); NULL
  }
)
if (!is.null(lmm_rs)) {
  anova_compare <- anova(lmm_full, lmm_rs)
  cat("\nLRT: random intercept vs random slope\n")
  print(anova_compare)
}

# Intra-class correlation (ICC)
icc_val <- performance::icc(lmm_full)
cat("\nICC:", round(icc_val$ICC_adjusted, 4),
    "(proportion of variance due to subjects)\n")

## emmeans and pairwise contrasts

# Arm/group differences at each time point
emm_arm_time<-emmeans::emmeans(lmm_full, ~ arm | time)
pairs(emm_arm_time, adjust="BH")

# Time differences within each arm/group
emm_time_arm <- emmeans::emmeans(lmm_full, ~ time | arm)
pairs(emm_time_arm, adjust="BH")

# Collecting emmeans for plotting
emm_plot_df <- as.data.frame(emm_arm_time) %>% rename(emmean_se = SE, lo = lower.CL, hi = upper.CL)

ggplot(emm_plot_df, aes(x=time, y=emmean, color=arm, group=arm,
                        ymin=lo, ymax=hi)) +
  geom_line(linewidth=1.1) +
  geom_point(size=3) +
  geom_errorbar(width=0.15) +
  labs(title="LMM: estimated marginal means ± 95% CI",
       x="Time", y="Estimated marginal mean") +
  theme_bw(base_size=14)

#- - - -
## Model 4 - Robust LMM (if LLM assumptions violated)

rlmm_fit<-robustlmm::rlmer(y ~ arm * time + age + (1 | ID),data = df)
summary(rlmm_fit)

#--------------------------------------------
## Step 5 - Model comparison
lmm_null<-lmerTest::lmer(y ~ time + age + (1|ID), data=df, REML=FALSE)
lrt<-anova(lmm_null, lmm_full)
lrt

cat("\nAIC comparison:\n")
cat("Null:", round(AIC(lmm_null),2), "| Full:", round(AIC(lmm_full),2), "\n")

#--------------------------------------------
## Step 6 - Reusable pipeline
run_lmm_pipeline<-function(data, outcome_col, fixed_effects,
                             random_intercept = "ID",
                             random_slope= NULL,
                             REML = FALSE,
                             run_robust= TRUE) {
  
  re_term <- if (!is.null(random_slope))
    paste0("(", random_slope, " | ", random_intercept, ")")
  else
    paste0("(1 | ", random_intercept, ")")
  
  form <- as.formula(paste(
    outcome_col, "~",
    paste(fixed_effects, collapse=" + "),
    "+", re_term))
  
  fit<-lmerTest::lmer(form, data=data, REML=REML)
  
  # Assumptions
  sw_p<-shapiro.test(residuals(fit))$p.value
  icc_v<-performance::icc(fit)$ICC_adjusted
  
  cat("\nLMM Pipeline Results\n")
  cat("Shapiro-Wilk p:", round(sw_p, 4),
      if(sw_p < 0.05) "*** NON-NORMAL ***" else "OK", "\n")
  cat("ICC:", round(icc_v, 4), "\n")
  print(anova(fit, type="III"))
  
  res<-list(fit=fit, icc=icc_v, sw_p=sw_p)
  
  if (run_robust) {
    rfit <- tryCatch(
      robustlmm::rlmer(form, data=data),
      error=function(e){message("Robust LMM failed."); NULL}
    )
    res$robust_fit <- rfit
  }
  
  return(res)
}

lmm_res<-run_lmm_pipeline(data= df,
  outcome_col= "y",
  fixed_effects = c("arm*time","age"),
  random_intercept= "ID",
  run_robust= TRUE)
