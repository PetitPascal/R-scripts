#-------------------------------------------------------------------------------
## Reproducible & generalisable G-computation script for a binary outcome
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup

rm(list=ls())
set.seed(123)

required_pkgs<-c("MASS", "dplyr", "tidyr", "ggplot2", "boot", "qgcomp")

is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(required_pkgs[!is_installed],repos="http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only=TRUE))

#--------------------------------------------
## Step 2: Simulating data (replace with your own)

n<-600
K<-6
rho<-0.45
Sigma<-rho^as.matrix(dist(1:K))
expo_df<-as.data.frame(exp(MASS::mvrnorm(n, mu=rep(0, K), Sigma=Sigma)/3))
colnames(expo_df)<-paste0("X", seq_len(K))
expo_names<-colnames(expo_df)

set.seed(123)
cov_df<-data.frame(age=rnorm(n, 50, 10),
                   sex=rbinom(n, 1, 0.5),
                   bmi=rnorm(n, 26, 4))

h_z<-as.numeric(as.matrix(log1p(expo_df)) %*% c(0.4, 0.3, 0.2, -0.1, 0.05, 0.05))
y_bin<-rbinom(n, 1, plogis(-1 + 0.015 * cov_df$age + 0.4 * cov_df$sex + 0.05 * cov_df$bmi + h_z))

sim_data<-cbind(data.frame(y=y_bin), cov_df, log1p(expo_df))

#--------------------------------------------
## Step 3: Assumption checks

# Positivity: ensuring overlap in exposure distribution across outcome groups
pos_check<-sim_data %>% group_by(y) %>% summarise(across(all_of(expo_names), mean))
pos_check

# Correlation among exposures
corr_expo<-cor(sim_data[, expo_names], method="spearman")
round(max(abs(corr_expo[upper.tri(corr_expo)])), 3)

# Checking for near-perfect separation
tab_check<-table(sim_data$y)
if(min(tab_check)<20) warning("Very few events/non-events — model may not converge reliably.")
min(tab_check)

#--------------------------------------------
## Step 4 - Fit logistic model + diagnostics

# Fitting outcome model
form<-as.formula(paste("y ~", paste(c(expo_names, "age", "sex", "bmi"), collapse=" + ")))
fit<-glm(form, data=sim_data, family=binomial())

# Checking model fit
summary(fit)

# Model characteristics
performance::performance(fit)

# Checking overdispersion
summary(fit)$deviance/summary(fit)$df.residual
check_overdispersion(fit)

# checking All assumptions at once
performance::check_model(fit) 

# Hosmer-Lemeshow test
hl_test<-function(model, g=10){
  obs<-model$y
  pred<-fitted(model)
  grp<-cut(pred, breaks=quantile(pred, probs=seq(0, 1, 1/g)), include.lowest=TRUE)
  tab<-tapply(obs, grp, function(x) c(obs=sum(x), exp=length(x) * mean(x), n=length(x)))
  chi2<-sum(sapply(tab, function(x) (x["obs"] - x["exp"])^2 / (x["exp"] * (1 - x["exp"]/x["n"]))))
  p<-pchisq(chi2, df=g - 2, lower.tail=FALSE)
  cat("\nHosmer-Lemeshow: chi2 =", round(chi2, 3), ", p =", round(p, 3), "(p>0.05=adequate fit)\n")
}
hl_test(fit)

#--------------------------------------------
## Step 5 - g-computation: Q25 vs. Q75

q25<-sapply(sim_data[, expo_names], quantile, 0.25)
q75<-sapply(sim_data[, expo_names], quantile, 0.75)

d_lo<-sim_data
d_lo[, expo_names]<-matrix(q25, n, K, byrow=TRUE)
d_hi<-sim_data
d_hi[, expo_names]<-matrix(q75, n, K, byrow=TRUE)

p_lo<-mean(predict(fit, d_lo, type="response"))
p_lo # P(Y=1 | Q25)
p_hi<-mean(predict(fit, d_hi, type="response"))
p_hi # P(Y=1 | Q75)
RD<-p_hi-p_lo # Risk difference
RD
RR<-p_hi/p_lo # Risk ratio
RR

#--------------------------------------------
## Step 6 - Bootstrap confidence intervals

boot_fn<-function(data, idx){
  d<-data[idx, ]
  fb<-glm(form, data=d, family=binomial())
  q25b<-sapply(d[, expo_names], quantile, 0.25)
  q75b<-sapply(d[, expo_names], quantile, 0.75)
  dl<-d
  dl[, expo_names]<-matrix(q25b, nrow(d), K, byrow=TRUE)
  dh<-d
  dh[, expo_names]<-matrix(q75b, nrow(d), K, byrow=TRUE)
  p_lo_b<-mean(predict(fb, dl, type="response"))
  p_hi_b<-mean(predict(fb, dh, type="response"))
  c(RD=p_hi_b-p_lo_b, RR=p_hi_b/p_lo_b)
}

set.seed(123)
boot_res<-boot::boot(sim_data, boot_fn, R=500)
ci_RD<-boot::boot.ci(boot_res, index=1, type="perc")
ci_RR<-boot::boot.ci(boot_res, index=2, type="perc")

cat("\nBootstrap 95% CIs:\n")
cat("Risk difference:", round(RD, 4), "[", round(ci_RD$percent[4], 4), ";", round(ci_RD$percent[5], 4), "]\n")
cat("Risk ratio:", round(RR, 4), "[", round(ci_RR$percent[4], 4), ";", round(ci_RR$percent[5], 4), "]\n")

#--------------------------------------------
## Step 7: Visualization

# Figure 1: Bootstrap distribution on risk difference
ggplot(data.frame(RD=boot_res$t[, 1]), aes(x=RD)) +
  geom_histogram(fill="#A6DDCE", color="black", bins=40) +
  geom_vline(xintercept=RD, color="red", linetype="dashed", linewidth=1) +
  geom_vline(xintercept=c(ci_RD$percent[4], ci_RD$percent[5]),
             color="navy", linetype="dotted", linewidth=0.8) +
  labs(title="Bootstrap distribution: risk difference",
       subtitle=paste0("RD=", round(RD, 3)," [95% CI: ", round(ci_RD$percent[4], 3), "; ", round(ci_RD$percent[5], 3), "]"),
       x="Risk difference", y="Count") +
  theme_bw(base_size=14)

# Figure 2: Predicted P(Y=1) across joint quantiles
q_seq<-seq(0.10, 0.90, by=0.10)
pred_seq<-sapply(q_seq, function(q) {
  d_q<-sim_data
  d_q[, expo_names]<-matrix(sapply(sim_data[, expo_names], quantile, q), n, K, byrow=TRUE)
  mean(predict(fit, d_q, type="response"))
})

ggplot(data.frame(quantile=q_seq, prob=pred_seq), aes(quantile, prob)) +
  geom_line(linewidth=1.2, color="#2166ac") +
  geom_point(size=3, color="#2166ac") +
  geom_hline(yintercept=mean(sim_data$y), linetype="dashed", color="gray50") +
  scale_x_continuous(labels=scales::percent) +
  labs(title="G-computation: P(Y=1) across joint exposure quantiles",
       x="Joint quantile", y="Predicted probability") +
  theme_bw(base_size=14)

#--------------------------------------------
## Step 8: - Exposure-specific weights

form_qgc<-as.formula(paste("y ~", paste(c(expo_names, "age", "sex", "bmi"), collapse=" + ")))
qgc_fit<-qgcomp::qgcomp.noboot(f=form_qgc, expnms=expo_names, data=sim_data, family=binomial(), q=4)

# Extracting weights
weights_df<-data.frame(exposure=names(qgc_fit$pos.weights),
  weight=qgc_fit$pos.weights,
  direction="Positive")

if(length(qgc_fit$neg.weights)>0){
  weights_df<-rbind(weights_df, data.frame(exposure=names(qgc_fit$neg.weights),
    weight=abs(qgc_fit$neg.weights),
    direction="Negative"))
}

# Plotting
ggplot(weights_df, aes(x=reorder(exposure, weight), y=weight, fill=direction)) +
  geom_col(color="black") +
  scale_fill_manual(values=c("Positive"="#E87D7D", "Negative"="#A6DDCE")) +
  coord_flip() +
  labs(title="Exposure-specific weights in mixture",
       x="", y="Absolute weight", fill="Direction") +
  theme_bw(base_size=14)

# With bootstrap
qgc_boot<-qgcomp::qgcomp.boot(f=form_qgc, expnms=expo_names,
                                data=sim_data, family=binomial(),
                                q=4, B=500, seed=123)
qgc_boot
