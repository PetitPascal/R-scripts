#-------------------------------------------------------------------------------
## Reproducible & generalisable G-computation script for a continuous outcome
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup

rm(list=ls())
set.seed(123)

required_pkgs<-c("MASS", "dplyr", "tidyr", "ggplot2", "boot","qgcomp", "margins", "broom","lmtest","zoo","performance")

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
y_cont<- -1 + 0.015 * cov_df$age + 0.4 * cov_df$sex + 0.05 * cov_df$bmi + h_z + rnorm(n, 0, 1.5)

sim_data<-cbind(data.frame(y=y_cont), cov_df,
                as.data.frame(log1p(expo_df)))  # log-transform exposures

#--------------------------------------------
## Step 3: Assumption checks

# Correlation among exposures
corr_expo<-cor(sim_data[, expo_names], method="spearman")
corr_expo

#--------------------------------------------
## Step 4: Fitting ourcome model + diagnostics

# Fitting outcome model
form<-as.formula(paste("y ~", paste(c(expo_names, "age", "sex", "bmi"), collapse=" + ")))
fit<-lm(form, data=sim_data)

# Checking model fit
summary(fit)

# Model characteristics
performance::performance(fit)

# Residual checks
par(mfrow=c(2, 2))
plot(fit, main="Linear model diagnostics")
par(mfrow=c(1, 1))
# performance::check_model(fit)

check_normality(fit) # Checking normality of residuals
lmtest::bptest(fit) # Checking homoscedasticity of residuals (Breusch-Pagan)
check_heteroscedasticity(fit) # Checking homoscedasticity of residuals
performance::check_model(fit) # checking all assumptions at once

#--------------------------------------------
## Step 5: g-computation: Q25 vs. Q75

q25<-sapply(sim_data[, expo_names], quantile, 0.25)
q75<-sapply(sim_data[, expo_names], quantile, 0.75)

d_lo<-sim_data
d_lo[, expo_names]<-matrix(q25, n, K, byrow=TRUE)
d_hi<-sim_data
d_hi[, expo_names]<-matrix(q75, n, K, byrow=TRUE)

# E[Y | Q25]
mu_lo<-mean(predict(fit, d_lo))
mu_lo

# E[Y | Q75]
mu_hi<-mean(predict(fit, d_hi))
mu_hi

# mean difference
diff<-mu_hi - mu_lo
diff

#--------------------------------------------
## Step 6: Bootstrap confidence intervals

boot_fn<-function(data, idx){
  d<-data[idx, ]
  fb<-lm(form, data=d)
  q25b<-sapply(d[, expo_names], quantile, 0.25)
  q75b<-sapply(d[, expo_names], quantile, 0.75)
  dl<-d
  dl[, expo_names]<-matrix(q25b, nrow(d), K, byrow=TRUE)
  dh<-d
  dh[, expo_names]<-matrix(q75b, nrow(d), K, byrow=TRUE)
  mean(predict(fb, dh)) - mean(predict(fb, dl))
}

set.seed(123)
boot_res<-boot::boot(sim_data, boot_fn, R=500)
ci<-boot::boot.ci(boot_res, type="perc")

cat("\nBootstrap 95% CI for mean difference:", round(diff, 4),
    "[", round(ci$percent[4], 4), ";", round(ci$percent[5], 4), "]\n")

#--------------------------------------------
## Step 7: Visualization

# Figure 1: Bootstrap distribution
ggplot(data.frame(diff=boot_res$t[, 1]), aes(x=diff)) +
  geom_histogram(fill="#A6DDCE", color="black", bins=40) +
  geom_vline(xintercept=diff, color="red", linetype="dashed", linewidth=1) +
  geom_vline(xintercept=c(ci$percent[4], ci$percent[5]),
             color="navy", linetype="dotted", linewidth=0.8) +
  labs(title="Bootstrap distribution: mean difference (g-computation)",
       subtitle=paste0("Estimate: ", round(diff, 4),
                         "  95% CI: [", round(ci$percent[4], 4), "; ", round(ci$percent[5], 4), "]"),
       x="Mean difference", y="Count") +
  theme_bw(base_size=14)

# Figure 2: Predicted E[Y] across joint quantile scenarios - marginal predictions plot across quantile scenarios
q_seq<-seq(0.10, 0.90, by=0.10)
pred_seq<-sapply(q_seq, function(q){
  d_q<-sim_data
  d_q[, expo_names]<-matrix(sapply(sim_data[, expo_names], quantile, q), n, K, byrow=TRUE)
  mean(predict(fit, d_q))
})

ggplot(data.frame(quantile=q_seq, mean_y=pred_seq), aes(quantile, mean_y)) +
  geom_line(linewidth=1.2, color="#2166ac") +
  geom_point(size=3, color="#2166ac") +
  scale_x_continuous(labels=scales::percent) +
  labs(title="G-computation: E[Y] across joint exposure quantiles",
       x="Joint quantile of all exposures", y="Predicted mean outcome") +
  theme_bw(base_size=14)

#--------------------------------------------
## Step 8: Exposure-specific weights

form_qgc<-as.formula(paste("y ~", paste(c(expo_names, "age", "sex", "bmi"), collapse=" + ")))
qgc_fit<-qgcomp::qgcomp.noboot(f=form_qgc, expnms=expo_names,data=sim_data, family=gaussian(), q=4)

# Extracting weights
weights_df<-data.frame(exposure=names(qgc_fit$pos.weights),
  weight=qgc_fit$pos.weights,
  direction="Positive")

if(length(qgc_fit$neg.weights) > 0){
  weights_df<-rbind(weights_df, data.frame(exposure =names(qgc_fit$neg.weights),
    weight=abs(qgc_fit$neg.weights),
    direction="Negative"))
}

# Plotting
ggplot(weights_df, aes(x=reorder(exposure, weight), y=weight, fill=direction)) +
  geom_col(color="black") +
  geom_text(aes(y=weight+0.05,label=signif(weight,3)))+
  scale_fill_manual(values=c("Positive"="#E87D7D", "Negative"="#A6DDCE")) +
  coord_flip() +
  labs(title="Exposure-specific weights in mixture",
       x="", y="Absolute weight", fill="Direction") +
  theme_bw(base_size=14)

# With bootstrap
qgc_boot<-qgcomp::qgcomp.boot(f=form_qgc, expnms=expo_names,
                                data=sim_data, family=gaussian(),
                                q=4, B=500, seed=123)
qgc_boot
