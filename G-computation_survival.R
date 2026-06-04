#-------------------------------------------------------------------------------
## Reproducible & generalisable G-computation script for a survival outcome (time-to-event)
#-------------------------------------------------------------------------------

#--------------------------------------------
## Step 1: Setup

rm(list=ls())
set.seed(123)

required_pkgs<-c("MASS", "dplyr", "tidyr", "ggplot2", "boot",  "survival", "flexsurv", "scales")

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
linpred<-0.02 * cov_df$age + 0.3 * cov_df$sex + 0.04 * cov_df$bmi + h_z

# Simulating Weibull event times + censoring
scale_par<-exp(-linpred/2)
t_event<-rweibull(n, shape=1.5, scale=scale_par)
t_censor<-runif(n, 5, 15)
time<-pmin(t_event, t_censor)
event<-as.integer(t_event<=t_censor)

sim_data<-cbind(data.frame(time=time, event=event),cov_df,as.data.frame(log1p(expo_df)))

#--------------------------------------------
## Step 3: Assumption checks

cat("Event rate:", round(mean(sim_data$event), 3), "\n")
if (mean(sim_data$event) < 0.05) warning("Very few events (<5%) — model estimates may be imprecise.")

cat("Median follow-up time:", round(median(sim_data$time), 2), "\n")
cat("Max pairwise correlation:", round(max(abs(cor(sim_data[, expo_names], method="spearman")[upper.tri(cor(sim_data[, expo_names]))])), 3), "\n")

#--------------------------------------------
## Step 4 - Parametric survival model (Weibull via flexsurv)

# Fitting the model
form_surv<-as.formula(paste("Surv(time, event) ~",paste(c(expo_names, "age", "sex", "bmi"), collapse=" + ")))
fit_weibull<-flexsurv::flexsurvreg(form_surv, data=sim_data, dist="weibull")
fit_weibull

# Optional: comparing with exponential and log-normal
fit_exp<-flexsurv::flexsurvreg(form_surv, data=sim_data, dist="exp")
fit_lnorm<-flexsurv::flexsurvreg(form_surv, data=sim_data, dist="lnorm")

aic_df<-data.frame(Distribution=c("Weibull", "Exponential", "Log-normal"), AIC=c(AIC(fit_weibull), AIC(fit_exp), AIC(fit_lnorm)))
aic_df
best_model<-fit_weibull # change if another fits better

#--------------------------------------------
## Step 5 - G-computation: cumulative risk at time t_eval (predicting S(t) for each person under Q25 and Q75 scenarios, then computing marginal cumulative risk=1 - mean(S(t)) )

t_eval<-median(sim_data$time) # time horizon for cumulative risk, adjust as needed
cat("\nTime horizon for g-computation:", round(t_eval, 2), "\n")

q25<-sapply(sim_data[, expo_names], quantile, 0.25)
q75<-sapply(sim_data[, expo_names], quantile, 0.75)

d_lo<-sim_data
d_lo[, expo_names]<-matrix(q25, n, K, byrow=TRUE)
d_hi<-sim_data
d_hi[, expo_names]<-matrix(q75, n, K, byrow=TRUE)

# Predicting survival probability S(t_eval) for each person under each scenario
get_survival_prob<-function(model, newdata, t){
  sapply(seq_len(nrow(newdata)), function(i){
    s<-summary(model, newdata=newdata[i, ], t=t, type="survival")
    s[[1]]$est
  })
}

S_lo<-get_survival_prob(best_model, d_lo, t_eval)
S_hi<-get_survival_prob(best_model, d_hi, t_eval)

risk_lo<-1-mean(S_lo)
risk_hi<-1-mean(S_hi)
RD_surv<-risk_hi-risk_lo
RR_surv<-risk_hi/risk_lo

cat("\nG-computation (survival) at t =", round(t_eval, 2), ":\n")
cat("Cumulative risk P(T<=t | Q25):", round(risk_lo, 4), "\n")
cat("Cumulative risk P(T<=t | Q75):", round(risk_hi, 4), "\n")
cat("Risk difference:", round(RD_surv, 4), "\n")
cat("Risk ratio:", round(RR_surv,  4), "\n")

#--------------------------------------------
## Step 6 - Bootstrap confidence intervals

boot_fn_surv<-function(data, idx){
  d<-data[idx, ]
  fb<-tryCatch(flexsurv::flexsurvreg(form_surv, data=d, dist="weibull"),error=function(e) NULL)
  if(is.null(fb)) return(c(NA, NA))
  
  q25b<-sapply(d[, expo_names], quantile, 0.25)
  q75b<-sapply(d[, expo_names], quantile, 0.75)
  dl<-d
  dl[, expo_names]<-matrix(q25b, nrow(d), K, byrow=TRUE)
  dh<-d
  dh[, expo_names]<-matrix(q75b, nrow(d), K, byrow=TRUE)
  
  S_lo_b<-get_survival_prob(fb, dl, t_eval)
  S_hi_b<-get_survival_prob(fb, dh, t_eval)
  r_lo_b<-1-mean(S_lo_b)
  r_hi_b<-1-mean(S_hi_b)
  c(RD=r_hi_b-r_lo_b, RR=r_hi_b/r_lo_b)
}

set.seed(123)
boot_res<-boot::boot(sim_data, boot_fn_surv, R=20) # use 500+ for publication
ci_RD<-boot::boot.ci(boot_res, index=1, type="perc")
ci_RR<-boot::boot.ci(boot_res, index=2, type="perc")

cat("\nBootstrap 95% CIs:\n")
cat("Risk difference:", round(RD_surv, 4), "[", round(ci_RD$percent[4], 4), ";", round(ci_RD$percent[5], 4), "]\n")
cat("Risk ratio:", round(RR_surv, 4), "[", round(ci_RR$percent[4], 4), ";", round(ci_RR$percent[5], 4), "]\n")

#--------------------------------------------
## Step 7: Visualization

# Figure 1: Bootstrap distribution of risk difference
boot_df<-data.frame(RD=boot_res$t[, 1]) %>% filter(!is.na(RD))
ggplot(boot_df, aes(x=RD)) +
  geom_histogram(fill="#A6DDCE", color="black", bins=35) +
  geom_vline(xintercept=RD_surv, color="red", linetype="dashed", linewidth=1) +
  geom_vline(xintercept=c(ci_RD$percent[4], ci_RD$percent[5]),
             color="navy", linetype="dotted", linewidth=0.8) +
  labs(title="Bootstrap distribution: risk difference (survival g-computation)",
       subtitle=paste0("t=", round(t_eval, 1),
                         "  RD=", round(RD_surv, 4),
                         "  95% CI: [", round(ci_RD$percent[4], 4),
                         "; ", round(ci_RD$percent[5], 4), "]"),
       x="Risk difference", y="Count") +
  theme_bw(base_size=14)

# Figure 2: Cumulative risk across joint quantile scenarios
q_seq<-seq(0.10, 0.90, by=0.10)
risk_seq<-sapply(q_seq, function(q){
  d_q<-sim_data
  d_q[, expo_names]<-matrix(sapply(sim_data[, expo_names], quantile, q), n, K, byrow=TRUE)
  S_q<-get_survival_prob(best_model, d_q, t_eval)
  1-mean(S_q)
})

ggplot(data.frame(quantile=q_seq, cum_risk=risk_seq), aes(quantile, cum_risk)) +
  geom_line(linewidth=1.2, color="#d6604d") +
  geom_point(size=3, color="#d6604d") +
  scale_x_continuous(labels=scales::percent) +
  scale_y_continuous(labels=scales::percent) +
  labs(title=paste0("G-computation: cumulative risk at t=", round(t_eval, 1),
                      " across joint exposure quantiles"),
       x="Joint quantile of all exposures", y="Cumulative risk") +
  theme_bw(base_size=14)

# Figure 3: Marginal survival curves under Q25 vs. Q75 scenarios
t_grid<-seq(0, max(sim_data$time), length.out=100)

get_marginal_survival<-function(model, scenario_data, t_grid){
  sapply(t_grid, function(t){
    s_vec<-get_survival_prob(model, scenario_data, t)
    mean(s_vec)
  })
}

S_curve_lo<-get_marginal_survival(best_model, d_lo, t_grid)
S_curve_hi<-get_marginal_survival(best_model, d_hi, t_grid)

surv_df<-data.frame(time=rep(t_grid, 2),
  survival=c(S_curve_lo, S_curve_hi),
  scenario=rep(c("Q25 (low exposure)", "Q75 (high exposure)"), each=length(t_grid)))

ggplot(surv_df, aes(x=time, y=survival, color=scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=c("Q25 (low exposure)"="#A6DDCE", "Q75 (high exposure)"="#d6604d")) +
  scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
  labs(title="Marginal survival curves: Q25 vs. Q75 (g-computation)",
       x="Time", y="Survival probability", color="") +
  theme_bw(base_size=14) +
  theme(legend.position="bottom")
