#------------------------------------------------------------------
#### Configurations  ####
#------------------------------------------------------------------

## Disabling memory torture
gctorture(FALSE)

## Installing packages if needed
required_packages <- c(
  "tidyverse", "lme4", "lmerTest", "emmeans", "nlme", "mice",
  "broom.mixed", "pwr", "janitor", "gt", "gtsummary", "performance","medicaldata",
  "ordinal", "geepack","conflicted", "patchwork","DataExplorer","fastDummies","flextable"
)

for (i in 1:length(required_packages)){
  if(required_packages[i]%in%.packages(all.available=TRUE)){
  }else{
    install.packages(required_packages[i])
  }
}

## Loading packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(nlme)
library(mice)
library(broom.mixed)
library(pwr)
library(janitor)
library(gt)
library(gtsummary)
library(performance)
library(ordinal)
library(geepack)
library(conflicted)
library(patchwork)
library(DataExplorer)
library(fastDummies)
library(flextable)

## Preventing package conflicts
conflict_prefer("select", "dplyr") 
conflict_prefer("filter", "dplyr") 
conflict_prefer("slice", "dplyr")
conflict_prefer("alpha", "scales")
conflict_prefer("intersect", "base")
conflict_prefer("lme4","lmer")
conflicted::conflicts_prefer(lme4::lmer)

## Setting work directory
here::here("RCT analysis")

#------------------------------------------------------------------
#### Data import ####
#------------------------------------------------------------------

### Creating a simulated RCT dataset
set.seed(123)

## Design
n_id<-200        # number of participants
n_time<-3        # time points, e.g. 0, 1, 2
arms<-c("A", "B", "C")

## Creating subject-level data with randomized arms
ID<-1:n_id
arm<-sample(arms, n_id, replace = TRUE)

## Expanding to long format: each ID appears at each time
dat <- expand.grid(ID=ID,time=0:(n_time - 1))

# Merging arm information
dat<-merge(dat, data.frame(ID=ID, arm=arm), by="ID")

## Simulating continuous variables

# variable 1
dat$cont2<-abs(rnorm(n_id)[match(dat$ID, unique(dat$ID))])  # subject-level covariate
dat$cont2<-dat$cont2[match(dat$ID, dat$ID)]                 # replicate across time

dat$cont1_mean <- with(dat, 10 +
                         ifelse(arm=="B", 1, 0) +
                         ifelse(arm=="C", 2, 0) +
                         0.5*time +
                         0.3*cont2)

dat$cont1 <- rnorm(nrow(dat), mean=dat$cont1_mean, sd=2)

# variable 2
dat$cont3 <- with(dat,rnorm(nrow(dat),
                            mean=5 + 0.2*time - 0.5*(arm=="C"),
                            sd=1.5))

## Simulating binary variables

# variable 1
logit_p1 <- with(dat,
                 -1 +
                   0.4*time +
                   0.3*(arm=="B") +
                   0.6*(arm=="C"))
p1 <- plogis(logit_p1)
dat$bin1 <- rbinom(nrow(dat), size = 1, prob = p1)

# variable 2
sex_id <- rbinom(n_id, 1, 0.5)
names(sex_id) <- ID
dat$bin2 <- sex_id[as.character(dat$ID)]

## Simulating ordinal variables
dat<-dat %>% mutate(ord1=case_when(cont1<11~1,
                                   cont1>=13.3~3,
                                   T~2),
                    ord2=case_when(cont3<4~1,
                                   cont3>6~3,
                                   T~2))

## Creating the final simulated dataset
sim_rct <- dat[, c("ID", "arm", "time",
                   "cont1", "cont3",
                   "bin1", "bin2",
                   "ord1","ord2")]

df<-as_tibble(sim_rct)

## Creating a dataset with the variable type
var_type<-as_tibble(data.frame(Var=colnames(df),type=c("nominal","binary","continuous","nominal","nominal",
                                                       "continuous","binary","ordinal","ordinal")))

#------------------------------------------------------------------
#### Descriptive analysis ####
#------------------------------------------------------------------

## Adding an overall group
df2<-df
df2$arm<-as.character(df2$arm)
df3<-df2 %>% mutate(arm="overall")
df2<-bind_rows(df2,df3)

## Continuous and ordinal variable characteristics
charact1<-df2 %>% 
          select(-ID) %>%
          select(all_of(c("arm","time",var_type %>% filter(type %in% c("continuous","ordinal")) %>% select(Var) %>% pull))) %>%
          pivot_longer(cols=-c(arm,time),names_to = "Variable",values_to="Value") %>%
          group_by(time,arm,Variable) %>%
          summarize(n_miss=length(Value)-length(na.omit(Value)),
                    mean=signif(mean(Value,na.rm=T),3),
                    sd=signif(sd(Value,na.rm=T),3),
                    median=signif(median(Value,na.rm=T),3),
                    Q1=signif(quantile(Value,0.25,na.rm=T),3),
                    Q3=signif(quantile(Value,0.75,na.rm=T),3),
                    min=signif(min(Value,na.rm=T),3),
                    max=signif(max(Value,na.rm=T),3)) %>%
          ungroup %>%
          mutate(Label=paste(mean," (",sd,"); ",median," [",Q1,"; ",Q3,"]; [",min,"; ",max,"]; n missing=",n_miss,sep="")) %>%
          select(time,arm,Variable,Label) %>%
          pivot_wider(names_from=arm, values_from=Label, names_prefix="arm_")

## Binary and ordinal variable characteristics
charact2<-df2 %>% select(all_of(c(var_type %>% filter(type!="continuous") %>% select(Var) %>% pull)))
charact2<-update_columns(data=charact2,
                   ind=colnames(charact2),
                   what=as.character)

charact2<-charact2 %>% 
  select(-ID) %>%
  pivot_longer(cols=-c(arm,time),names_to = "Variable",values_to="Value") %>%
  group_by(time,arm,Variable,Value) %>% 
  count %>%
  group_by(time,arm,Variable) %>%
  mutate(per=signif(n/sum(n)*100,3)) %>%
  ungroup %>%
  mutate(Label=paste(n," (",per,")",sep="")) %>%
  select(time,arm,Variable,Value,Label) %>%
  pivot_wider(names_from=arm, values_from=Label, names_prefix="arm_",values_fill = "0 (0%)")

## Exporting results
write.table(charact1,paste("Descriptive analysis_continuous var_",Sys.Date(),".csv",sep=""),sep=";",col.names=T,row.names=F)
write.table(charact2,paste("Descriptive analysis_categorical_",Sys.Date(),".csv",sep=""),sep=";",col.names=T,row.names=F)

#------------------------------------------------------------------
#### Data processing ####
#------------------------------------------------------------------

## Creating variable lists based on variable type
list_continuous_var<-var_type %>% filter(type %in% c("continuous")) %>% select(Var) %>% distinct %>% pull
list_binary_var<-var_type %>% filter(type %in% c("binary")) %>% select(Var) %>% distinct %>% pull
list_ordinal_var<-var_type %>% filter(type %in% c("ordinal")) %>% select(Var) %>% distinct %>% pull

## Ensuring variable type

# Ordered factors
df<-update_columns(data=df,
                   ind=list_ordinal_var,
                   what=function(x) factor(x,ordered = T,levels=sort(unique(x))))

# ID, time, arm
df$time <- factor(df$time, levels = c("0","1", "2"),labels=c("T0","T1","T2"))
df$ID <- factor(df$ID)
df$arm <- factor(df$arm, levels=sort(unique(df$arm)))

#------------------------------------------------------------------
#### Linear mixed effect model (continuous variables) ####
#------------------------------------------------------------------

Full_LLM_table<-c()
Full_robust_LLM_table<-c()
emmeans_LLM_table<-c()
robust_emmeans_LLM_table<-c()

for(i in 1:length(list_continuous_var)){
  
  # selecting variables
  df_tmp<-df %>% select(list_continuous_var[i],arm,time,ID)
  colnames(df_tmp)[1]<-"tmp"
  
  # creating the model
  model_lmm <- lmerTest::lmer(tmp ~ arm * time + (1 | ID), data = df_tmp, REML = FALSE)
  
  # full summary table
  lmm_results <- tidy(model_lmm, effects = "fixed", conf.int = TRUE)
  
  # testing model validity assumptions
  norm_check<-check_normality(model_lmm)[1] # normality of residuals
  var_check<-check_heteroscedasticity(model_lmm)[1] # homoscedasticity of residuals
  indep_check<-check_autocorrelation(model_lmm)[1] # independence of residuals
  
  # adding model test validity assumption
  lmm_results<-lmm_results %>% mutate(Variable=list_continuous_var[i],
                                        norm_check=norm_check,
                                        var_check=var_check,
                                        indep_check=indep_check)
  
  # saving results
  Full_LLM_table<-bind_rows(Full_LLM_table,lmm_results)
  
  #- - - - - - - - - - - - -
  ## Robust LMM in case of violation of one assumption
  robust_lmm <- robustlmm::rlmer(tmp ~ arm * time + (1 | ID), data = df_tmp)
  robust_lmm_results<-tidy(robust_lmm, effects = "fixed", conf.int = TRUE) %>% mutate(Variable=list_continuous_var[i])
  Full_robust_LLM_table<-bind_rows(Full_robust_LLM_table,robust_lmm_results)
  
  #- - - - - - - - - - - - -
  ### in case of interaction
  
  #- - - -
  ## Difference by time
  
  # Estimating marginal means (emmeans)
  emmeans_lmm <- emmeans(model_lmm, ~ arm | time)
  robust_emmeans_lmm <- emmeans(robust_lmm, ~ arm | time)
  
  # Pairwise contrasts between arms
  lmm_contrasts <- contrast(emmeans_lmm, method = "pairwise") %>% as_tibble() %>% 
                   mutate(Variable=list_continuous_var[i],
                          comparison="time")
  robust_lmm_contrasts <- contrast(robust_emmeans_lmm, method = "pairwise") %>% as_tibble() %>%
                          mutate(Variable=list_continuous_var[i],
                                 comparison="time")
  
  # saving results
  emmeans_lmm<-as_tibble(emmeans_lmm)
  colnames(emmeans_lmm)[c(4:7)]<-paste("emmean",colnames(emmeans_lmm)[c(4:7)],sep="_")
  robust_emmeans_lmm<-as_tibble(robust_emmeans_lmm)
  colnames(robust_emmeans_lmm)[c(4:7)]<-paste("emmean",colnames(robust_emmeans_lmm)[c(4:7)],sep="_")
  
  lmm_contrasts<-left_join(lmm_contrasts,emmeans_lmm,by="time")
  robust_lmm_contrasts<-left_join(robust_lmm_contrasts,robust_emmeans_lmm,by="time")
  
  emmeans_LLM_table<-bind_rows(emmeans_LLM_table,lmm_contrasts)
  robust_emmeans_LLM_table<-bind_rows(robust_emmeans_LLM_table,robust_lmm_contrasts)
  
  # Plot
  plot1<-df_tmp %>%
    group_by(arm, time) %>%
    summarise(mean_outcome = mean(tmp, na.rm = TRUE),
              se = sd(tmp, na.rm = TRUE)/sqrt(n())) %>%
    ggplot(aes(x = time, y = mean_outcome, color = arm, group = arm)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_outcome - se, ymax = mean_outcome + se), width = 0.1)+
    scale_y_continuous(list_continuous_var[i])+
    theme_bw()
  
  ggsave(plot1,
         file=paste("Interaction_",list_continuous_var[i],"_",
                    Sys.Date(),".pdf",sep=""),
         dpi=600,width=60,height=30,units = "cm",limitsize=F)
  
  #- - - -
  ##  Difference by arm
  
  # Estimating marginal means (emmeans)
  emmeans_lmm <- emmeans(model_lmm, ~ time | arm)
  robust_emmeans_lmm <- emmeans(robust_lmm, ~ time | arm)
  
  # Pairwise contrasts between arms
  lmm_contrasts <- contrast(emmeans_lmm, method = "pairwise") %>% as_tibble() %>% 
                   mutate(Variable=list_continuous_var[i],comparison="arm")
  robust_lmm_contrasts <- contrast(robust_emmeans_lmm, method = "pairwise") %>% as_tibble() %>%
    mutate(Variable=list_continuous_var[i],comparison="arm")
  
  # saving results
  emmeans_lmm<-as_tibble(emmeans_lmm)
  colnames(emmeans_lmm)[c(4:7)]<-paste("emmean",colnames(emmeans_lmm)[c(4:7)],sep="_")
  robust_emmeans_lmm<-as_tibble(robust_emmeans_lmm)
  colnames(robust_emmeans_lmm)[c(4:7)]<-paste("emmean",colnames(robust_emmeans_lmm)[c(4:7)],sep="_")
  
  lmm_contrasts<-left_join(lmm_contrasts,emmeans_lmm,by="arm")
  robust_lmm_contrasts<-left_join(robust_lmm_contrasts,robust_emmeans_lmm,by="arm")
  
  emmeans_LLM_table<-bind_rows(emmeans_LLM_table,lmm_contrasts)
  robust_emmeans_LLM_table<-bind_rows(robust_emmeans_LLM_table,robust_lmm_contrasts)

}

## Exporting results
write.table(Full_LLM_table,paste("LLM_",Sys.Date(),".csv",sep=""),sep=";",col.names=T,row.names=F)
write.table(Full_robust_LLM_table,paste("robut LLM_",Sys.Date(),".csv",sep=""),sep=";",col.names=T,row.names=F)
write.table(emmeans_LLM_table,paste("Pariwise LLM_",Sys.Date(),".csv",sep=""),sep=";",col.names=T,row.names=F)
write.table(robust_emmeans_LLM_table,paste("robust_pariwise LLM_",Sys.Date(),".csv",sep=""),sep=";",col.names=T,row.names=F)

#------------------------------------------------------------------
#### Logistic mixed effect model (binary variables) ####
#------------------------------------------------------------------

Full_GLM_table<-c()
Full_emmeans_GLM_table<-c()

for(i in 1:length(list_binary_var)){
  
  # selecting variables
  df_tmp<-df %>% select(list_binary_var[i],arm,time,ID)
  colnames(df_tmp)[1]<-"tmp"
  
  # creating the model
  model_glm <- glmer(tmp ~ arm * time + (1 | ID),family = binomial, data = df_tmp)
  
  # full summary table
  glm_results <- tidy(model_glm, effects = "fixed", conf.int = TRUE, exponentiate = T) %>%
                 mutate(Variable=list_binary_var[i])
  
  # saving results
  Full_GLM_table<-bind_rows(Full_GLM_table,glm_results)
  
  ## Pairwise comparison
  
  # Estimating marginal means (emmeans) - predicted probabilities
  emmeans_glm <- emmeans(model_glm, ~ arm * time, type = "response")
  
  # Pairwise contrasts between times
  time_compa<-pairs(emmeans(model_glm, ~ arm | time), adjust = "BH") %>% as_tibble() %>% 
              mutate(Variable=list_binary_var[i],compa="time")
  
  # Pairwise contrasts between arms
  arm_compa<-pairs(emmeans(model_glm, ~ time | arm), adjust = "BH") %>% as_tibble() %>% 
             mutate(Variable=list_binary_var[i],compa="arm")
  
  # saving results
  emmeans_glm<-as_tibble(emmeans_glm)
  colnames(emmeans_glm)[c(3:7)]<-paste("emmean",colnames(emmeans_glm)[c(3:7)],sep="_")
  
  time_compa <- left_join(time_compa, emmeans_glm, by = c("time"))
  arm_compa  <- left_join(arm_compa, emmeans_glm, by = c("arm"))
  
  Full_emmeans_GLM_table<-bind_rows(Full_emmeans_GLM_table,time_compa,arm_compa)
}

## Exporting results
write.table(Full_GLM_table,paste("GLM_",Sys.Date(),".csv",sep=""),sep=";",col.names=T,row.names=F)
write.table(Full_emmeans_GLM_table,paste("emmeans GLM_",Sys.Date(),".csv",sep=""),sep=";",col.names=T,row.names=F)

#------------------------------------------------------------------
#### Cumulative link mixed effect model (ordinal variables) ####
#------------------------------------------------------------------
Full_CLMM_table<-c()
Full_emmeans_CLMM_table<-c()

for(i in 1:length(list_ordinal_var)){
  
  # selecting variables
  df_tmp<-df %>% select(list_ordinal_var[i],arm,time,ID)
  colnames(df_tmp)[1]<-"tmp"
  
  # creating the model
  model_clmm <- clmm(tmp ~ arm * time + (1 | ID), data = df_tmp)
  
  # full summary table
  clmm_results <- tidy(model_clmm, effects = "fixed", conf.int = TRUE, exponentiate = T) %>%
                  mutate(Variable=list_ordinal_var[i]) %>%
                  filter(coef.type=="location")
  
  # saving results
  Full_CLMM_table<-bind_rows(Full_CLMM_table,clmm_results)
  
  ## Pairwise comparison
  
  # Estimating marginal means (emmeans) - predicted probabilities
  emmeans_clmm <- emmeans(model_clmm,~ arm * time | tmp,mode = "prob")
  
  # Pairwise contrasts between times
  time_compa<-pairs(emmeans(model_clmm,~ arm * time | tmp,mode = "prob"),
                    by = "time", adjust = "BH") %>% as_tibble() %>% 
    mutate(Variable=list_ordinal_var[i],compa="time") %>% select(-df)
  
  # Pairwise contrasts between arms
  arm_compa<-pairs(emmeans(model_clmm,~ arm * time | tmp,mode = "prob"),
                   by = "arm", adjust = "BH") %>% as_tibble() %>% 
    mutate(Variable=list_ordinal_var[i],compa="arm") %>% select(-df)
  
  # saving results
  emmeans_clmm<-as_tibble(emmeans_clmm) %>% select(-df)
  colnames(emmeans_clmm)[c(4:7)]<-paste("emmean",colnames(emmeans_clmm)[c(4:7)],sep="_")
  
  time_compa<-left_join(time_compa,emmeans_clmm,by="time")
  arm_compa<-left_join(arm_compa,emmeans_clmm,by="arm")
  
  Full_emmeans_CLMM_table<-bind_rows(Full_emmeans_CLMM_table,time_compa,arm_compa)
  
}

write.table(Full_CLMM_table, paste("CLMM_", Sys.Date(), ".csv", sep=""), sep=";", col.names=T, row.names=F)
write.table(Full_emmeans_CLMM_table, paste("emmeans CLMM_", Sys.Date(), ".csv", sep=""), sep=";", col.names=T, row.names=F)
