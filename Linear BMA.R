#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing packages if needed
pack_needed<-c("data.table","tidyverse","tibble","broom","rstanarm","BayesFactor","bayestestR","glmnet",
               "BMA","loo","bayesboot","grid","gridExtra","conflicted","usdm","fmsb","RColorBrewer","patchwork","here")

for (i in 1:length(pack_needed)){
  if(pack_needed[i]%in%.packages(all.available=TRUE)){
  }else{
    install.packages(pack_needed[i])
  }
}

## Package loading
library(tidyverse)
library(data.table)
library(broom)
library(rstanarm)
library(BayesFactor)
library(bayestestR)
library(glmnet)
library(BMA)
library(loo)
library(bayesboot)
library(grid)
library(gridExtra)
library(conflicted)
library(usdm)
library(fmsb)
library(RColorBrewer)
library(patchwork)
library(here)

## Preventing package conflicts
conflict_prefer("select", "dplyr") 
conflict_prefer("filter", "dplyr")
conflict_prefer("alpha", "scales")

## Setting the working directory
here::here("Linear BMA")

#----------------------------------------------------------------
#### Creating functions ####

## not including function (opposite function of %in%)
`%ni%`<-Negate('%in%')

#----------------------------------------------------------------
#### Data import ####

## Importing the clean dataset
clean_data<-as_tibble(na.omit(CardioDataSets::cardioRiskFactors_df))

## Renaming the outcome variable
colnames(clean_data)[which(colnames(clean_data)=="Uric")]<-"outcome"
clean_data_save<-clean_data

## Creating a table with the variable type
var_type_tab<-as_tibble(data.frame(feature=colnames(clean_data),
                                    Var_type=c("continuous","continuous","continuous","dichotomous","continuous","continuous","continuous","continuous",
                                               "continuous","continuous","continuous","dichotomous","continuous","continuous")))

#----------------------------------------------------------------
#### Feature selection ####

## Identifying and removing collinear variables
collinear<-usdm::vifstep(as.data.frame(clean_data %>% select(-outcome)),th=2.5)
clean_data<-clean_data %>% select(colnames(clean_data)[which(colnames(clean_data) %ni% c(collinear@excluded))])

## Regularization methods (LASSO: L1-penalized linear regression) to do feature selection and shrinkage simultaneously

# Preparing the data for LASSO
X <- as.matrix(clean_data %>% select(-outcome))
y <- clean_data$outcome

# Fitting LASSO linear regression with cross-validation to select lambda
set.seed(12345)
cvfit <- cv.glmnet(X, y, family = "gaussian", alpha = 1, nfolds = 5)

# Identifying lambda with minimum cross-validated error
best_lambda <- cvfit$lambda.min

# Fitting final model at best lambda
set.seed(12345)
lasso_model <- glmnet(X, y, family = "gaussian", alpha = 1, lambda = best_lambda)

# Extracting coefficients (nonzero are selected variables)
coef_selected <- coef(lasso_model)

# Identifying features to include in the BMA analysis
selected_vars <- rownames(coef_selected)[which(coef_selected[,1] != 0)]
selected_vars <- selected_vars[selected_vars != "(Intercept)"]

#----------------------------------------------------------------
#### BMA analysis ####

## Selecting variables of interest based on feature selection
feature_list<-clean_data_save %>% select(-outcome) %>% select(all_of(selected_vars))
outcome<-clean_data_save %>% select(outcome) %>% pull

## BMA
set.seed(12345)
bayes2<-bic.glm(feature_list,outcome,glm.family="gaussian")
summary(bayes2)

## Table with determinants identified

# Extracting full BMA coefficients
bma_full <- tibble(feature = names(bayes2$probne0),
                   mean_BMA = bayes2$postmean[2:(length(bayes2$probne0)+1)],
                   sd_BMA = bayes2$postsd[2:(length(bayes2$probne0)+1)],
                   post_prob = bayes2$probne0)

# Identifying the top 5 models by posterior probability
top5_models_idx <- order(bayes2$postprob, decreasing = TRUE)[1:5]
top5_models <- bayes2$which[top5_models_idx, ]
top5_probs <- bayes2$postprob[top5_models_idx]

# Converting top5_models to data frame with feature names
top5_models_df <- as.data.frame(top5_models)
colnames(top5_models_df) <- bma_full$feature

# Computing weighted average of coefficients for top 5 models
coef_top5 <- bma_full %>%
  rowwise() %>%
  mutate(mean_top5 = sum(top5_models_df[[feature]] * top5_probs) / sum(top5_probs)) %>%
  ungroup()

# Identifying the single top model by posterior probability
top_model_idx <- which.max(bayes2$postprob)
top_model <- bayes2$which[top_model_idx, ]

# Computing mean for top model
coef_top <- bma_full %>%
  rowwise() %>%
  mutate(mean_top = ifelse(top_model[which(bma_full$feature == feature)], mean_BMA, 0)) %>%
  ungroup()

# Combining into a comparison table
comparison_table <- bma_full %>%
  left_join(coef_top %>% select(feature, mean_top), by = "feature") %>%
  left_join(coef_top5 %>% select(feature, mean_top5), by = "feature") %>%
  select(feature, mean_top, mean_top5, mean_BMA, sd_BMA, post_prob)

comparison_table<-comparison_table %>% arrange(post_prob)

## Exporting the results
write.table(comparison_table,paste("Table_determinants_",Sys.Date(),".csv"),sep=";",col.names=T,row.names=F)

#- - -
## Model R2

# Extract posterior probabilities of models
posterior_probs <- bayes2$postprob

# Extract R2 for each model
r_squared_values <- 1 - (bayes2$deviance / sum((clean_data$outcome - mean(clean_data$outcome))^2))

# Compute Model-Averaged R² = variance explained by models
R2_BMA <- sum(posterior_probs * r_squared_values)
R2_BMA

# Compute Prediction-Based R2 = how well BMA predicts real values
y_pred <- predict(bayes2, type = "response",newdata=clean_data)
R2_pred <- 1 - sum((clean_data$outcome - y_pred)^2) / sum((clean_data$outcome - mean(clean_data$outcome))^2)
R2_pred

#- - -
## Model coefficient
coefficient<-bayes2$postmean
coefficient

#- - -
## Plotting the regression line
ggplot(data= clean_data, aes(x=y_pred, y=outcome)) +
  geom_point() +
  geom_smooth(method="glm",se=T,method.args = list(family = "gaussian"))+
  annotate(geom="text",label = paste("R^2","==",signif(R2_BMA,3),sep=""),parse = TRUE,x=median(clean_data$outcome,na.rm=T),y=max(clean_data$outcome,na.rm=T))+
  theme_bw()+
  scale_x_continuous('Predicted values (BMA)')+
  scale_y_continuous("Outcome")

#- - -
## Plotting BMA results
plot_df <- comparison_table %>%
  pivot_longer(cols = c(mean_top, mean_top5, mean_BMA),
               names_to = "ModelType",
               values_to = "Coefficient")

plot_df<-plot_df %>% mutate(ModelType=case_when(ModelType=="mean_top"~"mean top model",
                                                ModelType=="mean_top5"~"mean top 5 models",
                                                T~"mean BMA"))

plot_df$ModelType<-factor(plot_df$ModelType,levels=c("mean top model","mean BMA","mean top 5 models"))
plot_df$feature<-factor(plot_df$feature,levels=comparison_table$feature)

# Plot 1: Coefficient comparison plot
plot1<-ggplot(plot_df, aes(x = feature, y = Coefficient,
                           color = ModelType, shape = ModelType)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = Coefficient - sd_BMA, ymax = Coefficient + sd_BMA),
                width = 0.2,
                position = position_dodge(width = 0.6),
                color = "black",
                data = subset(plot_df, ModelType == "mean BMA")) +
  geom_hline(yintercept=0,linetype="dashed",color="grey62")+
  coord_flip() +
  labs(x = "",
       y = "Coefficient estimate") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=16,color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        legend.position='top',
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))

# Plot 2: Inclusion probability barplot
plot_df2<-plot_df %>% select(-c(ModelType,Coefficient)) %>% distinct
plot_df2$feature<-factor(plot_df2$feature,levels=comparison_table$feature)

plot2<-ggplot(plot_df2, aes(x = feature, y = post_prob)) +
  geom_col(fill = "steelblue",color="black") +
  geom_text(aes(x = feature, y = post_prob+4,label=post_prob),size=6)+
  scale_y_continuous(expand=c(0,0),limits=c(0,110))+
  coord_flip() +
  labs(x = "", y = "Posterior inclusion probability (%)") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        strip.background = element_rect(fill="#A6DDCE", colour="black", size=1),
        axis.text.x = element_text(size=16,color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title=element_text(size=16,face="bold",color="black"),
        legend.text = element_text(size = 16, face = "bold"),
        axis.line = element_line(color = "black",size = 0.1, linetype = "solid"))  

# Combining both plots
plot3<-plot1+plot2

# Exporting plot 3
ggsave(plot3,
       file=paste("BMA_results_",
                  Sys.Date(),".pdf",sep=""),
       dpi=600,width=60,height=30,units = "cm",limitsize=F)
