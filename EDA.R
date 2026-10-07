#-------------------------------------------------------------------------------
## Reproducible and generalizable exploratory data analysis (EDA) script
#
# covers:
#   - dataset-specific settings
#   - first look at the data (e.g., types, missing values, frequencies, descriptive tables)
#   - Table 1 by group (e.g., exposure status), then by each phenotype
#   - distribution of the key markers (e.g., density plots, histograms, summary statistics, values below the detection limit, extreme values), overall and by group
#   - pairwise plots of the key markers
#   - comparison of included vs. non-included individuals (standardised mean differences)
#   - flow of the study population
#   - correlations adapted to the type of each pair of variables
#   - collinearity (VIF) among exposures
#   - structure of the exposures: FAMD (i.e., PCA for mixed data)
#-------------------------------------------------------------------------------

#----------------------------------------------------------------
#### Configurations ####

## General R environment step-up

rm(list=ls()) # clearing the environment
gctorture(FALSE) # disabling memory torture
options(stringsAsFactors=FALSE) # no automatic conversion of characters into factors

## R packages

required_pkgs<-c("tidyverse","conflicted","gtsummary","flextable","DataExplorer","summarytools","psych","dlookr","labelled","car","FactoMineR","factoextra",
                 "GGally","here","survey")

is_installed<-required_pkgs %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed==FALSE)){
  install.packages(required_pkgs[!is_installed],repos="http://cran.us.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only=TRUE))

## Preventing package conflicts

conflict_prefer("select","dplyr")
conflict_prefer("filter","dplyr")
conflict_prefer("slice","dplyr")

## Setting the working directory

here::here("Exploratory data analysis")

## Setting seed

seed<-123

#----------------------------------------------------------------
#### Creating functions ####

#- - - - - -
## Custom ggplot theme

theme_Gaia<-function(){
  theme_bw() +
    theme(strip.text=element_text(size=14, colour="black", face="bold"),
          strip.background=element_rect(fill="#CAE1FF",colour="black"),
          axis.text=element_text(size=14, color="black"),
          axis.title=element_text(size=16, face="bold", color="black"),
          legend.text=element_text(size=14),
          legend.title=element_text(size=16, face="bold"),
          axis.line=element_line(color="black", linewidth=0.1))
}

#- - - - - -
## not-in operator

`%ni%`<-Negate(`%in%`)

#- - - - - -
## Formatting summary statistics for display

test_format<-function(x){
  x<-as.numeric(x)
  sign_x<-if_else(x<0,"neg","pos")
  x<-abs(x)
  x_raw<-x
  if(is.na(x)|is.infinite(x)) x_raw<-0
  if(x_raw>=100){
    virg_pos<-str_locate(as.character(x_raw),"[.]")[1]
    if(!is.na(virg_pos)&as.numeric(substr(x_raw,virg_pos+1,virg_pos+1))>=5){
      x<-x+1
      x<-as.numeric(substr(x,1,virg_pos-1))
    }
  }
  if(is.na(x)|is.infinite(x)){
    x<-""
  }else{
    if(x<0.05|x>=10000){
      x<-format(signif(x,3), scientific=TRUE)
    }else{
      x_save<-x
      x<-signif(x,3)
      if(nchar(x)==6) x<-as.numeric(substr(x,1,5))
      if(nchar(x)==5){
        if(as.numeric(substr(x,5,5))>=5){
          x<-x+0.01
          x<-substr(x,1,4)
        }else{
          x<-substr(x,1,4)
        }
      }else{
        if(x>=1000) x<-as.character(signif(x_save,4)) else x<-as.character(x)
      }
    }
  }
  if(sign_x=="neg"&x_raw!=0) x<-paste("-",x,sep="")
  if(x=="0e+00") x<-"0"
  return(x)
}

#- - - - - -
## Formatting all numeric columns of a table with test_format before export

format_tab<-function(tab){
  tab %>% mutate(across(where(is.numeric), ~ sapply(.x, test_format)))
}

#- - - - - -
## Exporting a table as a .csv file

export_tab<-function(tab, name){
  write.table(tab, paste(name, "_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)
}

#- - - - - -
## Flagging statistical outliers (> 3 IQR beyond quartiles), nothing is flagged if IQR = 0

flag_outliers<-function(x){
  q<-quantile(x, c(.25, .75), na.rm=TRUE)
  i<-diff(q)
  if(i==0) return(rep(FALSE, length(x)))
  !is.na(x) & (x < q[1] - 3 * i | x > q[2] + 3 * i)
}

#- - - - - -
## Guessing the type of a variable: continuous, binary, ordinal (ordered factor), nominal (factor with > 2 levels)

get_var_type<-function(x){
  if(is.factor(x) && is.ordered(x)) return("ordinal")
  if(is.factor(x) || is.character(x)) return(ifelse(length(unique(na.omit(x)))<=2, "binary", "nominal"))
  if(is.logical(x)) return("binary")
  if(length(unique(na.omit(x)))<=2) return("binary")
  "continuous"
}

#- - - - - -
## Standardised mean difference (SMD) between two groups (g=1 vs. g=0):
#
# - continuous variable: (mean1 - mean0)/pooled SD
# - categorical variable: one SMD per level, from proportions

smd_one<-function(x, g){
  if(is.numeric(x)){
    s<-sqrt((var(x[g==1], na.rm=TRUE) + var(x[g==0], na.rm=TRUE))/2)
    return(tibble(level="", smd=(mean(x[g==1], na.rm=TRUE) - mean(x[g==0], na.rm=TRUE))/s))
  }
  x<-factor(x)
  bind_rows(lapply(levels(x), function(l){
    p1<-mean(x[g==1]==l, na.rm=TRUE)
    p0<-mean(x[g==0]==l, na.rm=TRUE)
    tibble(level=l, smd=(p1 - p0)/sqrt((p1 * (1 - p1) + p0 * (1 - p0))/2))
  }))
}

#- - - - - -
## Summary statistics of a continuous marker

summarise_marker<-function(x, marker, group){
  lod_val<-ifelse(marker %in% names(lod_markers), lod_markers[[marker]], NA)
  tibble(marker=marker, group=group,
         n=sum(!is.na(x)), n_missing=sum(is.na(x)),
         median=median(x, na.rm=TRUE),
         q1=unname(quantile(x, .25, na.rm=TRUE)),
         q3=unname(quantile(x, .75, na.rm=TRUE)),
         min=min(x, na.rm=TRUE), max=max(x, na.rm=TRUE),
         n_below_lod=ifelse(is.na(lod_val), NA, sum(x<=lod_val, na.rm=TRUE)),
         per_below_lod=ifelse(is.na(lod_val), NA, mean(x<=lod_val, na.rm=TRUE)*100),
         n_extreme=sum(flag_outliers(x)))
}

#- - - - - -
## Description of continuous variables (one row per variable): n (% missing), mean +- SD, median [IQR], range

describe_continuous<-function(data, vars){
  data %>%
    select(all_of(vars)) %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    summarise(n=n(), n_NA=mean(is.na(value))*100,
              mean=mean(value, na.rm=TRUE), sd=sd(value, na.rm=TRUE),
              median=median(value, na.rm=TRUE), IQR=IQR(value, na.rm=TRUE),
              min=min(value, na.rm=TRUE), max=max(value, na.rm=TRUE), .groups="drop") %>%
    rowwise() %>%
    mutate(n_NA=paste("(", test_format(n_NA), "%)", sep=""),
           mean=paste(test_format(mean), " +- ", test_format(sd), sep=""),
           median=paste(test_format(median), " [", test_format(IQR), "]", sep=""),
           range=paste("[", test_format(min), "-", test_format(max), "]", sep="")) %>%
    ungroup() %>%
    mutate(n=paste(n, n_NA)) %>%
    select(variable=name, n, mean, median, range)
}

#- - - - - -
## Description of categorical variables (one row per variable and category): n and %

describe_categorical<-function(data, vars){
  data %>%
    select(all_of(vars)) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(everything()) %>%
    mutate(value=if_else(is.na(value), "NA", value)) %>%
    group_by(name, value) %>%
    count() %>%
    group_by(name) %>%
    mutate(per=n/sum(n)*100) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(per=test_format(per)) %>%
    ungroup() %>%
    mutate(n=as.character(n)) %>%
    select(variable=name, category=value, n, per)
}

#- - - - - -
## Test accounting for relatedness (cluster = family), used by add_p() in the tables
#
# - symmetric continuous variable: linear model with cluster-robust variance (svyglm + Wald test)
# - skewed continuous variable: cluster-robust rank test (Wilcoxon for 2 groups, Kruskal-Wallis for more)
# - categorical variable: Rao-Scott chi-squared test (chi-squared corrected for clustering)

test_family<-function(data, variable, by, ...){
  d<-data.frame(y=data[[variable]], g=factor(data[[by]]), fam=data[[family_var]])
  d<-d[complete.cases(d), ]
  des<-survey::svydesign(ids=~fam, data=d, nest=FALSE)
  if(is.numeric(d$y)){
    if(variable %in% skewed_vars){
      p<-survey::svyranktest(y ~ g, design=des, test=ifelse(nlevels(d$g)>2, "KruskalWallis", "wilcoxon"))$p.value
    }else{
      p<-survey::regTermTest(survey::svyglm(y ~ g, design=des), ~g)$p
    }
  }else{
    p<-survey::svychisq(~y + g, design=des, statistic="F")$p.value
  }
  tibble(p.value=as.numeric(p))
}

#- - - - - -
## Table 1 by a grouping variable
#
# - default: n (%) for categorical variables, mean (SD) for symmetric variables, median [Q1; Q3] for skewed ones 
# (sym_vars and skewed_vars are defined below from the data with a numerical skewness test)
#
# - detailed_table=TRUE: several lines per continuous variable (n non-missing, mean (SD), median (Q1; Q3), min; max)
#
# - show_pvalue=TRUE: p-value of a standard test (e.g., Wilcoxon/Kruskal-Wallis, chi-squared/Fisher). These tests ignore the
#   correlation between individuals (e.g., relatives): interpret with caution or turn off

make_table1<-function(data, by_var, vars=table_vars){
  data[[by_var]]<-factor(data[[by_var]])
  keep<-unique(c(vars, by_var, family_var)) # family_var is kept in the data (for the test) but not displayed (include=)
  data<-data %>% select(all_of(keep))
  
  if(detailed_table){
    tab<-data %>%
      tbl_summary(by=all_of(by_var),
                  include=all_of(vars),
                  statistic=list(all_continuous() ~ c("{N_nonmiss} ({p_miss})", "{mean} ({sd})", "{median} ({p25}; {p75})", "{min}; {max}"),
                                 all_categorical() ~ "{n} ({p}%)"),
                  type=all_continuous() ~ "continuous2",
                  missing_text="Missing",
                  missing="no")
  }else{
    stat<-list()
    if(length(sym_vars)>0) stat<-c(stat, list(all_of(sym_vars) ~ "{mean} ({sd})"))
    if(length(skewed_vars)>0) stat<-c(stat, list(all_of(skewed_vars) ~ "{median} [{p25}; {p75}]"))
    stat<-c(stat, list(all_categorical() ~ "{n} ({p}%)"))
    tab<-data %>%
      tbl_summary(by=all_of(by_var),
                  include=all_of(vars),
                  statistic=stat,
                  missing_text="Missing",
                  missing="ifany",
                  digits=list(all_continuous() ~ 2))
  }
  
  if(show_pvalue){
    if(is.null(family_var)){
      tab<-tab %>% add_p()
    }else{
      tab<-tab %>% add_p(test=list(all_continuous() ~ test_family, all_categorical() ~ test_family))
    }
  }
  
  tab %>% add_overall() %>% add_stat_label() %>% bold_labels()
}

#- - - - - -
## Association between two variables, with a coefficient adapted to their types:
#
# - continuous x continuous: Spearman
# - binary x binary: Phi coefficient (signed)
# - binary x continuous: point-biserial (Pearson with 0/1 coding)
# - ordinal x (ordinal, binary or continuous): Kendall's Tau-b
# - nominal (> 2 levels) x continuous: correlation ratio eta (unsigned, 0 to 1)
# - nominal x categorical: Cramer's V (unsigned, 0 to 1)

cor_pair<-function(x, y, tx, ty){
  d<-droplevels(na.omit(data.frame(x=x, y=y)))
  n<-nrow(d)
  if(n<10 || length(unique(d$x))<2 || length(unique(d$y))<2){
    return(tibble(r=NA_real_, p=NA_real_, method="not computable", n=n))
  }
  asnum<-function(v) if(is.factor(v)) as.integer(v) else as.numeric(v)
  types<-c(tx, ty)
  
  if("nominal" %in% types){
    if("continuous" %in% types){
      nom<-if(tx=="nominal") d$x else d$y
      con<-if(tx=="nominal") d$y else d$x
      fit<-summary(aov(asnum(con) ~ factor(nom)))[[1]]
      return(tibble(r=sqrt(fit[1, "Sum Sq"]/sum(fit[, "Sum Sq"])), p=fit[1, "Pr(>F)"], method="Correlation ratio eta (unsigned)", n=n))
    }
    tab<-table(d$x, d$y)
    ch<-suppressWarnings(chisq.test(tab, correct=FALSE))
    return(tibble(r=sqrt(unname(ch$statistic)/(n * (min(dim(tab)) - 1))), p=ch$p.value, method="Cramer's V (unsigned)", n=n))
  }
  
  if(all(types=="continuous")){
    ct<-suppressWarnings(cor.test(asnum(d$x), asnum(d$y), method="spearman", exact=FALSE))
    return(tibble(r=unname(ct$estimate), p=ct$p.value, method="Spearman", n=n))
  }
  
  if("ordinal" %in% types){
    ct<-suppressWarnings(cor.test(asnum(d$x), asnum(d$y), method="kendall", exact=FALSE))
    return(tibble(r=unname(ct$estimate), p=ct$p.value, method="Kendall's Tau-b", n=n))
  }
  
  if(all(types=="binary")){
    tab<-table(d$x, d$y)
    ch<-suppressWarnings(chisq.test(tab, correct=FALSE))
    sgn<-ifelse(all(dim(tab)==2), sign(tab[1, 1] * tab[2, 2] - tab[1, 2] * tab[2, 1]), 1)
    return(tibble(r=sgn * sqrt(unname(ch$statistic)/n), p=ch$p.value, method="Phi coefficient", n=n))
  }
  
  ct<-suppressWarnings(cor.test(asnum(d$x), asnum(d$y), method="pearson"))
  tibble(r=unname(ct$estimate), p=ct$p.value, method="Point-biserial", n=n)
}

#- - - - - -
## All pairwise associations among a set of variables (p-values adjusted by Benjamini-Hochberg over all pairs)

cor_mixed<-function(data, vars, alpha=0.05){
  types<-sapply(data[vars], get_var_type)
  pairs<-combn(vars, 2, simplify=FALSE)
  long<-bind_rows(lapply(pairs, function(p){
    bind_cols(tibble(var1=p[1], var2=p[2], type1=types[[p[1]]], type2=types[[p[2]]]),
              cor_pair(data[[p[1]]], data[[p[2]]], types[[p[1]]], types[[p[2]]]))
  }))
  long$p_adj<-p.adjust(long$p, method="BH")
  
  cm<-matrix(NA_real_, length(vars), length(vars), dimnames=list(vars, vars))
  sm<-matrix(FALSE, length(vars), length(vars), dimnames=list(vars, vars))
  diag(cm)<-1
  for(i in seq_len(nrow(long))){
    cm[long$var1[i], long$var2[i]]<-long$r[i]
    cm[long$var2[i], long$var1[i]]<-long$r[i]
    sm[long$var1[i], long$var2[i]]<-isTRUE(long$p_adj[i]<alpha)
    sm[long$var2[i], long$var1[i]]<-isTRUE(long$p_adj[i]<alpha)
  }
  list(long=long, r=cm, sig=sm, types=tibble(variable=vars, type=unname(types)))
}

#- - - - - -
## Heatmap of a correlation matrix (* = BH-adjusted p-value below alpha, if sig is provided)

plot_cor_matrix<-function(cm, legend_title, sig=NULL){
  cm_long<-as_tibble(as.data.frame(cm), rownames="var1") %>%
    pivot_longer(-var1, names_to="var2", values_to="r")
  if(!is.null(sig)){
    sig_long<-as_tibble(as.data.frame(sig), rownames="var1") %>%
      pivot_longer(-var1, names_to="var2", values_to="sig")
    cm_long<-left_join(cm_long, sig_long, by=c("var1","var2"))
  }else{
    cm_long$sig<-FALSE
  }
  cm_long %>%
    mutate(var1=factor(var1, levels=colnames(cm)), var2=factor(var2, levels=colnames(cm)),
           label=paste(round(r, 2), ifelse(sig, "*", ""), sep="")) %>%
    ggplot(aes(x=var1, y=var2, fill=r)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    geom_tile(color="black") +
    geom_text(aes(label=label), size=12/ggplot2::.pt) +
    scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", limits=c(-1, 1)) +
    labs(x="", y="", fill=legend_title) +
    theme_Gaia() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}

#- - - - - -
## Saving a plot (png)

#- - - - - -
## Saving a plot (png) with a size computed from its content

save_plot<-function(plot, name, width=NULL, height=NULL, panel_w=12, panel_h=8, cm_per_x=1.2, cm_per_y=1, matrix_cell=5,
                    margin_w=4, margin_h=3, width_max=30, height_max=60, dpi=300){
  if(is.null(width) | is.null(height)){
    if(inherits(plot, "ggmatrix")){
      w_auto<-plot$ncol * matrix_cell + margin_w
      h_auto<-plot$nrow * matrix_cell + margin_h
    }else{
      b<-ggplot2::ggplot_build(plot)
      lay<-b$layout$layout
      pp<-b$layout$panel_params[[1]]
      n_x<-tryCatch(length(pp$x$get_labels()), error=function(e) 0)
      n_y<-tryCatch(length(pp$y$get_labels()), error=function(e) 0)
      w_auto<-max(lay$COL) * max(panel_w, n_x * cm_per_x) + margin_w
      h_auto<-max(lay$ROW) * max(panel_h, n_y * cm_per_y) + margin_h
    }
    if(is.null(width)) width<-min(w_auto, width_max)
    if(is.null(height)) height<-min(h_auto, height_max)
  }
  ggsave(plot=plot, filename=paste(name, "_", Sys.Date(), ".png", sep=""), width=width, height=height, units="cm", dpi=dpi, bg="white", limitsize=FALSE)
  message("Figure saved: ", name, "_", Sys.Date(), ".png (", round(width, 1), " x ", round(height, 1), " cm)")
}

#- - - - - -
## Variance inflation factors for continuous and categorical variables
#
# - for a continuous or binary variable GVIF = VIF
# - for a factor with Df > 1, the comparable quantity is GVIF^(1/(2*Df)), squared: it is reported in the column "vif"
#
# - Note: the VIF does not depend on the outcome: a random dummy outcome is used to fit the model

compute_gvif<-function(d, vars){
  fit<-lm(rnorm(nrow(d)) ~ ., data=d[vars])
  v<-car::vif(fit)
  if(is.matrix(v)){
    tibble(variable=rownames(v), gvif=v[, 1], df=v[, 2], vif=v[, 3]^2)
  }else{
    tibble(variable=names(v), gvif=unname(v), df=1, vif=unname(v))
  }
}

#- - - - - -
## VIF stepwise exclusion on complete cases, in three stages:
#
# 1. constant variables (VIF undefined)
# 2. perfectly collinear variables (NA coefficient in lm): removed one at a time
# 3. the variable with the highest VIF (or a non-finite VIF) is removed, then all VIF are recomputed, until all VIF <= thr

select_vif<-function(data, vars, thr){
  removed<-tibble(variable=character(), vif=numeric(), reason=character())
  
  d<-data %>%
    select(all_of(vars)) %>%
    mutate(across(where(is.character), factor)) %>%
    mutate(across(where(is.ordered), ~ factor(.x, ordered=FALSE))) %>%
    na.omit() %>%
    droplevels()
  
  # 1. constant variables
  const<-vars[sapply(d[vars], function(x) length(unique(x))<=1)]
  if(length(const)>0){
    removed<-bind_rows(removed, tibble(variable=const, vif=NA_real_, reason="constant"))
    vars<-setdiff(vars, const)
  }
  
  # 2. perfectly collinear variables
  repeat{
    if(length(vars)<2) break
    fit<-lm(rnorm(nrow(d)) ~ ., data=d[vars])
    na_coef<-is.na(coef(fit))
    if(!any(na_coef)) break
    v<-attr(terms(fit), "term.labels")[fit$assign[na_coef]][1]
    removed<-bind_rows(removed, tibble(variable=v, vif=Inf, reason="perfect collinearity"))
    vars<-setdiff(vars, v)
  }
  
  # 3. stepwise exclusion
  repeat{
    if(length(vars)<3) break
    tab<-compute_gvif(d, vars)
    bad<-which(!is.finite(tab$vif))
    if(length(bad)>0){
      worst<-bad[1]
      reason<-"non-finite VIF"
    }else{
      worst<-which.max(tab$vif)
      if(tab$vif[worst]<=thr) break
      reason<-paste("VIF >", thr)
    }
    removed<-bind_rows(removed, tibble(variable=tab$variable[worst], vif=tab$vif[worst], reason=reason))
    vars<-setdiff(vars, tab$variable[worst])
  }
  
  list(kept=vars, removed=removed, final=if(length(vars)>=2) compute_gvif(d, vars) else tibble())
}

#----------------------------------------------------------------
#### User settings - edit with your own data ####

# Group variable used to split Table 1, if relevant (a factor or a 0/1 variable)
group_var<-"asthma_ever"

# Phenotype variables: one Table 1 is produced for each of them (factors with several levels), if relevant
phenotype_vars<-c("asthma_type","severity","control","persistence")

# Cluster variable (e.g., family identifier) for tests accounting for relatedness, NULL if individuals are independent
family_var<-"family_id"

# Variables shown in Table 1 (covariates and markers)
table_vars<-c("age","sex","center","ses","smoking","bmi","il6","crp","apwv","framingham","ntprobnp","troponin","sst2","cac_class")

# Key markers: continuous and categorical
marker_cont<-c("apwv","framingham","ntprobnp","troponin","sst2")
marker_cat<-c("cac_class")

# Limit of detection (LOD) of markers (in the analysis unit), empty if none
lod_markers<-c(troponin=2)

# Inclusion: indicator (1=included in the analysis, 0=not included) and optional reason for non-inclusion
inclusion_var<-"included"
reason_var<-"reason_non_inclusion" # NULL if not available

# Variables used to compare included vs. non-included individuals (must be available for both groups)
inclusion_compare_vars<-c("age","sex","ses","smoking","asthma_ever","fev1_prev","cv_event_prev","n_visits_prev","prs_cv")

# Variables for the correlation analysis
cor_vars<-c("age","sex","bmi","il6","crp","apwv","framingham","ntprobnp","troponin","sst2","cac_class","smoking","ses")

# Exposures (any type)
expo_vars<-c(sprintf("expo%02d", 1:9),"pets","rural","diet")

# Options
show_pvalue<-TRUE     # p-values in Table 1 (standard tests ignoring the correlation between individuals/subjects/samples)
detailed_table<-FALSE # TRUE: several lines per continuous variable in Table 1

# Thresholds
thr_skew<-1     # |skewness| above which a continuous variable is summarised by median [Q1; Q3] and plotted on a log scale
thr_smd<-0.1    # |SMD| above which a difference is flagged
thr_vif<-5      # VIF above which an exposure is considered collinear with the others
alpha_cor<-0.05 # threshold on BH-adjusted p-values for the stars of the correlation heatmaps
famd_ncp<-5     # number of dimensions kept by the FAMD

# Optional variable labels (name=label): variables absent from the dataset are ignored
var_labels<-list(age="Age (years)", sex="Sex (1=female)", center="Center", ses="Socio-economic status", smoking="Smoking status",
                 bmi="BMI (kg/m2)", il6="IL-6", crp="CRP", apwv="aPWV (m/s)", framingham="Framingham score (%)",
                 ntprobnp="NT-proBNP", troponin="Troponin", sst2="sST2", cac_class="CAC class")

#----------------------------------------------------------------
#### Simulating a dataset (replace with your own data) ####

# data_pop = all individuals of the source population
# markers are only available for included individuals (NA otherwise)

simulate_study<-function(n=2000){
  
  sex<-rbinom(n, 1, .5)
  age<-round(runif(n, 35, 80))
  center<-factor(sample(LETTERS[1:5], n, TRUE))
  family_id<-sample(seq_len(ceiling(n/2.5)), n, TRUE)
  ses<-factor(sample(c("low","mid","high"), n, TRUE, c(.3, .4, .3)), levels=c("low","mid","high"), ordered=TRUE)
  smoking<-factor(sample(c("never","former","current"), n, TRUE, c(.5, .3, .2)), levels=c("never","former","current"))
  bmi<-rnorm(n, 26, 4)
  expo<-as.data.frame(MASS::mvrnorm(n, rep(0, 8), 0.4^abs(outer(1:8, 1:8, "-"))))
  colnames(expo)<-sprintf("expo%02d", 1:8)
  expo$expo09<-expo$expo01 + rnorm(n, 0, .3) # nearly collinear with expo01 (illustrates the VIF)
  pets<-factor(rbinom(n, 1, .4))
  rural<-factor(rbinom(n, 1, .25))
  diet<-factor(sample(c("A","B","C"), n, TRUE)) # nominal exposure
  
  asthma_ever<-rbinom(n, 1, plogis(-1.2 + .5 * expo$expo01 + .3 * (smoking=="current")))
  onset_age<-ifelse(asthma_ever==1, pmax(1, pmin(age - 1, round(rgamma(n, 2, .06)))), NA)
  asthma_type<-factor(ifelse(asthma_ever==0, "none", ifelse(onset_age<16, "childhood", "adult")), levels=c("none","childhood","adult"))
  severity<-factor(ifelse(asthma_ever==0, "none", sample(c("mild","moderate","severe"), n, TRUE, c(.5, .3, .2))), levels=c("none","mild","moderate","severe"))
  control<-factor(ifelse(asthma_ever==0, "none", sample(c("controlled","uncontrolled"), n, TRUE, c(.6, .4))), levels=c("none","controlled","uncontrolled"))
  persistence<-factor(ifelse(asthma_ever==0, "none", sample(c("remitted","persistent"), n, TRUE, c(.4, .6))), levels=c("none","remitted","persistent"))
  
  il6<-exp(.3 + .25 * asthma_ever + .012 * (age - 55) + rnorm(n, 0, .5))
  crp<-exp(.6 + .25 * asthma_ever + .010 * (age - 55) + rnorm(n, 0, .6))
  apwv<-8 + .06 * (age - 55) + .4 * (smoking=="current") + .25 * asthma_ever + rnorm(n, 0, .8)
  framingham<-100 * plogis(-2.6 + .05 * (age - 55) - .6 * sex + .5 * (smoking=="current") + .2 * asthma_ever + rnorm(n, 0, .3))
  ntprobnp<-exp(4.5 + .02 * (age - 55) + .3 * asthma_ever + rnorm(n, 0, .6))
  troponin<-pmax(exp(1 + .01 * (age - 55) + .25 * asthma_ever + rnorm(n, 0, .7)), 2) # values below the LOD (2) recorded as the LOD
  sst2<-exp(3 + .01 * (age - 55) + .15 * asthma_ever + rnorm(n, 0, .4))
  cac_lat<--.5 + .05 * (age - 55) + .3 * asthma_ever + rnorm(n)
  cac_class<-cut(cac_lat, c(-Inf, .2, 1.2, Inf), labels=c("0","1-99",">=100"), ordered_result=TRUE)
  
  fev1_prev<-95 - 6 * asthma_ever + rnorm(n, 0, 12)
  cv_event_prev<-rbinom(n, 1, plogis(-3.5 + .03 * (age - 55)))
  n_visits_prev<-rpois(n, 3 + 2 * asthma_ever)
  prs_cv<-rnorm(n)
  
  included<-rbinom(n, 1, plogis(1.2 - .02 * (age - 55) - .5 * (smoking=="current") - .3 * asthma_ever))
  reason_non_inclusion<-ifelse(included==1, NA, sample(c("deceased","refused","lost to follow-up"), n, TRUE, c(.3, .3, .4)))
  
  d<-tibble(id=1:n, sex, age, center, family_id, ses, smoking, bmi, expo, pets, rural, diet,
            asthma_ever=factor(asthma_ever, levels=0:1, labels=c("No asthma","Ever asthma")),
            asthma_type, severity, control, persistence, il6, crp, apwv, framingham, ntprobnp, troponin, sst2, cac_class,
            fev1_prev, cv_event_prev, n_visits_prev, prs_cv, included, reason_non_inclusion)
  
  # markers and current variables are unknown for non-included individuals
  d<-d %>% mutate(across(c(bmi, il6, crp, apwv, framingham, ntprobnp, troponin, sst2, cac_class), ~ replace(., included==0, NA)))
  d
}

set.seed(seed)
data_pop<-simulate_study()

## With your own data (the clean dataset of the data management step, with variable types already set), replace the simulation by e.g.:
# data_pop<-as_tibble(read.csv("path/to/your_data.csv"))

# Analysis population = included individuals (characters converted into factors)
data_study<-data_pop %>%
  filter(.data[[inclusion_var]]==1) %>%
  mutate(across(where(is.character), factor))
data_study[[group_var]]<-factor(data_study[[group_var]])

# Applying labels (only to existing variables)
labels_ok<-var_labels[names(var_labels) %in% colnames(data_study)]
labelled::var_label(data_study)<-labels_ok

#----------------------------------------------------------------
#### First-look ####

data_study %>% glimpse()
data_study %>% summary()
data_study %>% psych::describe()
dlookr::diagnose(data_study)
plot_intro(data_study, ggtheme=theme_Gaia())

# frequencies of categorical variables and descriptive statistics of numerical variables
data_fix<-data_study %>%
  mutate(across(where(is.ordered), ~ factor(.x, ordered=FALSE))) %>%
  as.data.frame()

if(any(sapply(data_fix, is.factor))) summarytools::freq(data_fix %>% select(where(is.factor)))
data_fix %>% select(where(is.numeric)) %>% psych::describe(quant=c(.25, .75)) %>% as_tibble(rownames="variable")

# missing data
profile_missing(data_study) %>% arrange(desc(pct_missing))

#- - - - 
## Descriptive tables, one row per variable

all_cont<-colnames(data_study)[sapply(data_study, is.numeric)]
all_cat<-colnames(data_study)[sapply(data_study, function(x) is.factor(x) | is.character(x))]

table_conti<-describe_continuous(data_study, all_cont)
table_categ<-describe_categorical(data_study, all_cat)
print(table_conti, n=Inf)
print(table_categ, n=Inf)

#----------------------------------------------------------------
#### Table 1 ####

#- - - - 
## Choosing the summary statistic of each continuous variable from its numerical skewness

cont_vars<-table_vars[sapply(data_study[table_vars], is.numeric)]
cont_vars<-cont_vars[sapply(data_study[cont_vars], function(x) length(unique(na.omit(x)))>10)] # numeric variables with few values are summarized as categories
skew_vals<-sapply(data_study[cont_vars], psych::skew)
skewed_vars<-names(skew_vals)[abs(skew_vals)>thr_skew]
sym_vars<-setdiff(cont_vars, skewed_vars)

skew_tab<-tibble(variable=names(skew_vals), skewness=unname(skew_vals), statistic=ifelse(names(skew_vals) %in% skewed_vars, "median [Q1; Q3]", "mean (SD)"))
print(skew_tab, n=Inf)

#- - - - 
## Table 1 by group

table1_group<-make_table1(data_study, group_var)
table1_group

#- - - - 
## One Table 1 per phenotype

table1_phenotype<-lapply(phenotype_vars, function(v) make_table1(data_study, v))
names(table1_phenotype)<-phenotype_vars
table1_phenotype

#----------------------------------------------------------------
#### Distribution of the markers ####

#- - - - 
## Continuous markers: summary statistics, overall and by group

group_levels<-levels(data_study[[group_var]])

marker_summary<-bind_rows(lapply(marker_cont, function(v){
  bind_rows(summarise_marker(data_study[[v]], v, "All"),
            bind_rows(lapply(group_levels, function(g) summarise_marker(data_study[[v]][data_study[[group_var]] %in% g], v, g))))
}))
print(marker_summary, n=Inf)

#- - - - 
## Continuous markers: density plots and histograms by group
#
# - skewed markers (numerical test, |skewness| > thr_skew, positive values only) are shown on a log scale
# - dashed line = limit of detection (LOD)

marker_skew<-sapply(data_study[marker_cont], psych::skew)
marker_log<-names(marker_skew)[abs(marker_skew)>thr_skew & sapply(data_study[marker_cont], function(x) min(x, na.rm=TRUE)>0)]
tibble(marker=names(marker_skew), skewness=unname(marker_skew), log_scale=names(marker_skew) %in% marker_log)

data_long<-data_study %>%
  select(all_of(c(group_var, marker_cont))) %>%
  pivot_longer(-all_of(group_var), names_to="marker", values_to="value") %>%
  filter(!is.na(value)) %>%
  mutate(value_plot=if_else(marker %in% marker_log, log(value), value),
         marker_label=if_else(marker %in% marker_log, paste("log(", marker, ")", sep=""), marker))

lod_tab<-tibble(marker=names(lod_markers), lod=unname(lod_markers)) %>%
  mutate(lod_plot=if_else(marker %in% marker_log, log(lod), lod),
         marker_label=if_else(marker %in% marker_log, paste("log(", marker, ")", sep=""), marker))

fill_colors<-c("#A6DDCE","#F9CBC2","#F1CB0E","#6495ED")[seq_along(group_levels)]

plot_density<-ggplot(data_long, aes(x=value_plot, fill=.data[[group_var]], color=.data[[group_var]])) +
  geom_density(alpha=.4) +
  geom_vline(data=lod_tab, aes(xintercept=lod_plot), linetype="dashed") +
  facet_wrap(~marker_label, scales="free") +
  scale_fill_manual(group_var, values=fill_colors) +
  scale_color_manual(group_var, values=fill_colors) +
  labs(x="", y="Density") +
  theme_Gaia() +
  theme(legend.position="top")
plot_density

plot_hist<-ggplot(data_long, aes(x=value_plot, fill=.data[[group_var]])) +
  geom_histogram(bins=40, alpha=.6, position="identity", color="black") +
  geom_vline(data=lod_tab, aes(xintercept=lod_plot), linetype="dashed") +
  facet_wrap(~marker_label, scales="free") +
  scale_fill_manual(group_var, values=fill_colors) +
  labs(x="", y="Count") +
  theme_Gaia() +
  theme(legend.position="top")
plot_hist

#- - - - 
## Continuous markers: pairwise plots (same transformation as above)

data_pairs<-data_study %>%
  select(all_of(c(group_var, marker_cont))) %>%
  mutate(across(all_of(marker_log), log)) %>%
  rename_with(~ paste("log(", .x, ")", sep=""), all_of(marker_log)) %>%
  rename(group_pairs=all_of(group_var))

plot_pairs<-ggpairs(data_pairs, columns=2:ncol(data_pairs), mapping=aes(color=group_pairs, fill=group_pairs, alpha=.4)) + theme_bw()
plot_pairs

#- - - - 
## Categorical markers: frequencies, overall and by group

marker_cat_summary<-bind_rows(lapply(marker_cat, function(v){
  data_study %>%
    filter(!is.na(.data[[v]])) %>%
    count(group=.data[[group_var]], level=.data[[v]]) %>%
    group_by(group) %>%
    mutate(per=n/sum(n)*100) %>%
    ungroup() %>%
    mutate(marker=v)
}))
print(marker_cat_summary, n=Inf)

#----------------------------------------------------------------
#### Included vs. non-included ####

# SMDs (included minus non-included): only variables known for both groups can be compared

smd_tab<-bind_rows(lapply(inclusion_compare_vars, function(v){
  bind_cols(tibble(variable=v), smd_one(data_pop[[v]], data_pop[[inclusion_var]]))
})) %>%
  mutate(flag=abs(smd)>thr_smd)
print(smd_tab, n=Inf)

# Descriptive table
table_included<-data_pop %>%
  mutate(inclusion=factor(.data[[inclusion_var]], levels=c(0, 1), labels=c("Not included","Included"))) %>%
  select(all_of(unique(c(inclusion_compare_vars, "inclusion", family_var)))) %>%
  tbl_summary(by=inclusion, include=all_of(inclusion_compare_vars), missing_text="Missing", missing="ifany")
if(show_pvalue){
  if(is.null(family_var)){
    table_included<-table_included %>% add_p()
  }else{
    table_included<-table_included %>% add_p(test=list(all_continuous() ~ test_family, all_categorical() ~ test_family))
  }
}
table_included<-table_included %>% add_overall() %>% bold_labels()
table_included

# Flow of the study population
flow_tab<-tibble(step=c("Source population","Included in the analysis","Not included"),
                 n=c(nrow(data_pop), sum(data_pop[[inclusion_var]]==1), sum(data_pop[[inclusion_var]]==0)))
if(!is.null(reason_var)){
  reason_tab<-data_pop %>%
    filter(.data[[inclusion_var]]==0) %>%
    count(step=paste("Not included:", .data[[reason_var]]), name="n")
  flow_tab<-bind_rows(flow_tab, reason_tab)
}
print(flow_tab)

#----------------------------------------------------------------
#### Correlations ####

#- - - - 
## Markers and covariates (coefficient adapted to the types of each pair; see cor_pair)

cor_res<-cor_mixed(data_study, cor_vars, alpha=alpha_cor)
print(cor_res$types)
plot_cor_markers<-plot_cor_matrix(cor_res$r, "Association", sig=cor_res$sig)
plot_cor_markers

# strongest associations (no threshold: ranking only)
cor_top<-cor_res$long %>%
  mutate(abs_r=abs(r)) %>%
  arrange(desc(abs_r)) %>%
  select(-abs_r) %>%
  slice_head(n=10)
print(cor_top)

#- - - - 
## Exposures

cor_expo_res<-cor_mixed(data_study, expo_vars, alpha=alpha_cor)
plot_cor_expo<-plot_cor_matrix(cor_expo_res$r, "Association", sig=cor_expo_res$sig)
plot_cor_expo

cor_expo_top<-cor_expo_res$long %>%
  mutate(abs_r=abs(r)) %>%
  arrange(desc(abs_r)) %>%
  select(-abs_r) %>%
  slice_head(n=10)
cor_expo_top

#----------------------------------------------------------------
#### Collinearity of the exposures (VIF) ####

# VIF = 1/(1-R2) of each exposure regressed on all the others: it detects collinearity involving several variables,
# which pairwise correlations cannot. Computed on complete cases (n printed below)

set.seed(seed)
cat("Complete cases used for the VIF:", sum(complete.cases(data_study[expo_vars])), "\n")

vif_sel<-select_vif(data_study, expo_vars, thr=thr_vif)
vif_removed<-vif_sel$removed # exposures removed, in order of removal
vif_kept<-vif_sel$kept # exposures kept (to be used in the following steps)
vif_final<-vif_sel$final # VIF after exclusion
vif_removed
print(vif_final %>% arrange(desc(vif)), n=Inf)

# initial VIF (before any exclusion), when no variable is constant or perfectly collinear
vif_initial<-tryCatch(compute_gvif(data_study %>% select(all_of(expo_vars)) %>% mutate(across(where(is.character), factor)) %>%
                                     mutate(across(where(is.ordered), ~ factor(.x, ordered=FALSE))) %>% na.omit(), expo_vars),
                      error=function(e) tibble(note="initial VIF not computable (constant or perfectly collinear variable): see vif_removed"))
vif_initial

#----------------------------------------------------------------
#### Structure of the exposures: FAMD (PCA for mixed data) ####

# Unsupervised: no outcome is used. FAMD = PCA if all variables are continuous. Complete cases only (for many missing values, impute first, e.g., with missMDA::imputeFAMD)

complete_id<-complete.cases(data_study[expo_vars])

data_famd<-data_study[complete_id, ] %>%
  select(all_of(expo_vars)) %>%
  mutate(across(where(is.character), factor)) %>%
  mutate(across(where(is.ordered), ~ factor(.x, ordered=FALSE))) %>%
  as.data.frame()

group_famd<-data_study[[group_var]][complete_id]

famd<-FactoMineR::FAMD(data_famd, ncp=famd_ncp, graph=FALSE)

# eigenvalues
famd_eig<-tibble(dimension=seq_len(nrow(famd$eig)), eigenvalue=famd$eig[, 1], per_variance=famd$eig[, 2], cum_variance=famd$eig[, 3])
print(famd_eig, n=Inf)

# contribution of each exposure to the first two dimensions
famd_contrib<-tibble(variable=rownames(famd$var$contrib), dim1=famd$var$contrib[, 1], dim2=famd$var$contrib[, 2])

plot_contrib<-famd_contrib %>%
  pivot_longer(-variable, names_to="dimension", values_to="contribution") %>%
  ggplot(aes(x=reorder(variable, contribution), y=contribution, fill=dimension)) +
  geom_col(position="dodge", color="black") +
  coord_flip() +
  scale_fill_manual("", values=c("#A6DDCE","#F9CBC2")) +
  labs(x="", y="Contribution (%)") +
  theme_Gaia() +
  theme(legend.position="top")
plot_contrib

# scree plot
plot_scree<-factoextra::fviz_screeplot(famd, addlabels=TRUE, ncp=10, barfill="#A6DDCE", barcolor="black") + theme_Gaia()
plot_scree

# all variables (quantitative and categorical) on dimensions 1-2, coloured by cos2 (quality of representation)
plot_var<-factoextra::fviz_famd_var(famd, "var", repel=TRUE, col.var="cos2", gradient.cols=c("#2166ac","#F1CB0E","#b2182b")) + theme_Gaia()
plot_var

# correlation circle (quantitative variables only)
plot_circle<-factoextra::fviz_famd_var(famd, "quanti.var", repel=TRUE, col.var="cos2", gradient.cols=c("#2166ac","#F1CB0E","#b2182b")) + theme_Gaia()
plot_circle

# categories of the categorical variables
plot_quali<-factoextra::fviz_famd_var(famd, "quali.var", repel=TRUE, col.var="cos2", gradient.cols=c("#2166ac","#F1CB0E","#b2182b")) + theme_Gaia()
plot_quali

# cos2 of the variables on dimensions 1-2, and contributions to dimensions 1 and 2
plot_cos2<-factoextra::fviz_cos2(famd, choice="var", axes=1:2, fill="#A6DDCE", color="black") + theme_Gaia()
plot_contrib1<-factoextra::fviz_contrib(famd, choice="var", axes=1, fill="#A6DDCE", color="black") + theme_Gaia()
plot_contrib2<-factoextra::fviz_contrib(famd, choice="var", axes=2, fill="#F9CBC2", color="black") + theme_Gaia()

# individuals coloured by group, with confidence ellipses (descriptive only)
plot_ind<-factoextra::fviz_famd_ind(famd, geom="point", habillage=group_famd, addEllipses=TRUE, alpha.ind=.4) + theme_Gaia()
plot_ind

# biplot: individuals + quantitative variables as arrows (not provided by factoextra for FAMD: built by hand)
famd_ind<-tibble(dim1=famd$ind$coord[, 1], dim2=famd$ind$coord[, 2], group=group_famd)
quanti_coord<-tibble(variable=rownames(famd$quanti.var$coord), dim1=famd$quanti.var$coord[, 1], dim2=famd$quanti.var$coord[, 2])
scale_arrow<-0.8 * min(max(abs(famd_ind$dim1)), max(abs(famd_ind$dim2)))

plot_biplot<-ggplot(famd_ind, aes(x=dim1, y=dim2, color=group)) +
  geom_point(alpha=.3) +
  geom_segment(data=quanti_coord, aes(x=0, y=0, xend=dim1 * scale_arrow, yend=dim2 * scale_arrow), arrow=arrow(length=unit(.2, "cm")), color="black", inherit.aes=FALSE) +
  geom_text(data=quanti_coord, aes(x=dim1 * scale_arrow * 1.1, y=dim2 * scale_arrow * 1.1, label=variable), color="black", inherit.aes=FALSE) +
  scale_color_manual(group_var, values=fill_colors) +
  labs(x=paste("Dim 1 (", round(famd_eig$per_variance[1], 1), "%)", sep=""), y=paste("Dim 2 (", round(famd_eig$per_variance[2], 1), "%)", sep="")) +
  theme_Gaia() +
  theme(legend.position="top")
plot_biplot

#----------------------------------------------------------------
#### Exporting and saving results ####

# Tables 1 (Word)
table1_group %>% as_flex_table() %>% save_as_docx(path=paste("Table 1_", group_var, "_", Sys.Date(), ".docx", sep=""))
for(v in phenotype_vars){
  table1_phenotype[[v]] %>% as_flex_table() %>% save_as_docx(path=paste("Table 1_", v, "_", Sys.Date(), ".docx", sep=""))
}
table_included %>% as_flex_table() %>% save_as_docx(path=paste("Table_included vs. non-included_", Sys.Date(), ".docx", sep=""))

# Descriptive tables
export_tab(table_conti, "Description_continuous variables")
export_tab(table_categ, "Description_categorical variables")

# Statistic chosen for each continuous variable
export_tab(format_tab(skew_tab), "Table 1_statistic choice")

# Markers
export_tab(format_tab(marker_summary), "Marker summary")
export_tab(format_tab(marker_cat_summary), "Categorical marker summary")

# Included vs. non-included and flow
export_tab(format_tab(smd_tab), "Included vs. non-included_SMD")
export_tab(format_tab(flow_tab), "Flow")

# Associations (long format: coefficient, method, p-value, adjusted p-value)
export_tab(format_tab(cor_res$long), "Associations_markers and covariates")
export_tab(format_tab(cor_expo_res$long), "Associations_exposures")

# Collinearity and FAMD
export_tab(format_tab(vif_initial), "VIF_initial")
export_tab(format_tab(vif_final), "VIF_after exclusion")
export_tab(format_tab(vif_removed), "VIF_excluded exposures")
export_tab(format_tab(famd_eig), "FAMD_eigenvalues")
export_tab(format_tab(tibble(variable=rownames(famd$var$contrib), dim1=famd$var$contrib[, 1], dim2=famd$var$contrib[, 2])), "FAMD_contributions")

# Figures
save_plot(plot_density, "Marker density", dpi=300)
save_plot(plot_hist, "Marker histograms", dpi=300)
save_plot(plot_pairs, "Marker pairs", dpi=300)
save_plot(plot_cor_markers, "Associations_markers and covariates", cm_per_x=1.6, cm_per_y=1.4)
save_plot(plot_cor_expo, "Associations_exposures", cm_per_x=1.6, cm_per_y=1.4)
save_plot(plot_scree, "FAMD_scree plot", dpi=300)
save_plot(plot_var, "FAMD_variables", dpi=300)
save_plot(plot_circle, "FAMD_correlation circle", dpi=300)
save_plot(plot_quali, "FAMD_categories", dpi=300)
save_plot(plot_cos2, "FAMD_cos2", dpi=300)
save_plot(plot_contrib1, "FAMD_contributions dim 1", dpi=300)
save_plot(plot_contrib2, "FAMD_contributions dim 2", dpi=300)
save_plot(plot_contrib, "FAMD_contributions", dpi=300)
save_plot(plot_ind, "FAMD_individuals", dpi=300)
save_plot(plot_biplot, "FAMD_biplot", dpi=300)

# Session information (reproducibility)
writeLines(capture.output(sessionInfo()), paste("Session info_", Sys.Date(), ".txt", sep=""))
