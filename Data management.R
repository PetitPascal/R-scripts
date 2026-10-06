#-------------------------------------------------------------------------------
## Reproducible and generalizable data management script
#
# covers:
#   - dataset-specific settings
#   - cleaning: duplicates, harmonization, units, plausibility ranges, outlier flagging
#   - variable dictionary and automatic type recoding
#   - derived variables
#   - missing data description, exclusion rules, structural missing values
#   - single imputation and multiple imputation
#   - numerical skewness test and log-transformation of skewed variables only
#   - one-hot encoding
#   - export of clean datasets and logs
#-------------------------------------------------------------------------------

#----------------------------------------------------------------
#### Configurations ####

## General R environment step-up

rm(list=ls()) # clearing the environment
gctorture(FALSE) # disabling memory torture
options(stringsAsFactors=FALSE) # no automatic conversion of characters into factors

## R packages

required_pkgs<-c("tidyr","ggplot2","DataExplorer","missRanger","mice","labelled","dplyr","tidyverse","conflicted","psych","dlookr","here","fastDummies","qs2","nanoparquet")

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

here::here("Data management")

## Setting seed

seed<-123

#----------------------------------------------------------------
#### Creating functions ####

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
## Custom ggplot theme

theme_Gaia<-function(){
  theme_bw() +
    theme(strip.text=element_text(size=12, colour="black", face="bold"),
          strip.background=element_rect(fill="#CAE1FF",colour="black"),
          axis.text=element_text(size=12, color="black"),
          axis.title=element_text(size=12, face="bold", color="black"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12, face="bold"),
          axis.line=element_line(color="black", linewidth=0.1))
}

#- - - - - -
## not-in operator

`%ni%`<-Negate(`%in%`)

#- - - - - -
## Converting factor/character to numeric safely (e.g., factor "0"/"1" -> 0/1)

num<-function(x) as.numeric(as.character(x))

#- - - - - -
## Flagging statistical outliers (> 3 IQR beyond quartiles), nothing is flagged if IQR = 0

flag_outliers<-function(x){
  q<-quantile(x, c(.25, .75), na.rm=TRUE)
  i<-diff(q)
  if(i==0) return(rep(FALSE, length(x)))
  !is.na(x) & (x < q[1] - 3 * i | x > q[2] + 3 * i)
}

#- - - - - -
## Adding a line to the data management log

add_log<-function(log, step, n, total){
  bind_rows(log, tibble(step=step, n=n, per=n/total*100))
}

#- - - - - -
## Guessing the nature of a variable

guess_var_type<-function(x){
  if(inherits(x, c("Date","POSIXct"))) return("date")
  if(is.logical(x)) return("binary")
  if(is.factor(x) && is.ordered(x)) return("ordinal")
  if(is.character(x) || is.factor(x)){
    return(ifelse(length(unique(na.omit(x)))<=2, "binary", "nominal"))
  }
  if(length(unique(na.omit(x)))<=2) return("binary")
  "continuous"
}

#- - - - - -
## Building a variable dictionary (variable, type): id / date / continuous / binary / nominal / ordinal

build_var_dict<-function(data, id_vars, type_override=c()){
  dict<-tibble(variable=colnames(data), type=sapply(data, guess_var_type))
  dict$type[dict$variable %in% id_vars]<-"id"
  if(length(type_override)>0){
    ov<-type_override[names(type_override) %in% dict$variable]
    dict$type[match(names(ov), dict$variable)]<-unname(ov)
  }
  dict
}

#- - - - - -
## Applying the types of the dictionary with update_columns()

apply_var_types<-function(data, dict){
  vars<-function(t) dict$variable[dict$type==t & dict$variable %in% colnames(data)]
  if(length(vars("id"))>0) data<-update_columns(data=data, ind=vars("id"), what=as.character)
  if(length(vars("date"))>0) data<-update_columns(data=data, ind=vars("date"), what=as.Date)
  if(length(vars("continuous"))>0) data<-update_columns(data=data, ind=vars("continuous"), what=num)
  if(length(vars("binary"))>0) data<-update_columns(data=data, ind=vars("binary"), what=function(x) factor(x, levels=sort(unique(x))))
  if(length(vars("nominal"))>0) data<-update_columns(data=data, ind=vars("nominal"), what=factor)
  if(length(vars("ordinal"))>0) data<-update_columns(data=data, ind=vars("ordinal"), what=function(x) factor(x, ordered=TRUE, levels=sort(unique(x))))
  data
}

#- - - - - -
## One-hot encoding: nominal -> dummy columns, binary factor -> 0/1, ordinal -> integer codes

one_hot_encode<-function(data, dict, remove_first_dummy=FALSE){
  vars<-function(t) dict$variable[dict$type==t & dict$variable %in% colnames(data)]
  for(v in vars("binary")) if(is.factor(data[[v]])) data[[v]]<-as.integer(data[[v]])-1L
  for(v in vars("ordinal")) data[[v]]<-as.integer(data[[v]])
  if(length(vars("nominal"))>0){
    data<-fastDummies::dummy_cols(data, select_columns=vars("nominal"), remove_first_dummy=remove_first_dummy,
                                  remove_selected_columns=TRUE, ignore_na=TRUE)
  }
  as_tibble(data)
}

#- - - - - -
## Risk score

compute_framingham<-function(age, female, sbp, chol, smoker, diabetes){
  100 * plogis(-9 + .07 * age - .5 * female + .018 * sbp + .15 * chol + .6 * smoker + .7 * diabetes)
}

#- - - - - -
## Variables that depend on (possibly imputed) measurements: recomputed in every imputed dataset

derive_after_imputation<-function(data){
  if("cac" %in% colnames(data)){
    data$cac_class<-cut(num(data$cac), c(-Inf, 0, 99, Inf), labels=c("0","1-99",">=100"), ordered_result=TRUE)
  }
  for(v in names(lod)){
    if(v %in% colnames(data)) data[[paste0(v, "_detectable")]]<-as.integer(data[[v]] > lod[[v]])
  }
  if(all(c("age_years","sex","sbp","chol_mmol","smoking","diabetes") %in% colnames(data))){
    data$framingham<-compute_framingham(age=data$age_years, female=num(data$sex), sbp=data$sbp, chol=data$chol_mmol,
                                        smoker=as.integer(data$smoking=="current"), diabetes=num(data$diabetes))
    data$score2<-data$framingham * .8
  }
  snps<-intersect(names(prs_weights), colnames(data))
  if(length(snps)>0) data$prs_cv<-as.numeric(scale(as.matrix(data[snps]) %*% prs_weights[snps]))
  data
}

#- - - - - -
## Log-transforming only the variables flagged "log" by the skewness test

transform_skewed<-function(data, skew_tab){
  todo<-skew_tab %>% filter(action=="log")
  for(i in seq_len(nrow(todo))){
    v<-todo$variable[i]
    if(v %in% colnames(data)) data[[paste0("log_", v)]]<-log(data[[v]] + todo$offset[i])
  }
  data
}

#----------------------------------------------------------------
#### User settings - edit with your own data ####

## Identifiers
id_var<-"id" # individual/subject/sample identifier
cluster_var<-"family_id" # cluster identifier (e.g., family, center), NULL if none

## Dates: all columns whose name matches this pattern are treated as dates
date_pattern<-"^date_"

## Variables whose NA are structural (e.g., age at onset in people without a disease) are excluded from the missingness rules and from imputation
structural_vars<-c("onset_age")

## Variables for which the automatic type guess must be overridden
type_override<-c(cac="continuous", education="ordinal")

## Limit of detection (LOD) of biomarkers, in the unit used after harmonisation
lod<-c(troponin=2)

## Plausibility ranges: values outside are set to NA
ranges<-list(bmi=c(12, 70),
             sbp=c(70, 260),
             chol_mmol=c(2, 12),
             apwv=c(3, 20),
             age_years=c(10, 100))

## Exclusion thresholds
thr_var<-30   # excluding variables with > thr_var % missing values
thr_indiv<-50 # excluding individuals with > thr_indiv % missing values
# thr_rare<-10  # excluding binary variables with fewer than thr_rare individuals in the smaller category

## Transformation thresholds
thr_skew<-1 # |skewness| above which a variable is considered skewed
thr_zero<-25 # variables with more than thr_zero % of zeros are not log-transformed (zero-inflated)
no_transform_pattern<-"^snp|^prs|^age|^onset" # variables never transformed

## Imputation
m_imp<-5      # number of imputed datasets (mice)
maxit_imp<-10 # number of iterations (mice)

## Polygenic risk score weights (can be replace with published weights)
set.seed(seed)
prs_weights<-setNames(round(runif(10, -.3, .3), 3), sprintf("snp%02d", 1:10))

#----------------------------------------------------------------
#### Simulating a messy dataset (replace with your own data) ####

simulate_raw<-function(n=800){
  
  id<-sprintf("P%04d", 1:n)
  family_id<-sample(seq_len(ceiling(n/2.5)), n, TRUE)
  sex<-sample(c("F", "M", "female", "male", NA), n, TRUE, c(.38, .38, .09, .09, .06))
  date_birth<-as.Date("1945-01-01") + sample(0:(365 * 35), n, TRUE)
  date_1<-date_birth + round(365.25 * runif(n, 18, 50))
  date_4<-date_1 + round(365.25 * 30)
  asthma<-rbinom(n, 1, .3)
  date_asthma_onset<-as.Date(ifelse(asthma==1, date_birth + round(365.25 * runif(n, 1, 60)), NA), origin="1970-01-01")
  bad<-sample(which(asthma==1), 5)
  date_asthma_onset[bad]<-date_birth[bad] - 100
  date_ics_start<-as.Date(ifelse(asthma==1 & runif(n) < .7, date_asthma_onset + 365 * 2, NA), origin="1970-01-01")
  ics_daily_dose_ug<-ifelse(is.na(date_ics_start), NA, sample(c(200, 400, 800), n, TRUE))
  smoking<-sample(c("never", "Never", "former", "ex", "current", "1", NA), n, TRUE, c(.3, .1, .2, .1, .2, .05, .05))
  age_years<-as.numeric(date_4 - date_birth)/365.25
  age_years[sample(n, 3)]<-c(-5, 250, 999)
  bmi<-rnorm(n, 26, 4)
  bmi[sample(n, 3)]<-c(250, 3, 90)
  sbp<-rnorm(n, 125, 15)
  chol_mmol<-rnorm(n, 5.2, 1)
  diabetes<-rbinom(n, 1, .08)
  il6<-exp(rnorm(n, .5 + .3 * asthma, .6))
  crp<-exp(rnorm(n, .8 + .3 * asthma, .8))
  ntprobnp<-exp(rnorm(n, 4.5 + .3 * asthma, .7))
  troponin<-pmax(exp(rnorm(n, 1 + .2 * asthma, .7)), 2) # values below the LOD (2 ng/L) are recorded as the LOD
  troponin_unit<-sample(c("ng/L", "ug/L"), n, TRUE, c(.9, .1))
  troponin[troponin_unit=="ug/L"]<-troponin[troponin_unit=="ug/L"]/1000
  apwv<-rnorm(n, 8 + .05 * (age_years - 55) + .3 * asthma, 1)
  cac<-ifelse(runif(n) < .5, 0, round(rexp(n, 1/100)))
  expo<-as.data.frame(matrix(rnorm(n * 6), n, 6, dimnames=list(NULL, sprintf("expo%02d", 1:6))))
  expo$expo_const<-1
  expo$expo_sparse<-ifelse(runif(n) < .6, NA, rnorm(n))
  expo$expo_rare<-ifelse(runif(n) < .01, 1, 0)
  G<-matrix(rbinom(n * 10, 2, .3), n, 10, dimnames=list(NULL, sprintf("snp%02d", 1:10)))
  d<-data.frame(id, family_id, sex, date_birth, date_1, date_4, date_asthma_onset,
                date_ics_start, ics_daily_dose_ug, smoking, age_years, bmi, sbp, chol_mmol, diabetes,
                il6, crp, ntprobnp, troponin, troponin_unit, apwv, cac, expo, G)
  
  for(v in c("bmi", "sbp", "il6", "apwv", "expo01", "expo02", "chol_mmol"))
    d[sample(n, round(.05 * n)), v]<-NA
  d[sample(n, 1), -(1:2)]<-NA
  rbind(d, d[sample(n, 5), ])
}

set.seed(seed)
data_raw<-as_tibble(simulate_raw())

## With your own data, replace the simulation by e.g.:
# data_raw<-as_tibble(read.csv("path/to/your_data.csv"))

#----------------------------------------------------------------
#### First-look ####

data_raw %>% glimpse()
data_raw %>% summary()
data_raw %>% str()
data_raw %>% psych::describe()
data_raw %>% diagnose
data_raw %>% dlookr::describe()
summary(overview(data_raw))
plot_intro(data_raw, ggtheme=theme_Gaia())

#----------------------------------------------------------------
#### Cleaning ####

#- - - - 
## Removing duplicated rows

nb_row<-nrow(data_raw)
data_raw<-data_raw %>% distinct
n_duplicates<-nb_row-nrow(data_raw)

log_data<-tibble(step="Duplicated rows removed", n=n_duplicates, per=n_duplicates/nb_row*100)

# Same identifier with different content -> real conflict that must be solved manually
n_conflict<-sum(duplicated(data_raw[[id_var]]))
if(n_conflict>0) warning(n_conflict, " identifier(s) appear several times with different values: check manually")

#- - - -
## Harmonizing codings and units

# sex (1=female, 0=male)
data_raw$sex<-dplyr::case_when(data_raw$sex %in% c("F", "female") ~ 1L, data_raw$sex %in% c("M", "male") ~ 0L, TRUE ~ NA_integer_)

# smoking
data_raw$smoking<-dplyr::case_when(tolower(data_raw$smoking)=="never" ~ "never",
                                   tolower(data_raw$smoking) %in% c("former", "ex") ~ "former",
                                   tolower(data_raw$smoking) %in% c("current", "1") ~ "current", TRUE ~ NA_character_)

# converting troponin values from ug/L to ng/L, then dropping the unit column
data_raw$troponin<-ifelse(data_raw$troponin_unit=="ug/L", data_raw$troponin * 1000, data_raw$troponin)
data_raw$troponin_unit<-NULL

# recomputing age from dates
data_raw$age_years<-as.numeric(data_raw$date_4 - data_raw$date_birth)/365.25

# impossible/inconsistent dates -> setting to NA
bad_onset<-!is.na(data_raw$date_asthma_onset) & data_raw$date_asthma_onset < data_raw$date_birth
data_raw$date_asthma_onset[bad_onset]<-NA

log_data<-add_log(log_data, "Asthma onset before birth (set to NA)", sum(bad_onset), nrow(data_raw))

# plausibility ranges -> setting out-of-range values to NA
for(v in names(ranges)){
  o<-!is.na(data_raw[[v]]) & (data_raw[[v]] < ranges[[v]][1] | data_raw[[v]] > ranges[[v]][2])
  data_raw[[v]][o]<-NA
  log_data<-add_log(log_data, paste("Out-of-range ", v, " set to NA", sep=""), sum(o), nrow(data_raw))
}

# structural zeros: no treatment start date = not treated -> dose is 0, not missing
data_raw$ics_daily_dose_ug[is.na(data_raw$date_ics_start)]<-0

#- - - -
## Recoding variable types: variable dictionary + update_columns

date_cols<-grep(date_pattern, colnames(data_raw), value=TRUE)
data_raw<-update_columns(data=data_raw, ind=date_cols, what=as.Date) # dates first, so that the type guess is correct

var_dict<-build_var_dict(data_raw, id_vars=c(id_var, cluster_var), type_override=type_override)
print(var_dict, n=Inf) # Check this table: correct it with `type_override` if a type is wrong
data_raw<-apply_var_types(data_raw, var_dict)

#- - - -
## Derived variables that only depend on dates

# asthma phenotypes
data_raw$asthma_ever<-as.integer(!is.na(data_raw$date_asthma_onset))
data_raw$onset_age<-as.numeric(data_raw$date_asthma_onset - data_raw$date_birth)/365.25
data_raw$asthma_type<-factor(ifelse(data_raw$asthma_ever==0, "none", ifelse(data_raw$onset_age < 16, "childhood", "adult")),
                             levels=c("none", "childhood", "adult"))
data_raw$asthma_incident<-as.integer(data_raw$asthma_ever==1 & data_raw$date_asthma_onset > data_raw$date_1)

# cumulative inhaled corticosteroid exposure (duration x dose)
data_raw$ics_years<-ifelse(is.na(data_raw$date_ics_start), 0, pmax(0, as.numeric(data_raw$date_4 - data_raw$date_ics_start)/365.25))
data_raw$ics_cum_dose_mg<-data_raw$ics_years * 365.25 * data_raw$ics_daily_dose_ug/1000

#----------------------------------------------------------------
#### Outliers (flagged, not removed) ####

analysis_cols<-setdiff(colnames(data_raw), c(id_var, cluster_var, date_cols, structural_vars))
cont_vars<-var_dict$variable[var_dict$type=="continuous" & var_dict$variable %in% analysis_cols]

outliers_tab<-tibble(variable=cont_vars, n_outliers=sapply(data_raw[cont_vars], function(x) sum(flag_outliers(x)))) %>%
  filter(n_outliers>0) %>%
  arrange(desc(n_outliers)) %>%
  mutate(per_outliers=n_outliers/nrow(data_raw)*100)

# data_raw %>% select(find_outliers(.)) %>% diagnose

#----------------------------------------------------------------
#### Missing data: description ####

# by variable (dates, identifiers and structural variables are not concerned)
miss_var<-data_raw %>%
  select(all_of(analysis_cols)) %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(), names_to="variable", values_to="per_missing") %>%
  arrange(desc(per_missing))

plot_na_pareto(data_raw %>% select(all_of(analysis_cols)), only_na=TRUE, grade=list(High=0.1, Middle=0.3, Low=1))

# by individual
miss_indiv<-tibble(id=data_raw[[id_var]], per_missing=rowMeans(is.na(data_raw[analysis_cols])) * 100)

#----------------------------------------------------------------
#### Variable and individual exclusion ####

# variables: constant, > thr_var % missing, binary with < thr_rare individuals in the smaller category
excl_var_const<-analysis_cols[sapply(data_raw[analysis_cols], function(x) length(unique(na.omit(x)))<=1)]
excl_var_missing<-miss_var$variable[miss_var$per_missing > thr_var]
binary_vars<-var_dict$variable[var_dict$type=="binary" & var_dict$variable %in% analysis_cols]
# excl_var_rare<-binary_vars[sapply(data_raw[binary_vars], function(x) length(unique(na.omit(x)))==2 && min(table(x))<thr_rare)]
excl_var<-unique(c(excl_var_const, 
                   # excl_var_rare,
                   excl_var_missing))

log_data<-add_log(log_data, "Constant variables excluded", length(excl_var_const), ncol(data_raw))
log_data<-add_log(log_data, paste("Variables excluded (> ", thr_var, "% missing)", sep=""), length(excl_var_missing), ncol(data_raw))
# log_data<-add_log(log_data, paste("Binary variables excluded (< ", thr_rare, " in the smaller category)", sep=""), length(excl_var_rare), ncol(data_raw))

# individuals: > thr_indiv % missing
analysis_cols<-setdiff(analysis_cols, excl_var)
miss_indiv<-tibble(id=data_raw[[id_var]], per_missing=rowMeans(is.na(data_raw[analysis_cols])) * 100)
excl_indiv_miss<-miss_indiv$id[miss_indiv$per_missing > thr_indiv]
log_data<-add_log(log_data, paste("Individuals excluded (> ", thr_indiv, "% missing)", sep=""), length(excl_indiv_miss), nrow(data_raw))

data_clean<-data_raw %>%
  select(-any_of(excl_var)) %>%
  filter(.data[[id_var]] %ni% excl_indiv_miss)

#----------------------------------------------------------------
#### Single imputation: random forests + predictive mean matching ####

# Dates, identifiers and structural variables are kept aside and put back after imputation
static_vars<-c(id_var, cluster_var, intersect(date_cols, colnames(data_clean)), intersect(structural_vars, colnames(data_clean)))
imp_data<-data_clean %>% select(-all_of(static_vars))

imp_data<-update_columns(data=imp_data, ind=colnames(imp_data)[sapply(imp_data, is.character)], what=factor)
bad_class<-colnames(imp_data)[!sapply(imp_data, function(x) is.numeric(x) | is.factor(x))]
if(length(bad_class)>0) stop("Variables with a class not accepted by the imputation: ", paste(bad_class, collapse=", "))

set.seed(seed)
data_single<-missRanger::missRanger(imp_data, pmm.k=3, num.trees=100, maxiter=5, seed=seed, verbose=0)
data_single<-bind_cols(data_clean %>% select(all_of(static_vars)), data_single)

# variables depending on imputed measurements are (re)computed after imputation
data_single<-derive_after_imputation(data_single)

#----------------------------------------------------------------
#### Skewness: numerical test, then log-transformation of skewed variables only ####

cand_vars<-colnames(data_single)[sapply(data_single, is.numeric)]
cand_vars<-setdiff(cand_vars, c(id_var, cluster_var, structural_vars))
cand_vars<-cand_vars[!grepl(no_transform_pattern, cand_vars)]

skew_tab<-bind_rows(lapply(cand_vars, function(v){
  x<-data_single[[v]]
  tibble(variable=v,
         n_unique=length(unique(x)),
         skewness=psych::skew(x),
         per_zero=mean(x==0)*100,
         min_value=min(x),
         offset=ifelse(min(x)>0, 0, ifelse(any(x>0), min(x[x>0])/2, 1))) # offset only needed if zeros are present
})) %>%
  mutate(action=case_when(n_unique<=10 ~ "none (few distinct values)",
                          abs(skewness)<=thr_skew ~ "none (not skewed)",
                          min_value<0 ~ "none (negative values)",
                          per_zero>thr_zero ~ "none (zero-inflated)",
                          TRUE ~ "log"))

data_single<-transform_skewed(data_single, skew_tab)

# skewness after transformation (control)
skew_tab$skewness_after<-sapply(seq_len(nrow(skew_tab)), function(i){
  v<-ifelse(skew_tab$action[i]=="log", paste0("log_", skew_tab$variable[i]), skew_tab$variable[i])
  psych::skew(data_single[[v]])
})
print(skew_tab, n=Inf)

#----------------------------------------------------------------
#### Multiple imputation (mice) ####

imp<-mice::mice(imp_data, m=m_imp, maxit=maxit_imp, seed=seed, printFlag=FALSE) # default methods: pmm (numeric), logreg (binary), polyreg (nominal), polr (ordinal)

# long format: one block of rows per imputed dataset (.imp = dataset number, .id = row number of imp_data)
data_multiple<-mice::complete(imp, action="long", include=FALSE)

# putting identifiers, dates and structural variables back (same row order as imp_data)
data_static<-data_clean %>% select(all_of(static_vars)) %>% mutate(.id=row_number())
data_multiple<-left_join(data_multiple, data_static, by=".id")

# list of m datasets, with the variables depending on imputed measurements recomputed in each of them
data_multiple_list<-split(data_multiple, data_multiple$.imp)
data_multiple_list<-lapply(data_multiple_list, derive_after_imputation)
data_multiple_list<-lapply(data_multiple_list, transform_skewed, skew_tab=skew_tab)
data_multiple<-bind_rows(data_multiple_list)

# # Rubin's rules: example of the generic principle (the model is only illustrative)
# fit_list<-lapply(data_multiple_list, function(d) lm(apwv ~ asthma_ever + age_years + sex, data=d))
# summary(mice::pool(mice::as.mira(fit_list)))

#----------------------------------------------------------------
#### One-hot encoding ####

var_dict_final<-build_var_dict(data_single, id_vars=c(id_var, cluster_var), type_override=type_override)
data_onehot<-one_hot_encode(data_single, var_dict_final, remove_first_dummy=FALSE)

#----------------------------------------------------------------
#### Exporting and saving results ####

# Data management log
log_data<-log_data %>% rowwise %>% mutate_if(is.numeric,test_format) %>% ungroup
write.table(log_data, paste("Data management log_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)

# Variable dictionary
write.table(var_dict, paste("Variable dictionary_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)

# Outlier table
outliers_tab<-outliers_tab %>% rowwise %>% mutate_if(is.numeric,test_format) %>% ungroup
write.table(outliers_tab, paste("Outliers table_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)

# Missing data tables
miss_var<-miss_var %>% rowwise %>% mutate_if(is.numeric,test_format) %>% ungroup
write.table(miss_var, paste("Missing data by variable_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)

# Skewness table
skew_tab<-skew_tab %>% rowwise %>% mutate_if(is.numeric,test_format) %>% ungroup
write.table(skew_tab, paste("Skewness table_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)

# Clean dataset without imputation
write.table(data_clean, paste("Clean dataset_no imputation_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)

# Single imputation clean dataset
write.table(data_single, paste("Clean dataset_single imputation_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)

# Single imputation clean dataset, one-hot encoded
write.table(data_onehot, paste("Clean dataset_single imputation_one hot encoded_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)

# Multiple imputation clean dataset
write.table(data_multiple, paste("Clean dataset_multiple imputation_", Sys.Date(), ".csv", sep=""), sep=";", dec=".", row.names=FALSE, col.names=TRUE)

# Mice object (to reuse it without re-running the imputation)
qs2::qs_save(imp, paste("mice dataset_", Sys.Date(), ".qs2", sep=""))

# Session information (reproducibility)
writeLines(capture.output(sessionInfo()), paste("Session info_", Sys.Date(), ".txt", sep=""))
