**Collection of reproducible data analysis scripts in R**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17752875.svg)](https://doi.org/10.5281/zenodo.17752875)

This repository provides a curated collection of R scripts designed to support reproducible data analysis.

The scripts have been developed and used in several of my publications and cover a range of analytical approaches including survival analysis, Bayesian model averaging, and machine learning methods such as Extreme Gradient Boosting (XGBoost).

---

***Overview***

The primary goal of this repository is to promote transparency, reproducibility, and reusability of statistical analyses.

Each script is fully documented and includes, where applicable:

•	Example datasets for demonstration

•	Usage instructions

•	References to associated publications

Researchers and users are encouraged to adapt these scripts for their own data and workflows.

---

***Contents***

•	*Cox proportional hazards regression.R*

Implementation of a comprehensive Cox proportional hazards regression analysis, suitable for time-to-event data.

•	*Survival XGBoost Cox analysis.R*

Implementation of Extreme Gradient Boosting (XGBoost) adapted for Cox proportional hazards modeling, suitable for time-to-event data.

•	*XGBoost_binary classification_nested CV.R*

Implements XGBoost for binary classification tasks with nested cross-validation, allowing for robust model evaluation and hyperparameter tuning.

•	*XGBoost_binary classification_k-fold nested CV.R*

Implements XGBoost for binary classification tasks with 50x repeated 5-fold nested cross-validation, allowing for robust model evaluation and hyperparameter tuning.

•	*XGBoost_linear regression_nested CV.R*

Implements XGBoost for linear regression task with nested cross-validation, allowing for robust model evaluation and hyperparameter tuning.

•	*XGBoost_ordinal logistic classification.R*

Implements XGBoost for logistic ordinal classification tasks with nested cross-validation, allowing for robust model evaluation and hyperparameter tuning.

•	*Multi-class XGBoost.R*

Implements a multi-class prediction using XGBoost. Three reproducible and generalizable examples are provided:

-Example 1: Predicting epidemic cluster based on epidemic curves
  
-Example 2: Predicting chemical exposure cluster based on age, BMI, physical activity
  
-Example 3: Predicting cardiovascular risk cluster based on age, BMI, physical activity, BP, cholesterol

•	*Continuous_Multioutput_XGBoost.R*

Implements a generalizable R and Python workflow for continuous multi-output regression using XGBoost.Three reproducible and generalizable examples are provided:

-Example 1: Predict epidemic parameters (beta, mu, S0) from simulated epidemic curves (SIR model).
  
-Example 2: Predict chemical exposure (exposure_A, exposure_B) from age, BMI, and physical activity level.
  
-Example 3: Predict blood pressure and cholesterol from age, BMI, and physical activity level..

•	*Linear BMA.R*

Performs Bayesian Model Averaging (BMA) for probabilistic model selection in linear regression.

•	*Logistic BMA.R*

Performs BMA for probabilistic model selection in logistic regression.

•	*SHAP_direction.R*

Function for determining feature direction based on SHAP values using Rcpp pairwise counting.

•	*BKMR.R*

Performs a Bayesian kernel machine regression (BKMR) analysis.

•	*WQS.R*

Performs a Weighted Quantile Sum regression (WQS) analysis.

•	*PSM.R*

Performs a propensity score matching (PSM) analysis.

•	*RCS.R*

Performs a Restricted Cubic Spline (RCS) analysis.

•	*RCT.R*

Performs a descriptive analysis, linear mixed effect models (continuous variables), logistic mixed effect models (binary variables), and cumulative link mixed effect models (ordinal variables) for analyzing randomized controlled trial (RCT) data.

•	*Convert_semi transparent color to opaque.R*

This function converts a semi-transparent color to its opaque equivalent, considering the visual effect on a white background.

•	*Loading and using a model.R*

This script loads a saved model, prepares new data, and makes predictions based on the model.

•	*Carbon_footprint.R*

Initial attempt to estimate energy consumption and carbon emissions from running R code. This function attempts to monitor CPU, RAM, and GPU usage of R code, estimate energy consumption, and calculate associated carbon emissions (gCO2).

•	*Randomization.R*

Creates random assignments for experiments (e.g., clinical trial).

•	*gmean.R*

Function computing the geometric mean of a numeric vector.

•	*GSD.R*

Function computing the geometric standard deviation of a numeric vector.

•	*Limit_value_exceedance_plot.R*

Assesses chemical concentration exceedance against a specified limit value and visualizes the results.

•	*Continuous_variable_distribution_fig.R*

Script for visualizing continuous variable distributions with multiple plot types.

---

***Requirements***

•	R version ≥ 4.5.0

•	Additional dependencies are listed in the header of each script.

---

***Usage***
1.	Clone or download this repository:

                                 git clone https://github.com/PetitPascal/R-scripts.git
  	

3.	Open R or RStudio.

4.	Load the desired script and follow the inline documentation and example usage.
   
5.	Ensure all required packages are installed before running the scripts.

---

***Reproducibility and citation***

All scripts are provided to ensure transparent and reproducible analyses.

If you use any part of this repository in your research or publications, please cite the corresponding DOI of the release used.

•	Initial release: Collection of R scripts for reproducible data analysis.

•	Version v1.0.1: Includes additional logistic BMA and binary classification with XGBoost and nested CV scripts.

Refer to the DOI listed in the Releases.

All scripts are provided to promote transparency and reproducibility in data analysis.

If you use or adapt these scripts in your own work, please cite the corresponding release DOI:

Cite as (examples):

•	Petit, P. (2025). Collection of reproducible data analysis scripts in R (Version v1.0.1). [R]. Zenodo.
https://doi.org/10.5281/zenodo.17047856.

•	Petit, P. (2025). Weighted survival XGBoot with R [R script]. In R-scripts (Version 1.0.0) [R]. Zenodo. https://zenodo.org/records/17752875.

•	Petit, P. (2025). Logistic Bayesian Model Averaging (BMA) analysis for probabilistic model selection in logistic regression [R script]. In R-scripts (Version 1.0.1) [R]. Zenodo. https://zenodo.org/records/17752875.

•	Petit, P. (2025). Extreme Gradient Boosting (XGBoost) for binary classification with nested cross-validation [R script]. In R-scripts (Version 1.0.1) [R]. Zenodo. https://zenodo.org/records/17752875.

•	Petit, P. (2025). Bayesian kernel machine regression (BKMR) analysis [R script]. In R-scripts (Version 1.0.2) [R]. Zenodo. https://zenodo.org/records/17752875.

---

***Version History***

•	v1.0.0:	Initial release -	XGBoost for Cox proportional hazards modeling

•	v1.0.1:	Update release	- Added Logistic BMA and XGBoost (binary classification with nested CV)

•	v1.0.2:	Update release	- Added several new scripts

---

***Related Publications***

These scripts have been developed to support analyses in the following works:

•	Petit P, Nübel J, Walter J, Butter C, Heinze M, Ignatyev J, Haase-Fielitz A, Vuillerme N, Muehlensiepen F. Telemedicine adoption in cardiology: determinants and predictors identified using Bayesian model averaging and machine learning. PLOS Digit Health. 2026;5(4):e0001359. doi: 10.1371/journal.pdig.0001359.

•	Petit P, Vuillerme N, Gehrmann J, Stephan J, Muehlensiepen F, Nübel J, Hahn F, Martens E. Determinants and predictors of telemedicine use among physicians in the German outpatient sector: a secondary analysis of a web-based survey. Submitted.

•	Muehlensiepen F, Petit P, Knitza J, Welcker M, Vuillerme N. Prediction of the acceptance of telemedicine among rheumatic patients. A machine learning powered secondary analysis of German survey data. Rheumatol Int. 2024;44(3):523-534. doi: 10.1007/s00296-023-05518-9.

•	Muehlensiepen F, Petit P, Knitza J, Welcker M, Vuillerme N. Identification of Motivational Determinants for Telemedicine Use Among Patients With Rheumatoid Arthritis in Germany: Secondary Analysis of Data From a Nationwide Cross-Sectional Survey Study. J Med Internet Res. 2024;26:e47733. doi: 10.2196/47733.

•	Muehlensiepen F, Petit P, Knitza J, Welcker M, Vuillerme N. Factors Associated With Telemedicine Use Among Patients With Rheumatic and Musculoskeletal Disease: Secondary Analysis of Data From a German Nationwide Survey. J Med Internet Res. 2023;25:e40912. doi: 10.2196/40912.

•	Muehlensiepen F, Petit P, Knitza J, Welcker M, Vuillerme N. Factors associated with Telemedicine Usage among German General Practitioners and Rheumatologists: Secondary Analysis of Data from a Nationwide Survey. J Med Internet Res. 2022;24(11):e40304. doi: 10.2196/40304.

•	Muehlensiepen F, May S, Seifert F, Wengemuth E, Johannsen O, Middeke M, Heinze M, Petit P, Vuillerme N, Spethmann S, Bruch D. User Profiles and Engagement in a Hypertension Self-Management App: Cross-Sectional Survey. J Med Internet Res. 2026;28:e83075. doi: 10.2196/83075.

•	May S, Darkow R, Knitza J, Boy K, Klemm P, Heinze M, Vuillerme N, Petit P, Steffens-Korbanka P, Kladny H, Hornig J, Aries P, Welcker M, Muehlensiepen F. Digital Transformation of Rheumatology Care in Germany: Cross-Sectional National Survey. J Med Internet Res. 2025;27:e52601. doi: 10.2196/52601.

•	Grube L, Petit P, Vuillerme N, Nitschke M, Nwosu OB, Knitza J, Krusche M, Seifer AK, Eskofier BM, Schett G, Morf H. Complementary App-based Yoga Home Exercise Therapy for patients with axial Spondyloarthritis: A Usability Study. JMIR Form Res. 2024;8:e57185. doi: 10.2196/57185. 

•	Blaskowitz PPVA et al. Impact of the digital health application ViViRA on spinal mobility, physical function, quality of life and pain perception in spondyloarthritides patients: a randomized controlled trial. Arthritis Res Ther. 2024. 2024;26(1):208. doi: 10.1186/s13075-024-03443-1.

•	Petit P, Berger F, Bonneterre V, Vuillerme N. Exploring Alzheimer's Disease Risk Factors in Farmers with Explainable Machine Learning and Administrative Health Data. Submitted.

•	Petit P, Berger F, Bonneterre V, Vuillerme N. Leveraging Machine Learning with Real-World Data to Identify exposomic Risk Factors in Parkinson’s Disease among Farmers. Submitted.

•	Petit P, Bonneterre V, Vuillerme N. Using Machine Learning and Nationwide Population-based Data to Unravel Predictors of Treated Depression in Farmers. Ment Illn. 2025;2025:17. doi: 10.1155/mij/5570491.

•	Petit P, Berger F, Bonneterre V, Vuillerme N. Investigating Parkinson’s disease risk across farming activities using data mining and large-scale administrative health data. npj Parkinsons Dis. 2025;11:13. doi: 10.1038/s41531-024-00864-2.

•	Petit P, Leroyer A, Chamot S, Fumery M, Bonneterre V. Farming activities and risk of inflammatory bowel disease: a French nationwide population-based cohort study. J Crohns Colitis. 2024;jjae050. doi: 10.1093/ecco-jcc/jjae050.

•	Petit P, Chamot S, Al-Salameh A, Cancé C, Desailloud R, Bonneterre V. Farming activity and risk of treated thyroid disorders: insights from the TRACTOR project, a nationwide cohort study. Environ Res. 2024;249:118458. doi: 10.1016/j.envres.2024.118458.

•	Petit P, Gondard E, Gandon G, Moreaud O, Sauvée M, Bonneterre V. Agricultural activities and risk of Alzheimer’s disease: the TRACTOR project, a nationwide retrospective cohort study. Eur J Epidemiol. 2024;39(3):271-287. doi: 10.1007/s10654-023-01079-0. 

•	Petit P, Gandon G, Dubuc M, Vuillerme N, Bonneterre V. Agricultural activities and risk of treatment for depressive disorders among the entire French agricultural workforce: the TRACTOR project, a nationwide retrospective cohort study. Lancet Reg Health Eur. 2023;31:100674. doi: 10.1016/j.lanepe.2023.100674.

•	Petit P, Gandon G, Chabardes S, Bonneterre V. Agricultural activities and risk of central nervous system tumors among French farm managers: results from the TRACTOR project. Int J Cancer. 2022;151(10):1737-1749. doi: 10.1002/ijc.34197.

• Petit P, Bosson-Rieutort D, Maugard C, Gondard E, Ozenfant D, Joubert N, François O, Bonneterre V. The TRACTOR Project: TRACking and MoniToring Occupational Risks in Agriculture Using French Insurance Health Data (MSA). Ann Work Expo Health. 2022;66(3);402-411. doi: 10.1093/annweh/wxab083.

•	Chamot S, Petit P, Al-Salameh A, Bonneterre V, Cancé C, Decocq G, Desailloud R. Environmental Pollution and the Risk of Congenital Hypothyroidism: Insights from a French Nationwide Retrospective Ecological Cohort Study. J Hazard Mater Adv. 2025;17:100560. doi : 10.1016/j.hazadv.2024.100560.

•	Chamot S, Al-Salameh A, Balcaen T, Petit P, Bonneterre V, Cancé C, Desailloud R. Congenital and acquired hypothyroidism: temporal and spatial trends in France from 2014 to 2019. Ann Epidemiol. 2024;98:18-24. doi: 10.1016/j.annepidem.2024.07.091.

•	Chamot S, Al-Salameh A, Petit P, Bonneterre V, Cancé C, Decocq G, Boullier A, Braun K, Desailloud R. Does prenatal exposure to multiple airborne and tap-water pollutants increase neonatal thyroid-stimulating hormone concentrations? Data from the Picardy region, France. Sci Tot Environ. 2023;905:167089. doi: 10.1016/j.scitotenv.2023.167089.

•	Aix ML, Petit P, Bicout DJ. Air pollution and health impacts during the COVID-19 lockdowns in Grenoble, France. Environ Pollut. 2022;303:119134. doi: 10.1016/j.envpol.2022.119134.

•	Clauzel A, Persoons R, Maître A, Balducci F, Petit P. Review of environmental airborne pyrene/benzo[a]pyrene levels from industrial emissions for the improvement of 1-hydroxypyrene biomonitoring interpretation. J Toxicol Environ Health B Crit Rev. 2024;27(5-6):1-21. doi: 10.1080/10937404.2024.2362632.

• Valière M, Petit P, Persoons R, Demeilliers C, Maître A. Consistency between air and biological monitoring for assessing polycyclic aromatic hydrocarbon exposure and cancer risk of workers. Environ Res. 2022;112268. doi: 10.1016/j.envres.2021.112268. 

•	Petit P, Bicout DJ. Health risk assessment with multiple reference indices. Sci Total Environ. 2022;804:149971. doi: 10.1016/j.scitotenv.2021.149971.

•	Petit P, Maître A, Bicout D. A consensus approach for estimating health risk: application to inhalation cancer risks. Environ Res. 2021;196:110436. doi: 10.1016/j.envres.2020.110436.

•	Petit P, Maître A, Persoons R, Bicout DJ. Lung cancer risk assessment for workers exposed to polycyclic aromatic hydrocarbons in various industries. Environ Int. 2019;124:109-120. doi: 10.1016/j.envint.2018.12.058.

•	Petit P, Maître A, Persoons R, Bicout DJ. Modeling the exposure functions of atmospheric polycyclic aromatic hydrocarbon mixtures in occupational environments. Sci Tot Environ. 2017;584-585:1185-1197. doi: 10.1016/j.scitotenv.2017.01.182.

•	Petit P, Bicout DJ, Persoons R, Bonneterre V, Barbeau D, Maître A. Constructing a database of similar exposure groups: the application of the Exporisq-HAP database from 1995 to 2015. Ann Work Expo Health. 2017;61(4):440-456. doi: 10.1093/annweh/wxx017.

•	Maître A, Petit P, Marques M, Hervé C, Montlevier S, Persoons R, Bicout DJ. Exporisq-HAP database: 20 years of monitoring French occupational exposure to polycyclic aromatic hydrocarbon mixtures and identification of exposure determinants. Int J Hyg Environ Health. 2018;221(2):334-346. doi: 10.1016/j.ijheh.2017.12.008.

•	Petit P. Toxicological and Exposure Database Inventory: a review. Int J Hyg Environ Health. 2022;246:114055. doi: 10.1016/j.ijheh.2022.114055.

•	Choueiri J, Petit P, Balducci F, Bicout DJ, Demeilliers C. Literature-based inventory of chemical substance concentrations measured in organic food consumed in Europe. Data. 2024;9(7):89. doi: 10.3390/data9070089.

•	Milane T, Vuillerme N, Petit P, Warmerdam E, Romijnders R, Bianchini E, Maetzler W, Hansen C. Impact of forward and backward walking on gait parameters across Parkinson’s disease stages and severity: a prospective observational study. BMC Neurol. 2025;29(1):379. doi : 10.1186/s12883-025-04321-2.

•	Chardon M, Barbieri FA, Hansen C, Petit P, Vuillerme N. Impact of Overweight on Spatial-Temporal Gait Parameters During Obstacle Crossing in Young Adults: A Cross-Sectional Study. Sensors. 2024;24(23):7867. doi: 10.3390/s24237867.

•	Chardon M, Barbieri FA, Petit P, Vuillerme N. Reliability of obstacle-crossing parameters during overground walking in young adults. Sensors. 2024;24(11):3387. doi: 10.3390/s24113387.

---

***Contact***

For questions regarding the scripts or associated studies:

**Pascal Petit**

email: pascal.petit@univ-grenoble-alpes.fr

•	*ORCID*: https://orcid.org/0000-0001-9015-5230

•	*ResearchGate*: https://www.researchgate.net/profile/Pascal-Petit-3

•	*Google Scholar*: https://scholar.google.fr/citations?user=ja8PT6MAAAAJ&hl=fr

•	*Web Of Science*: https://www.webofscience.com/wos/author/record/M-4351-2017

•	*HAL*: https://hal.science/search/index/q/*/authIdHal_s/pascal-petit

•	*Thèse.fr*: https://theses.fr/223750166

***Current affiliation***: Univ. Grenoble Alpes, CNRS, Grenoble INP*, LIG, 38000 Grenoble, France

*Institute of Engineering Univ. Grenoble Alpes

***Former affiliations***:

•	Univ. Grenoble Alpes, AGEIS, 38000 Grenoble, France

•	Univ. Grenoble Alpes, CNRS, UMR 5525, VetAgro Sup, Grenoble INP, TIMC, 38000 Grenoble, France
                      
•	CHU Grenoble Alpes, Centre Régional de Pathologies Professionnelles et Environnementales, 38000 Grenoble, France

---
***Funding***

My research has been partially supported by:

•	The French government, through the National Research Agency (ANR - Agence Nationale de la Recherche), under the *France 2030* program (MIAI Cluster), grant **ANR-23-IACL-0006** (February 2025 – present).

•	The French government, through the ANR, under the *Investissements d’avenir* program, grants **ANR-10-AIRT-0005** and **ANR-15-IDEX-0002** (September 2022 – April 2026).

•	The *Agence nationale de sécurité sanitaire de l’alimentation, de l’environnement et du travail* (ANSES), grants **2016-CRD-03_PPV16/534B** and **2018-CRD-14_PPV18** (October 2018 – December 2020).

•	*Mutualité Sociale Agricole* (MSA), grant **MSA-2020-STOP** (January 2021 – December 2022).

•	*Fondation pour la Recherche sur Alzheimer*, grant **2020-A-01** (January 2021 – December 2021).

•	*MIAI@Grenoble Alpes*, grant **ANR-19-P3IA-0003** (October 2018 – January 2025).

---
If you find these scripts useful, please star this repository and cite the DOI in your research!
