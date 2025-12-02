**Collection of reproducible data analysis scripts in R**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17752875.svg)](https://doi.org/10.5281/zenodo.17752875)

This repository provides a curated collection of R scripts designed to support reproducible data analysis.

The scripts have been developed and used in several of my publications and cover a range of analytical approaches including survival analysis, Bayesian model averaging, and machine learning methods such as Extreme Gradient Boosting (XGBoost).

***Overview***

The primary goal of this repository is to promote transparency, reproducibility, and reusability of statistical analyses.

Each script is fully documented and includes, where applicable:

•	Example datasets for demonstration

•	Usage instructions

•	References to associated publications

Researchers and users are encouraged to adapt these scripts for their own data and workflows.


***Contents***

•	*Survival XGBoost Cox analysis.R*

Implementation of Extreme Gradient Boosting (XGBoost) adapted for Cox proportional hazards modeling, suitable for time-to-event data.

•	*XGBoost_binary classification_nested CV.R*

Implements XGBoost for binary classification tasks with nested cross-validation, allowing for robust model evaluation and hyperparameter tuning.

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

***Requirements***

•	R version ≥ 4.5.0

•	Additional dependencies are listed in the header of each script.

***Usage***
1.	Clone or download this repository:

                                 git clone https://github.com/PetitPascal/R-scripts.git
  	

3.	Open R or RStudio.

4.	Load the desired script and follow the inline documentation and example usage.
   
5.	Ensure all required packages are installed before running the scripts.

***Reproducibility and citation***

All scripts are provided to ensure transparent and reproducible analyses.

If you use any part of this repository in your research or publications, please cite the corresponding DOI of the release used.

•	Initial release: Collection of R scripts for reproducible data analysis.

•	Version v1.0.1: Includes additional logistic BMA and binary classification with XGBoost and nested CV scripts.

Refer to the DOI listed in the Releases.

All scripts are provided to promote transparency and reproducibility in data analysis.

If you use or adapt these scripts in your own work, please cite the corresponding release DOI:

Cite as:

•	Petit, P. (2025). Collection of reproducible data analysis scripts in R (Version v1.0.1). [R]. Zenodo.
https://doi.org/10.5281/zenodo.17047856.

•	Petit, P. (2025). Weighted survival XGBoot with R [R script]. In R-scripts (Version 1.0.0) [R]. Zenodo. https://doi.org/10.5281/zenodo.17047856. 

•	Petit, P. (2025). Logistic Bayesian Model Averaging (BMA) analysis for probabilistic model selection in logistic regression [R script]. In R-scripts (Version 1.0.1) [R]. Zenodo. https://doi.org/10.5281/zenodo.17047856.

•	Petit, P. (2025). Extreme Gradient Boosting (XGBoost) for binary classification with nested cross-validation [R script]. In R-scripts (Version 1.0.1) [R]. Zenodo. https://doi.org/10.5281/zenodo.17047856.

***Version History***

•	v1.0.0:	Initial release -	XGBoost for Cox proportional hazards modeling

•	v1.0.1:	Update release	- Added Logistic BMA and XGBoost (binary classification with nested CV)

•	v1.0.2:	Update release	- Added several new scripts

***Related Publications***

These scripts have been developed to support analyses in the following works:

•	Petit P, Nübel J, Walter J, Butter C, Heinze M, Ignatyev J, Haase-Fielitz A, Vuillerme N, Muehlensiepen F. Predicting telemedicine use among healthcare professionals in cardiology: a Bayesian model averaging and machine learning analysis of German survey data. Submitted.

•	Petit P, Berger F, Bonneterre V, Vuillerme N. Exploring Alzheimer's Disease Risk Factors in Farmers with Explainable Machine Learning and Administrative Health Data. Submitted.

•	Petit P, Berger F, Bonneterre V, Vuillerme N. Leveraging Machine Learning with Real-World Data to Identify Risk Factors in Parkinson’s Disease among Farmers. Submitted.

•	Petit P, Bonneterre V, Vuillerme N. Using Machine Learning and Nationwide Population-based Data to Unravel Predictors of Treated Depression in Farmers. Ment Illn. 2025;2025:17. doi: 10.1155/mij/5570491.

•	Muehlensiepen F, Petit P, Knitza J, Welcker M, Vuillerme N. Prediction of the acceptance of telemedicine among rheumatic patients. A machine learning powered secondary analysis of German survey data. Rheumatol Int. 2024;44(3):523-534. doi: 10.1007/s00296-023-05518-9.

•	Muehlensiepen F, Petit P, Knitza J, Welcker M, Vuillerme N. Identification of Motivational Determinants for Telemedicine Use Among Patients With Rheumatoid Arthritis in Germany: Secondary Analysis of Data From a Nationwide Cross-Sectional Survey Study. J Med Internet Res. 2024;26:e47733. doi: 10.2196/47733.

•	Muehlensiepen F, Petit P, Knitza J, Welcker M, Vuillerme N. Factors Associated With Telemedicine Use Among Patients With Rheumatic and Musculoskeletal Disease: Secondary Analysis of Data From a German Nationwide Survey. J Med Internet Res. 2023;25:e40912. doi: 10.2196/40912.

•	Muehlensiepen F, Petit P, Knitza J, Welcker M, Vuillerme N. Factors associated with Telemedicine Usage among German General Practitioners and Rheumatologists: Secondary Analysis of Data from a Nationwide Survey. J Med Internet Res. 2022;24(11):e40304. doi: 10.2196/40304.

***Contact***

For questions, suggestions, feedback, or collaborations, please contact:

**Pascal Petit**

email: pascal.petit@univ-grenoble-alpes.fr

🌐 ORCID: https://orcid.org/0000-0001-9015-5230

🌐 ResearchGate: https://www.researchgate.net/profile/Pascal-Petit-3

🌐 Google Scholar: https://scholar.google.fr/citations?user=ja8PT6MAAAAJ&hl=fr

🌐 Web Of Science: https://www.webofscience.com/wos/author/record/M-4351-2017

🌐 HAL: https://hal.science/search/index/q/*/authIdHal_s/pascal-petit

🌐 Thèse.fr: https://theses.fr/223750166

🏛️ Current affiliation: Univ. Grenoble Alpes, AGEIS, 38000 Grenoble, France

🏛️Former affiliations:

•	Univ. Grenoble Alpes, CNRS, UMR 5525, VetAgro Sup, Grenoble INP, TIMC, 38000 Grenoble, France
                      
•	CHU Grenoble Alpes, Centre Régional de Pathologies Professionnelles et Environnementales, 38000 Grenoble, France

⭐ If you find these scripts useful, please star this repository and cite the DOI in your research!
