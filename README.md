# Comparing-Simulated-and-Synthetic-Data-Types-Crohns




## Description

This contains R code for creating a simulation trial of Ustekinumab data and corresponding synthetic datasets including analysis and graphs. 


Simulated clinical trial data was generated for three different types for trials (randomized control trial, observational study, and external data) based on probabilities from Ustekinumab data. The outcome being looked at is reduction of Crohn's Disease symptoms, and the treatment/exposure variable of interest is the Ustekinumab drug.  References for original datasets can be found below.

Synthetic data was generated using the R library Synthpop. Synthetic data is created by subsetting the original data to include only variables of interest, and then creating a cookbook of values.  Three methods of synthetic data generation were used, each corresponding to its own R script: CART (Categorical and Regression Tree modelling), random sampling, and linear/logistic regression. More information on how Synthpop generates synthetic data can be found here:
https://www.synthpop.org.uk/about-synthpop.html#methodology

The simulated and synthetic data are compared between data types and sample sizes by measuring the treatment effect using the chi-squared analysis.  The standard mean difference was also calulcated and graphed.

Sim_Synth_Ustek.R contains the pipeline for synthetic data generated using the CART method.
Sim_Synth_Ustek_Rand.R contains the pipeline for synthetic data generated using the random sampling method.
Sim_Synth_Ustek_LL.R contains the pipeline for synthetic data generated using the linear/logisitc regression method.


## Installation

The script can be downloaded and run in R Studio.

## Libraries

The R Libraries used include the following:


- tidyverse
- sandwich
- stargazer
- truncnorm
- tableone
- survival
- ggplot2
- ggmap
- table1
- lubridate
- plyr
- reshape
- MASS
- reshape2
- synthpop
- dplyr
- stddiff
- Matching
- survey
- ggstance
- hrbrthemes
- viridis
- rgp
- rsimsum



## Usage

This code is for comparing the quality of synthetic controls between different data types.  The simulated probabilties can be altered to fit new datasets for further testing.

## Citations

### More information on Synthpop Synthesis Methods

- Nowok, B., Raab, G.M. and Dibben, C. (2016) ‘synthpop: Bespoke Creation of Synthetic Data in R’, Journal of Statistical Software, 74, pp. 1–26. Available at: https://doi.org/10.18637/jss.v074.i11.


### Randomized Control Trial Original Data

- Feagan, B.G. et al. (2016) ‘Ustekinumab as Induction and Maintenance Therapy for Crohn’s Disease’, New England Journal of Medicine, 375(20), pp. 1946–1960. Available at: https://doi.org/10.1056/NEJMoa1602773.


### Observational Study Original Data

- Biemans, V.B.C. et al. (2020) ‘Ustekinumab for Crohn’s Disease: Results of the ICC Registry, a Nationwide Prospective Observational Cohort Study’, Journal of Crohn’s and Colitis, 14(1), pp. 33–45. Available at: https://doi.org/10.1093/ecco-jcc/jjz119.

### External Original Data

- Bello, F. et al. (2024) ‘Long-term real-world data of ustekinumab in Crohn’s disease: the Stockholm ustekinumab study’, Therapeutic Advances in Gastroenterology, 17, p. 17562848241242700. Available at: https://doi.org/10.1177/17562848241242700.

