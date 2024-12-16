
# import packages
library(tidyverse)
library(sandwich)
library(stargazer)

library(truncnorm)
library(tableone)
library(survival)
library(ggplot2)
library(ggmap)
library(table1)
library(lubridate)
library(plyr)
library(reshape)
library(MASS)
library(reshape2)
library(synthpop)

library(dplyr)

library(stddiff)
library(Matching)
library(survey)
library(ggstance)

library(hrbrthemes)
library(viridis)

library(rgp)

library(rsimsum)

# creating scenarios for number of observations
# halve the n due to pre and post being added together
nsize <- function (scenario, datatype, simsyn) { 
  if (scenario == 1) { # scenario one is n = 100 to n = 100
    N <- 50
  }   else if (scenario == 2 ) { # scenario two is n = 20000 to n = 20000
    N <- 10000
  }   else if (scenario == 3 & datatype == "rct") { # scenario three is realistic sample sizes
    N <- 628
  }   else if (scenario == 3 & datatype == "ob"){
    N <- 221
  }   else if (scenario == 3 & datatype == "ex"){
    N <- 322
  }   else if (scenario == 4 & simsyn == "sim"){ # scenario four is n = 100 to n = 20000
    N <- 50
  }   else if (scenario == 4 & simsyn == "syn"){
    N <- 20000
  }
  return(N)
} 








# for loop for sim and analysis

nSim <- 40000 # number of simulations
set.seed(9999)


# creating results lists: 

sim_results_rct <- vector("list", length = nSim) # empty list to store df
sim_results_ob <- vector("list", length = nSim) 
sim_results_ex <- vector("list", length = nSim) 
syn_results_rct <- vector("list", length = nSim)
syn_results_ob <- vector("list", length = nSim)
syn_results_ex <- vector("list", length = nSim)

chi_results_rct <- vector("list", length = nSim) # empty list to store chi square results
chi_results_ob <- vector("list", length = nSim)
chi_results_ex <- vector("list", length = nSim)
chi_results_rct_syn <- vector("list", length = nSim)
chi_results_ob_syn <- vector("list", length = nSim)
chi_results_ex_syn <- vector("list", length = nSim)


table_rct_data <- vector("list", length = nSim) # create empty list for tables 
table_ob_data <- vector("list", length = nSim) 
table_ex_data <- vector("list", length = nSim) 


rct_age_smd <- vector("list", length = nSim) # create empty list for smd for all variables
#rct_race_smd <- vector("list", length = nSim)
rct_sex_smd <- vector("list", length = nSim)
rct_treat_smd <- vector("list", length = nSim)
rct_outcome_smd <- vector("list", length = nSim)

ob_age_smd <- vector("list", length = nSim) 
#ob_race_smd <- vector("list", length = nSim)
ob_sex_smd <- vector("list", length = nSim)
ob_exposed_smd <- vector("list", length = nSim)
ob_case_smd <- vector("list", length = nSim)

ex_age_smd <- vector("list", length = nSim) 
#ex_race_smd <- vector("list", length = nSim)
ex_sex_smd <- vector("list", length = nSim)
ex_treatment_smd <- vector("list", length = nSim)
ex_case_smd <- vector("list", length = nSim)



# for smd syn and sim comps:

rct_ex_age_smd <- vector("list", length = nSim) # create empty list for smd for all variables
#rct_race_smd <- vector("list", length = nSim)
rct_ex_sex_smd <- vector("list", length = nSim)
rct_ex_treat_smd <- vector("list", length = nSim)
rct_ex_case_smd <- vector("list", length = nSim)

ob_ex_age_smd <- vector("list", length = nSim) 
#ob_race_smd <- vector("list", length = nSim)
ob_ex_sex_smd <- vector("list", length = nSim)
ob_ex_treat_smd <- vector("list", length = nSim)
ob_ex_case_smd <- vector("list", length = nSim)

rct_ob_age_smd <- vector("list", length = nSim) 
#ex_race_smd <- vector("list", length = nSim)
rct_ob_sex_smd <- vector("list", length = nSim)
rct_ob_treat_smd <- vector("list", length = nSim)
rct_ob_case_smd <- vector("list", length = nSim)





rct_ex_age_smd_sim <- vector("list", length = nSim) # create empty list for smd for all variables
#rct_race_smd <- vector("list", length = nSim)
rct_ex_sex_smd_sim <- vector("list", length = nSim)
rct_ex_treat_smd_sim <- vector("list", length = nSim)
rct_ex_case_smd_sim <- vector("list", length = nSim)

ob_ex_age_smd_sim <- vector("list", length = nSim) 
#ob_race_smd <- vector("list", length = nSim)
ob_ex_sex_smd_sim <- vector("list", length = nSim)
ob_ex_treat_smd_sim <- vector("list", length = nSim)
ob_ex_case_smd_sim <- vector("list", length = nSim)

rct_ob_age_smd_sim <- vector("list", length = nSim) 
#ex_race_smd <- vector("list", length = nSim)
rct_ob_sex_smd_sim <- vector("list", length = nSim)
rct_ob_treat_smd_sim <- vector("list", length = nSim)
rct_ob_case_smd_sim <- vector("list", length = nSim)



# loop for gen simulations
for (i in 1:nSim) {
  
  print('Start')
  
  # rct
  # treatment
  
  scenario <- 1 
  if (i-1 %% 100 == 0) {
    scenario <- scenario + 1
    # changes scenario every 10000 iterations of loop
  }
  
  N <- nsize(scenario, "rct", "sim") # selecting scenario
  
  
  treat <- sample(c(1, 0), N, replace = TRUE, prob = c(0.5, 0.5))
  
  # compliance to treatment
  comply <- ifelse(treat == 1, sample(c(1, 0), N, replace = TRUE, prob = c(1, 0)), sample(c(1, 0), N, replace = TRUE, prob = c(0, 1)))
  
  # other variables
  sex <- rbinom(N, 1, 0.522)
  #race <- sample(c(1, 2, 3, 4), N, replace = TRUE, prob = c(0.079, 0.028, 0.018, 0.689))
  age <- rtruncnorm(N, a = 18, b = 72, mean = 37.3, sd = 11.8) # b is guessed
  age <- round(age, 0) # round age to whole numbers
  random_var <- rbinom(N, 1, prob = 0.5)
  post <- rep(0, N)  # pre and post assignment
  
  # create df combining variables
  df <- data.frame(treat, comply, sex, age, random_var, post)
  
  # add id
  df <- df %>% dplyr::mutate(id = row_number())
  
  # create post df
  df2 <- df %>% mutate(post = 1)
  
  # combine pre and post
  df_combined <- rbind(df, df2)
  
  # define treatment effect, only for post group
  df_combined <- df_combined %>%
    mutate(treat_effect =
             ifelse(treat == 1 & comply == 1, 0.488, # use remission rates, chose week 12 vs placebo
                    ifelse(treat == 1 & comply == 0, 0.359,
                           ifelse(treat == 0 & comply == 1, 0.359,
                                  ifelse(treat == 0 & comply == 0, 0.359, NA))))) %>%
    mutate(treat_effect = ifelse(post == 0, NA, treat_effect))  # Filter treatment effect to apply only to post-treatment group
  
  # calculate sex_effect 
  df_combined <- df_combined %>%
    mutate(sex_effect = ifelse(sex == 1 | sex == 0, 0, NA)) # Sex effect not recorded or reported / only baseline / assumed as zero
  
  # create overall treatment effect including demographic influences
  df_combined <- df_combined %>%
    mutate(overall_treat_effect = treat_effect + sex_effect)
  
  
  # create pre and post outcomes (baseline vs treatment)
  df_combined <- df_combined %>%
    mutate(outcome = ifelse(post == 0, rbinom(sum(post == 0), 1, 0.359), # placebo is used instead of baseline since it is soc
                            rbinom(sum(post == 1), 1, overall_treat_effect)))
  # store dfs in results
  sim_results_rct[[i]] <- df_combined
  
  # remove df2
  rm(df2)
  
  
  
  
  
  # observational
  # treatment
  
  N <- nsize(scenario, "ob", "sim") # selecting scenario
  
  exposed=sample(c(1,0), N,replace=TRUE, prob=c(.76,.14)) # probability of being exposed to the treatment or not
  
  
  # other variables
  sex <- rbinom(N, 1, 0.602)
  #race <- sample(c(1, 2, 3, 4), N, replace = TRUE, prob = c(0.079, 0.028, 0.018, 0.689))
  age <- rtruncnorm(N, a = 16, b = 65, mean = 38.2, sd = 16.96) # values estimated since only median, iqr, and minimum given
  age <- round(age, 0) # round age to whole numbers
  random_var <- rbinom(N, 1, prob = 0.5)
  post <- rep(0, N)  # pre and post assignment
  
  # create df combining variables
  ob_df <- data.frame(exposed, sex, age, random_var, post)
  
  # add id
  ob_df <- ob_df %>% dplyr::mutate(id = row_number())
  
  # create post df
  ob_df2 <- ob_df %>% mutate(post = 1)
  
  # combine pre and post
  ob_df_combined <- rbind(ob_df, ob_df2)
  
  # define treatment effect, only for post group
  # define exposure effect
  ob_df_combined <- ob_df_combined %>% mutate(exposure_effect=
                                                ifelse(exposed==1, .387,
                                                       ifelse(exposed==0, .242,       
                                                              NA)))%>%
    mutate(exposure_effect = ifelse(post == 0, NA, exposure_effect))
  
  
  
  # calculate sex_effect and race_effect
  ob_df_combined <- ob_df_combined %>%
    mutate(sex_effect = ifelse(sex == 1 | sex == 0, 0, NA))
  
  # create overall treatment effect including demographic influences
  ob_df_combined <- ob_df_combined %>%
    mutate(overall_treat_effect = exposure_effect + sex_effect)
  
  
  # create pre and post outcomes (testing for COVID)
  ob_df_combined <- ob_df_combined %>%
    mutate(case = ifelse(post == 0, rbinom(sum(post == 0), 1, 0.242), 
                         rbinom(sum(post == 1), 1, overall_treat_effect)))
  # store dfs in results
  sim_results_ob[[i]] <- ob_df_combined
  
  # remove df2
  rm(ob_df2)
  
  
  
  
  # external
  # treatment
  N <- nsize(scenario, "ex", "sim") # selecting scenario
  
  ex_case=sample(c(1,0), N*2,replace=TRUE, prob=c(.64,.36))
  
  # other variables
  treatment=sample(c(1,0), N*2,replace=TRUE, prob=c(.88,.12))
  sex <- rbinom(N, 1, 0.48)
  #race <- sample(c(1, 2, 3, 4), N*2, replace = TRUE, prob = c(0.079, 0.028, 0.018, 0.689))
  age <- rtruncnorm(N*2, a = 15, b = 70, mean = 38, sd = 17.78) # none provided, all estimated
  age <- round(age, 0) # round age to whole numbers
  random_var <- rbinom(N*2, 1, prob = 0.5)
  
  
  
  # create df combining variables
  ex_df <- data.frame(ex_case, treatment, sex, age, random_var)
  
  # add id
  ex_df <- ex_df %>% dplyr::mutate(id = row_number())
  
  # store dfs in results
  sim_results_ex[[i]] <- ex_df
  
  
  
  
  
  
  
  
  
  
  
  
  # begin synthetic data creation:
  
  
  
  # MODEL - rct
  rct_mydata_con = subset(df_combined, select = c("id", "age"))
  
  # logreg requires no NA values
  rct_mydata_con <- na.omit(rct_mydata_con)
  
  
  
  
  
  codebook.syn(rct_mydata_con)
  
  
  
  N <- nsize(scenario, "rct", "syn") # selecting scenario
  
  # minimum number of observations needed is 10
  rct_mysyn_con <- syn(rct_mydata_con, method = "norm", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)
  
  summary(rct_mysyn_con)
  
  
  
  
  
  
  rct_mydata_cat = subset(df_combined, select = c( "sex", "treat", "outcome"))
  # logreg requires no NA values
  rct_mydata_cat <- na.omit(rct_mydata_cat)
  
  
  codebook.syn(rct_mydata_cat)
  
  
  
  N <- nsize(scenario, "rct", "syn") # selecting scenario
  
  # minimum number of observations needed is 10
  rct_mysyn_cat <- syn(rct_mydata_cat, method = "logreg", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)
  
  summary(rct_mysyn_cat)
  
  
  # combining the two dataframes
  
  rct_mysyn_con <- rct_mysyn_con$syn
  rct_mysyn_cat <- rct_mysyn_cat$syn
  
  # Create a common key column in each dataframe
  rct_mysyn_con$Key <- 1:nrow(rct_mysyn_con)
  rct_mysyn_cat$Key <- 1:nrow(rct_mysyn_cat)
  
  # Merge the dataframes based on the common key column
  rct_mysyn <- merge(rct_mysyn_con, rct_mysyn_cat, by = "Key", all = TRUE)
  
  
  # drop key
  rct_mysyn = subset(rct_mysyn, select = c("id", "age",  "sex", "treat", "outcome"))
  
  
  
  # store dfs in results
  syn_results_rct[[i]] <- rct_mysyn
  
  
  
  
  
  
  
  
  
  
  # MODEL - ob
  ob_mydata_con = subset(ob_df_combined, select = c("id", "age"))
  # logreg requires no NA values
  ob_mydata_con <- na.omit(ob_mydata_con)
  
  
  codebook.syn(ob_mydata_con)
  
  
  
  N <- nsize(scenario, "ob", "syn") # selecting scenario
  
  # minimum number of observations needed is 10
  ob_mysyn_con <- syn(ob_mydata_con, method = "norm", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)
  
  summary(ob_mysyn_con)
  
  
  
  
  
  
  ob_mydata_cat = subset(ob_df_combined, select = c( "sex", "exposed", "case"))
  # logreg requires no NA values
  ob_mydata_cat <- na.omit(ob_mydata_cat)
  
  
  codebook.syn(ob_mydata_cat)
  
  
  
  N <- nsize(scenario, "ob", "syn") # selecting scenario
  
  # minimum number of observations needed is 10
  ob_mysyn_cat <- syn(ob_mydata_cat, method = "logreg", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)
  
  summary(ob_mysyn_cat)
  
  
  # combining the two dataframes
  
  ob_mysyn_con <- ob_mysyn_con$syn
  ob_mysyn_cat <- ob_mysyn_cat$syn
  
  # Create a common key column in each dataframe
  ob_mysyn_con$Key <- 1:nrow(ob_mysyn_con)
  ob_mysyn_cat$Key <- 1:nrow(ob_mysyn_cat)
  
  # Merge the dataframes based on the common key column
  ob_mysyn <- merge(ob_mysyn_con, ob_mysyn_cat, by = "Key", all = TRUE)
  
  
  # drop key
  ob_mysyn = subset(ob_mysyn, select = c("id", "age",  "sex", "exposed", "case"))
  
  
  
  # store dfs in results
  syn_results_ob[[i]] <- ob_mysyn
  
  
  
  
  
  
  # MODEL - ex
  ex_mydata_con = subset(ex_df, select = c("id", "age"))
  # logreg requires no NA values
  ex_mydata_con <- na.omit(ex_mydata_con)
  
  
  codebook.syn(ex_mydata_con)
  
  
  
  N <- nsize(scenario, "ex", "syn") # selecting scenario
  
  # minimum number of observations needed is 10
  ex_mysyn_con <- syn(ex_mydata_con, method = "norm", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)
  
  summary(ex_mysyn_con)
  
  
  
  
  
  
  ex_mydata_cat = subset(ex_df, select = c( "sex", "treatment", "ex_case"))
  # logreg requires no NA values
  ex_mydata_cat <- na.omit(ex_mydata_cat)
  
  
  codebook.syn(ex_mydata_cat)
  
  
  
  N <- nsize(scenario, "ob", "syn") # selecting scenario
  
  # minimum number of observations needed is 10
  ex_mysyn_cat <- syn(ex_mydata_cat, method = "logreg", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)
  
  summary(ex_mysyn_cat)
  
  
  # combining the two dataframes
  
  ex_mysyn_con <- ex_mysyn_con$syn
  ex_mysyn_cat <- ex_mysyn_cat$syn
  
  # Create a common key column in each dataframe
  ex_mysyn_con$Key <- 1:nrow(ex_mysyn_con)
  ex_mysyn_cat$Key <- 1:nrow(ex_mysyn_cat)
  
  # Merge the dataframes based on the common key column
  ex_mysyn <- merge(ex_mysyn_con, ex_mysyn_cat, by = "Key", all = TRUE)
  
  
  # drop key
  ex_mysyn = subset(ex_mysyn, select = c("id", "age",  "sex", "treatment", "ex_case"))
  
  
  
  # store dfs in results
  syn_results_ex[[i]] <- ex_mysyn
  
  
  
  
  # testing treatment effect
  try(test_rct <- chisq.test(table(df_combined$treat, df_combined$outcome)))
  try(chi_results_rct[[i]] <- test_rct)
  try(test_rct_new <- chisq.test(table(rct_mysyn$treat, rct_mysyn$outcome)))
  try(chi_results_rct_syn[[i]] <- test_rct_new)
  
  try(test_ob <- chisq.test(table(ob_df_combined$exposed, ob_df_combined$case)))
  try(chi_results_ob[[i]] <- test_ob)
  try(test_ob_new <- chisq.test(table(ob_mysyn$exposed, ob_mysyn$case)))
  try(chi_results_ob_syn[[i]] <- test_ob_new)
  
  try(test_ex <- chisq.test(table(ex_df$treatment, ex_df$ex_case)))
  try(chi_results_ex[[i]] <- test_ex)
  try(test_ex_new <- chisq.test(table(ex_mysyn$treatment, ex_mysyn$ex_case)))
  try(chi_results_ex_syn[[i]] <- test_ex_new)
  
  
  
  
  
  
  
  
  # calculating smd
  
  
  
  
  # create a common key column in each dataframe
  df_combined$Key <- 1:nrow(df_combined)
  rct_mysyn$Key <- 1:nrow(rct_mysyn)
  
  # merge the dataframes based on the common key column
  rct <- merge(df_combined, rct_mysyn, by = "Key", all = TRUE)
  
  
  # drop key
  rct = rct[,-1]
  
  rct_age <- data.frame(age = c(rct[,"age.x"], rct[,"age.y"]))
  #rct_race<- data.frame(race = c(rct[,"race.x"], rct[,"race.y"]))
  rct_sex<- data.frame(sex = c(rct[,"sex.x"], rct[,"sex.y"]))
  rct_treat <- data.frame(treat = c(rct[,"treat.x"], rct[,"treat.y"]))
  rct_outcome <- data.frame(outcome = c(rct[,"outcome.x"], rct[,"outcome.y"]))
  
  rct_all <- data.frame(rct_age$age,  rct_sex$sex, rct_treat$treat, rct_outcome$outcome)
  
  # table for sd
  ## Covariates
  vars <- c("rct_age.age", "rct_race.race", "rct_sex.sex", "rct_treat.treat", "rct_outcome.outcome"
  )
  
  ## Construct a table
  table_rct_data[[i]] <- CreateTableOne(vars = vars, data = rct_all, test = FALSE)
  
  
  
  
  
  
  
  
  # smd 
  
  rct_age_sd <- sd(rct_all$rct_age.age, na.rm = T)
  rct_age_md <- mean(rct$age.x, na.rm = T)-mean(rct$age.y, na.rm = T)
  rct_age_smd[[i]] <- rct_age_md/rct_age_sd
  

  rct_sex_sd <- sd(rct_all$rct_sex.sex, na.rm = T)
  rct_sex_md <- mean(rct$sex.x, na.rm = T)-mean(rct$sex.y, na.rm = T)
  rct_sex_smd[[i]] <- rct_sex_md/rct_sex_sd
  
  rct_treat_sd <- sd(rct_all$rct_treat.treat, na.rm = T)
  rct_treat_md <- mean(rct$treat.x, na.rm = T)-mean(rct$treat.y, na.rm = T)
  rct_treat_smd[[i]] <- rct_treat_md/rct_treat_sd
  
  rct_outcome_sd <- sd(rct_all$rct_outcome.outcome, na.rm = T)
  rct_outcome_md <- mean(rct$outcome.x, na.rm = T)-mean(rct$outcome.y, na.rm = T)
  rct_outcome_smd[[i]] <- rct_outcome_md/rct_outcome_sd
  
  
  
  
  
  # obs
  
  
  
  # Create a common key column in each dataframe
  ob_df_combined$Key <- 1:nrow(ob_df_combined)
  ob_mysyn$Key <- 1:nrow(ob_mysyn)
  
  # Merge the dataframes based on the common key column
  ob <- merge(ob_df_combined, ob_mysyn, by = "Key", all = TRUE)
  
  
  # drop key
  ob = ob[,-1]
  
  ob_age <- data.frame(age = c(ob[,"age.x"], ob[,"age.y"]))
 # ob_race<- data.frame(race = c(ob[,"race.x"], ob[,"race.y"]))
  ob_sex<- data.frame(sex = c(ob[,"sex.x"], ob[,"sex.y"]))
  ob_exposed <- data.frame(exposed = c(ob[,"exposed.x"], ob[,"exposed.y"]))
  ob_case <- data.frame(case = c(ob[,"case.x"], ob[,"case.y"]))
  
  ob_all <- data.frame(ob_age$age,  ob_sex$sex, ob_exposed$exposed, ob_case$case)
  
  # table for sd
  ## Covariates
  vars <- c("ob_age.age", "ob_race.race", "ob_sex.sex", "ob_exposed.exposed", "ob_case.case"
  )
  
  ## Construct a table
  table_ob_data[[i]] <- CreateTableOne(vars = vars, data = ob_all, test = FALSE)
  
  
  
  
  
  
  
  # smd 
  
  ob_age_sd <- sd(ob_all$ob_age.age, na.rm = T)
  ob_age_md <- mean(ob$age.x, na.rm = T)-mean(ob$age.y, na.rm = T)
  ob_age_smd[[i]] <- ob_age_md/ob_age_sd
  

  ob_sex_sd <- sd(ob_all$ob_sex.sex, na.rm = T)
  ob_sex_md <- mean(ob$sex.x, na.rm = T)-mean(ob$sex.y, na.rm = T)
  ob_sex_smd[[i]] <- ob_sex_md/ob_sex_sd
  
  ob_exposed_sd <- sd(ob_all$ob_exposed.exposed, na.rm = T)
  ob_exposed_md <- mean(ob$exposed.x, na.rm = T)-mean(ob$exposed.y, na.rm = T)
  ob_exposed_smd[[i]] <- ob_exposed_md/ob_exposed_sd
  
  ob_case_sd <- sd(ob_all$ob_case.case, na.rm = T)
  ob_case_md <- mean(ob$case.x, na.rm = T)-mean(ob$case.y, na.rm = T)
  ob_case_smd[[i]] <- ob_case_md/ob_case_sd
  
  
  
  # ext
  
  
  
  # Create a common key column in each dataframe
  ex_df$Key <- 1:nrow(ex_df)
  ex_mysyn$Key <- 1:nrow(ex_mysyn)
  
  # Merge the dataframes based on the common key column
  ex <- merge(ex_df, ex_mysyn, by = "Key", all = TRUE)
  
  
  # drop key
  ex = ex[,-1]
  
  ex_age <- data.frame(age = c(ex[,"age.x"], ex[,"age.y"]))
 # ex_race<- data.frame(race = c(ex[,"race.x"], ex[,"race.y"]))
  ex_sex<- data.frame(sex = c(ex[,"sex.x"], ex[,"sex.y"]))
  ex_treatment <- data.frame(treatment = c(ex[,"treatment.x"], ex[,"treatment.y"]))
  ex_case <- data.frame(case = c(ex[,"ex_case.x"], ex[,"ex_case.y"]))
  
  ex_all <- data.frame(ex_age$age,  ex_sex$sex, ex_treatment$treatment, ex_case$case)
  
  # table for sd
  ## Covariates
  vars <- c("ex_age.age", "ex_race.race", "ex_sex.sex", "ex_treatment.treatment", "ex_case.case"
  )
  
  ## Construct a table
  table_ex_data[[i]] <- CreateTableOne(vars = vars, data = ex_all, test = FALSE)
  
  
  
  
  
  
  
  # smd 
  
  ex_age_sd <- sd(ex_all$ex_age.age, na.rm = T)
  ex_age_md <- mean(ex$age.x, na.rm = T)-mean(ex$age.y, na.rm = T)
  ex_age_smd[[i]] <- ex_age_md/ex_age_sd
  

  ex_sex_sd <- sd(ex_all$ex_sex.sex, na.rm = T)
  ex_sex_md <- mean(ex$sex.x, na.rm = T)-mean(ex$sex.y, na.rm = T)
  ex_sex_smd[[i]] <- ex_sex_md/ex_sex_sd
  
  ex_treatment_sd <- sd(ex_all$ex_treatment.treatment, na.rm = T)
  ex_treatment_md <- mean(ex$treatment.x, na.rm = T)-mean(ex$treatment.y, na.rm = T)
  ex_treatment_smd[[i]] <- ex_treatment_md/ex_treatment_sd
  
  ex_case_sd <- sd(ex_all$ex_case.case, na.rm = T)
  ex_case_md <- mean(ex$ex_case.x, na.rm = T)-mean(ex$ex_case.y, na.rm = T)
  ex_case_smd[[i]] <- ex_case_md/ex_case_sd
  
  
  
  
  
  
  
  
  
  
  
  
  
  # change column names to make sense:
  
  ex_all <- ex_all %>% 
    dplyr::rename(
      treat = ex_treatment.treatment,
      sex = ex_sex.sex,
      age = ex_age.age,
      case = ex_case.case
    )
  
  rct_all <- rct_all %>% 
    dplyr::rename(
      treat = rct_treat.treat,
      sex = rct_sex.sex,
      age = rct_age.age,
      case = rct_outcome.outcome
    )
  
  ob_all <- ob_all %>% 
    dplyr::rename(
      treat = ob_exposed.exposed,
      sex = ob_sex.sex,
      age = ob_age.age,
      case = ob_case.case
    )
  
  
  # smd synth compared to synth
  # x is sim, y is synth
  
  
  
  rct_ex_age <- data.frame(age = c(ex_all$age, rct_all$age))
  rct_ex_age_sd_syn <- sd(rct_ex_age$age, na.rm = T)
  rct_ex_age_md_syn <- mean(rct$age.y, na.rm = T)-mean(ex$age.y, na.rm = T)
  rct_ex_age_smd[[i]] <- rct_ex_age_md_syn/rct_ex_age_sd_syn
  
  
  
  rct_ex_sex  <- data.frame(sex = c(ex_all$sex, rct_all$sex))
  rct_ex_sex_sd_syn <- sd(rct_ex_sex$sex, na.rm = T)
  rct_ex_sex_md_syn <- mean(rct$sex.y, na.rm = T)-mean(ex$sex.y, na.rm = T)
  rct_ex_sex_smd[[i]] <- rct_ex_sex_md_syn/rct_ex_sex_sd_syn
  
  
  rct_ex_treat <- data.frame(treat = c(ex_all$treat, rct_all$treat))
  rct_ex_treat_sd_syn <- sd(rct_ex_treat$treat, na.rm = T)
  rct_ex_treat_md_syn <- mean(rct$treat.y, na.rm = T)-mean(rct$treat.y, na.rm = T)
  rct_ex_treat_smd[[i]] <- rct_ex_treat_md_syn/rct_ex_treat_sd_syn
  
  rct_ex_case  <- data.frame(case = c(ex_all$case, rct_all$case))
  rct_ex_case_sd_syn <- sd(rct_ex_case$case, na.rm = T)
  rct_ex_case_md_syn <- mean(rct$outcome.y, na.rm = T)-mean(ex$ex_case.y, na.rm = T)
  rct_ex_case_smd[[i]] <- rct_ex_case_md_syn/rct_ex_case_sd_syn
  
  
  
  
  
  
  
  rct_ob_age <-  data.frame(age = c(ob_all$age, rct_all$age))
  rct_ob_age_sd_syn <- sd(rct_ob_age$age, na.rm = T)
  rct_ob_age_md_syn <- mean(rct$age.y, na.rm = T)-mean(ob$age.y, na.rm = T)
  rct_ob_age_smd[[i]] <- rct_ob_age_md_syn/rct_ob_age_sd_syn
  
  rct_ob_sex <-data.frame(sex = c(ob_all$sex, rct_all$sex))
  rct_ob_sex_sd_syn <- sd(rct_ob_sex$sex, na.rm = T)
  rct_ob_sex_md_syn <- mean(rct$sex.y, na.rm = T)-mean(ob$sex.y, na.rm = T)
  rct_ob_sex_smd[[i]] <- rct_ob_sex_md_syn/rct_ob_sex_sd_syn
  
  rct_ob_treat <- data.frame(treat = c(ob_all$treat, rct_all$treat))
  rct_ob_treat_sd_syn <- sd(rct_ob_treat$treat, na.rm = T)
  rct_ob_treat_md_syn <- mean(rct$treat.y, na.rm = T)-mean(ob$exposed.y, na.rm = T)
  rct_ob_treat_smd[[i]] <- rct_ob_treat_md_syn/rct_ob_treat_sd_syn
  
  
  rct_ob_case <- data.frame(case = c(ob_all$case, rct_all$case))
  rct_ob_case_sd_syn <- sd(rct_ob_case$case, na.rm = T)
  rct_ob_case_md_syn <- mean(rct$outcome.y, na.rm = T)-mean(ob$case.y, na.rm = T)
  rct_ob_case_smd[[i]] <- rct_ob_case_md_syn/rct_ob_case_sd_syn
  
  
  
  ob_ex_age <- data.frame(age = c(ob_all$age, ex_all$age))
  ob_ex_age_sd_syn <- sd(ob_ex_age$age, na.rm = T)
  ob_ex_age_md_syn <- mean(ob$age.y, na.rm = T)-mean(ex$age.y, na.rm = T)
  ob_ex_age_smd[[i]] <- ob_ex_age_md_syn/ob_ex_age_sd_syn
  
  ob_ex_sex <- data.frame(sex = c(ob_all$sex, ex_all$sex))
  ob_ex_sex_sd_syn <- sd(ob_ex_sex$sex, na.rm = T)
  ob_ex_sex_md_syn <- mean(ob$sex.y, na.rm = T)-mean(ex$sex.y, na.rm = T)
  ob_ex_sex_smd[[i]] <- ob_ex_sex_md_syn/ob_ex_sex_sd_syn
  
  ob_ex_treat <- data.frame(treat = c(ob_all$treat, ex_all$treat))
  ob_ex_treat_sd_syn <- sd(ob_ex_treat$treat, na.rm = T)
  ob_ex_treat_md_syn <- mean(ob$exposed.y, na.rm = T)-mean(ex$treatment.y, na.rm = T)
  ob_ex_treat_smd[[i]] <- ob_ex_treat_md_syn/ob_ex_treat_sd_syn
  
  ob_ex_case <- data.frame(case = c(ob_all$case, ex_all$case))
  ob_ex_case_sd_syn <- sd(ob_ex_case$case, na.rm = T)
  ob_ex_case_md_syn <- mean(ob$case.y, na.rm = T)-mean(ex$ex_case.y, na.rm = T)
  ob_ex_case_smd[[i]] <- ob_ex_case_md_syn/ob_ex_case_sd_syn
  
  
  # repeat for sim
  
  # x is sim, y is synth
  
  rct_ex_age <- data.frame(age = c(rct_all$age, ex_all$age))
  rct_ex_age_sd_sim <- sd(rct_ex_age$age, na.rm = T)
  rct_ex_age_md_sim <- mean(rct$age.x, na.rm = T)-mean(ex$age.x, na.rm = T)
  rct_ex_age_smd_sim[[i]] <- rct_ex_age_md_sim/rct_ex_age_sd_sim
  
  
  
  rct_ex_sex <- data.frame(sex = c(rct_all$sex, ex_all$sex))
  rct_ex_sex_sd_sim <- sd(rct_ex_sex$sex, na.rm = T)
  rct_ex_sex_md_sim <- mean(rct$sex.x, na.rm = T)-mean(ex$sex.x, na.rm = T)
  rct_ex_sex_smd_sim[[i]] <- rct_ex_sex_md_sim/rct_ex_sex_sd_sim
  
  
  rct_ex_treat <- data.frame(treat = c(rct_all$treat, ex_all$treat))
  rct_ex_treat_sd_sim <- sd(rct_ex_treat$treat, na.rm = T)
  rct_ex_treat_md_sim <- mean(rct$treat.x, na.rm = T)-mean(rct$treat.x, na.rm = T)
  rct_ex_treat_smd_sim[[i]] <- rct_ex_treat_md_sim/rct_ex_treat_sd_sim
  
  rct_ex_case <- data.frame(case = c(rct_all$case, ex_all$case))
  rct_ex_case_sd_sim <- sd(rct_ex_case$case, na.rm = T)
  rct_ex_case_md_sim <- mean(rct$outcome.x, na.rm = T)-mean(ex$ex_case.x, na.rm = T)
  rct_ex_case_smd_sim[[i]] <- rct_ex_case_md_sim/rct_ex_case_sd_sim
  
  
  
  
  
  rct_ob_age <-data.frame(age = c(rct_all$age, ob_all$age))
  rct_ob_age_sd_sim <- sd(rct_ob_age$age, na.rm = T)
  rct_ob_age_md_sim <- mean(rct$age.x, na.rm = T)-mean(ob$age.x, na.rm = T)
  rct_ob_age_smd_sim[[i]] <- rct_ob_age_md_sim/rct_ob_age_sd_sim
  
  rct_ob_sex <- data.frame(sex = c(rct_all$sex, ob_all$sex))
  rct_ob_sex_sd_sim <- sd(rct_ob_sex$sex, na.rm = T)
  rct_ob_sex_md_sim <- mean(rct$sex.x, na.rm = T)-mean(ob$sex.x, na.rm = T)
  rct_ob_sex_smd_sim[[i]] <- rct_ob_sex_md_sim/rct_ob_sex_sd_sim
  
  rct_ob_treat <- data.frame(treat = c(rct_all$treat, ob_all$treat))
  rct_ob_treat_sd_sim <- sd(rct_ob_treat$treat, na.rm = T)
  rct_ob_treat_md_sim <- mean(rct$treat.x, na.rm = T)-mean(ob$exposed.x, na.rm = T)
  rct_ob_treat_smd_sim[[i]] <- rct_ob_treat_md_sim/rct_ob_treat_sd_sim
  
  
  rct_ob_case <- data.frame(case = c(rct_all$case, ob_all$case))
  rct_ob_case_sd_sim <- sd(rct_ob_case$case, na.rm = T)
  rct_ob_case_md_sim <- mean(rct$outcome.x, na.rm = T)-mean(ob$case.x, na.rm = T)
  rct_ob_case_smd_sim[[i]] <- rct_ob_case_md_sim/rct_ob_case_sd_sim
  
  
  
  ob_ex_age <- data.frame(age = c(ex_all$age, ob_all$age))
  ob_ex_age_sd_sim <- sd(ob_ex_age$age, na.rm = T)
  ob_ex_age_md_sim <- mean(ob$age.x, na.rm = T)-mean(ex$age.x, na.rm = T)
  ob_ex_age_smd_sim[[i]] <- ob_ex_age_md_sim/ob_ex_age_sd_sim
  
  ob_ex_sex <- data.frame(sex = c(ex_all$sex, ob_all$sex))
  ob_ex_sex_sd_sim <- sd(ob_ex_sex$sex, na.rm = T)
  ob_ex_sex_md_sim <- mean(ob$sex.x, na.rm = T)-mean(ex$sex.x, na.rm = T)
  ob_ex_sex_smd_sim[[i]] <- ob_ex_sex_md_sim/ob_ex_sex_sd_sim
  
  ob_ex_treat <- data.frame(treat = c(ex_all$treat, ob_all$treat))
  ob_ex_treat_sd_sim <- sd(ob_ex_treat$treat, na.rm = T)
  ob_ex_treat_md_sim <- mean(ob$exposed.x, na.rm = T)-mean(ex$treatment.x, na.rm = T)
  ob_ex_treat_smd_sim[[i]] <- ob_ex_treat_md_sim/ob_ex_treat_sd_sim
  
  ob_ex_case <- data.frame(case = c(ex_all$case, ob_all$case))
  ob_ex_case_sd_sim <- sd(ob_ex_case$case, na.rm = T)
  ob_ex_case_md_sim <- mean(ob$case.x, na.rm = T)-mean(ex$ex_case.x, na.rm = T)
  ob_ex_case_smd_sim[[i]] <- ob_ex_case_md_sim/ob_ex_case_sd_sim
  
  
  
  
  
}









# ANALYSIS AND GRAPHING




# function to split lists
distribute_into_lists <- function(input_list, sizes) {
  # Check if the sum of sizes is equal to the length of the input list
  if (sum(sizes) != length(input_list)) {
    stop("The sum of sizes must equal the length of the input list")
  }
  
  # Initialize variables
  result <- vector("list", length(sizes))
  start_index <- 1
  
  # Distribute elements based on specified sizes
  for (i in seq_along(sizes)) {
    end_index <- start_index + sizes[i] - 1
    result[[i]] <- input_list[start_index:end_index]
    start_index <- end_index + 1
  }
  
  return(result)
  
}




# lists contain all scenarios combined, in chronological order
# need to splice them into four to display results

list_size <- c(10000, 10000, 10000, 10000)

sim_results_rct <- distribute_into_lists(sim_results_rct, list_size)
sim_results_ob  <- distribute_into_lists(sim_results_ob, list_size)
sim_results_ex <- distribute_into_lists(sim_results_ex, list_size)
syn_results_rct <- distribute_into_lists(syn_results_rct, list_size)
syn_results_ob <- distribute_into_lists(syn_results_ob, list_size)
syn_results_ex <- distribute_into_lists(syn_results_ex, list_size)

chi_results_rct <- distribute_into_lists(chi_results_rct, list_size)
chi_results_ob <- distribute_into_lists(chi_results_ob, list_size)
chi_results_ex <- distribute_into_lists(chi_results_ex, list_size)
chi_results_rct_syn <- distribute_into_lists(chi_results_rct_syn, list_size)
chi_results_ob_syn <- distribute_into_lists(chi_results_ob_syn, list_size)
chi_results_ex_syn <- distribute_into_lists(chi_results_ex_syn, list_size)


table_rct_data <- distribute_into_lists(table_rct_data, list_size)
table_ob_data <- distribute_into_lists(table_ob_data, list_size)
table_ex_data <- distribute_into_lists(table_ex_data, list_size)

rct_age_smd <- distribute_into_lists(rct_age_smd, list_size)
#rct_race_smd <- distribute_into_lists(rct_race_smd, list_size)
rct_sex_smd <- distribute_into_lists(rct_sex_smd, list_size)
rct_treat_smd <- distribute_into_lists(rct_treat_smd, list_size)
rct_outcome_smd <- distribute_into_lists(rct_outcome_smd, list_size)

ob_age_smd <- distribute_into_lists(ob_age_smd, list_size)
#ob_race_smd <- distribute_into_lists(ob_race_smd, list_size)
ob_sex_smd <- distribute_into_lists(ob_sex_smd, list_size)
ob_exposed_smd <- distribute_into_lists(ob_exposed_smd, list_size)
ob_case_smd <- distribute_into_lists(ob_case_smd, list_size)

ex_age_smd  <- distribute_into_lists(ex_age_smd, list_size)
#ex_race_smd <- distribute_into_lists(ex_race_smd, list_size)
ex_sex_smd <- distribute_into_lists(ex_sex_smd, list_size)
ex_treatment_smd <- distribute_into_lists(ex_treatment_smd, list_size)
ex_case_smd <- distribute_into_lists(ex_case_smd, list_size)






rct_ex_age_smd <- distribute_into_lists(rct_ex_age_smd, list_size)
rct_ex_sex_smd <- distribute_into_lists(rct_ex_sex_smd, list_size)
rct_ex_treat_smd <- distribute_into_lists(rct_ex_treat_smd, list_size)
rct_ex_case_smd <- distribute_into_lists(rct_ex_case_smd, list_size)

rct_ob_age_smd <- distribute_into_lists(rct_ob_age_smd, list_size)
rct_ob_sex_smd <- distribute_into_lists(rct_ob_sex_smd, list_size)
rct_ob_treat_smd <- distribute_into_lists(rct_ob_treat_smd, list_size)
rct_ob_case_smd <- distribute_into_lists(rct_ob_case_smd, list_size)

ob_ex_age_smd  <- distribute_into_lists(ob_ex_age_smd, list_size)
ob_ex_sex_smd <- distribute_into_lists(ob_ex_sex_smd, list_size)
ob_ex_treat_smd <- distribute_into_lists(ob_ex_treat_smd, list_size)
ob_ex_case_smd <- distribute_into_lists(ob_ex_case_smd, list_size)



rct_ex_age_smd_sim <- distribute_into_lists(rct_ex_age_smd_sim, list_size)
rct_ex_sex_smd_sim <- distribute_into_lists(rct_ex_sex_smd_sim, list_size)
rct_ex_treat_smd_sim <- distribute_into_lists(rct_ex_treat_smd_sim, list_size)
rct_ex_case_smd_sim <- distribute_into_lists(rct_ex_case_smd_sim, list_size)

rct_ob_age_smd_sim <- distribute_into_lists(rct_ob_age_smd_sim, list_size)
rct_ob_sex_smd_sim <- distribute_into_lists(rct_ob_sex_smd_sim, list_size)
rct_ob_treat_smd_sim <- distribute_into_lists(rct_ob_treat_smd_sim, list_size)
rct_ob_case_smd_sim <- distribute_into_lists(rct_ob_case_smd_sim, list_size)

ob_ex_age_smd_sim  <- distribute_into_lists(ob_ex_age_smd_sim, list_size)
ob_ex_sex_smd_sim <- distribute_into_lists(ob_ex_sex_smd_sim, list_size)
ob_ex_treat_smd_sim <- distribute_into_lists(ob_ex_treat_smd_sim, list_size)
ob_ex_case_smd_sim <- distribute_into_lists(ob_ex_case_smd_sim, list_size)








# check p values of chi square tests

chi_rct_p <- (sapply(chi_results_rct[[1]],"[[",3)) 
chi_ob_p <- (sapply(chi_results_ob[[1]],"[[",3)) 
chi_ex_p <- (sapply(chi_results_ex[[1]],"[[",3)) 
chi_rct_p_syn <- (sapply(chi_results_rct_syn[[1]],"[[",3)) 
chi_ob_p_syn <- (sapply(chi_results_ob_syn[[1]],"[[",3)) 
chi_ex_p_syn <- (sapply(chi_results_ex_syn[[1]],"[[",3)) 


counter <- 0
for (i in na.omit(chi_rct_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ob_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ex_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)

counter <- 0
for (i in na.omit(chi_rct_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ob_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ex_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)





chi_rct_p <- (sapply(chi_results_rct[[2]],"[[",3)) 
chi_ob_p <- (sapply(chi_results_ob[[2]],"[[",3)) 
chi_ex_p <- (sapply(chi_results_ex[[2]],"[[",3)) 
chi_rct_p_syn <- (sapply(chi_results_rct_syn[[2]],"[[",3)) 
chi_ob_p_syn <- (sapply(chi_results_ob_syn[[2]],"[[",3)) 
chi_ex_p_syn <- (sapply(chi_results_ex_syn[[2]],"[[",3)) 


counter <- 0
for (i in na.omit(chi_rct_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ob_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ex_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)

counter <- 0
for (i in na.omit(chi_rct_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ob_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ex_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)


chi_rct_p <- (sapply(chi_results_rct[[3]],"[[",3)) 
chi_ob_p <- (sapply(chi_results_ob[[3]],"[[",3)) 
chi_ex_p <- (sapply(chi_results_ex[[3]],"[[",3)) 
chi_rct_p_syn <- (sapply(chi_results_rct_syn[[3]],"[[",3)) 
chi_ob_p_syn <- (sapply(chi_results_ob_syn[[3]],"[[",3)) 
chi_ex_p_syn <- (sapply(chi_results_ex_syn[[3]],"[[",3)) 


counter <- 0
for (i in na.omit(chi_rct_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ob_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ex_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)

counter <- 0
for (i in na.omit(chi_rct_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ob_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ex_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)


chi_rct_p <- (sapply(chi_results_rct[[4]],"[[",3)) 
chi_ob_p <- (sapply(chi_results_ob[[4]],"[[",3)) 
chi_ex_p <- (sapply(chi_results_ex[[4]],"[[",3)) 
chi_rct_p_syn <- (sapply(chi_results_rct_syn[[4]],"[[",3)) 
chi_ob_p_syn <- (sapply(chi_results_ob_syn[[4]],"[[",3)) 
chi_ex_p_syn <- (sapply(chi_results_ex_syn[[4]],"[[",3)) 


counter <- 0
for (i in na.omit(chi_rct_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ob_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ex_p)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)

counter <- 0
for (i in na.omit(chi_rct_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ob_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)
counter <- 0
for (i in na.omit(chi_ex_p_syn)){
  if(i > 0.05){
    print(i)
    counter <- counter+1
  } 
  else if (i < 0.05){
    counter <- counter+0
  }
}
print(counter)









# smd plot



# calculate mean smd for each variable - for scenario 1
rct_age_smd_av <- mean(sapply(rct_age_smd[[1]], mean))
#rct_race_smd_av <- mean(sapply(rct_race_smd[[1]], mean))
rct_sex_smd_av <- mean(sapply(rct_sex_smd[[1]], mean))
rct_treat_smd_av <- mean(sapply(rct_treat_smd[[1]], mean))
rct_outcome_smd_av <- mean(sapply(rct_outcome_smd[[1]], mean))

ob_age_smd_av <- mean(sapply(ob_age_smd[[1]], mean))
#ob_race_smd_av <- mean(sapply(ob_race_smd[[1]], mean))
ob_sex_smd_av <- mean(sapply(ob_sex_smd[[1]], mean))
ob_exposed_smd_av <- mean(sapply(ob_exposed_smd[[1]], mean))
ob_case_smd_av <- mean(sapply(ob_case_smd[[1]], mean))


ex_age_smd_av <- mean(sapply(ex_age_smd[[1]], mean))
#ex_race_smd_av <- mean(sapply(ex_race_smd[[1]], mean))
ex_sex_smd_av <- mean(sapply(ex_sex_smd[[1]], mean))
ex_treatment_smd_av <- mean(sapply(ex_treatment_smd[[1]], mean))
ex_case_smd_av <- mean(sapply(ex_case_smd[[1]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age", "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT","RCT","RCT","RCT","Observational",
           "Observational","Observational","Observational",
           "Exernal","Exernal","Exernal","Exernal"),
  SMD= c(rct_age_smd_av,  rct_sex_smd_av, rct_treat_smd_av, rct_outcome_smd_av, 
         ob_age_smd_av,  ob_sex_smd_av, ob_exposed_smd_av, ob_case_smd_av,
         ex_age_smd_av,  ex_sex_smd_av, ex_treatment_smd_av, ex_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen1_plot_ll_ub.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 1
rct_age_smd_sap <- sapply(rct_age_smd[[1]], mean)
#rct_race_smd_sap <- sapply(rct_race_smd[[1]], mean)
rct_sex_smd_sap <-sapply(rct_sex_smd[[1]], mean)
rct_treat_smd_sap <- sapply(rct_treat_smd[[1]], mean)
rct_outcome_smd_sap <- sapply(rct_outcome_smd[[1]], mean)

ob_age_smd_sap <- sapply(ob_age_smd[[1]], mean)
#ob_race_smd_sap <- sapply(ob_race_smd[[1]], mean)
ob_sex_smd_sap <- sapply(ob_sex_smd[[1]], mean)
ob_exposed_smd_sap <- sapply(ob_exposed_smd[[1]], mean)
ob_case_smd_sap <- sapply(ob_case_smd[[1]], mean)


ex_age_smd_sap <- sapply(ex_age_smd[[1]], mean)
#ex_race_smd_sap <- sapply(ex_race_smd[[1]], mean)
ex_sex_smd_sap <- sapply(ex_sex_smd[[1]], mean)
ex_treatment_smd_sap <- sapply(ex_treatment_smd[[1]], mean)
ex_case_smd_sap <- sapply(ex_case_smd[[1]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_age_smd_sap,  rct_sex_smd_sap, rct_treat_smd_sap, rct_outcome_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age", "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT","RCT","RCT","RCT","Observational",
           "Observational","Observational","Observational",
           "Exernal","Exernal","Exernal","Exernal"),
  SMD= c(rct_age_smd_sap, rct_sex_smd_sap, rct_treat_smd_sap, rct_outcome_smd_sap, 
         ob_age_smd_sap,  ob_sex_smd_sap, ob_exposed_smd_sap, ob_case_smd_sap,
         ex_age_smd_sap,  ex_sex_smd_sap, ex_treatment_smd_sap, ex_case_smd_sap)
  
)



# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen1_box_ll_ub.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 1)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen1_box_gray_ll_ub.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1) +
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 1)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()









# calculate mean smd for each variable - for scenario 2
rct_age_smd_av <- mean(sapply(rct_age_smd[[2]], mean))
#rct_race_smd_av <- mean(sapply(rct_race_smd[[2]], mean))
rct_sex_smd_av <- mean(sapply(rct_sex_smd[[2]], mean))
rct_treat_smd_av <- mean(sapply(rct_treat_smd[[2]], mean))
rct_outcome_smd_av <- mean(sapply(rct_outcome_smd[[2]], mean))

ob_age_smd_av <- mean(sapply(ob_age_smd[[2]], mean))
#ob_race_smd_av <- mean(sapply(ob_race_smd[[2]], mean))
ob_sex_smd_av <- mean(sapply(ob_sex_smd[[2]], mean))
ob_exposed_smd_av <- mean(sapply(ob_exposed_smd[[2]], mean))
ob_case_smd_av <- mean(sapply(ob_case_smd[[2]], mean))


ex_age_smd_av <- mean(sapply(ex_age_smd[[2]], mean))
#ex_race_smd_av <- mean(sapply(ex_race_smd[[2]], mean))
ex_sex_smd_av <- mean(sapply(ex_sex_smd[[2]], mean))
ex_treatment_smd_av <- mean(sapply(ex_treatment_smd[[2]], mean))
ex_case_smd_av <- mean(sapply(ex_case_smd[[2]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age", "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT","RCT","RCT","RCT","Observational",
           "Observational","Observational","Observational",
           "Exernal","Exernal","Exernal","Exernal"),
  SMD= c(rct_age_smd_av, rct_sex_smd_av, rct_treat_smd_av, rct_outcome_smd_av, 
         ob_age_smd_av,  ob_sex_smd_av, ob_exposed_smd_av, ob_case_smd_av,
         ex_age_smd_av,  ex_sex_smd_av, ex_treatment_smd_av, ex_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen2_plot_ll_ub.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 2
rct_age_smd_sap <- sapply(rct_age_smd[[2]], mean)
#rct_race_smd_sap <- sapply(rct_race_smd[[2]], mean)
rct_sex_smd_sap <-sapply(rct_sex_smd[[2]], mean)
rct_treat_smd_sap <- sapply(rct_treat_smd[[2]], mean)
rct_outcome_smd_sap <- sapply(rct_outcome_smd[[2]], mean)

ob_age_smd_sap <- sapply(ob_age_smd[[2]], mean)
#ob_race_smd_sap <- sapply(ob_race_smd[[2]], mean)
ob_sex_smd_sap <- sapply(ob_sex_smd[[2]], mean)
ob_exposed_smd_sap <- sapply(ob_exposed_smd[[2]], mean)
ob_case_smd_sap <- sapply(ob_case_smd[[2]], mean)


ex_age_smd_sap <- sapply(ex_age_smd[[2]], mean)
#ex_race_smd_sap <- sapply(ex_race_smd[[2]], mean)
ex_sex_smd_sap <- sapply(ex_sex_smd[[2]], mean)
ex_treatment_smd_sap <- sapply(ex_treatment_smd[[2]], mean)
ex_case_smd_sap <- sapply(ex_case_smd[[2]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_age_smd_sap,  rct_sex_smd_sap, rct_treat_smd_sap, rct_outcome_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age", "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT","RCT","RCT","RCT","Observational",
           "Observational","Observational","Observational",
           "Exernal","Exernal","Exernal","Exernal"),
  SMD= c(rct_age_smd_sap, rct_sex_smd_sap, rct_treat_smd_sap, rct_outcome_smd_sap, 
         ob_age_smd_sap,  ob_sex_smd_sap, ob_exposed_smd_sap, ob_case_smd_sap,
         ex_age_smd_sap,  ex_sex_smd_sap, ex_treatment_smd_sap, ex_case_smd_sap)
  
)



# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen2_box_ll_ub.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() +ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 2)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen2_box_gray_ll_ub.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1) +
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 2)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()









# calculate mean smd for each variable - for scenario 3
rct_age_smd_av <- mean(sapply(rct_age_smd[[3]], mean))
#rct_race_smd_av <- mean(sapply(rct_race_smd[[3]], mean))
rct_sex_smd_av <- mean(sapply(rct_sex_smd[[3]], mean))
rct_treat_smd_av <- mean(sapply(rct_treat_smd[[3]], mean))
rct_outcome_smd_av <- mean(sapply(rct_outcome_smd[[3]], mean))

ob_age_smd_av <- mean(sapply(ob_age_smd[[3]], mean))
#ob_race_smd_av <- mean(sapply(ob_race_smd[[3]], mean))
ob_sex_smd_av <- mean(sapply(ob_sex_smd[[3]], mean))
ob_exposed_smd_av <- mean(sapply(ob_exposed_smd[[3]], mean))
ob_case_smd_av <- mean(sapply(ob_case_smd[[3]], mean))


ex_age_smd_av <- mean(sapply(ex_age_smd[[3]], mean))
#ex_race_smd_av <- mean(sapply(ex_race_smd[[3]], mean))
ex_sex_smd_av <- mean(sapply(ex_sex_smd[[3]], mean))
ex_treatment_smd_av <- mean(sapply(ex_treatment_smd[[3]], mean))
ex_case_smd_av <- mean(sapply(ex_case_smd[[3]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT","RCT","RCT","RCT","Observational",
           "Observational","Observational","Observational",
           "Exernal","Exernal","Exernal","Exernal"),
  SMD= c(rct_age_smd_av,rct_sex_smd_av, rct_treat_smd_av, rct_outcome_smd_av, 
         ob_age_smd_av,  ob_sex_smd_av, ob_exposed_smd_av, ob_case_smd_av,
         ex_age_smd_av,  ex_sex_smd_av, ex_treatment_smd_av, ex_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen3_plot_ll_ub.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 3
rct_age_smd_sap <- sapply(rct_age_smd[[3]], mean)
#rct_race_smd_sap <- sapply(rct_race_smd[[3]], mean)
rct_sex_smd_sap <-sapply(rct_sex_smd[[3]], mean)
rct_treat_smd_sap <- sapply(rct_treat_smd[[3]], mean)
rct_outcome_smd_sap <- sapply(rct_outcome_smd[[3]], mean)

ob_age_smd_sap <- sapply(ob_age_smd[[3]], mean)
#ob_race_smd_sap <- sapply(ob_race_smd[[3]], mean)
ob_sex_smd_sap <- sapply(ob_sex_smd[[3]], mean)
ob_exposed_smd_sap <- sapply(ob_exposed_smd[[3]], mean)
ob_case_smd_sap <- sapply(ob_case_smd[[3]], mean)


ex_age_smd_sap <- sapply(ex_age_smd[[3]], mean)
#ex_race_smd_sap <- sapply(ex_race_smd[[3]], mean)
ex_sex_smd_sap <- sapply(ex_sex_smd[[3]], mean)
ex_treatment_smd_sap <- sapply(ex_treatment_smd[[3]], mean)
ex_case_smd_sap <- sapply(ex_case_smd[[3]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_age_smd_sap,  rct_sex_smd_sap, rct_treat_smd_sap, rct_outcome_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT","RCT","RCT","RCT","Observational",
           "Observational","Observational","Observational",
           "Exernal","Exernal","Exernal","Exernal"),
  SMD= c(rct_age_smd_sap,  rct_sex_smd_sap, rct_treat_smd_sap, rct_outcome_smd_sap, 
         ob_age_smd_sap,  ob_sex_smd_sap, ob_exposed_smd_sap, ob_case_smd_sap,
         ex_age_smd_sap,  ex_sex_smd_sap, ex_treatment_smd_sap, ex_case_smd_sap)
  
)



# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen3_box_ll_ub.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 3)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen3_box_gray_ll_ub.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1) +
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 3)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()













# calculate mean smd for each variable - for scenario 4
rct_age_smd_av <- mean(sapply(rct_age_smd[[4]], mean))
#rct_race_smd_av <- mean(sapply(rct_race_smd[[4]], mean))
rct_sex_smd_av <- mean(sapply(rct_sex_smd[[4]], mean))
rct_treat_smd_av <- mean(sapply(rct_treat_smd[[4]], mean))
rct_outcome_smd_av <- mean(sapply(rct_outcome_smd[[4]], mean))

ob_age_smd_av <- mean(sapply(ob_age_smd[[4]], mean))
#ob_race_smd_av <- mean(sapply(ob_race_smd[[4]], mean))
ob_sex_smd_av <- mean(sapply(ob_sex_smd[[4]], mean))
ob_exposed_smd_av <- mean(sapply(ob_exposed_smd[[4]], mean))
ob_case_smd_av <- mean(sapply(ob_case_smd[[4]], mean))


ex_age_smd_av <- mean(sapply(ex_age_smd[[4]], mean))
#ex_race_smd_av <- mean(sapply(ex_race_smd[[4]], mean))
ex_sex_smd_av <- mean(sapply(ex_sex_smd[[4]], mean))
ex_treatment_smd_av <- mean(sapply(ex_treatment_smd[[4]], mean))
ex_case_smd_av <- mean(sapply(ex_case_smd[[4]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age", "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT","RCT","RCT","RCT","Observational",
           "Observational","Observational","Observational",
           "Exernal","Exernal","Exernal","Exernal"),
  SMD= c(rct_age_smd_av,  rct_sex_smd_av, rct_treat_smd_av, rct_outcome_smd_av, 
         ob_age_smd_av,  ob_sex_smd_av, ob_exposed_smd_av, ob_case_smd_av,
         ex_age_smd_av,  ex_sex_smd_av, ex_treatment_smd_av, ex_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen4_plot_ll_ub.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 4
rct_age_smd_sap <- sapply(rct_age_smd[[4]], mean)
#rct_race_smd_sap <- sapply(rct_race_smd[[4]], mean)
rct_sex_smd_sap <-sapply(rct_sex_smd[[4]], mean)
rct_treat_smd_sap <- sapply(rct_treat_smd[[4]], mean)
rct_outcome_smd_sap <- sapply(rct_outcome_smd[[4]], mean)

ob_age_smd_sap <- sapply(ob_age_smd[[4]], mean)
#ob_race_smd_sap <- sapply(ob_race_smd[[4]], mean)
ob_sex_smd_sap <- sapply(ob_sex_smd[[4]], mean)
ob_exposed_smd_sap <- sapply(ob_exposed_smd[[4]], mean)
ob_case_smd_sap <- sapply(ob_case_smd[[4]], mean)


ex_age_smd_sap <- sapply(ex_age_smd[[4]], mean)
#ex_race_smd_sap <- sapply(ex_race_smd[[4]], mean)
ex_sex_smd_sap <- sapply(ex_sex_smd[[4]], mean)
ex_treatment_smd_sap <- sapply(ex_treatment_smd[[4]], mean)
ex_case_smd_sap <- sapply(ex_case_smd[[4]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_age_smd_sap, rct_sex_smd_sap, rct_treat_smd_sap, rct_outcome_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT","RCT","RCT","RCT","Observational",
           "Observational","Observational","Observational",
           "Exernal","Exernal","Exernal","Exernal"),
  SMD= c(rct_age_smd_sap,  rct_sex_smd_sap, rct_treat_smd_sap, rct_outcome_smd_sap, 
         ob_age_smd_sap,  ob_sex_smd_sap, ob_exposed_smd_sap, ob_case_smd_sap,
         ex_age_smd_sap,  ex_sex_smd_sap, ex_treatment_smd_sap, ex_case_smd_sap)
  
)



# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen4_box_ll_ub.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + 
  ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen4_box_gray_ll_ub.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1)+
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()


















# graph to show comparison across synthetic data sets


# take the standard mean difference of the synthetic data compared, before and after. show that on graph instead.

# NEED TO RECALC NEW SMD INSIDE LOOPS AND RERUN !!!







# calculate mean smd for each variable - for scenario 1, SYNTH COMPARISON
rct_ex_age_smd_av <- mean(sapply(rct_ex_age_smd[[1]], mean))

rct_ex_sex_smd_av <- mean(sapply(rct_ex_sex_smd[[1]], mean))
rct_ex_treat_smd_av <- mean(sapply(rct_ex_treat_smd[[1]], mean))
rct_ex_case_smd_av <- mean(sapply(rct_ex_case_smd[[1]], mean))

ob_ex_age_smd_av <- mean(sapply(ob_ex_age_smd[[1]], mean))

ob_ex_sex_smd_av <- mean(sapply(ob_ex_sex_smd[[1]], mean))
ob_ex_treat_smd_av <- mean(sapply(ob_ex_treat_smd[[1]], mean))
ob_ex_case_smd_av <- mean(sapply(ob_ex_case_smd[[1]], mean))


rct_ob_age_smd_av <- mean(sapply(rct_ob_age_smd[[1]], mean))

rct_ob_sex_smd_av <- mean(sapply(rct_ob_sex_smd[[1]], mean))
rct_ob_treat_smd_av <- mean(sapply(rct_ob_treat_smd[[1]], mean))
rct_ob_case_smd_av <- mean(sapply(rct_ob_case_smd[[1]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT/External","RCT/External","RCT/External","RCT/External","RCT/Observational",
           "RCT/Observational","RCT/Observational","RCT/Observational",
           "Observational/Exernal","Observational/Exernal","Observational/Exernal","Observational/Exernal"),
  SMD= c(rct_ex_age_smd_av,  rct_ex_sex_smd_av, rct_ex_treat_smd_av, rct_ex_case_smd_av, 
         ob_ex_age_smd_av, ob_ex_sex_smd_av, ob_ex_treat_smd_av, ob_ex_case_smd_av,
         rct_ob_age_smd_av,  rct_ob_sex_smd_av, rct_ob_treat_smd_av, rct_ob_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen1_plot_ub_ll_syn.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 4
rct_ex_age_smd_sap <- sapply(rct_ex_age_smd[[1]], mean)

rct_ex_sex_smd_sap <-sapply(rct_ex_sex_smd[[1]], mean)
rct_ex_treat_smd_sap <- sapply(rct_ex_treat_smd[[1]], mean)
rct_ex_case_smd_sap <- sapply(rct_ex_case_smd[[1]], mean)

rct_ob_age_smd_sap <- sapply(rct_ob_age_smd[[1]], mean)

rct_ob_sex_smd_sap <- sapply(rct_ob_sex_smd[[1]], mean)
rct_ob_treat_smd_sap <- sapply(rct_ob_treat_smd[[1]], mean)
rct_ob_case_smd_sap <- sapply(rct_ob_case_smd[[1]], mean)


ob_ex_age_smd_sap <- sapply(ob_ex_age_smd[[1]], mean)

ob_ex_sex_smd_sap <- sapply(ob_ex_sex_smd[[1]], mean)
ob_ex_treat_smd_sap <- sapply(ob_ex_treat_smd[[1]], mean)
ob_ex_case_smd_sap <- sapply(ob_ex_case_smd[[1]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_ex_age_smd_sap,  rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT/Ex","RCT/Ex","RCT/Ex","RCT/Ex",
           "Obs/Ex","Obs/Ex","Obs/Ex","Obs/Ex",
           "RCT/Obs","RCT/Obs","RCT/Obs","RCT/Obs"),
  SMD= c(rct_ex_age_smd_sap, rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap, 
         ob_ex_age_smd_sap, ob_ex_sex_smd_sap, ob_ex_treat_smd_sap, ob_ex_case_smd_sap,
         rct_ob_age_smd_sap, rct_ob_sex_smd_sap, rct_ob_treat_smd_sap, rct_ob_case_smd_sap)
  
)



# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen1_box_ub_ll_syn.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + 
  ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen1_box_gray_ub_ll_syn.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1)+
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()












# test and repeat for other scenarios



# calculate mean smd for each variable - for scenario 1, SYNTH COMPARISON
rct_ex_age_smd_av <- mean(sapply(rct_ex_age_smd[[2]], mean))

rct_ex_sex_smd_av <- mean(sapply(rct_ex_sex_smd[[2]], mean))
rct_ex_treat_smd_av <- mean(sapply(rct_ex_treat_smd[[2]], mean))
rct_ex_case_smd_av <- mean(sapply(rct_ex_case_smd[[2]], mean))

ob_ex_age_smd_av <- mean(sapply(ob_ex_age_smd[[2]], mean))

ob_ex_sex_smd_av <- mean(sapply(ob_ex_sex_smd[[2]], mean))
ob_ex_treat_smd_av <- mean(sapply(ob_ex_treat_smd[[2]], mean))
ob_ex_case_smd_av <- mean(sapply(ob_ex_case_smd[[2]], mean))


rct_ob_age_smd_av <- mean(sapply(rct_ob_age_smd[[2]], mean))

rct_ob_sex_smd_av <- mean(sapply(rct_ob_sex_smd[[2]], mean))
rct_ob_treat_smd_av <- mean(sapply(rct_ob_treat_smd[[2]], mean))
rct_ob_case_smd_av <- mean(sapply(rct_ob_case_smd[[2]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT/External","RCT/External","RCT/External","RCT/External","RCT/Observational",
           "RCT/Observational","RCT/Observational","RCT/Observational",
           "Observational/Exernal","Observational/Exernal","Observational/Exernal","Observational/Exernal"),
  SMD= c(rct_ex_age_smd_av,  rct_ex_sex_smd_av, rct_ex_treat_smd_av, rct_ex_case_smd_av, 
         ob_ex_age_smd_av, ob_ex_sex_smd_av, ob_ex_treat_smd_av, ob_ex_case_smd_av,
         rct_ob_age_smd_av,  rct_ob_sex_smd_av, rct_ob_treat_smd_av, rct_ob_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen2_plot_ub_ll_syn.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 4
rct_ex_age_smd_sap <- sapply(rct_ex_age_smd[[2]], mean)

rct_ex_sex_smd_sap <-sapply(rct_ex_sex_smd[[2]], mean)
rct_ex_treat_smd_sap <- sapply(rct_ex_treat_smd[[2]], mean)
rct_ex_case_smd_sap <- sapply(rct_ex_case_smd[[2]], mean)

rct_ob_age_smd_sap <- sapply(rct_ob_age_smd[[2]], mean)

rct_ob_sex_smd_sap <- sapply(rct_ob_sex_smd[[2]], mean)
rct_ob_treat_smd_sap <- sapply(rct_ob_treat_smd[[2]], mean)
rct_ob_case_smd_sap <- sapply(rct_ob_case_smd[[2]], mean)


ob_ex_age_smd_sap <- sapply(ob_ex_age_smd[[2]], mean)

ob_ex_sex_smd_sap <- sapply(ob_ex_sex_smd[[2]], mean)
ob_ex_treat_smd_sap <- sapply(ob_ex_treat_smd[[2]], mean)
ob_ex_case_smd_sap <- sapply(ob_ex_case_smd[[2]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_ex_age_smd_sap,  rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT/Ex","RCT/Ex","RCT/Ex","RCT/Ex",
           "Obs/Ex","Obs/Ex","Obs/Ex","Obs/Ex",
           "RCT/Obs","RCT/Obs","RCT/Obs","RCT/Obs"),
  SMD= c(rct_ex_age_smd_sap, rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap, 
         ob_ex_age_smd_sap, ob_ex_sex_smd_sap, ob_ex_treat_smd_sap, ob_ex_case_smd_sap,
         rct_ob_age_smd_sap, rct_ob_sex_smd_sap, rct_ob_treat_smd_sap, rct_ob_case_smd_sap)
  
)



# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen2_box_ub_ll_syn.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + 
  ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen2_box_gray_ub_ll_syn.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1)+
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()













# calculate mean smd for each variable - for scenario 1, SYNTH COMPARISON
rct_ex_age_smd_av <- mean(sapply(rct_ex_age_smd[[3]], mean))

rct_ex_sex_smd_av <- mean(sapply(rct_ex_sex_smd[[3]], mean))
rct_ex_treat_smd_av <- mean(sapply(rct_ex_treat_smd[[3]], mean))
rct_ex_case_smd_av <- mean(sapply(rct_ex_case_smd[[3]], mean))

ob_ex_age_smd_av <- mean(sapply(ob_ex_age_smd[[3]], mean))

ob_ex_sex_smd_av <- mean(sapply(ob_ex_sex_smd[[3]], mean))
ob_ex_treat_smd_av <- mean(sapply(ob_ex_treat_smd[[3]], mean))
ob_ex_case_smd_av <- mean(sapply(ob_ex_case_smd[[3]], mean))


rct_ob_age_smd_av <- mean(sapply(rct_ob_age_smd[[3]], mean))

rct_ob_sex_smd_av <- mean(sapply(rct_ob_sex_smd[[3]], mean))
rct_ob_treat_smd_av <- mean(sapply(rct_ob_treat_smd[[3]], mean))
rct_ob_case_smd_av <- mean(sapply(rct_ob_case_smd[[3]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT/External","RCT/External","RCT/External","RCT/External","RCT/Observational",
           "RCT/Observational","RCT/Observational","RCT/Observational",
           "Observational/Exernal","Observational/Exernal","Observational/Exernal","Observational/Exernal"),
  SMD= c(rct_ex_age_smd_av,  rct_ex_sex_smd_av, rct_ex_treat_smd_av, rct_ex_case_smd_av, 
         ob_ex_age_smd_av, ob_ex_sex_smd_av, ob_ex_treat_smd_av, ob_ex_case_smd_av,
         rct_ob_age_smd_av,  rct_ob_sex_smd_av, rct_ob_treat_smd_av, rct_ob_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen3_plot_ub_ll_syn.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 4
rct_ex_age_smd_sap <- sapply(rct_ex_age_smd[[3]], mean)

rct_ex_sex_smd_sap <-sapply(rct_ex_sex_smd[[3]], mean)
rct_ex_treat_smd_sap <- sapply(rct_ex_treat_smd[[3]], mean)
rct_ex_case_smd_sap <- sapply(rct_ex_case_smd[[3]], mean)

rct_ob_age_smd_sap <- sapply(rct_ob_age_smd[[3]], mean)

rct_ob_sex_smd_sap <- sapply(rct_ob_sex_smd[[3]], mean)
rct_ob_treat_smd_sap <- sapply(rct_ob_treat_smd[[3]], mean)
rct_ob_case_smd_sap <- sapply(rct_ob_case_smd[[3]], mean)


ob_ex_age_smd_sap <- sapply(ob_ex_age_smd[[3]], mean)

ob_ex_sex_smd_sap <- sapply(ob_ex_sex_smd[[3]], mean)
ob_ex_treat_smd_sap <- sapply(ob_ex_treat_smd[[3]], mean)
ob_ex_case_smd_sap <- sapply(ob_ex_case_smd[[3]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_ex_age_smd_sap,  rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap )
)
graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT/Ex","RCT/Ex","RCT/Ex","RCT/Ex",
           "Obs/Ex","Obs/Ex","Obs/Ex","Obs/Ex",
           "RCT/Obs","RCT/Obs","RCT/Obs","RCT/Obs"),
  SMD= c(rct_ex_age_smd_sap, rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap, 
         ob_ex_age_smd_sap, ob_ex_sex_smd_sap, ob_ex_treat_smd_sap, ob_ex_case_smd_sap,
         rct_ob_age_smd_sap, rct_ob_sex_smd_sap, rct_ob_treat_smd_sap, rct_ob_case_smd_sap)
  
)




# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen3_box_ub_ll_syn.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + 
  ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen3_box_gray_ub_ll_syn.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1)+
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()













# calculate mean smd for each variable - for scenario 1, SYNTH COMPARISON
rct_ex_age_smd_av <- mean(sapply(rct_ex_age_smd[[4]], mean))

rct_ex_sex_smd_av <- mean(sapply(rct_ex_sex_smd[[4]], mean))
rct_ex_treat_smd_av <- mean(sapply(rct_ex_treat_smd[[4]], mean))
rct_ex_case_smd_av <- mean(sapply(rct_ex_case_smd[[4]], mean))

ob_ex_age_smd_av <- mean(sapply(ob_ex_age_smd[[4]], mean))

ob_ex_sex_smd_av <- mean(sapply(ob_ex_sex_smd[[4]], mean))
ob_ex_treat_smd_av <- mean(sapply(ob_ex_treat_smd[[4]], mean))
ob_ex_case_smd_av <- mean(sapply(ob_ex_case_smd[[4]], mean))


rct_ob_age_smd_av <- mean(sapply(rct_ob_age_smd[[4]], mean))

rct_ob_sex_smd_av <- mean(sapply(rct_ob_sex_smd[[4]], mean))
rct_ob_treat_smd_av <- mean(sapply(rct_ob_treat_smd[[4]], mean))
rct_ob_case_smd_av <- mean(sapply(rct_ob_case_smd[[4]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT/External","RCT/External","RCT/External","RCT/External","RCT/Observational",
           "RCT/Observational","RCT/Observational","RCT/Observational",
           "Observational/Exernal","Observational/Exernal","Observational/Exernal","Observational/Exernal"),
  SMD= c(rct_ex_age_smd_av,  rct_ex_sex_smd_av, rct_ex_treat_smd_av, rct_ex_case_smd_av, 
         ob_ex_age_smd_av, ob_ex_sex_smd_av, ob_ex_treat_smd_av, ob_ex_case_smd_av,
         rct_ob_age_smd_av,  rct_ob_sex_smd_av, rct_ob_treat_smd_av, rct_ob_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen4_plot_ub_ll_syn.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 4
rct_ex_age_smd_sap <- sapply(rct_ex_age_smd[[4]], mean)

rct_ex_sex_smd_sap <-sapply(rct_ex_sex_smd[[4]], mean)
rct_ex_treat_smd_sap <- sapply(rct_ex_treat_smd[[4]], mean)
rct_ex_case_smd_sap <- sapply(rct_ex_case_smd[[4]], mean)

rct_ob_age_smd_sap <- sapply(rct_ob_age_smd[[4]], mean)

rct_ob_sex_smd_sap <- sapply(rct_ob_sex_smd[[4]], mean)
rct_ob_treat_smd_sap <- sapply(rct_ob_treat_smd[[4]], mean)
rct_ob_case_smd_sap <- sapply(rct_ob_case_smd[[4]], mean)


ob_ex_age_smd_sap <- sapply(ob_ex_age_smd[[4]], mean)

ob_ex_sex_smd_sap <- sapply(ob_ex_sex_smd[[4]], mean)
ob_ex_treat_smd_sap <- sapply(ob_ex_treat_smd[[4]], mean)
ob_ex_case_smd_sap <- sapply(ob_ex_case_smd[[4]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_ex_age_smd_sap,  rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT/Ex","RCT/Ex","RCT/Ex","RCT/Ex",
           "Obs/Ex","Obs/Ex","Obs/Ex","Obs/Ex",
           "RCT/Obs","RCT/Obs","RCT/Obs","RCT/Obs"),
  SMD= c(rct_ex_age_smd_sap, rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap, 
         ob_ex_age_smd_sap, ob_ex_sex_smd_sap, ob_ex_treat_smd_sap, ob_ex_case_smd_sap,
         rct_ob_age_smd_sap, rct_ob_sex_smd_sap, rct_ob_treat_smd_sap, rct_ob_case_smd_sap)
  
)



# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen4_box_ub_ll_syn.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + 
  ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen4_box_gray_ub_ll_syn.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1)+
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()

























# REP for sim





# graph to show comparison across sim data sets


# take the standard mean difference of the synthetic data compared, before and after. show that on graph instead.

# NEED TO RECALC NEW SMD INSIDE LOOPS AND RERUN !!!







# calculate mean smd for each variable - for scenario 1, SYNTH COMPARISON
rct_ex_age_smd_av <- mean(sapply(rct_ex_age_smd_sim[[1]], mean))

rct_ex_sex_smd_av <- mean(sapply(rct_ex_sex_smd_sim[[1]], mean))
rct_ex_treat_smd_av <- mean(sapply(rct_ex_treat_smd_sim[[1]], mean))
rct_ex_case_smd_av <- mean(sapply(rct_ex_case_smd_sim[[1]], mean))

ob_ex_age_smd_av <- mean(sapply(ob_ex_age_smd_sim[[1]], mean))

ob_ex_sex_smd_av <- mean(sapply(ob_ex_sex_smd_sim[[1]], mean))
ob_ex_treat_smd_av <- mean(sapply(ob_ex_treat_smd_sim[[1]], mean))
ob_ex_case_smd_av <- mean(sapply(ob_ex_case_smd_sim[[1]], mean))


rct_ob_age_smd_av <- mean(sapply(rct_ob_age_smd_sim[[1]], mean))

rct_ob_sex_smd_av <- mean(sapply(rct_ob_sex_smd_sim[[1]], mean))
rct_ob_treat_smd_av <- mean(sapply(rct_ob_treat_smd_sim[[1]], mean))
rct_ob_case_smd_av <- mean(sapply(rct_ob_case_smd_sim[[1]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT/External","RCT/External","RCT/External","RCT/External","RCT/Observational",
           "RCT/Observational","RCT/Observational","RCT/Observational",
           "Observational/Exernal","Observational/Exernal","Observational/Exernal","Observational/Exernal"),
  SMD= c(rct_ex_age_smd_av,  rct_ex_sex_smd_av, rct_ex_treat_smd_av, rct_ex_case_smd_av, 
         ob_ex_age_smd_av, ob_ex_sex_smd_av, ob_ex_treat_smd_av, ob_ex_case_smd_av,
         rct_ob_age_smd_av,  rct_ob_sex_smd_av, rct_ob_treat_smd_av, rct_ob_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen1_plot_ub_ll_sim.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 4
rct_ex_age_smd_sap <- sapply(rct_ex_age_smd_sim[[1]], mean)

rct_ex_sex_smd_sap <-sapply(rct_ex_sex_smd_sim[[1]], mean)
rct_ex_treat_smd_sap <- sapply(rct_ex_treat_smd_sim[[1]], mean)
rct_ex_case_smd_sap <- sapply(rct_ex_case_smd_sim[[1]], mean)

rct_ob_age_smd_sap <- sapply(rct_ob_age_smd_sim[[1]], mean)

rct_ob_sex_smd_sap <- sapply(rct_ob_sex_smd_sim[[1]], mean)
rct_ob_treat_smd_sap <- sapply(rct_ob_treat_smd_sim[[1]], mean)
rct_ob_case_smd_sap <- sapply(rct_ob_case_smd_sim[[1]], mean)


ob_ex_age_smd_sap <- sapply(ob_ex_age_smd_sim[[1]], mean)

ob_ex_sex_smd_sap <- sapply(ob_ex_sex_smd_sim[[1]], mean)
ob_ex_treat_smd_sap <- sapply(ob_ex_treat_smd_sim[[1]], mean)
ob_ex_case_smd_sap <- sapply(ob_ex_case_smd_sim[[1]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_ex_age_smd_sap,  rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT/Ex","RCT/Ex","RCT/Ex","RCT/Ex",
           "Obs/Ex","Obs/Ex","Obs/Ex","Obs/Ex",
           "RCT/Obs","RCT/Obs","RCT/Obs","RCT/Obs"),
  SMD= c(rct_ex_age_smd_sap, rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap, 
         ob_ex_age_smd_sap, ob_ex_sex_smd_sap, ob_ex_treat_smd_sap, ob_ex_case_smd_sap,
         rct_ob_age_smd_sap, rct_ob_sex_smd_sap, rct_ob_treat_smd_sap, rct_ob_case_smd_sap)
  
)





# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen1_box_ub_ll_sim.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + 
  ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen1_box_gray_ub_ll_sim.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1)+
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()












# test and repeat for other scenarios



# calculate mean smd for each variable - for scenario 1, SYNTH COMPARISON
rct_ex_age_smd_av <- mean(sapply(rct_ex_age_smd_sim[[2]], mean))

rct_ex_sex_smd_av <- mean(sapply(rct_ex_sex_smd_sim[[2]], mean))
rct_ex_treat_smd_av <- mean(sapply(rct_ex_treat_smd_sim[[2]], mean))
rct_ex_case_smd_av <- mean(sapply(rct_ex_case_smd_sim[[2]], mean))

ob_ex_age_smd_av <- mean(sapply(ob_ex_age_smd_sim[[2]], mean))

ob_ex_sex_smd_av <- mean(sapply(ob_ex_sex_smd_sim[[2]], mean))
ob_ex_treat_smd_av <- mean(sapply(ob_ex_treat_smd_sim[[2]], mean))
ob_ex_case_smd_av <- mean(sapply(ob_ex_case_smd_sim[[2]], mean))


rct_ob_age_smd_av <- mean(sapply(rct_ob_age_smd_sim[[2]], mean))

rct_ob_sex_smd_av <- mean(sapply(rct_ob_sex_smd_sim[[2]], mean))
rct_ob_treat_smd_av <- mean(sapply(rct_ob_treat_smd_sim[[2]], mean))
rct_ob_case_smd_av <- mean(sapply(rct_ob_case_smd_sim[[2]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT/External","RCT/External","RCT/External","RCT/External","RCT/Observational",
           "RCT/Observational","RCT/Observational","RCT/Observational",
           "Observational/Exernal","Observational/Exernal","Observational/Exernal","Observational/Exernal"),
  SMD= c(rct_ex_age_smd_av,  rct_ex_sex_smd_av, rct_ex_treat_smd_av, rct_ex_case_smd_av, 
         ob_ex_age_smd_av, ob_ex_sex_smd_av, ob_ex_treat_smd_av, ob_ex_case_smd_av,
         rct_ob_age_smd_av,  rct_ob_sex_smd_av, rct_ob_treat_smd_av, rct_ob_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen2_plot_ub_ll_sim.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 4
rct_ex_age_smd_sap <- sapply(rct_ex_age_smd_sim[[2]], mean)

rct_ex_sex_smd_sap <-sapply(rct_ex_sex_smd_sim[[2]], mean)
rct_ex_treat_smd_sap <- sapply(rct_ex_treat_smd_sim[[2]], mean)
rct_ex_case_smd_sap <- sapply(rct_ex_case_smd_sim[[2]], mean)

rct_ob_age_smd_sap <- sapply(rct_ob_age_smd_sim[[2]], mean)

rct_ob_sex_smd_sap <- sapply(rct_ob_sex_smd_sim[[2]], mean)
rct_ob_treat_smd_sap <- sapply(rct_ob_treat_smd_sim[[2]], mean)
rct_ob_case_smd_sap <- sapply(rct_ob_case_smd_sim[[2]], mean)


ob_ex_age_smd_sap <- sapply(ob_ex_age_smd_sim[[2]], mean)

ob_ex_sex_smd_sap <- sapply(ob_ex_sex_smd_sim[[2]], mean)
ob_ex_treat_smd_sap <- sapply(ob_ex_treat_smd_sim[[2]], mean)
ob_ex_case_smd_sap <- sapply(ob_ex_case_smd_sim[[2]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_ex_age_smd_sap,  rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT/Ex","RCT/Ex","RCT/Ex","RCT/Ex",
           "Obs/Ex","Obs/Ex","Obs/Ex","Obs/Ex",
           "RCT/Obs","RCT/Obs","RCT/Obs","RCT/Obs"),
  SMD= c(rct_ex_age_smd_sap, rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap, 
         ob_ex_age_smd_sap, ob_ex_sex_smd_sap, ob_ex_treat_smd_sap, ob_ex_case_smd_sap,
         rct_ob_age_smd_sap, rct_ob_sex_smd_sap, rct_ob_treat_smd_sap, rct_ob_case_smd_sap)
  
)





# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen2_box_ub_ll_sim.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + 
  ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen2_box_gray_ub_ll_sim.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1)+
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()













# calculate mean smd for each variable - for scenario 1, SYNTH COMPARISON
rct_ex_age_smd_av <- mean(sapply(rct_ex_age_smd_sim[[3]], mean))

rct_ex_sex_smd_av <- mean(sapply(rct_ex_sex_smd_sim[[3]], mean))
rct_ex_treat_smd_av <- mean(sapply(rct_ex_treat_smd_sim[[3]], mean))
rct_ex_case_smd_av <- mean(sapply(rct_ex_case_smd_sim[[3]], mean))

ob_ex_age_smd_av <- mean(sapply(ob_ex_age_smd_sim[[3]], mean))

ob_ex_sex_smd_av <- mean(sapply(ob_ex_sex_smd_sim[[3]], mean))
ob_ex_treat_smd_av <- mean(sapply(ob_ex_treat_smd_sim[[3]], mean))
ob_ex_case_smd_av <- mean(sapply(ob_ex_case_smd_sim[[3]], mean))


rct_ob_age_smd_av <- mean(sapply(rct_ob_age_smd_sim[[3]], mean))

rct_ob_sex_smd_av <- mean(sapply(rct_ob_sex_smd_sim[[3]], mean))
rct_ob_treat_smd_av <- mean(sapply(rct_ob_treat_smd_sim[[3]], mean))
rct_ob_case_smd_av <- mean(sapply(rct_ob_case_smd_sim[[3]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT/External","RCT/External","RCT/External","RCT/External","RCT/Observational",
           "RCT/Observational","RCT/Observational","RCT/Observational",
           "Observational/Exernal","Observational/Exernal","Observational/Exernal","Observational/Exernal"),
  SMD= c(rct_ex_age_smd_av,  rct_ex_sex_smd_av, rct_ex_treat_smd_av, rct_ex_case_smd_av, 
         ob_ex_age_smd_av, ob_ex_sex_smd_av, ob_ex_treat_smd_av, ob_ex_case_smd_av,
         rct_ob_age_smd_av,  rct_ob_sex_smd_av, rct_ob_treat_smd_av, rct_ob_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen3_plot_ub_ll_sim.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 4
rct_ex_age_smd_sap <- sapply(rct_ex_age_smd_sim[[3]], mean)

rct_ex_sex_smd_sap <-sapply(rct_ex_sex_smd_sim[[3]], mean)
rct_ex_treat_smd_sap <- sapply(rct_ex_treat_smd_sim[[3]], mean)
rct_ex_case_smd_sap <- sapply(rct_ex_case_smd_sim[[3]], mean)

rct_ob_age_smd_sap <- sapply(rct_ob_age_smd_sim[[3]], mean)

rct_ob_sex_smd_sap <- sapply(rct_ob_sex_smd_sim[[3]], mean)
rct_ob_treat_smd_sap <- sapply(rct_ob_treat_smd_sim[[3]], mean)
rct_ob_case_smd_sap <- sapply(rct_ob_case_smd_sim[[3]], mean)


ob_ex_age_smd_sap <- sapply(ob_ex_age_smd_sim[[3]], mean)

ob_ex_sex_smd_sap <- sapply(ob_ex_sex_smd_sim[[3]], mean)
ob_ex_treat_smd_sap <- sapply(ob_ex_treat_smd_sim[[3]], mean)
ob_ex_case_smd_sap <- sapply(ob_ex_case_smd_sim[[3]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_ex_age_smd_sap,  rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT/Ex","RCT/Ex","RCT/Ex","RCT/Ex",
           "Obs/Ex","Obs/Ex","Obs/Ex","Obs/Ex",
           "RCT/Obs","RCT/Obs","RCT/Obs","RCT/Obs"),
  SMD= c(rct_ex_age_smd_sap, rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap, 
         ob_ex_age_smd_sap, ob_ex_sex_smd_sap, ob_ex_treat_smd_sap, ob_ex_case_smd_sap,
         rct_ob_age_smd_sap, rct_ob_sex_smd_sap, rct_ob_treat_smd_sap, rct_ob_case_smd_sap)
  
)





# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen3_box_ub_ll_sim.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + 
  ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen3_box_gray_ub_ll_sim.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1)+
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()













# calculate mean smd for each variable - for scenario 1, SYNTH COMPARISON
rct_ex_age_smd_av <- mean(sapply(rct_ex_age_smd_sim[[4]], mean))

rct_ex_sex_smd_av <- mean(sapply(rct_ex_sex_smd_sim[[4]], mean))
rct_ex_treat_smd_av <- mean(sapply(rct_ex_treat_smd_sim[[4]], mean))
rct_ex_case_smd_av <- mean(sapply(rct_ex_case_smd_sim[[4]], mean))

ob_ex_age_smd_av <- mean(sapply(ob_ex_age_smd_sim[[4]], mean))

ob_ex_sex_smd_av <- mean(sapply(ob_ex_sex_smd_sim[[4]], mean))
ob_ex_treat_smd_av <- mean(sapply(ob_ex_treat_smd_sim[[4]], mean))
ob_ex_case_smd_av <- mean(sapply(ob_ex_case_smd_sim[[4]], mean))


rct_ob_age_smd_av <- mean(sapply(rct_ob_age_smd_sim[[4]], mean))

rct_ob_sex_smd_av <- mean(sapply(rct_ob_sex_smd_sim[[4]], mean))
rct_ob_treat_smd_av <- mean(sapply(rct_ob_treat_smd_sim[[4]], mean))
rct_ob_case_smd_av <- mean(sapply(rct_ob_case_smd_sim[[4]], mean))






graph_data <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome"),
  type = c("RCT/External","RCT/External","RCT/External","RCT/External","RCT/Observational",
           "RCT/Observational","RCT/Observational","RCT/Observational",
           "Observational/Exernal","Observational/Exernal","Observational/Exernal","Observational/Exernal"),
  SMD= c(rct_ex_age_smd_av,  rct_ex_sex_smd_av, rct_ex_treat_smd_av, rct_ex_case_smd_av, 
         ob_ex_age_smd_av, ob_ex_sex_smd_av, ob_ex_treat_smd_av, ob_ex_case_smd_av,
         rct_ob_age_smd_av,  rct_ob_sex_smd_av, rct_ob_treat_smd_av, rct_ob_case_smd_av)
  
)



# plot
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen4_plot_ub_ll_sim.jpg", width = 600, height = 400)

ggplot(graph_data, aes(y = label, x = SMD, shape=type)) +
  geom_point(size = 3, position = position_dodgev(height = 0.5)) +
  ggtitle("SMD for RCT Real and Synthetic Data, n= xxxxx")+
  xlab("Average Standardised Mean Difference (SMD)") +
  scale_shape_manual(values=c(0,5,2))+
  ylab("Variable") + xlim(-0.5, 0.5) +  # Set the x-axis limits to -0.5 to 0.5 +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))

dev.off()









# boxplot 


# sapply each value into numbers, not list of value 4
rct_ex_age_smd_sap <- sapply(rct_ex_age_smd_sim[[4]], mean)

rct_ex_sex_smd_sap <-sapply(rct_ex_sex_smd_sim[[4]], mean)
rct_ex_treat_smd_sap <- sapply(rct_ex_treat_smd_sim[[4]], mean)
rct_ex_case_smd_sap <- sapply(rct_ex_case_smd_sim[[4]], mean)

rct_ob_age_smd_sap <- sapply(rct_ob_age_smd_sim[[4]], mean)

rct_ob_sex_smd_sap <- sapply(rct_ob_sex_smd_sim[[4]], mean)
rct_ob_treat_smd_sap <- sapply(rct_ob_treat_smd_sim[[4]], mean)
rct_ob_case_smd_sap <- sapply(rct_ob_case_smd_sim[[4]], mean)


ob_ex_age_smd_sap <- sapply(ob_ex_age_smd_sim[[4]], mean)

ob_ex_sex_smd_sap <- sapply(ob_ex_sex_smd_sim[[4]], mean)
ob_ex_treat_smd_sap <- sapply(ob_ex_treat_smd_sim[[4]], mean)
ob_ex_case_smd_sap <- sapply(ob_ex_case_smd_sim[[4]], mean)



data <- data.frame(
  name=c(rep("Age",100),  rep("Sex",100), rep("Treatment",100), rep('Outcome', 100)),
  value=c(rct_ex_age_smd_sap,  rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap )
)

graph_data_sap <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ## This provides an order to the data
  label = c("Age",  "Sex", "Treatment", "Outcome",
            "Age",  "Sex", "Treatment", "Outcome",
            "Age", "Sex", "Treatment", "Outcome"),
  type = c("RCT/Ex","RCT/Ex","RCT/Ex","RCT/Ex",
           "Obs/Ex","Obs/Ex","Obs/Ex","Obs/Ex",
           "RCT/Obs","RCT/Obs","RCT/Obs","RCT/Obs"),
  SMD= c(rct_ex_age_smd_sap, rct_ex_sex_smd_sap, rct_ex_treat_smd_sap, rct_ex_case_smd_sap, 
         ob_ex_age_smd_sap, ob_ex_sex_smd_sap, ob_ex_treat_smd_sap, ob_ex_case_smd_sap,
         rct_ob_age_smd_sap, rct_ob_sex_smd_sap, rct_ob_treat_smd_sap, rct_ob_case_smd_sap)
  
)



# basic
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen4_box_ub_ll_sim.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot( aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + 
  ylim(-1, 1) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))

dev.off()


# basic grayscale
jpeg("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures/Scen4_box_gray_ub_ll_sim.jpg", width = 600, height = 400)

graph_data_sap %>%
  ggplot(aes(x=label, y=SMD, fill=type)) +
  geom_boxplot() + ylim(-1, 1)+
  scale_fill_manual(values = c("white", "lightgray", "darkgray")) +  # Custom grayscale values
  theme_ipsum() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black"),  # Ensure axis titles are visible
    
  ) +
  ggtitle("Standard Mean Difference by Data Type (Scenario 4)") +
  xlab("Variables") + 
  guides(fill = guide_legend(title = "Data Type"))


dev.off()



















# summary statistics for smd:

rct_age_smd_sap <- sapply(rct_age_smd[[1]], mean)
#rct_race_smd_sap <- sapply(rct_race_smd[[1]], mean)
rct_sex_smd_sap <-sapply(rct_sex_smd[[1]], mean)
rct_treat_smd_sap <- sapply(rct_treat_smd[[1]], mean)
rct_outcome_smd_sap <- sapply(rct_outcome_smd[[1]], mean)

ob_age_smd_sap <- sapply(ob_age_smd[[1]], mean)
#ob_race_smd_sap <- sapply(ob_race_smd[[1]], mean)
ob_sex_smd_sap <- sapply(ob_sex_smd[[1]], mean)
ob_exposed_smd_sap <- sapply(ob_exposed_smd[[1]], mean)
ob_case_smd_sap <- sapply(ob_case_smd[[1]], mean)


ex_age_smd_sap <- sapply(ex_age_smd[[1]], mean)
#ex_race_smd_sap <- sapply(ex_race_smd[[1]], mean)
ex_sex_smd_sap <- sapply(ex_sex_smd[[1]], mean)
ex_treatment_smd_sap <- sapply(ex_treatment_smd[[1]], mean)
ex_case_smd_sap <- sapply(ex_case_smd[[1]], mean)

print(summary(rct_age_smd_sap))
print(summary(rct_race_smd_sap))
print(summary(rct_sex_smd_sap))
print(summary(rct_treat_smd_sap))
print(summary(rct_outcome_smd_sap))

print(summary(ob_age_smd_sap))
print(summary(ob_race_smd_sap))
print(summary(ob_sex_smd_sap))
print(summary(ob_exposed_smd_sap))
print(summary(ob_case_smd_sap))

print(summary(ex_age_smd_sap))
print(summary(ex_race_smd_sap))
print(summary(ex_sex_smd_sap))
print(summary(ex_treatment_smd_sap))
print(summary(ex_case_smd_sap))




rct_age_smd_sap <- sapply(rct_age_smd[[2]], mean)
rct_race_smd_sap <- sapply(rct_race_smd[[2]], mean)
rct_sex_smd_sap <-sapply(rct_sex_smd[[2]], mean)
rct_treat_smd_sap <- sapply(rct_treat_smd[[2]], mean)
rct_outcome_smd_sap <- sapply(rct_outcome_smd[[2]], mean)

ob_age_smd_sap <- sapply(ob_age_smd[[2]], mean)
ob_race_smd_sap <- sapply(ob_race_smd[[2]], mean)
ob_sex_smd_sap <- sapply(ob_sex_smd[[2]], mean)
ob_exposed_smd_sap <- sapply(ob_exposed_smd[[2]], mean)
ob_case_smd_sap <- sapply(ob_case_smd[[2]], mean)


ex_age_smd_sap <- sapply(ex_age_smd[[2]], mean)
ex_race_smd_sap <- sapply(ex_race_smd[[2]], mean)
ex_sex_smd_sap <- sapply(ex_sex_smd[[2]], mean)
ex_treatment_smd_sap <- sapply(ex_treatment_smd[[2]], mean)
ex_case_smd_sap <- sapply(ex_case_smd[[2]], mean)

print(summary(rct_age_smd_sap))
print(summary(rct_race_smd_sap))
print(summary(rct_sex_smd_sap))
print(summary(rct_treat_smd_sap))
print(summary(rct_outcome_smd_sap))

print(summary(ob_age_smd_sap))
print(summary(ob_race_smd_sap))
print(summary(ob_sex_smd_sap))
print(summary(ob_exposed_smd_sap))
print(summary(ob_case_smd_sap))

print(summary(ex_age_smd_sap))
print(summary(ex_race_smd_sap))
print(summary(ex_sex_smd_sap))
print(summary(ex_treatment_smd_sap))
print(summary(ex_case_smd_sap))



rct_age_smd_sap <- sapply(rct_age_smd[[3]], mean)
rct_race_smd_sap <- sapply(rct_race_smd[[3]], mean)
rct_sex_smd_sap <-sapply(rct_sex_smd[[3]], mean)
rct_treat_smd_sap <- sapply(rct_treat_smd[[3]], mean)
rct_outcome_smd_sap <- sapply(rct_outcome_smd[[3]], mean)

ob_age_smd_sap <- sapply(ob_age_smd[[3]], mean)
ob_race_smd_sap <- sapply(ob_race_smd[[3]], mean)
ob_sex_smd_sap <- sapply(ob_sex_smd[[3]], mean)
ob_exposed_smd_sap <- sapply(ob_exposed_smd[[3]], mean)
ob_case_smd_sap <- sapply(ob_case_smd[[3]], mean)


ex_age_smd_sap <- sapply(ex_age_smd[[3]], mean)
ex_race_smd_sap <- sapply(ex_race_smd[[3]], mean)
ex_sex_smd_sap <- sapply(ex_sex_smd[[3]], mean)
ex_treatment_smd_sap <- sapply(ex_treatment_smd[[3]], mean)
ex_case_smd_sap <- sapply(ex_case_smd[[3]], mean)


print(summary(rct_age_smd_sap))
print(summary(rct_race_smd_sap))
print(summary(rct_sex_smd_sap))
print(summary(rct_treat_smd_sap))
print(summary(rct_outcome_smd_sap))

print(summary(ob_age_smd_sap))
print(summary(ob_race_smd_sap))
print(summary(ob_sex_smd_sap))
print(summary(ob_exposed_smd_sap))
print(summary(ob_case_smd_sap))

print(summary(ex_age_smd_sap))
print(summary(ex_race_smd_sap))
print(summary(ex_sex_smd_sap))
print(summary(ex_treatment_smd_sap))
print(summary(ex_case_smd_sap))



rct_age_smd_sap <- sapply(rct_age_smd[[4]], mean)
rct_race_smd_sap <- sapply(rct_race_smd[[4]], mean)
rct_sex_smd_sap <-sapply(rct_sex_smd[[4]], mean)
rct_treat_smd_sap <- sapply(rct_treat_smd[[4]], mean)
rct_outcome_smd_sap <- sapply(rct_outcome_smd[[4]], mean)

ob_age_smd_sap <- sapply(ob_age_smd[[4]], mean)
ob_race_smd_sap <- sapply(ob_race_smd[[4]], mean)
ob_sex_smd_sap <- sapply(ob_sex_smd[[4]], mean)
ob_exposed_smd_sap <- sapply(ob_exposed_smd[[4]], mean)
ob_case_smd_sap <- sapply(ob_case_smd[[4]], mean)


ex_age_smd_sap <- sapply(ex_age_smd[[4]], mean)
ex_race_smd_sap <- sapply(ex_race_smd[[4]], mean)
ex_sex_smd_sap <- sapply(ex_sex_smd[[4]], mean)
ex_treatment_smd_sap <- sapply(ex_treatment_smd[[4]], mean)
ex_case_smd_sap <- sapply(ex_case_smd[[4]], mean)

print(summary(rct_age_smd_sap))
print(summary(rct_race_smd_sap))
print(summary(rct_sex_smd_sap))
print(summary(rct_treat_smd_sap))
print(summary(rct_outcome_smd_sap))

print(summary(ob_age_smd_sap))
print(summary(ob_race_smd_sap))
print(summary(ob_sex_smd_sap))
print(summary(ob_exposed_smd_sap))
print(summary(ob_case_smd_sap))

print(summary(ex_age_smd_sap))
print(summary(ex_race_smd_sap))
print(summary(ex_sex_smd_sap))
print(summary(ex_treatment_smd_sap))
print(summary(ex_case_smd_sap))

