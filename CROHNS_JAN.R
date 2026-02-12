# Crohns




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
library(gridExtra)
library(rsimsum)
library(parallel)
library(BSDA)
library(summarytools)
# updated sample sizes for better runtime:
# 100 moved to 150 to avoid issues with synthpop
# 20,000 shrunk to 500 to speed up runtime

# creating scenarios for number of observations
# halve the n due to pre and post being added together later
nsize <- function (scenario, datatype, simsyn) {
  if (scenario == 1) { # scenario one is n = 150 to n = 150
    N <- 75
  }   else if (scenario == 2 ) { # scenario two is n = 500 to n = 500
    N <- 250
  }   else if (scenario == 3 & datatype == "rct") { # scenario three is realistic sample sizes
    N <- 628
  }   else if (scenario == 3 & datatype == "ob"){
    N <- 221
  }   else if (scenario == 3 & datatype == "ex"){
    N <- 322
  }   else if (scenario == 4 & simsyn == "sim"){ # scenario four is n = 150 to n = 500
    N <- 75
  }   else if (scenario == 4 & simsyn == "syn"){
    N <- 250
  }
  return(N)
}



# function for smd calculation:
process_data <- function(df1, df2, var_names, table_list_name, smd_list_prefix, i) {
  # merge df
  df1$Key <- 1:nrow(df1)
  df2$Key <- 1:nrow(df2)
  merged_df <- merge(df1, df2, by = "Key", all = TRUE)
  merged_df <- merged_df[, -which(names(merged_df) == "Key")]

  # make a single df
  combined_df <- data.frame(
    lapply(var_names, function(var) {
      c(merged_df[[paste0(var, ".x")]], merged_df[[paste0(var, ".y")]])
    })
  )
  colnames(combined_df) <- var_names

  # save table
  temp <- get(table_list_name, envir = .GlobalEnv)
  temp[[i]] <- CreateTableOne(vars = var_names, data = combined_df, test = FALSE)
  assign(table_list_name, temp, envir = .GlobalEnv)

  # calculate smd
  for (var in var_names) {
    x <- merged_df[[paste0(var, ".x")]]
    y <- merged_df[[paste0(var, ".y")]]

    pooled_sd <- sqrt((sd(x, na.rm = TRUE)^2 + sd(y, na.rm = TRUE)^2) / 2)
    mean_diff <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
    smd_value <- mean_diff / pooled_sd

    smd_list_name <- paste0(smd_list_prefix, "_", var, "_smd_all")
    smd_list <- get(smd_list_name, envir = .GlobalEnv)
    smd_list[[i]] <- smd_value
    assign(smd_list_name, smd_list, envir = .GlobalEnv)
  }
}


# function for smd calculation across data types
compute_smd_across <- function(df1, df2, vars, prefix, suffix, i) {
  for (var in vars) {
    x1 <- df1[[paste0(var, ".", suffix)]]
    x2 <- df2[[paste0(var, ".", suffix)]]

    sd1 <- sd(x1, na.rm = TRUE)
    sd2 <- sd(x2, na.rm = TRUE)
    pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
    mean_diff <- mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)

    # by summix
    smd_list_name <- if (suffix == "x") {
      paste0(prefix, "_", var, "_smd_sim")
    } else if (suffix == "y") {
      paste0(prefix, "_", var, "_smd_syn")
    } else {
      stop("compute_smd_across: suffix must be 'x' or 'y'")
    }

    smd_list <- get(smd_list_name, envir = .GlobalEnv)
    smd_list[[i]] <- mean_diff / pooled_sd
    assign(smd_list_name, smd_list, envir = .GlobalEnv)
  }
}
# function for making empty variable lists
create_smd_lists <- function(non_combo_prefixes, combo_prefixes, variables, n) {

  for (pref in non_combo_prefixes) {
    for (var in variables) {
      list_name <- paste0(pref, "_", var, "_smd_all")
      assign(list_name, vector("list", length = n), envir = .GlobalEnv)
    }
  }
  # combos for sim/syn
  for (pref in combo_prefixes) {
    for (var in variables) {
      assign(paste0(pref, "_", var, "_smd_sim"), vector("list", length = n), envir = .GlobalEnv)
      assign(paste0(pref, "_", var, "_smd_syn"), vector("list", length = n), envir = .GlobalEnv)
    }
  }
}
# function for making empty table lists
create_table_lists <- function(names, n) {
  for (name in names) {
    assign(name, vector("list", length = n), envir = .GlobalEnv)
  }
}


# function for control response rate

control_response <- function(df, list_name, i) {
  num_responders <- df %>%
    filter(treat == 0, outcome == 1) %>%
    nrow()

  num_control <- df %>%
    filter(treat == 0) %>%
    nrow()

  rate <- ifelse(num_control > 0, num_responders / num_control, NA)

  tmp <- get(list_name, envir = .GlobalEnv)
  tmp[[i]] <- rate
  assign(list_name, tmp, envir = .GlobalEnv)
}


# for loop for sim and analysis

nSim <- 400 # number of simulations # small for testing
set.seed(9999)


#  make empty lists

vars <- c("age", "sex", "treat", "outcome")
non_combo_prefixes <- c("rct", "ob", "ex")                   # these get *_smd_all
combo_prefixes     <- c("rct_ob", "ob_ex", "rct_ex")         # these get *_smd_sim and *_smd_syn


create_smd_lists(non_combo_prefixes, combo_prefixes, vars, nSim)

create_table_lists(c("table_rct_data", "table_ob_data", "table_ex_data"), nSim)

sim_results_rct <- vector("list", length = nSim)
sim_results_ob  <- vector("list", length = nSim)
sim_results_ex  <- vector("list", length = nSim)
syn_results_rct <- vector("list", length = nSim)
syn_results_ob <- vector("list", length = nSim)
syn_results_ex <- vector("list", length = nSim)

sim_resp_rct <- vector("list", length = nSim)
sim_resp_ob  <- vector("list", length = nSim)
sim_resp_ex  <- vector("list", length = nSim)
syn_resp_rct <- vector("list", length = nSim)
syn_resp_ob <- vector("list", length = nSim)
syn_resp_ex <- vector("list", length = nSim)

# list for sum stats

sim_sum_rct_mean <- vector("list", length = nSim)
sim_sum_ob_mean  <- vector("list", length = nSim)
sim_sum_ex_mean  <- vector("list", length = nSim)
syn_sum_rct_mean <- vector("list", length = nSim)
syn_sum_ob_mean <- vector("list", length = nSim)
syn_sum_ex_mean <- vector("list", length = nSim)


sim_sum_rct_med <- vector("list", length = nSim)
sim_sum_ob_med  <- vector("list", length = nSim)
sim_sum_ex_med  <- vector("list", length = nSim)
syn_sum_rct_med <- vector("list", length = nSim)
syn_sum_ob_med <- vector("list", length = nSim)
syn_sum_ex_med <- vector("list", length = nSim)

sim_sum_rct_range <- vector("list", length = nSim)
sim_sum_ob_range  <- vector("list", length = nSim)
sim_sum_ex_range  <- vector("list", length = nSim)
syn_sum_rct_range <- vector("list", length = nSim)
syn_sum_ob_range <- vector("list", length = nSim)
syn_sum_ex_range <- vector("list", length = nSim)

sim_sum_rct_sd <- vector("list", length = nSim)
sim_sum_ob_sd  <- vector("list", length = nSim)
sim_sum_ex_sd  <- vector("list", length = nSim)
syn_sum_rct_sd <- vector("list", length = nSim)
syn_sum_ob_sd <- vector("list", length = nSim)
syn_sum_ex_sd <- vector("list", length = nSim)

# list for dist
sim_kurt_rct <- vector("list", length = nSim)
sim_skew_rct <- vector("list", length = nSim)

sim_kurt_ob <- vector("list", length = nSim)
sim_skew_ob <- vector("list", length = nSim)

sim_kurt_ex <- vector("list", length = nSim)
sim_skew_ex <- vector("list", length = nSim)

syn_kurt_rct <- vector("list", length = nSim)
syn_skew_rct <- vector("list", length = nSim)

syn_kurt_ob <- vector("list", length = nSim)
syn_skew_ob <- vector("list", length = nSim)

syn_kurt_ex <- vector("list", length = nSim)
syn_skew_ex <- vector("list", length = nSim)

scenario <- 1  #  scenario before the loop


# loop for gen simulations


for (i in 1:nSim) {

  if ((i - 1) %% (nSim/4) == 0 & i > 1) {  # change scenario every 10,000 iterations
    scenario <- scenario + 1
  }



  N <- nsize(scenario, "rct", "sim") # selecting scenario


  treat <- sample(c(1, 0), N, replace = TRUE, prob = c(0.63, 0.37))

  # compliance to treatment
  comply <- ifelse(treat == 1, sample(c(1, 0), N, replace = TRUE, prob = c(1, 0)), sample(c(1, 0), N, replace = TRUE, prob = c(0, 1)))

  # other variables
  sex <- rbinom(N, 1, 0.522)
  age <- rtruncnorm(N, a = 18, b = 72, mean = 37.3, sd = 11.8)
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
    mutate(treat_effect = ifelse(post == 0, NA, treat_effect))  # filter treatment effect to apply only to post-treatment group

  # calculate sex_effect
  df_combined <- df_combined %>%
    mutate(sex_effect = ifelse(sex == 1 | sex == 0, 0, NA)) # sex effect not recorded or reported / only baseline / assumed as zero

  # create overall treatment effect including demographic influences
  df_combined <- df_combined %>%
    mutate(overall_treat_effect = ifelse(post == 0, 0, treat_effect + sex_effect))

  # create pre and post outcomes (baseline vs treatment)
  df_combined <- df_combined %>%
    mutate(outcome = ifelse(post == 0, rbinom(nrow(df_combined[df_combined$post == 0, ]), 1, 0.359),
                            sapply(overall_treat_effect, function(prob) rbinom(1, 1, prob))))
  # store dfs in results
  sim_results_rct[[i]] <- df_combined

  # remove df2
  rm(df2)





  # observational
  # treatment

  N <- nsize(scenario, "ob", "sim") # selecting scenario

  treat=sample(c(1,0), N,replace=TRUE, prob=c(.66,.34)) # probability of being exposed to the treatment or not


  # other variables
  sex <- rbinom(N, 1, 0.602)
  age <- rtruncnorm(N, a = 16, b = 65, mean = 38.2, sd = 16.96)
  age <- round(age, 0) # round age to whole numbers
  random_var <- rbinom(N, 1, prob = 0.5)
  post <- rep(0, N)  # pre and post assignment

  # create df combining variables
  ob_df <- data.frame(treat, sex, age, random_var, post)

  # add id
  ob_df <- ob_df %>% dplyr::mutate(id = row_number())

  # create post df
  ob_df2 <- ob_df %>% mutate(post = 1)

  # combine pre and post
  ob_df_combined <- rbind(ob_df, ob_df2)

  # define treatment effect, only for post group
  # define treat
  ob_df_combined <- ob_df_combined %>% mutate(treat_effect=
                                                ifelse(treat==1, .387,
                                                       ifelse(treat==0, .242,
                                                              NA)))%>%
    mutate(treat_effect = ifelse(post == 0, NA, treat_effect))



  # calculate sex_effect and race_effect
  ob_df_combined <- ob_df_combined %>%
    mutate(sex_effect = ifelse(sex == 1 | sex == 0, 0, NA))

  # create overall treatment effect including demographic influences
  ob_df_combined <- ob_df_combined %>%
    mutate(overall_treat_effect = ifelse(post == 0, 0, treat_effect + sex_effect))


  # create pre and post outcomes
  ob_df_combined <- ob_df_combined %>%
    mutate(outcome = ifelse(post == 0, rbinom(sum(post == 0), 1, 0.242),
                            rbinom(sum(post == 1), 1, overall_treat_effect)))
  # store dfs in results
  sim_results_ob[[i]] <- ob_df_combined

  # remove df2
  rm(ob_df2)




  # external
  # treatment
  N <- nsize(scenario, "ex", "sim") # selecting scenario

  outcome=sample(c(1,0), N*2,replace=TRUE, prob=c(.64,.36))

  # other variables
  treat=sample(c(1,0), N*2,replace=TRUE, prob=c(.88,.12))
  sex <- rbinom(N, 1, 0.48)
  age <- rtruncnorm(N*2, a = 15, b = 70, mean = 38, sd = 17.78)
  age <- round(age, 0) # round age to whole numbers
  random_var <- rbinom(N*2, 1, prob = 0.5)



  # create df combining variables
  ex_df <- data.frame(outcome, treat, sex, age, random_var)

  # add id
  ex_df <- ex_df %>% dplyr::mutate(id = row_number())

  # store dfs in results
  sim_results_ex[[i]] <- ex_df











  # begin synthetic data creation:
  # MODEL - rct
  rct_mydata = subset(df_combined, select = c("id", "age", "sex", "treat", "outcome"))

  N <- nsize(scenario, "rct", "syn") # selecting scenario

  # minimum number of observations needed is 10
  rct_mysyn <- syn(rct_mydata, k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)

  summary(rct_mysyn)

  # store dfs in results
  syn_results_rct[[i]] <- rct_mysyn$syn




  # MODEL - obs
  ob_mydata = subset(ob_df_combined, select = c("id", "age", "sex", "treat", "outcome"))



  N <- nsize(scenario, "ob", "syn") # selecting scenario

  # minimum number of observations needed is 10
  ob_mysyn <- syn(ob_mydata, k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)

  summary(ob_mysyn)

  # store dfs in results
  syn_results_ob[[i]] <- ob_mysyn$syn





  # MODEL - ex
  ex_mydata = subset(ex_df, select = c("id", "age", "sex", "treat", "outcome"))


  N <- nsize(scenario, "ex", "syn") # selecting scenario


  print(i) # sim number

  # minimum number of observations needed is 10
  ex_mysyn <- syn(ex_mydata, k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)

  summary(ex_mysyn)

  # store dfs in results
  syn_results_ex[[i]] <- ex_mysyn$syn





  # recreate synthetic data in df format
  rct_mysyn <- rct_mysyn$syn
  ob_mysyn <- ob_mysyn$syn
  ex_mysyn <- ex_mysyn$syn



  # dist stats
  sim_kurt_rct[[i]] <- data.frame(age = kurtosis(df_combined$age),
                                  sex = kurtosis(df_combined$sex),
                                  treat = kurtosis(df_combined$treat),
                                  outcome = kurtosis(df_combined$outcome))


  sim_kurt_ob[[i]] <- data.frame(age = kurtosis(ob_df_combined$age),
                                 sex = kurtosis(ob_df_combined$sex),
                                 treat = kurtosis(ob_df_combined$treat),
                                 outcome = kurtosis(ob_df_combined$outcome))



  sim_kurt_ex[[i]] <- data.frame(age = kurtosis(ex_df$age),
                                 sex = kurtosis(ex_df$sex),
                                 treat = kurtosis(ex_df$treat),
                                 outcome = kurtosis(ex_df$outcome))






  sim_skew_rct[[i]] <- data.frame(age = skewness(df_combined$age),
                                  sex = skewness(df_combined$sex),
                                  treat = skewness(df_combined$treat),
                                  outcome = skewness(df_combined$outcome))


  sim_skew_ob[[i]] <- data.frame(age = skewness(ob_df_combined$age),
                                 sex = skewness(ob_df_combined$sex),
                                 treat = skewness(ob_df_combined$treat),
                                 outcome = skewness(ob_df_combined$outcome))



  sim_skew_ex[[i]] <- data.frame(age = skewness(ex_df$age),
                                 sex = skewness(ex_df$sex),
                                 treat = skewness(ex_df$treat),
                                 outcome = skewness(ex_df$outcome))





  syn_kurt_rct[[i]] <- data.frame(age = kurtosis(rct_mysyn$age),
                                  sex = kurtosis(rct_mysyn$sex),
                                  treat = kurtosis(rct_mysyn$treat),
                                  outcome = kurtosis(rct_mysyn$outcome))


  syn_kurt_ob[[i]] <- data.frame(age = kurtosis(ob_mysyn$age),
                                 sex = kurtosis(ob_mysyn$sex),
                                 treat = kurtosis(ob_mysyn$treat),
                                 outcome = kurtosis(ob_mysyn$outcome))



  syn_kurt_ex[[i]] <- data.frame(age = kurtosis(ex_mysyn$age),
                                 sex = kurtosis(ex_mysyn$sex),
                                 treat = kurtosis(ex_mysyn$treat),
                                 outcome = kurtosis(ex_mysyn$outcome))






  syn_skew_rct[[i]] <- data.frame(age = skewness(rct_mysyn$age),
                                  sex = skewness(rct_mysyn$sex),
                                  treat = skewness(rct_mysyn$treat),
                                  outcome = skewness(rct_mysyn$outcome))


  syn_skew_ob[[i]] <- data.frame(age = skewness(ob_mysyn$age),
                                 sex = skewness(ob_mysyn$sex),
                                 treat = skewness(ob_mysyn$treat),
                                 outcome = skewness(ob_mysyn$outcome))



  syn_skew_ex[[i]] <- data.frame(age = skewness(ex_mysyn$age),
                                 sex = skewness(ex_mysyn$sex),
                                 treat = skewness(ex_mysyn$treat),
                                 outcome = skewness(ex_mysyn$outcome))



  # summary stats

  # mean, median, range, sd,

  sim_sum_rct_mean[[i]] <- data.frame(age = mean(df_combined$age),
                                      sex = mean(df_combined$sex),
                                      treat = mean(df_combined$treat),
                                      outcome = mean(df_combined$outcome))


  sim_sum_rct_med[[i]] <- data.frame(age = median(df_combined$age),
                                     sex = median(df_combined$sex),
                                     treat = median(df_combined$treat),
                                     outcome = median(df_combined$outcome))

  sim_sum_rct_range[[i]] <- data.frame(age = range(df_combined$age),
                                       sex = range(df_combined$sex),
                                       treat = range(df_combined$treat),
                                       outcome = range(df_combined$outcome))


  sim_sum_rct_sd[[i]] <- data.frame(age = sd(df_combined$age),
                                    sex = sd(df_combined$sex),
                                    treat = sd(df_combined$treat),
                                    outcome = sd(df_combined$outcome))






  sim_sum_ob_mean[[i]]  <- data.frame(age = mean(ob_df_combined$age),
                                      sex = mean(ob_df_combined$sex),
                                      treat = mean(ob_df_combined$treat),
                                      outcome = mean(ob_df_combined$outcome))


  sim_sum_ob_med[[i]]  <- data.frame(age = median(ob_df_combined$age),
                                     sex = median(ob_df_combined$sex),
                                     treat = median(ob_df_combined$treat),
                                     outcome = median(ob_df_combined$outcome))


  sim_sum_ob_range[[i]]  <- data.frame(age = range(ob_df_combined$age),
                                       sex = range(ob_df_combined$sex),
                                       treat = range(ob_df_combined$treat),
                                       outcome = range(ob_df_combined$outcome))



  sim_sum_ob_sd[[i]]  <- data.frame(age = sd(ob_df_combined$age),
                                    sex = sd(ob_df_combined$sex),
                                    treat = sd(ob_df_combined$treat),
                                    outcome = sd(ob_df_combined$outcome))





  sim_sum_ex_mean[[i]]  <- data.frame(age = mean(ex_df$age),
                                      sex = mean(ex_df$sex),
                                      treat = mean(ex_df$treat),
                                      outcome = mean(ex_df$outcome))


  sim_sum_ex_med[[i]]  <- data.frame(age = median(ex_df$age),
                                     sex = median(ex_df$sex),
                                     treat = median(ex_df$treat),
                                     outcome = median(ex_df$outcome))

  sim_sum_ex_range[[i]]  <- data.frame(age = range(ex_df$age),
                                       sex = range(ex_df$sex),
                                       treat = range(ex_df$treat),
                                       outcome = range(ex_df$outcome))

  sim_sum_ex_sd[[i]]  <- data.frame(age = sd(ex_df$age),
                                    sex = sd(ex_df$sex),
                                    treat = sd(ex_df$treat),
                                    outcome = sd(ex_df$outcome))











  syn_sum_rct_mean[[i]] <- data.frame(age = mean(rct_mysyn$age),
                                      sex = mean(rct_mysyn$sex),
                                      treat = mean(rct_mysyn$treat),
                                      outcome = mean(rct_mysyn$outcome))


  syn_sum_rct_med[[i]] <- data.frame(age = median(rct_mysyn$age),
                                     sex = median(rct_mysyn$sex),
                                     treat = median(rct_mysyn$treat),
                                     outcome = median(rct_mysyn$outcome))

  syn_sum_rct_range[[i]] <- data.frame(age = range(rct_mysyn$age),
                                       sex = range(rct_mysyn$sex),
                                       treat = range(rct_mysyn$treat),
                                       outcome = range(rct_mysyn$outcome))



  syn_sum_rct_sd[[i]] <- data.frame(age = sd(rct_mysyn$age),
                                    sex = sd(rct_mysyn$sex),
                                    treat = sd(rct_mysyn$treat),
                                    outcome = sd(rct_mysyn$outcome))






  syn_sum_ob_mean[[i]] <- data.frame(age = mean(ob_mysyn$age),
                                     sex = mean(ob_mysyn$sex),
                                     treat = mean(ob_mysyn$treat),
                                     outcome = mean(ob_mysyn$outcome))


  syn_sum_ob_med[[i]] <- data.frame(age = median(ob_mysyn$age),
                                    sex = median(ob_mysyn$sex),
                                    treat = median(ob_mysyn$treat),
                                    outcome = median(ob_mysyn$outcome))

  syn_sum_ob_range[[i]] <- data.frame(age = range(ob_mysyn$age),
                                      sex = range(ob_mysyn$sex),
                                      treat = range(ob_mysyn$treat),
                                      outcome = range(ob_mysyn$outcome))



  syn_sum_ob_sd[[i]] <- data.frame(age = sd(ob_mysyn$age),
                                   sex = sd(ob_mysyn$sex),
                                   treat = sd(ob_mysyn$treat),
                                   outcome = sd(ob_mysyn$outcome))








  syn_sum_ex_mean[[i]] <- data.frame(age = mean(ex_mysyn$age),
                                     sex = mean(ex_mysyn$sex),
                                     treat = mean(ex_mysyn$treat),
                                     outcome = mean(ex_mysyn$outcome))


  syn_sum_ex_med[[i]] <- data.frame(age = median(ex_mysyn$age),
                                    sex = median(ex_mysyn$sex),
                                    treat = median(ex_mysyn$treat),
                                    outcome = median(ex_mysyn$outcome))

  syn_sum_ex_range[[i]] <- data.frame(age = range(ex_mysyn$age),
                                      sex = range(ex_mysyn$sex),
                                      treat = range(ex_mysyn$treat),
                                      outcome = range(ex_mysyn$outcome))



  syn_sum_ex_sd[[i]] <- data.frame(age = sd(ex_mysyn$age),
                                   sex = sd(ex_mysyn$sex),
                                   treat = sd(ex_mysyn$treat),
                                   outcome = sd(ex_mysyn$outcome))




  # calculating smd
  process_data(df_combined, rct_mysyn, vars, "table_rct_data", "rct", i)
  process_data(ob_df_combined, ob_mysyn, vars, "table_ob_data", "ob", i)
  process_data(ex_df, ex_mysyn, vars, "table_ex_data", "ex", i)

  # Synthetic vs Synthetic
  rct_ex <- merge(rct_mysyn, ex_mysyn, by = "id", suffixes = c(".x", ".y"))
  rct_ob <- merge(rct_mysyn, ob_mysyn, by = "id", suffixes = c(".x", ".y"))
  ob_ex  <- merge(ob_mysyn, ex_mysyn, by = "id", suffixes = c(".x", ".y"))

  compute_smd_across(rct_ex, rct_ex, vars, "rct_ex", "y", i)
  compute_smd_across(rct_ob, rct_ob, vars, "rct_ob", "y", i)
  compute_smd_across(ob_ex,  ob_ex,  vars, "ob_ex",  "y", i)

  compute_smd_across(rct_ex, rct_ex, vars, "rct_ex", "x", i)
  compute_smd_across(rct_ob, rct_ob, vars, "rct_ob", "x", i)
  compute_smd_across(ob_ex,  ob_ex,  vars, "ob_ex",  "x", i)



  # calculate control response rate:

  control_response(rct_mysyn, "syn_resp_rct", i)
  control_response(ob_mysyn,  "syn_resp_ob",  i)
  control_response(ex_mysyn,  "syn_resp_ex",  i)

  control_response(df_combined,     "sim_resp_rct", i)
  control_response(ob_df_combined,  "sim_resp_ob",  i)
  control_response(ex_df,           "sim_resp_ex",  i)
}









# ANALYSIS AND GRAPHING



# function to split lists
distribute_into_lists <- function(input_list, sizes) {
  if (sum(sizes) != length(input_list)) {
    stop("The sum of sizes must equal the length of the input list")
  }

  result <- vector("list", length(sizes))
  start_index <- 1

  for (i in seq_along(sizes)) {
    end_index <- start_index + sizes[i] - 1
    result[[i]] <- input_list[start_index:end_index]
    start_index <- end_index + 1
  }

  return(result)
}


list_size <- c(nSim/4, nSim/4, nSim/4, nSim/4)
save_dir <- "C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Saved_Rdata/August/CR_CART/"


# save data:

# non-combos
non_combos <- c("rct", "ob", "ex")
variables <- c("age", "sex", "treat", "outcome")
non_combo_names <- unlist(lapply(non_combos, function(pref) {
  paste0(pref, "_", variables, "_smd_all")
}))

# combos
combos <- c("rct_ex", "rct_ob", "ob_ex")
combo_names <- unlist(lapply(combos, function(pref) {
  paste0(rep(pref, each = length(variables)), "_", variables, "_smd_", rep(c("sim", "syn"), each = length(variables) * 1))
}))

# misc
extra_objects <- c(
  "sim_results_rct", "sim_results_ob", "sim_results_ex",
  "syn_results_rct", "syn_results_ob", "syn_results_ex",
  "sim_resp_rct", "sim_resp_ob", "sim_resp_ex",
  "syn_resp_rct", "syn_resp_ob", "syn_resp_ex",
  "table_rct_data", "table_ob_data", "table_ex_data"
)

# final list of all object names
obj_names <- c(extra_objects, non_combo_names, combo_names)

# save all
for (obj_name in obj_names) {

  if (!exists(obj_name)) {
    warning(paste("Object not found:", obj_name))
    next
  }

  obj <- get(obj_name)

  # check for list
  if (!is.list(obj)) {
    warning(paste("Object is not a list:", obj_name))
    next
  }

  obj_split <- distribute_into_lists(obj, list_size)

  for (scenario_index in seq_along(obj_split)) {
    saveRDS(
      obj_split[[scenario_index]],
      file = paste0(save_dir, obj_name, "_scenario", scenario_index, ".Rdata")
    )
  }
}

# dist stats

sim_kurt_rct = do.call(rbind, sim_kurt_rct)
sim_kurt_ob = do.call(rbind,sim_kurt_ob)
sim_kurt_ex = do.call(rbind, sim_kurt_ex)
sim_skew_rct = do.call(rbind, sim_skew_rct)
sim_skew_ob = do.call(rbind,sim_skew_ob)
sim_skew_ex = do.call(rbind, sim_skew_ex)

syn_kurt_rct = do.call(rbind, syn_kurt_rct)
syn_kurt_ob = do.call(rbind,syn_kurt_ob)
syn_kurt_ex = do.call(rbind, syn_kurt_ex)
syn_skew_rct = do.call(rbind, syn_skew_rct)
syn_skew_ob = do.call(rbind,syn_skew_ob)
syn_skew_ex = do.call(rbind, syn_skew_ex)

write.csv(sim_kurt_rct, file = "cr_sim_kurt_rct_cart.csv", row.names = FALSE)
write.csv(sim_kurt_ob, file = "cr_sim_kurt_ob_cart.csv", row.names = FALSE)
write.csv(sim_kurt_ex, file = "cr_sim_kurt_ex_cart.csv", row.names = FALSE)
write.csv(sim_skew_rct, file = "cr_sim_skew_rct_cart.csv", row.names = FALSE)
write.csv(sim_skew_ob, file = "cr_sim_skew_ob_cart.csv", row.names = FALSE)
write.csv(sim_skew_ex, file = "cr_sim_skew_ex_cart.csv", row.names = FALSE)

write.csv(syn_kurt_rct, file = "cr_syn_kurt_rct_cart.csv", row.names = FALSE)
write.csv(syn_kurt_ob, file = "cr_syn_kurt_ob_cart.csv", row.names = FALSE)
write.csv(syn_kurt_ex, file = "cr_syn_kurt_ex_cart.csv", row.names = FALSE)
write.csv(syn_skew_rct, file = "cr_syn_skew_rct_cart.csv", row.names = FALSE)
write.csv(syn_skew_ob, file = "cr_syn_skew_ob_cart.csv", row.names = FALSE)
write.csv(syn_skew_ex, file = "cr_syn_skew_ex_cart.csv", row.names = FALSE)



# dealing with and saving summary stats csv

sim_sum_rct_mean = do.call(rbind, sim_sum_rct_mean)
sim_sum_rct_med = do.call(rbind, sim_sum_rct_med)
sim_sum_rct_range = do.call(rbind, sim_sum_rct_range)
sim_sum_rct_sd = do.call(rbind, sim_sum_rct_sd)

sim_sum_ob_mean = do.call(rbind, sim_sum_ob_mean)
sim_sum_ob_med = do.call(rbind, sim_sum_ob_med)
sim_sum_ob_range = do.call(rbind, sim_sum_ob_range)
sim_sum_ob_sd = do.call(rbind, sim_sum_ob_sd)


sim_sum_ex_mean = do.call(rbind, sim_sum_ex_mean)
sim_sum_ex_med = do.call(rbind, sim_sum_ex_med)
sim_sum_ex_range = do.call(rbind, sim_sum_ex_range)
sim_sum_ex_sd = do.call(rbind, sim_sum_ex_sd)



syn_sum_rct_mean = do.call(rbind, syn_sum_rct_mean)
syn_sum_rct_med = do.call(rbind, syn_sum_rct_med)
syn_sum_rct_range = do.call(rbind, syn_sum_rct_range)
syn_sum_rct_sd = do.call(rbind, syn_sum_rct_sd)

syn_sum_ob_mean = do.call(rbind, syn_sum_ob_mean)
syn_sum_ob_med = do.call(rbind, syn_sum_ob_med)
syn_sum_ob_range = do.call(rbind, syn_sum_ob_range)
syn_sum_ob_sd = do.call(rbind, syn_sum_ob_sd)


syn_sum_ex_mean = do.call(rbind, syn_sum_ex_mean)
syn_sum_ex_med = do.call(rbind, syn_sum_ex_med)
syn_sum_ex_range = do.call(rbind, syn_sum_ex_range)
syn_sum_ex_sd = do.call(rbind, syn_sum_ex_sd)








write.csv(sim_sum_rct_mean, file = "cr_sim_sum_rct_cart_mean.csv", row.names = FALSE)
write.csv(sim_sum_rct_med, file = "cr_sim_sum_rct_cart_med.csv", row.names = FALSE)
write.csv(sim_sum_rct_range, file = "cr_sim_sum_rct_cart_range.csv", row.names = FALSE)
write.csv(sim_sum_rct_sd, file = "cr_sim_sum_rct_cart_sd.csv", row.names = FALSE)

write.csv(sim_sum_ob_mean, file = "cr_sim_sum_ob_cart_mean.csv", row.names = FALSE)
write.csv(sim_sum_ob_med, file = "cr_sim_sum_ob_cart_med.csv", row.names = FALSE)
write.csv(sim_sum_ob_range, file = "cr_sim_sum_ob_cart_range.csv", row.names = FALSE)
write.csv(sim_sum_ob_sd, file = "cr_sim_sum_ob_cart_sd.csv", row.names = FALSE)

write.csv(sim_sum_ex_mean, file = "cr_sim_sum_ex_cart_mean.csv", row.names = FALSE)
write.csv(sim_sum_ex_med, file = "cr_sim_sum_ex_cart_med.csv", row.names = FALSE)
write.csv(sim_sum_ex_range, file = "cr_sim_sum_ex_cart_range.csv", row.names = FALSE)
write.csv(sim_sum_ex_sd, file = "cr_sim_sum_ex_cart_sd.csv", row.names = FALSE)


write.csv(syn_sum_rct_mean, file = "cr_syn_sum_rct_cart_mean.csv", row.names = FALSE)
write.csv(syn_sum_rct_med, file = "cr_syn_sum_rct_cart_med.csv", row.names = FALSE)
write.csv(syn_sum_rct_range, file = "cr_syn_sum_rct_cart_range.csv", row.names = FALSE)
write.csv(syn_sum_rct_sd, file = "cr_syn_sum_rct_cart_sd.csv", row.names = FALSE)

write.csv(syn_sum_ob_mean, file = "cr_syn_sum_ob_cart_mean.csv", row.names = FALSE)
write.csv(syn_sum_ob_med, file = "cr_syn_sum_ob_cart_med.csv", row.names = FALSE)
write.csv(syn_sum_ob_range, file = "cr_syn_sum_ob_cart_range.csv", row.names = FALSE)
write.csv(syn_sum_ob_sd, file = "cr_syn_sum_ob_cart_sd.csv", row.names = FALSE)

write.csv(syn_sum_ex_mean, file = "cr_syn_sum_ex_cart_mean.csv", row.names = FALSE)
write.csv(syn_sum_ex_med, file = "cr_syn_sum_ex_cart_med.csv", row.names = FALSE)
write.csv(syn_sum_ex_range, file = "cr_syn_sum_ex_cart_range.csv", row.names = FALSE)
write.csv(syn_sum_ex_sd, file = "cr_syn_sum_ex_cart_sd.csv", row.names = FALSE)





# plotting smd figures

variables <- c("Age", "Sex", "Treatment", "Outcome")
types <- c("RCT", "Observational", "External")

# build variable names
make_varname <- function(prefix, var, suffix = "") {
  paste0(prefix, "_", tolower(var), "_smd", suffix)
}

# get relevant lists for each suffix and scenario
get_smd_lists <- function(scenario_num, suffix) {
  get_slice <- function(varname) {
    full_list <- get(varname)
    total_len <- length(full_list)
    num_scenarios <- 4
    num_sims_per_scenario <- total_len / num_scenarios

    if (num_sims_per_scenario != floor(num_sims_per_scenario)) {
      stop("Number of simulations per scenario is not an integer. Check data length.")
    }

    start_idx <- (scenario_num - 1) * num_sims_per_scenario + 1
    end_idx <- scenario_num * num_sims_per_scenario
    full_list[start_idx:end_idx]
  }

  if (suffix == "_all") {
    rct_list <- list(
      get_slice(make_varname("rct", "age", suffix)),
      get_slice(make_varname("rct", "sex", suffix)),
      get_slice(make_varname("rct", "treat", suffix)),
      get_slice(make_varname("rct", "outcome", suffix))
    )
    ob_list <- list(
      get_slice(make_varname("ob", "age", suffix)),
      get_slice(make_varname("ob", "sex", suffix)),
      get_slice(make_varname("ob", "treat", suffix)),
      get_slice(make_varname("ob", "outcome", suffix))
    )
    ex_list <- list(
      get_slice(make_varname("ex", "age", suffix)),
      get_slice(make_varname("ex", "sex", suffix)),
      get_slice(make_varname("ex", "treat", suffix)),
      get_slice(make_varname("ex", "outcome", suffix))
    )
  } else if (suffix %in% c("_sim", "_syn")) {
    rct_list <- list(
      get_slice(make_varname("rct_ob", "age", suffix)),
      get_slice(make_varname("rct_ob", "sex", suffix)),
      get_slice(make_varname("rct_ob", "treat", suffix)),
      get_slice(make_varname("rct_ob", "outcome", suffix))
    )
    ob_list <- list(
      get_slice(make_varname("ob_ex", "age", suffix)),
      get_slice(make_varname("ob_ex", "sex", suffix)),
      get_slice(make_varname("ob_ex", "treat", suffix)),
      get_slice(make_varname("ob_ex", "outcome", suffix))
    )
    ex_list <- list(
      get_slice(make_varname("rct_ex", "age", suffix)),
      get_slice(make_varname("rct_ex", "sex", suffix)),
      get_slice(make_varname("rct_ex", "treat", suffix)),
      get_slice(make_varname("rct_ex", "outcome", suffix))
    )
  } else {
    stop("Suffix not recognized.")
  }

  list(RCT = rct_list, Observational = ob_list, External = ex_list)
}

# mean smd plot
create_mean_smd_plot <- function(scenario_num, suffix) {
  all_lists <- get_smd_lists(scenario_num, suffix)
  smd_mean <- do.call(rbind, lapply(names(all_lists), function(type_name) {
    vars_list <- all_lists[[type_name]]
    means <- sapply(vars_list, function(x) {
      v <- unlist(x)
      if (length(v) == 0) return(NA_real_)
      mean(v[is.finite(v)], na.rm = TRUE)
    })
    data.frame(label = variables, type = type_name, SMD = means, stringsAsFactors = FALSE)
  }))

  # drop rows with NA SMD
  smd_mean <- smd_mean[is.finite(smd_mean$SMD), ]
  smd_mean$type <- factor(smd_mean$type, levels = types)

  ggplot(smd_mean, aes(y = label, x = SMD, shape = type)) +
    geom_point(size = 3, position = position_dodgev(height = 0.5)) +
    ggtitle(paste0("SMD - ",
                   ifelse(suffix == "_all", "Original Data",
                          ifelse(suffix == "_sim", "Simulated Combos", "Synthetic Combos")),
                   " - Scenario ", scenario_num)) +
    xlab("Average Standardised Mean Difference (SMD)") +
    scale_shape_manual(values = c(0, 5, 2)) +
    ylab("Variable") +
    xlim(-0.5, 0.5) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"))
}

# color boxplot
create_boxplot_color <- function(scenario_num, suffix) {
  all_lists <- get_smd_lists(scenario_num, suffix)

  rows <- list()
  for (type_name in names(all_lists)) {
    vars_list <- all_lists[[type_name]]
    for (j in seq_along(vars_list)) {
      raw_vals <- unlist(vars_list[[j]])
      raw_vals <- raw_vals[is.finite(raw_vals)]
      if (length(raw_vals) == 0) next
      df_j <- data.frame(
        label = variables[j],
        type = type_name,
        SMD = raw_vals,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- df_j
    }
  }

  if (length(rows) == 0) {
    warning("No finite SMD values found for scenario ", scenario_num, " suffix ", suffix)
    return(ggplot() + ggtitle("No data"))
  }

  boxplot_df <- do.call(rbind, rows)
  boxplot_df$type <- factor(boxplot_df$type, levels = types)
  boxplot_df$label <- factor(boxplot_df$label, levels = variables)

  ggplot(boxplot_df, aes(x = label, y = SMD, fill = type)) +
    geom_boxplot(color = "black", alpha = 0.9, outlier.size = 1, coef = 1.5) +
    coord_cartesian(ylim = c(-1, 1)) +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.background = element_blank()
    ) +
    ggtitle(paste0("Standard Mean Difference by Data Type - ",
                   ifelse(suffix == "_all", "Original Data",
                          ifelse(suffix == "_sim", "Simulated Combos", "Synthetic Combos")),
                   " - Scenario ", scenario_num)) +
    xlab("Variables") +
    guides(fill = guide_legend(title = "Data Type"))
}

# run and save plots
scenarios <- 1:4
suffixes <- c("_all", "_sim", "_syn")
fig_dir <- "C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Aug/CR_CART"

for (suffix in suffixes) {
  mean_plots <- lapply(scenarios, create_mean_smd_plot, suffix = suffix)
  boxplots_color <- lapply(scenarios, create_boxplot_color, suffix = suffix)

  jpeg(file.path(fig_dir, paste0("AllScenarios_mean_smd_plot", suffix, ".jpg")), width = 1200, height = 800)
  grid.arrange(grobs = mean_plots, ncol = 2, nrow = 2)
  dev.off()

  jpeg(file.path(fig_dir, paste0("AllScenarios_smdbox_color", suffix, ".jpg")), width = 1200, height = 800)
  grid.arrange(grobs = boxplots_color, ncol = 2, nrow = 2)
  dev.off()
}












# create df from response list
make_resp_df <- function(resp_list, data_type, method) {
  data.frame(
    rate = unlist(resp_list),
    iteration = seq_along(unlist(resp_list)),
    data_type = data_type,
    method = method
  )
}

# combine all response data into one data frame
resp_df <- dplyr::bind_rows(
  make_resp_df(sim_resp_rct, "RCT", "Sim"),
  make_resp_df(syn_resp_rct, "RCT", "Syn"),
  make_resp_df(sim_resp_ob,  "OB",  "Sim"),
  make_resp_df(syn_resp_ob,  "OB",  "Syn"),
  make_resp_df(sim_resp_ex,  "EX",  "Sim"),
  make_resp_df(syn_resp_ex,  "EX",  "Syn")
)

# assign scenarios to resp df
resp_df <- resp_df %>%
  dplyr::mutate(
    scenario = ceiling(iteration / (nSim / 4))
  )

# summarise mean and se
summary_df <- resp_df %>%
  dplyr::group_by(data_type, method, scenario) %>%
  dplyr::summarise(
    mean_rate = mean(rate, na.rm = TRUE),
    se = sd(rate, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    scenario_label = factor(paste("Scenario", scenario),
                            levels = paste("Scenario", 1:4))
  )

# plot
p <- ggplot(summary_df, aes(x = method, y = mean_rate, color = data_type)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = mean_rate - se, ymax = mean_rate + se),
                position = position_dodge(width = 0.5), width = 0.2) +
  facet_wrap(~scenario_label, ncol = 1) +
  labs(title = "Control Response Rate by Method, Data Type, and Scenario",
       x = "Method",
       y = "Mean Control Response Rate",
       color = "Data Type") +
  scale_color_brewer(palette = "Dark2") +

  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    panel.spacing = unit(1.2, "lines")
  )

ggsave("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Aug/CR_CART/control_response_rate_means.jpg",
       p, width = 6, height = 12)














# get exact values for reporting:

# function to split response rate values into sep dfs:
split_to_dfs <- function(values, k = 4, prefix = "df") {
  # k is number of groups
  split_values <- split(values, cut(seq_along(values), k, labels = FALSE))

  # loops
  for (i in seq_along(split_values)) {
    df_name <- paste0(prefix, i)
    assign(df_name, data.frame(value = split_values[[i]]), envir = .GlobalEnv)
  }
}

split_to_dfs(sim_resp_rct, k=4)
sim_resp_rct_1 <- data.frame(value = as.vector(t(df1)))
sim_resp_rct_2 <- data.frame(value = as.vector(t(df2)))
sim_resp_rct_3 <- data.frame(value = as.vector(t(df3)))
sim_resp_rct_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(sim_resp_ob, k=4)
sim_resp_ob_1 <- data.frame(value = as.vector(t(df1)))
sim_resp_ob_2 <- data.frame(value = as.vector(t(df2)))
sim_resp_ob_3 <- data.frame(value = as.vector(t(df3)))
sim_resp_ob_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(sim_resp_ex, k=4)
sim_resp_ex_1 <- data.frame(value = as.vector(t(df1)))
sim_resp_ex_2 <- data.frame(value = as.vector(t(df2)))
sim_resp_ex_3 <- data.frame(value = as.vector(t(df3)))
sim_resp_ex_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(syn_resp_rct, k=4)
syn_resp_rct_1 <- data.frame(value = as.vector(t(df1)))
syn_resp_rct_2 <- data.frame(value = as.vector(t(df2)))
syn_resp_rct_3 <- data.frame(value = as.vector(t(df3)))
syn_resp_rct_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(syn_resp_ob, k=4)
syn_resp_ob_1 <- data.frame(value = as.vector(t(df1)))
syn_resp_ob_2 <- data.frame(value = as.vector(t(df2)))
syn_resp_ob_3 <- data.frame(value = as.vector(t(df3)))
syn_resp_ob_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(syn_resp_ex, k=4)
syn_resp_ex_1 <- data.frame(value = as.vector(t(df1)))
syn_resp_ex_2 <- data.frame(value = as.vector(t(df2)))
syn_resp_ex_3 <- data.frame(value = as.vector(t(df3)))
syn_resp_ex_4 <- data.frame(value = as.vector(t(df4)))



mean(sim_resp_rct_1$value)
mean(sim_resp_rct_2$value)
mean(sim_resp_rct_3$value)
mean(sim_resp_rct_4$value)
sd(sim_resp_rct_1$value)
sd(sim_resp_rct_2$value)
sd(sim_resp_rct_3$value)
sd(sim_resp_rct_4$value)
range(sim_resp_rct_1$value)
range(sim_resp_rct_2$value)
range(sim_resp_rct_3$value)
range(sim_resp_rct_4$value)

mean(sim_resp_ob_1$value)
mean(sim_resp_ob_2$value)
mean(sim_resp_ob_3$value)
mean(sim_resp_ob_4$value)
sd(sim_resp_ob_1$value)
sd(sim_resp_ob_2$value)
sd(sim_resp_ob_3$value)
sd(sim_resp_ob_4$value)
range(sim_resp_ob_1$value)
range(sim_resp_ob_2$value)
range(sim_resp_ob_3$value)
range(sim_resp_ob_4$value)

mean(sim_resp_ex_1$value)
mean(sim_resp_ex_2$value)
mean(sim_resp_ex_3$value)
mean(sim_resp_ex_4$value)
sd(sim_resp_ex_1$value)
sd(sim_resp_ex_2$value)
sd(sim_resp_ex_3$value)
sd(sim_resp_ex_4$value)
range(sim_resp_ex_1$value)
range(sim_resp_ex_2$value)
range(sim_resp_ex_3$value)
range(sim_resp_ex_4$value)




mean(syn_resp_rct_1$value)
mean(syn_resp_rct_2$value)
mean(syn_resp_rct_3$value)
mean(syn_resp_rct_4$value)
sd(syn_resp_rct_1$value)
sd(syn_resp_rct_2$value)
sd(syn_resp_rct_3$value)
sd(syn_resp_rct_4$value)
range(syn_resp_rct_1$value)
range(syn_resp_rct_2$value)
range(syn_resp_rct_3$value)
range(syn_resp_rct_4$value)

mean(syn_resp_ob_1$value)
mean(syn_resp_ob_2$value)
mean(syn_resp_ob_3$value)
mean(syn_resp_ob_4$value)
sd(syn_resp_ob_1$value)
sd(syn_resp_ob_2$value)
sd(syn_resp_ob_3$value)
sd(syn_resp_ob_4$value)
range(syn_resp_ob_1$value)
range(syn_resp_ob_2$value)
range(syn_resp_ob_3$value)
range(syn_resp_ob_4$value)

mean(syn_resp_ex_1$value)
mean(syn_resp_ex_2$value)
mean(syn_resp_ex_3$value)
mean(syn_resp_ex_4$value)
sd(syn_resp_ex_1$value)
sd(syn_resp_ex_2$value)
sd(syn_resp_ex_3$value)
sd(syn_resp_ex_4$value)
range(syn_resp_ex_1$value)
range(syn_resp_ex_2$value)
range(syn_resp_ex_3$value)
range(syn_resp_ex_4$value)


test_resp_rct_1 <- tsum.test(
  mean.x = mean(sim_resp_rct_1$value), s.x = sd(sim_resp_rct_1$value), n.x = length(sim_resp_rct_1$value),
  mean.y = mean(syn_resp_rct_1$value), s.y = sd(syn_resp_rct_1$value), n.y = length(syn_resp_rct_1$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_rct_2 <- tsum.test(
  mean.x = mean(sim_resp_rct_2$value), s.x = sd(sim_resp_rct_2$value), n.x = length(sim_resp_rct_2$value),
  mean.y = mean(syn_resp_rct_2$value), s.y = sd(syn_resp_rct_2$value), n.y = length(syn_resp_rct_2$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_rct_3 <- tsum.test(
  mean.x = mean(sim_resp_rct_3$value), s.x = sd(sim_resp_rct_3$value), n.x = length(sim_resp_rct_3$value),
  mean.y = mean(syn_resp_rct_3$value), s.y = sd(syn_resp_rct_3$value), n.y = length(syn_resp_rct_3$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_rct_4 <- tsum.test(
  mean.x = mean(sim_resp_rct_4$value), s.x = sd(sim_resp_rct_4$value), n.x = length(sim_resp_rct_4$value),
  mean.y = mean(syn_resp_rct_4$value), s.y = sd(syn_resp_rct_4$value), n.y = length(syn_resp_rct_4$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)




test_resp_ob_1 <- tsum.test(
  mean.x = mean(sim_resp_ob_1$value), s.x = sd(sim_resp_ob_1$value), n.x = length(sim_resp_ob_1$value),
  mean.y = mean(syn_resp_ob_1$value), s.y = sd(syn_resp_ob_1$value), n.y = length(syn_resp_ob_1$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ob_2 <- tsum.test(
  mean.x = mean(sim_resp_ob_2$value), s.x = sd(sim_resp_ob_2$value), n.x = length(sim_resp_ob_2$value),
  mean.y = mean(syn_resp_ob_2$value), s.y = sd(syn_resp_ob_2$value), n.y = length(syn_resp_ob_2$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ob_3 <- tsum.test(
  mean.x = mean(sim_resp_ob_3$value), s.x = sd(sim_resp_ob_3$value), n.x = length(sim_resp_ob_3$value),
  mean.y = mean(syn_resp_ob_3$value), s.y = sd(syn_resp_ob_3$value), n.y = length(syn_resp_ob_3$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ob_4 <- tsum.test(
  mean.x = mean(sim_resp_ob_4$value), s.x = sd(sim_resp_ob_4$value), n.x = length(sim_resp_ob_4$value),
  mean.y = mean(syn_resp_ob_4$value), s.y = sd(syn_resp_ob_4$value), n.y = length(syn_resp_ob_4$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)




test_resp_ex_1 <- tsum.test(
  mean.x = mean(sim_resp_ex_1$value), s.x = sd(sim_resp_ex_1$value), n.x = length(sim_resp_ex_1$value),
  mean.y = mean(syn_resp_ex_1$value), s.y = sd(syn_resp_ex_1$value), n.y = length(syn_resp_ex_1$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ex_2 <- tsum.test(
  mean.x = mean(sim_resp_ex_2$value), s.x = sd(sim_resp_ex_2$value), n.x = length(sim_resp_ex_2$value),
  mean.y = mean(syn_resp_ex_2$value), s.y = sd(syn_resp_ex_2$value), n.y = length(syn_resp_ex_2$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ex_3 <- tsum.test(
  mean.x = mean(sim_resp_ex_3$value), s.x = sd(sim_resp_ex_3$value), n.x = length(sim_resp_ex_3$value),
  mean.y = mean(syn_resp_ex_3$value), s.y = sd(syn_resp_ex_3$value), n.y = length(syn_resp_ex_3$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ex_4 <- tsum.test(
  mean.x = mean(sim_resp_ex_4$value), s.x = sd(sim_resp_ex_4$value), n.x = length(sim_resp_ex_4$value),
  mean.y = mean(syn_resp_ex_4$value), s.y = sd(syn_resp_ex_4$value), n.y = length(syn_resp_ex_4$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)






# Extract t-test results with p-values, confidence intervals, and effect sizes
ttest_pvalues <- data.frame(
  Test = c("RCT_1", "RCT_2", "RCT_3", "RCT_4", "OB_1", "OB_2", "OB_3", "OB_4", "EX_1", "EX_2", "EX_3", "EX_4"),
  P_Value = c(test_resp_rct_1$p.value, test_resp_rct_2$p.value, test_resp_rct_3$p.value, test_resp_rct_4$p.value,
              test_resp_ob_1$p.value, test_resp_ob_2$p.value, test_resp_ob_3$p.value, test_resp_ob_4$p.value,
              test_resp_ex_1$p.value, test_resp_ex_2$p.value, test_resp_ex_3$p.value, test_resp_ex_4$p.value),
  CI_Lower = c(test_resp_rct_1$conf.int[1], test_resp_rct_2$conf.int[1], test_resp_rct_3$conf.int[1], test_resp_rct_4$conf.int[1],
               test_resp_ob_1$conf.int[1], test_resp_ob_2$conf.int[1], test_resp_ob_3$conf.int[1], test_resp_ob_4$conf.int[1],
               test_resp_ex_1$conf.int[1], test_resp_ex_2$conf.int[1], test_resp_ex_3$conf.int[1], test_resp_ex_4$conf.int[1]),
  CI_Upper = c(test_resp_rct_1$conf.int[2], test_resp_rct_2$conf.int[2], test_resp_rct_3$conf.int[2], test_resp_rct_4$conf.int[2],
               test_resp_ob_1$conf.int[2], test_resp_ob_2$conf.int[2], test_resp_ob_3$conf.int[2], test_resp_ob_4$conf.int[2],
               test_resp_ex_1$conf.int[2], test_resp_ex_2$conf.int[2], test_resp_ex_3$conf.int[2], test_resp_ex_4$conf.int[2]),
  T_Statistic = c(test_resp_rct_1$statistic, test_resp_rct_2$statistic, test_resp_rct_3$statistic, test_resp_rct_4$statistic,
                  test_resp_ob_1$statistic, test_resp_ob_2$statistic, test_resp_ob_3$statistic, test_resp_ob_4$statistic,
                  test_resp_ex_1$statistic, test_resp_ex_2$statistic, test_resp_ex_3$statistic, test_resp_ex_4$statistic)
)
write.csv(ttest_pvalues, file = "CROHNS_ttest_pvalues_cart.csv", row.names = FALSE)




# for smd:
split_to_dfs(rct_age_smd_all, k=4)
rct_age_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_age_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_age_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_age_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(rct_sex_smd_all, k=4)
rct_sex_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_sex_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_sex_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_sex_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(rct_treat_smd_all, k=4)
rct_treat_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_treat_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_treat_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_treat_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(rct_outcome_smd_all, k=4)
rct_outcome_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_outcome_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_outcome_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_outcome_smd_4 <- data.frame(value = as.vector(t(df4)))


split_to_dfs(ob_age_smd_all, k=4)
ob_age_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_age_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_age_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_age_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ob_sex_smd_all, k=4)
ob_sex_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_sex_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_sex_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_sex_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ob_treat_smd_all, k=4)
ob_treat_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_treat_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_treat_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_treat_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ob_outcome_smd_all, k=4)
ob_outcome_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_outcome_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_outcome_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_outcome_smd_4 <- data.frame(value = as.vector(t(df4)))



split_to_dfs(ex_age_smd_all, k=4)
ex_age_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_age_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_age_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_age_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ex_sex_smd_all, k=4)
ex_sex_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_sex_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_sex_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_sex_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ex_treat_smd_all, k=4)
ex_treat_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_treat_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_treat_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_treat_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ex_outcome_smd_all, k=4)
ex_outcome_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_outcome_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_outcome_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_outcome_smd_4 <- data.frame(value = as.vector(t(df4)))



# just median, range, sd

range(rct_age_smd_1$value)
median(rct_age_smd_1$value)
sd(rct_age_smd_1$value)
range(rct_age_smd_2$value)
median(rct_age_smd_2$value)
sd(rct_age_smd_2$value)
range(rct_age_smd_3$value)
median(rct_age_smd_3$value)
sd(rct_age_smd_3$value)
range(rct_age_smd_4$value)
median(rct_age_smd_4$value)
sd(rct_age_smd_4$value)

range(rct_sex_smd_1$value)
median(rct_sex_smd_1$value)
sd(rct_sex_smd_1$value)
range(rct_sex_smd_2$value)
median(rct_sex_smd_2$value)
sd(rct_sex_smd_2$value)
range(rct_sex_smd_3$value)
median(rct_sex_smd_3$value)
sd(rct_sex_smd_3$value)
range(rct_sex_smd_4$value)
median(rct_sex_smd_4$value)
sd(rct_sex_smd_4$value)


range(rct_treat_smd_1$value)
median(rct_treat_smd_1$value)
sd(rct_treat_smd_1$value)
range(rct_treat_smd_2$value)
median(rct_treat_smd_2$value)
sd(rct_treat_smd_2$value)
range(rct_treat_smd_3$value)
median(rct_treat_smd_3$value)
sd(rct_treat_smd_3$value)
range(rct_treat_smd_4$value)
median(rct_treat_smd_4$value)
sd(rct_treat_smd_4$value)



range(rct_outcome_smd_1$value)
median(rct_outcome_smd_1$value)
sd(rct_outcome_smd_1$value)
range(rct_outcome_smd_2$value)
median(rct_outcome_smd_2$value)
sd(rct_outcome_smd_2$value)

range(rct_outcome_smd_3$value)
median(rct_outcome_smd_3$value)
sd(rct_outcome_smd_3$value)
range(rct_outcome_smd_4$value)
median(rct_outcome_smd_4$value)
sd(rct_outcome_smd_4$value)


range(ob_age_smd_1$value)
median(ob_age_smd_1$value)
sd(ob_age_smd_1$value)
range(ob_age_smd_2$value)
median(ob_age_smd_2$value)
sd(ob_age_smd_2$value)
range(ob_age_smd_3$value)
median(ob_age_smd_3$value)
sd(ob_age_smd_3$value)
range(ob_age_smd_4$value)
median(ob_age_smd_4$value)
sd(ob_age_smd_4$value)

range(ob_sex_smd_1$value)
median(ob_sex_smd_1$value)
sd(ob_sex_smd_1$value)
range(ob_sex_smd_2$value)
median(ob_sex_smd_2$value)
sd(ob_sex_smd_2$value)
range(ob_sex_smd_3$value)
median(ob_sex_smd_3$value)
sd(ob_sex_smd_3$value)
range(ob_sex_smd_4$value)
median(ob_sex_smd_4$value)
sd(ob_sex_smd_4$value)


range(ob_treat_smd_1$value)
median(ob_treat_smd_1$value)
sd(ob_treat_smd_1$value)
range(ob_treat_smd_2$value)
median(ob_treat_smd_2$value)
sd(ob_treat_smd_2$value)
range(ob_treat_smd_3$value)
median(ob_treat_smd_3$value)
sd(ob_treat_smd_3$value)
range(ob_treat_smd_4$value)
median(ob_treat_smd_4$value)
sd(ob_treat_smd_4$value)


range(ob_outcome_smd_1$value)
median(ob_outcome_smd_1$value)
sd(ob_outcome_smd_1$value)
range(ob_outcome_smd_2$value)
median(ob_outcome_smd_2$value)
sd(ob_outcome_smd_2$value)
range(ob_outcome_smd_3$value)
median(ob_outcome_smd_3$value)
sd(ob_outcome_smd_3$value)

range(ob_outcome_smd_4$value)
median(ob_outcome_smd_4$value)
sd(ob_outcome_smd_4$value)



range(ex_age_smd_1$value)
median(ex_age_smd_1$value)
sd(ex_age_smd_1$value)
range(ex_age_smd_2$value)
median(ex_age_smd_2$value)
sd(ex_age_smd_2$value)
range(ex_age_smd_3$value)
median(ex_age_smd_3$value)
sd(ex_age_smd_3$value)
range(ex_age_smd_4$value)
median(ex_age_smd_4$value)
sd(ex_age_smd_4$value)

range(ex_sex_smd_1$value)
median(ex_sex_smd_1$value)
sd(ex_sex_smd_1$value)
range(ex_sex_smd_2$value)
median(ex_sex_smd_2$value)
sd(ex_sex_smd_2$value)
range(ex_sex_smd_3$value)
median(ex_sex_smd_3$value)
sd(ex_sex_smd_3$value)
range(ex_sex_smd_4$value)
median(ex_sex_smd_4$value)
sd(ex_sex_smd_4$value)


range(ex_treat_smd_1$value)
median(ex_treat_smd_1$value)
sd(ex_treat_smd_1$value)
range(ex_treat_smd_2$value)
median(ex_treat_smd_2$value)
sd(ex_treat_smd_2$value)
range(ex_treat_smd_3$value)
median(ex_treat_smd_3$value)
sd(ex_treat_smd_3$value)
range(ex_treat_smd_4$value)
median(ex_treat_smd_4$value)
sd(ex_treat_smd_4$value)


range(ex_outcome_smd_1$value)
median(ex_outcome_smd_1$value)
sd(ex_outcome_smd_1$value)
range(ex_outcome_smd_2$value)
median(ex_outcome_smd_2$value)
sd(ex_outcome_smd_2$value)
range(ex_outcome_smd_3$value)
median(ex_outcome_smd_3$value)
sd(ex_outcome_smd_3$value)
range(ex_outcome_smd_4$value)
median(ex_outcome_smd_4$value)
sd(ex_outcome_smd_4$value)










































































# repeat for ll:

# updated sample sizes for better runtime:
# 100 moved to 150 to avoid issues with synthpop
# 20,000 shrunk to 500 to speed up runtime

# creating scenarios for number of observations
# halve the n due to pre and post being added together later
nsize <- function (scenario, datatype, simsyn) {
  if (scenario == 1) { # scenario one is n = 150 to n = 150
    N <- 75
  }   else if (scenario == 2 ) { # scenario two is n = 500 to n = 500
    N <- 250
  }   else if (scenario == 3 & datatype == "rct") { # scenario three is realistic sample sizes
    N <- 628
  }   else if (scenario == 3 & datatype == "ob"){
    N <- 221
  }   else if (scenario == 3 & datatype == "ex"){
    N <- 322
  }   else if (scenario == 4 & simsyn == "sim"){ # scenario four is n = 150 to n = 500
    N <- 75
  }   else if (scenario == 4 & simsyn == "syn"){
    N <- 250
  }
  return(N)
}



# function for smd calculation:
process_data <- function(df1, df2, var_names, table_list_name, smd_list_prefix, i) {
  # merge df
  df1$Key <- 1:nrow(df1)
  df2$Key <- 1:nrow(df2)
  merged_df <- merge(df1, df2, by = "Key", all = TRUE)
  merged_df <- merged_df[, -which(names(merged_df) == "Key")]

  # make a single df
  combined_df <- data.frame(
    lapply(var_names, function(var) {
      c(merged_df[[paste0(var, ".x")]], merged_df[[paste0(var, ".y")]])
    })
  )
  colnames(combined_df) <- var_names

  # save table
  temp <- get(table_list_name, envir = .GlobalEnv)
  temp[[i]] <- CreateTableOne(vars = var_names, data = combined_df, test = FALSE)
  assign(table_list_name, temp, envir = .GlobalEnv)

  # calculate smd
  for (var in var_names) {
    x <- merged_df[[paste0(var, ".x")]]
    y <- merged_df[[paste0(var, ".y")]]

    pooled_sd <- sqrt((sd(x, na.rm = TRUE)^2 + sd(y, na.rm = TRUE)^2) / 2)
    mean_diff <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
    smd_value <- mean_diff / pooled_sd

    smd_list_name <- paste0(smd_list_prefix, "_", var, "_smd_all")
    smd_list <- get(smd_list_name, envir = .GlobalEnv)
    smd_list[[i]] <- smd_value
    assign(smd_list_name, smd_list, envir = .GlobalEnv)
  }
}


# function for smd calculation across data types
compute_smd_across <- function(df1, df2, vars, prefix, suffix, i) {
  for (var in vars) {
    x1 <- df1[[paste0(var, ".", suffix)]]
    x2 <- df2[[paste0(var, ".", suffix)]]

    sd1 <- sd(x1, na.rm = TRUE)
    sd2 <- sd(x2, na.rm = TRUE)
    pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
    mean_diff <- mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)

    # by summix
    smd_list_name <- if (suffix == "x") {
      paste0(prefix, "_", var, "_smd_sim")
    } else if (suffix == "y") {
      paste0(prefix, "_", var, "_smd_syn")
    } else {
      stop("compute_smd_across: suffix must be 'x' or 'y'")
    }

    smd_list <- get(smd_list_name, envir = .GlobalEnv)
    smd_list[[i]] <- mean_diff / pooled_sd
    assign(smd_list_name, smd_list, envir = .GlobalEnv)
  }
}
# function for making empty variable lists
create_smd_lists <- function(non_combo_prefixes, combo_prefixes, variables, n) {

  for (pref in non_combo_prefixes) {
    for (var in variables) {
      list_name <- paste0(pref, "_", var, "_smd_all")
      assign(list_name, vector("list", length = n), envir = .GlobalEnv)
    }
  }
  # combos for sim/syn
  for (pref in combo_prefixes) {
    for (var in variables) {
      assign(paste0(pref, "_", var, "_smd_sim"), vector("list", length = n), envir = .GlobalEnv)
      assign(paste0(pref, "_", var, "_smd_syn"), vector("list", length = n), envir = .GlobalEnv)
    }
  }
}
# function for making empty table lists
create_table_lists <- function(names, n) {
  for (name in names) {
    assign(name, vector("list", length = n), envir = .GlobalEnv)
  }
}


# function for control response rate

control_response <- function(df, list_name, i) {
  num_responders <- df %>%
    filter(treat == 0, outcome == 1) %>%
    nrow()

  num_control <- df %>%
    filter(treat == 0) %>%
    nrow()

  rate <- ifelse(num_control > 0, num_responders / num_control, NA)

  tmp <- get(list_name, envir = .GlobalEnv)
  tmp[[i]] <- rate
  assign(list_name, tmp, envir = .GlobalEnv)
}


# for loop for sim and analysis

nSim <- 40 # number of simulations # small for testing
set.seed(9999)


#  make empty lists

vars <- c("age", "sex", "treat", "outcome")
non_combo_prefixes <- c("rct", "ob", "ex")                   # these get *_smd_all
combo_prefixes     <- c("rct_ob", "ob_ex", "rct_ex")         # these get *_smd_sim and *_smd_syn


create_smd_lists(non_combo_prefixes, combo_prefixes, vars, nSim)

create_table_lists(c("table_rct_data", "table_ob_data", "table_ex_data"), nSim)

sim_results_rct <- vector("list", length = nSim)
sim_results_ob  <- vector("list", length = nSim)
sim_results_ex  <- vector("list", length = nSim)
syn_results_rct <- vector("list", length = nSim)
syn_results_ob <- vector("list", length = nSim)
syn_results_ex <- vector("list", length = nSim)

sim_resp_rct <- vector("list", length = nSim)
sim_resp_ob  <- vector("list", length = nSim)
sim_resp_ex  <- vector("list", length = nSim)
syn_resp_rct <- vector("list", length = nSim)
syn_resp_ob <- vector("list", length = nSim)
syn_resp_ex <- vector("list", length = nSim)



# list for sum stats

# list for sum stats

sim_sum_rct_mean <- vector("list", length = nSim)
sim_sum_ob_mean  <- vector("list", length = nSim)
sim_sum_ex_mean  <- vector("list", length = nSim)
syn_sum_rct_mean <- vector("list", length = nSim)
syn_sum_ob_mean <- vector("list", length = nSim)
syn_sum_ex_mean <- vector("list", length = nSim)


sim_sum_rct_med <- vector("list", length = nSim)
sim_sum_ob_med  <- vector("list", length = nSim)
sim_sum_ex_med  <- vector("list", length = nSim)
syn_sum_rct_med <- vector("list", length = nSim)
syn_sum_ob_med <- vector("list", length = nSim)
syn_sum_ex_med <- vector("list", length = nSim)

sim_sum_rct_range <- vector("list", length = nSim)
sim_sum_ob_range  <- vector("list", length = nSim)
sim_sum_ex_range  <- vector("list", length = nSim)
syn_sum_rct_range <- vector("list", length = nSim)
syn_sum_ob_range <- vector("list", length = nSim)
syn_sum_ex_range <- vector("list", length = nSim)

sim_sum_rct_sd <- vector("list", length = nSim)
sim_sum_ob_sd  <- vector("list", length = nSim)
sim_sum_ex_sd  <- vector("list", length = nSim)
syn_sum_rct_sd <- vector("list", length = nSim)
syn_sum_ob_sd <- vector("list", length = nSim)
syn_sum_ex_sd <- vector("list", length = nSim)

# list for dist
sim_kurt_rct <- vector("list", length = nSim)
sim_skew_rct <- vector("list", length = nSim)

sim_kurt_ob <- vector("list", length = nSim)
sim_skew_ob <- vector("list", length = nSim)

sim_kurt_ex <- vector("list", length = nSim)
sim_skew_ex <- vector("list", length = nSim)

syn_kurt_rct <- vector("list", length = nSim)
syn_skew_rct <- vector("list", length = nSim)

syn_kurt_ob <- vector("list", length = nSim)
syn_skew_ob <- vector("list", length = nSim)

syn_kurt_ex <- vector("list", length = nSim)
syn_skew_ex <- vector("list", length = nSim)

scenario <- 1  #  scenario before the loop


# loop for gen simulations


for (i in 1:nSim) {

  if ((i - 1) %% (nSim/4) == 0 & i > 1) {  # change scenario every 10,000 iterations
    scenario <- scenario + 1
  }



  N <- nsize(scenario, "rct", "sim") # selecting scenario


  treat <- sample(c(1, 0), N, replace = TRUE, prob = c(0.63, 0.37))

  # compliance to treatment
  comply <- ifelse(treat == 1, sample(c(1, 0), N, replace = TRUE, prob = c(1, 0)), sample(c(1, 0), N, replace = TRUE, prob = c(0, 1)))

  # other variables
  sex <- rbinom(N, 1, 0.522)
  age <- rtruncnorm(N, a = 18, b = 72, mean = 37.3, sd = 11.8)
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
    mutate(treat_effect = ifelse(post == 0, NA, treat_effect))  # filter treatment effect to apply only to post-treatment group

  # calculate sex_effect
  df_combined <- df_combined %>%
    mutate(sex_effect = ifelse(sex == 1 | sex == 0, 0, NA)) # sex effect not recorded or reported / only baseline / assumed as zero

  # create overall treatment effect including demographic influences
  df_combined <- df_combined %>%
    mutate(overall_treat_effect = ifelse(post == 0, 0, treat_effect + sex_effect))

  # create pre and post outcomes (baseline vs treatment)
  df_combined <- df_combined %>%
    mutate(outcome = ifelse(post == 0, rbinom(nrow(df_combined[df_combined$post == 0, ]), 1, 0.359),
                            sapply(overall_treat_effect, function(prob) rbinom(1, 1, prob))))
  # store dfs in results
  sim_results_rct[[i]] <- df_combined

  # remove df2
  rm(df2)





  # observational
  # treatment

  N <- nsize(scenario, "ob", "sim") # selecting scenario

  treat=sample(c(1,0), N,replace=TRUE, prob=c(.66,.34)) # probability of being exposed to the treatment or not


  # other variables
  sex <- rbinom(N, 1, 0.602)
  age <- rtruncnorm(N, a = 16, b = 65, mean = 38.2, sd = 16.96)
  age <- round(age, 0) # round age to whole numbers
  random_var <- rbinom(N, 1, prob = 0.5)
  post <- rep(0, N)  # pre and post assignment

  # create df combining variables
  ob_df <- data.frame(treat, sex, age, random_var, post)

  # add id
  ob_df <- ob_df %>% dplyr::mutate(id = row_number())

  # create post df
  ob_df2 <- ob_df %>% mutate(post = 1)

  # combine pre and post
  ob_df_combined <- rbind(ob_df, ob_df2)

  # define treatment effect, only for post group
  # define treat
  ob_df_combined <- ob_df_combined %>% mutate(treat_effect=
                                                ifelse(treat==1, .387,
                                                       ifelse(treat==0, .242,
                                                              NA)))%>%
    mutate(treat_effect = ifelse(post == 0, NA, treat_effect))



  # calculate sex_effect and race_effect
  ob_df_combined <- ob_df_combined %>%
    mutate(sex_effect = ifelse(sex == 1 | sex == 0, 0, NA))

  # create overall treatment effect including demographic influences
  ob_df_combined <- ob_df_combined %>%
    mutate(overall_treat_effect = ifelse(post == 0, 0, treat_effect + sex_effect))


  # create pre and post outcomes
  ob_df_combined <- ob_df_combined %>%
    mutate(outcome = ifelse(post == 0, rbinom(sum(post == 0), 1, 0.242),
                            rbinom(sum(post == 1), 1, overall_treat_effect)))
  # store dfs in results
  sim_results_ob[[i]] <- ob_df_combined

  # remove df2
  rm(ob_df2)




  # external
  # treatment
  N <- nsize(scenario, "ex", "sim") # selecting scenario

  outcome=sample(c(1,0), N*2,replace=TRUE, prob=c(.64,.36))

  # other variables
  treat=sample(c(1,0), N*2,replace=TRUE, prob=c(.88,.12))
  sex <- rbinom(N, 1, 0.48)
  age <- rtruncnorm(N*2, a = 15, b = 70, mean = 38, sd = 17.78)
  age <- round(age, 0) # round age to whole numbers
  random_var <- rbinom(N*2, 1, prob = 0.5)



  # create df combining variables
  ex_df <- data.frame(outcome, treat, sex, age, random_var)

  # add id
  ex_df <- ex_df %>% dplyr::mutate(id = row_number())

  # store dfs in results
  sim_results_ex[[i]] <- ex_df











  # begin synthetic data creation:
  # MODEL - rct


  rct_mydata_con1 = subset(df_combined, select = c("id", "age"))

  # logreg requires no NA values
  rct_mydata_con1 <- na.omit(rct_mydata_con1)

  N <- nsize(scenario, "rct", "syn") # selecting scenario

  # minimum number of observations needed is 10
  rct_mysyn_con1 <- syn(rct_mydata_con1, method = "norm", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)


  rct_mydata_cat2 = subset(df_combined, select = c( "sex", "treat", "outcome"))
  rct_mydata_cat2 <- na.omit(rct_mydata_cat2)

  N <- nsize(scenario, "rct", "syn") # selecting scenario

  # minimum number of observations needed is 10
  rct_mysyn_cat2 <- syn(rct_mydata_cat2, method = "logreg", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)
  summary(rct_mysyn_cat2)

  # combining the two dataframes
  rct_mysyn_con <- rct_mysyn_con1$syn
  rct_mysyn_cat <- rct_mysyn_cat2$syn
  # create a common key column in each dataframe
  rct_mysyn_con$Key <- 1:nrow(rct_mysyn_con)
  rct_mysyn_cat$Key <- 1:nrow(rct_mysyn_cat)
  # merge the dataframes based on the common key column
  rct_mysyn <- merge(rct_mysyn_con, rct_mysyn_cat, by = "Key", all = TRUE)
  # drop key

  rct_mysyn = subset(rct_mysyn, select = c("id", "age",  "sex", "treat", "outcome"))


  # store dfs in results
  syn_results_rct[[i]] <- rct_mysyn






  # MODEL - ob
  ob_mydata_con = subset(ob_df_combined, select = c("id", "age"))
  # logreg requires no NA values
  ob_mydata_con <- na.omit(ob_mydata_con)


  N <- nsize(scenario, "ob", "syn") # selecting scenario

  # minimum number of observations needed is 10
  ob_mysyn_con <- syn(ob_mydata_con, method = "norm", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)

  ob_mydata_cat = subset(ob_df_combined, select = c( "sex", "treat", "outcome"))
  # logreg requires no NA values
  ob_mydata_cat <- na.omit(ob_mydata_cat)




  N <- nsize(scenario, "ob", "syn") # selecting scenario

  # minimum number of observations needed is 10
  ob_mysyn_cat <- syn(ob_mydata_cat, method = "logreg", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)


  # combining the two dataframes

  ob_mysyn_con <- ob_mysyn_con$syn
  ob_mysyn_cat <- ob_mysyn_cat$syn

  # create a common key column in each dataframe
  ob_mysyn_con$Key <- 1:nrow(ob_mysyn_con)
  ob_mysyn_cat$Key <- 1:nrow(ob_mysyn_cat)

  # merge the dataframes based on the common key column
  ob_mysyn <- merge(ob_mysyn_con, ob_mysyn_cat, by = "Key", all = TRUE)

  # drop key
  ob_mysyn = subset(ob_mysyn, select = c("id", "age",  "sex", "treat", "outcome"))

  # store dfs in results
  syn_results_ob[[i]] <- ob_mysyn






  # MODEL - ex
  ex_mydata_con = subset(ex_df, select = c("id", "age"))
  # logreg requires no NA values
  ex_mydata_con <- na.omit(ex_mydata_con)



  N <- nsize(scenario, "ex", "syn") # selecting scenario

  # minimum number of observations needed is 10
  ex_mysyn_con <- syn(ex_mydata_con, method = "norm", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)

  summary(ex_mysyn_con)



  ex_mydata_cat = subset(ex_df, select = c( "sex", "treat", "outcome"))
  # logreg requires no NA values
  ex_mydata_cat <- na.omit(ex_mydata_cat)


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
  ex_mysyn = subset(ex_mysyn, select = c("id", "age",  "sex", "treat", "outcome"))



  # store dfs in results
  syn_results_ex[[i]] <- ex_mysyn




  # dist stats
  sim_kurt_rct[[i]] <- data.frame(age = kurtosis(df_combined$age),
                                  sex = kurtosis(df_combined$sex),
                                  treat = kurtosis(df_combined$treat),
                                  outcome = kurtosis(df_combined$outcome))


  sim_kurt_ob[[i]] <- data.frame(age = kurtosis(ob_df_combined$age),
                                 sex = kurtosis(ob_df_combined$sex),
                                 treat = kurtosis(ob_df_combined$treat),
                                 outcome = kurtosis(ob_df_combined$outcome))



  sim_kurt_ex[[i]] <- data.frame(age = kurtosis(ex_df$age),
                                 sex = kurtosis(ex_df$sex),
                                 treat = kurtosis(ex_df$treat),
                                 outcome = kurtosis(ex_df$outcome))






  sim_skew_rct[[i]] <- data.frame(age = skewness(df_combined$age),
                                  sex = skewness(df_combined$sex),
                                  treat = skewness(df_combined$treat),
                                  outcome = skewness(df_combined$outcome))


  sim_skew_ob[[i]] <- data.frame(age = skewness(ob_df_combined$age),
                                 sex = skewness(ob_df_combined$sex),
                                 treat = skewness(ob_df_combined$treat),
                                 outcome = skewness(ob_df_combined$outcome))



  sim_skew_ex[[i]] <- data.frame(age = skewness(ex_df$age),
                                 sex = skewness(ex_df$sex),
                                 treat = skewness(ex_df$treat),
                                 outcome = skewness(ex_df$outcome))





  syn_kurt_rct[[i]] <- data.frame(age = kurtosis(rct_mysyn$age),
                                  sex = kurtosis(rct_mysyn$sex),
                                  treat = kurtosis(rct_mysyn$treat),
                                  outcome = kurtosis(rct_mysyn$outcome))


  syn_kurt_ob[[i]] <- data.frame(age = kurtosis(ob_mysyn$age),
                                 sex = kurtosis(ob_mysyn$sex),
                                 treat = kurtosis(ob_mysyn$treat),
                                 outcome = kurtosis(ob_mysyn$outcome))



  syn_kurt_ex[[i]] <- data.frame(age = kurtosis(ex_mysyn$age),
                                 sex = kurtosis(ex_mysyn$sex),
                                 treat = kurtosis(ex_mysyn$treat),
                                 outcome = kurtosis(ex_mysyn$outcome))






  syn_skew_rct[[i]] <- data.frame(age = skewness(rct_mysyn$age),
                                  sex = skewness(rct_mysyn$sex),
                                  treat = skewness(rct_mysyn$treat),
                                  outcome = skewness(rct_mysyn$outcome))


  syn_skew_ob[[i]] <- data.frame(age = skewness(ob_mysyn$age),
                                 sex = skewness(ob_mysyn$sex),
                                 treat = skewness(ob_mysyn$treat),
                                 outcome = skewness(ob_mysyn$outcome))



  syn_skew_ex[[i]] <- data.frame(age = skewness(ex_mysyn$age),
                                 sex = skewness(ex_mysyn$sex),
                                 treat = skewness(ex_mysyn$treat),
                                 outcome = skewness(ex_mysyn$outcome))





  # summary stats

  # mean, median, range, sd,

  sim_sum_rct_mean[[i]] <- data.frame(age = mean(df_combined$age),
                                      sex = mean(df_combined$sex),
                                      treat = mean(df_combined$treat),
                                      outcome = mean(df_combined$outcome))


  sim_sum_rct_med[[i]] <- data.frame(age = median(df_combined$age),
                                     sex = median(df_combined$sex),
                                     treat = median(df_combined$treat),
                                     outcome = median(df_combined$outcome))

  sim_sum_rct_range[[i]] <- data.frame(age = range(df_combined$age),
                                       sex = range(df_combined$sex),
                                       treat = range(df_combined$treat),
                                       outcome = range(df_combined$outcome))


  sim_sum_rct_sd[[i]] <- data.frame(age = sd(df_combined$age),
                                    sex = sd(df_combined$sex),
                                    treat = sd(df_combined$treat),
                                    outcome = sd(df_combined$outcome))






  sim_sum_ob_mean[[i]]  <- data.frame(age = mean(ob_df_combined$age),
                                      sex = mean(ob_df_combined$sex),
                                      treat = mean(ob_df_combined$treat),
                                      outcome = mean(ob_df_combined$outcome))


  sim_sum_ob_med[[i]]  <- data.frame(age = median(ob_df_combined$age),
                                     sex = median(ob_df_combined$sex),
                                     treat = median(ob_df_combined$treat),
                                     outcome = median(ob_df_combined$outcome))


  sim_sum_ob_range[[i]]  <- data.frame(age = range(ob_df_combined$age),
                                       sex = range(ob_df_combined$sex),
                                       treat = range(ob_df_combined$treat),
                                       outcome = range(ob_df_combined$outcome))



  sim_sum_ob_sd[[i]]  <- data.frame(age = sd(ob_df_combined$age),
                                    sex = sd(ob_df_combined$sex),
                                    treat = sd(ob_df_combined$treat),
                                    outcome = sd(ob_df_combined$outcome))





  sim_sum_ex_mean[[i]]  <- data.frame(age = mean(ex_df$age),
                                      sex = mean(ex_df$sex),
                                      treat = mean(ex_df$treat),
                                      outcome = mean(ex_df$outcome))


  sim_sum_ex_med[[i]]  <- data.frame(age = median(ex_df$age),
                                     sex = median(ex_df$sex),
                                     treat = median(ex_df$treat),
                                     outcome = median(ex_df$outcome))

  sim_sum_ex_range[[i]]  <- data.frame(age = range(ex_df$age),
                                       sex = range(ex_df$sex),
                                       treat = range(ex_df$treat),
                                       outcome = range(ex_df$outcome))

  sim_sum_ex_sd[[i]]  <- data.frame(age = sd(ex_df$age),
                                    sex = sd(ex_df$sex),
                                    treat = sd(ex_df$treat),
                                    outcome = sd(ex_df$outcome))











  syn_sum_rct_mean[[i]] <- data.frame(age = mean(rct_mysyn$age),
                                      sex = mean(rct_mysyn$sex),
                                      treat = mean(rct_mysyn$treat),
                                      outcome = mean(rct_mysyn$outcome))


  syn_sum_rct_med[[i]] <- data.frame(age = median(rct_mysyn$age),
                                     sex = median(rct_mysyn$sex),
                                     treat = median(rct_mysyn$treat),
                                     outcome = median(rct_mysyn$outcome))

  syn_sum_rct_range[[i]] <- data.frame(age = range(rct_mysyn$age),
                                       sex = range(rct_mysyn$sex),
                                       treat = range(rct_mysyn$treat),
                                       outcome = range(rct_mysyn$outcome))



  syn_sum_rct_sd[[i]] <- data.frame(age = sd(rct_mysyn$age),
                                    sex = sd(rct_mysyn$sex),
                                    treat = sd(rct_mysyn$treat),
                                    outcome = sd(rct_mysyn$outcome))






  syn_sum_ob_mean[[i]] <- data.frame(age = mean(ob_mysyn$age),
                                     sex = mean(ob_mysyn$sex),
                                     treat = mean(ob_mysyn$treat),
                                     outcome = mean(ob_mysyn$outcome))


  syn_sum_ob_med[[i]] <- data.frame(age = median(ob_mysyn$age),
                                    sex = median(ob_mysyn$sex),
                                    treat = median(ob_mysyn$treat),
                                    outcome = median(ob_mysyn$outcome))

  syn_sum_ob_range[[i]] <- data.frame(age = range(ob_mysyn$age),
                                      sex = range(ob_mysyn$sex),
                                      treat = range(ob_mysyn$treat),
                                      outcome = range(ob_mysyn$outcome))



  syn_sum_ob_sd[[i]] <- data.frame(age = sd(ob_mysyn$age),
                                   sex = sd(ob_mysyn$sex),
                                   treat = sd(ob_mysyn$treat),
                                   outcome = sd(ob_mysyn$outcome))








  syn_sum_ex_mean[[i]] <- data.frame(age = mean(ex_mysyn$age),
                                     sex = mean(ex_mysyn$sex),
                                     treat = mean(ex_mysyn$treat),
                                     outcome = mean(ex_mysyn$outcome))


  syn_sum_ex_med[[i]] <- data.frame(age = median(ex_mysyn$age),
                                    sex = median(ex_mysyn$sex),
                                    treat = median(ex_mysyn$treat),
                                    outcome = median(ex_mysyn$outcome))

  syn_sum_ex_range[[i]] <- data.frame(age = range(ex_mysyn$age),
                                      sex = range(ex_mysyn$sex),
                                      treat = range(ex_mysyn$treat),
                                      outcome = range(ex_mysyn$outcome))



  syn_sum_ex_sd[[i]] <- data.frame(age = sd(ex_mysyn$age),
                                   sex = sd(ex_mysyn$sex),
                                   treat = sd(ex_mysyn$treat),
                                   outcome = sd(ex_mysyn$outcome))








  # calculating smd
  process_data(df_combined, rct_mysyn, vars, "table_rct_data", "rct", i)
  process_data(ob_df_combined, ob_mysyn, vars, "table_ob_data", "ob", i)
  process_data(ex_df, ex_mysyn, vars, "table_ex_data", "ex", i)

  # Synthetic vs Synthetic
  rct_ex <- merge(rct_mysyn, ex_mysyn, by = "id", suffixes = c(".x", ".y"))
  rct_ob <- merge(rct_mysyn, ob_mysyn, by = "id", suffixes = c(".x", ".y"))
  ob_ex  <- merge(ob_mysyn, ex_mysyn, by = "id", suffixes = c(".x", ".y"))

  compute_smd_across(rct_ex, rct_ex, vars, "rct_ex", "y", i)
  compute_smd_across(rct_ob, rct_ob, vars, "rct_ob", "y", i)
  compute_smd_across(ob_ex,  ob_ex,  vars, "ob_ex",  "y", i)

  compute_smd_across(rct_ex, rct_ex, vars, "rct_ex", "x", i)
  compute_smd_across(rct_ob, rct_ob, vars, "rct_ob", "x", i)
  compute_smd_across(ob_ex,  ob_ex,  vars, "ob_ex",  "x", i)



  # calculate control response rate:

  control_response(rct_mysyn, "syn_resp_rct", i)
  control_response(ob_mysyn,  "syn_resp_ob",  i)
  control_response(ex_mysyn,  "syn_resp_ex",  i)

  control_response(df_combined,     "sim_resp_rct", i)
  control_response(ob_df_combined,  "sim_resp_ob",  i)
  control_response(ex_df,           "sim_resp_ex",  i)
}









# ANALYSIS AND GRAPHING



# function to split lists
distribute_into_lists <- function(input_list, sizes) {
  if (sum(sizes) != length(input_list)) {
    stop("The sum of sizes must equal the length of the input list")
  }

  result <- vector("list", length(sizes))
  start_index <- 1

  for (i in seq_along(sizes)) {
    end_index <- start_index + sizes[i] - 1
    result[[i]] <- input_list[start_index:end_index]
    start_index <- end_index + 1
  }

  return(result)
}


list_size <- c(nSim/4, nSim/4, nSim/4, nSim/4)
save_dir <- "C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Saved_Rdata/August/CR_LL/"


# save data:

# non-combos
non_combos <- c("rct", "ob", "ex")
variables <- c("age", "sex", "treat", "outcome")
non_combo_names <- unlist(lapply(non_combos, function(pref) {
  paste0(pref, "_", variables, "_smd_all")
}))

# combos
combos <- c("rct_ex", "rct_ob", "ob_ex")
combo_names <- unlist(lapply(combos, function(pref) {
  paste0(rep(pref, each = length(variables)), "_", variables, "_smd_", rep(c("sim", "syn"), each = length(variables) * 1))
}))

# misc
extra_objects <- c(
  "sim_results_rct", "sim_results_ob", "sim_results_ex",
  "syn_results_rct", "syn_results_ob", "syn_results_ex",
  "sim_resp_rct", "sim_resp_ob", "sim_resp_ex",
  "syn_resp_rct", "syn_resp_ob", "syn_resp_ex",
  "table_rct_data", "table_ob_data", "table_ex_data"
)

# final list of all object names
obj_names <- c(extra_objects, non_combo_names, combo_names)

# save all
for (obj_name in obj_names) {

  if (!exists(obj_name)) {
    warning(paste("Object not found:", obj_name))
    next
  }

  obj <- get(obj_name)

  # check for list
  if (!is.list(obj)) {
    warning(paste("Object is not a list:", obj_name))
    next
  }

  obj_split <- distribute_into_lists(obj, list_size)

  for (scenario_index in seq_along(obj_split)) {
    saveRDS(
      obj_split[[scenario_index]],
      file = paste0(save_dir, obj_name, "_scenario", scenario_index, ".Rdata")
    )
  }
}
# dist stats

sim_kurt_rct = do.call(rbind, sim_kurt_rct)
sim_kurt_ob = do.call(rbind,sim_kurt_ob)
sim_kurt_ex = do.call(rbind, sim_kurt_ex)
sim_skew_rct = do.call(rbind, sim_skew_rct)
sim_skew_ob = do.call(rbind,sim_skew_ob)
sim_skew_ex = do.call(rbind, sim_skew_ex)

syn_kurt_rct = do.call(rbind, syn_kurt_rct)
syn_kurt_ob = do.call(rbind,syn_kurt_ob)
syn_kurt_ex = do.call(rbind, syn_kurt_ex)
syn_skew_rct = do.call(rbind, syn_skew_rct)
syn_skew_ob = do.call(rbind,syn_skew_ob)
syn_skew_ex = do.call(rbind, syn_skew_ex)

write.csv(sim_kurt_rct, file = "cr_sim_kurt_rct_ll.csv", row.names = FALSE)
write.csv(sim_kurt_ob, file = "cr_sim_kurt_ob_ll.csv", row.names = FALSE)
write.csv(sim_kurt_ex, file = "cr_sim_kurt_ex_ll.csv", row.names = FALSE)
write.csv(sim_skew_rct, file = "cr_sim_skew_rct_ll.csv", row.names = FALSE)
write.csv(sim_skew_ob, file = "cr_sim_skew_ob_ll.csv", row.names = FALSE)
write.csv(sim_skew_ex, file = "cr_sim_skew_ex_ll.csv", row.names = FALSE)

write.csv(syn_kurt_rct, file = "cr_syn_kurt_rct_ll.csv", row.names = FALSE)
write.csv(syn_kurt_ob, file = "cr_syn_kurt_ob_ll.csv", row.names = FALSE)
write.csv(syn_kurt_ex, file = "cr_syn_kurt_ex_ll.csv", row.names = FALSE)
write.csv(syn_skew_rct, file = "cr_syn_skew_rct_ll.csv", row.names = FALSE)
write.csv(syn_skew_ob, file = "cr_syn_skew_ob_ll.csv", row.names = FALSE)
write.csv(syn_skew_ex, file = "cr_syn_skew_ex_ll.csv", row.names = FALSE)


sim_sum_rct_mean = do.call(rbind, sim_sum_rct_mean)
sim_sum_rct_med = do.call(rbind, sim_sum_rct_med)
sim_sum_rct_range = do.call(rbind, sim_sum_rct_range)
sim_sum_rct_sd = do.call(rbind, sim_sum_rct_sd)

sim_sum_ob_mean = do.call(rbind, sim_sum_ob_mean)
sim_sum_ob_med = do.call(rbind, sim_sum_ob_med)
sim_sum_ob_range = do.call(rbind, sim_sum_ob_range)
sim_sum_ob_sd = do.call(rbind, sim_sum_ob_sd)


sim_sum_ex_mean = do.call(rbind, sim_sum_ex_mean)
sim_sum_ex_med = do.call(rbind, sim_sum_ex_med)
sim_sum_ex_range = do.call(rbind, sim_sum_ex_range)
sim_sum_ex_sd = do.call(rbind, sim_sum_ex_sd)



syn_sum_rct_mean = do.call(rbind, syn_sum_rct_mean)
syn_sum_rct_med = do.call(rbind, syn_sum_rct_med)
syn_sum_rct_range = do.call(rbind, syn_sum_rct_range)
syn_sum_rct_sd = do.call(rbind, syn_sum_rct_sd)

syn_sum_ob_mean = do.call(rbind, syn_sum_ob_mean)
syn_sum_ob_med = do.call(rbind, syn_sum_ob_med)
syn_sum_ob_range = do.call(rbind, syn_sum_ob_range)
syn_sum_ob_sd = do.call(rbind, syn_sum_ob_sd)


syn_sum_ex_mean = do.call(rbind, syn_sum_ex_mean)
syn_sum_ex_med = do.call(rbind, syn_sum_ex_med)
syn_sum_ex_range = do.call(rbind, syn_sum_ex_range)
syn_sum_ex_sd = do.call(rbind, syn_sum_ex_sd)



write.csv(sim_sum_rct_mean, file = "cr_sim_sum_rct_ll_mean.csv", row.names = FALSE)
write.csv(sim_sum_rct_med, file = "cr_sim_sum_rct_ll_med.csv", row.names = FALSE)
write.csv(sim_sum_rct_range, file = "cr_sim_sum_rct_ll_range.csv", row.names = FALSE)
write.csv(sim_sum_rct_sd, file = "cr_sim_sum_rct_ll_sd.csv", row.names = FALSE)

write.csv(sim_sum_ob_mean, file = "cr_sim_sum_ob_ll_mean.csv", row.names = FALSE)
write.csv(sim_sum_ob_med, file = "cr_sim_sum_ob_ll_med.csv", row.names = FALSE)
write.csv(sim_sum_ob_range, file = "cr_sim_sum_ob_ll_range.csv", row.names = FALSE)
write.csv(sim_sum_ob_sd, file = "cr_sim_sum_ob_ll_sd.csv", row.names = FALSE)

write.csv(sim_sum_ex_mean, file = "cr_sim_sum_ex_ll_mean.csv", row.names = FALSE)
write.csv(sim_sum_ex_med, file = "cr_sim_sum_ex_ll_med.csv", row.names = FALSE)
write.csv(sim_sum_ex_range, file = "cr_sim_sum_ex_ll_range.csv", row.names = FALSE)
write.csv(sim_sum_ex_sd, file = "cr_sim_sum_ex_ll_sd.csv", row.names = FALSE)


write.csv(syn_sum_rct_mean, file = "cr_syn_sum_rct_ll_mean.csv", row.names = FALSE)
write.csv(syn_sum_rct_med, file = "cr_syn_sum_rct_ll_med.csv", row.names = FALSE)
write.csv(syn_sum_rct_range, file = "cr_syn_sum_rct_ll_range.csv", row.names = FALSE)
write.csv(syn_sum_rct_sd, file = "cr_syn_sum_rct_ll_sd.csv", row.names = FALSE)

write.csv(syn_sum_ob_mean, file = "cr_syn_sum_ob_ll_mean.csv", row.names = FALSE)
write.csv(syn_sum_ob_med, file = "cr_syn_sum_ob_ll_med.csv", row.names = FALSE)
write.csv(syn_sum_ob_range, file = "cr_syn_sum_ob_ll_range.csv", row.names = FALSE)
write.csv(syn_sum_ob_sd, file = "cr_syn_sum_ob_ll_sd.csv", row.names = FALSE)

write.csv(syn_sum_ex_mean, file = "cr_syn_sum_ex_ll_mean.csv", row.names = FALSE)
write.csv(syn_sum_ex_med, file = "cr_syn_sum_ex_ll_med.csv", row.names = FALSE)
write.csv(syn_sum_ex_range, file = "cr_syn_sum_ex_ll_range.csv", row.names = FALSE)
write.csv(syn_sum_ex_sd, file = "cr_syn_sum_ex_ll_sd.csv", row.names = FALSE)





# plotting smd figures

variables <- c("Age", "Sex", "Treatment", "Outcome")
types <- c("RCT", "Observational", "External")

# build variable names
make_varname <- function(prefix, var, suffix = "") {
  paste0(prefix, "_", tolower(var), "_smd", suffix)
}

# get relevant lists for each suffix and scenario
get_smd_lists <- function(scenario_num, suffix) {
  get_slice <- function(varname) {
    full_list <- get(varname)
    total_len <- length(full_list)
    num_scenarios <- 4
    num_sims_per_scenario <- total_len / num_scenarios

    if (num_sims_per_scenario != floor(num_sims_per_scenario)) {
      stop("Number of simulations per scenario is not an integer. Check data length.")
    }

    start_idx <- (scenario_num - 1) * num_sims_per_scenario + 1
    end_idx <- scenario_num * num_sims_per_scenario
    full_list[start_idx:end_idx]
  }

  if (suffix == "_all") {
    rct_list <- list(
      get_slice(make_varname("rct", "age", suffix)),
      get_slice(make_varname("rct", "sex", suffix)),
      get_slice(make_varname("rct", "treat", suffix)),
      get_slice(make_varname("rct", "outcome", suffix))
    )
    ob_list <- list(
      get_slice(make_varname("ob", "age", suffix)),
      get_slice(make_varname("ob", "sex", suffix)),
      get_slice(make_varname("ob", "treat", suffix)),
      get_slice(make_varname("ob", "outcome", suffix))
    )
    ex_list <- list(
      get_slice(make_varname("ex", "age", suffix)),
      get_slice(make_varname("ex", "sex", suffix)),
      get_slice(make_varname("ex", "treat", suffix)),
      get_slice(make_varname("ex", "outcome", suffix))
    )
  } else if (suffix %in% c("_sim", "_syn")) {
    rct_list <- list(
      get_slice(make_varname("rct_ob", "age", suffix)),
      get_slice(make_varname("rct_ob", "sex", suffix)),
      get_slice(make_varname("rct_ob", "treat", suffix)),
      get_slice(make_varname("rct_ob", "outcome", suffix))
    )
    ob_list <- list(
      get_slice(make_varname("ob_ex", "age", suffix)),
      get_slice(make_varname("ob_ex", "sex", suffix)),
      get_slice(make_varname("ob_ex", "treat", suffix)),
      get_slice(make_varname("ob_ex", "outcome", suffix))
    )
    ex_list <- list(
      get_slice(make_varname("rct_ex", "age", suffix)),
      get_slice(make_varname("rct_ex", "sex", suffix)),
      get_slice(make_varname("rct_ex", "treat", suffix)),
      get_slice(make_varname("rct_ex", "outcome", suffix))
    )
  } else {
    stop("Suffix not recognized.")
  }

  list(RCT = rct_list, Observational = ob_list, External = ex_list)
}

# mean smd plot
create_mean_smd_plot <- function(scenario_num, suffix) {
  all_lists <- get_smd_lists(scenario_num, suffix)
  smd_mean <- do.call(rbind, lapply(names(all_lists), function(type_name) {
    vars_list <- all_lists[[type_name]]
    means <- sapply(vars_list, function(x) {
      v <- unlist(x)
      if (length(v) == 0) return(NA_real_)
      mean(v[is.finite(v)], na.rm = TRUE)
    })
    data.frame(label = variables, type = type_name, SMD = means, stringsAsFactors = FALSE)
  }))

  # drop rows with NA SMD
  smd_mean <- smd_mean[is.finite(smd_mean$SMD), ]
  smd_mean$type <- factor(smd_mean$type, levels = types)

  ggplot(smd_mean, aes(y = label, x = SMD, shape = type)) +
    geom_point(size = 3, position = position_dodgev(height = 0.5)) +
    ggtitle(paste0("SMD - ",
                   ifelse(suffix == "_all", "Original Data",
                          ifelse(suffix == "_sim", "Simulated Combos", "Synthetic Combos")),
                   " - Scenario ", scenario_num)) +
    xlab("Average Standardised Mean Difference (SMD)") +
    scale_shape_manual(values = c(0, 5, 2)) +
    ylab("Variable") +
    xlim(-0.5, 0.5) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"))
}

# color boxplot
create_boxplot_color <- function(scenario_num, suffix) {
  all_lists <- get_smd_lists(scenario_num, suffix)

  rows <- list()
  for (type_name in names(all_lists)) {
    vars_list <- all_lists[[type_name]]
    for (j in seq_along(vars_list)) {
      raw_vals <- unlist(vars_list[[j]])
      raw_vals <- raw_vals[is.finite(raw_vals)]
      if (length(raw_vals) == 0) next
      df_j <- data.frame(
        label = variables[j],
        type = type_name,
        SMD = raw_vals,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- df_j
    }
  }

  if (length(rows) == 0) {
    warning("No finite SMD values found for scenario ", scenario_num, " suffix ", suffix)
    return(ggplot() + ggtitle("No data"))
  }

  boxplot_df <- do.call(rbind, rows)
  boxplot_df$type <- factor(boxplot_df$type, levels = types)
  boxplot_df$label <- factor(boxplot_df$label, levels = variables)

  ggplot(boxplot_df, aes(x = label, y = SMD, fill = type)) +
    geom_boxplot(color = "black", alpha = 0.9, outlier.size = 1, coef = 1.5) +
    coord_cartesian(ylim = c(-1, 1)) +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.background = element_blank()
    ) +
    ggtitle(paste0("Standard Mean Difference by Data Type - ",
                   ifelse(suffix == "_all", "Original Data",
                          ifelse(suffix == "_sim", "Simulated Combos", "Synthetic Combos")),
                   " - Scenario ", scenario_num)) +
    xlab("Variables") +
    guides(fill = guide_legend(title = "Data Type"))
}

# run and save plots
scenarios <- 1:4
suffixes <- c("_all", "_sim", "_syn")
fig_dir <- "C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Aug/CR_LL"

for (suffix in suffixes) {
  mean_plots <- lapply(scenarios, create_mean_smd_plot, suffix = suffix)
  boxplots_color <- lapply(scenarios, create_boxplot_color, suffix = suffix)

  jpeg(file.path(fig_dir, paste0("AllScenarios_mean_smd_plot", suffix, ".jpg")), width = 1200, height = 800)
  grid.arrange(grobs = mean_plots, ncol = 2, nrow = 2)
  dev.off()

  jpeg(file.path(fig_dir, paste0("AllScenarios_smdbox_color", suffix, ".jpg")), width = 1200, height = 800)
  grid.arrange(grobs = boxplots_color, ncol = 2, nrow = 2)
  dev.off()
}












# create df from response list
make_resp_df <- function(resp_list, data_type, method) {
  data.frame(
    rate = unlist(resp_list),
    iteration = seq_along(unlist(resp_list)),
    data_type = data_type,
    method = method
  )
}

# combine all response data into one data frame
resp_df <- dplyr::bind_rows(
  make_resp_df(sim_resp_rct, "RCT", "Sim"),
  make_resp_df(syn_resp_rct, "RCT", "Syn"),
  make_resp_df(sim_resp_ob,  "OB",  "Sim"),
  make_resp_df(syn_resp_ob,  "OB",  "Syn"),
  make_resp_df(sim_resp_ex,  "EX",  "Sim"),
  make_resp_df(syn_resp_ex,  "EX",  "Syn")
)

# assign scenarios to resp df
resp_df <- resp_df %>%
  dplyr::mutate(
    scenario = ceiling(iteration / (nSim / 4))
  )

# summarise mean and se
summary_df <- resp_df %>%
  dplyr::group_by(data_type, method, scenario) %>%
  dplyr::summarise(
    mean_rate = mean(rate, na.rm = TRUE),
    se = sd(rate, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    scenario_label = factor(paste("Scenario", scenario),
                            levels = paste("Scenario", 1:4))
  )

# plot
p <- ggplot(summary_df, aes(x = method, y = mean_rate, color = data_type)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = mean_rate - se, ymax = mean_rate + se),
                position = position_dodge(width = 0.5), width = 0.2) +
  facet_wrap(~scenario_label, ncol = 1) +
  labs(title = "Control Response Rate by Method, Data Type, and Scenario",
       x = "Method",
       y = "Mean Control Response Rate",
       color = "Data Type") +
  scale_color_brewer(palette = "Dark2") +

  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    panel.spacing = unit(1.2, "lines")
  )

ggsave("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Aug/CR_LL/control_response_rate_means.jpg",
       p, width = 6, height = 12)





# get exact values for reporting:

# function to split response rate values into sep dfs:
split_to_dfs <- function(values, k = 4, prefix = "df") {
  # k is number of groups
  split_values <- split(values, cut(seq_along(values), k, labels = FALSE))

  # loops
  for (i in seq_along(split_values)) {
    df_name <- paste0(prefix, i)
    assign(df_name, data.frame(value = split_values[[i]]), envir = .GlobalEnv)
  }
}

split_to_dfs(sim_resp_rct, k=4)
sim_resp_rct_1 <- data.frame(value = as.vector(t(df1)))
sim_resp_rct_2 <- data.frame(value = as.vector(t(df2)))
sim_resp_rct_3 <- data.frame(value = as.vector(t(df3)))
sim_resp_rct_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(sim_resp_ob, k=4)
sim_resp_ob_1 <- data.frame(value = as.vector(t(df1)))
sim_resp_ob_2 <- data.frame(value = as.vector(t(df2)))
sim_resp_ob_3 <- data.frame(value = as.vector(t(df3)))
sim_resp_ob_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(sim_resp_ex, k=4)
sim_resp_ex_1 <- data.frame(value = as.vector(t(df1)))
sim_resp_ex_2 <- data.frame(value = as.vector(t(df2)))
sim_resp_ex_3 <- data.frame(value = as.vector(t(df3)))
sim_resp_ex_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(syn_resp_rct, k=4)
syn_resp_rct_1 <- data.frame(value = as.vector(t(df1)))
syn_resp_rct_2 <- data.frame(value = as.vector(t(df2)))
syn_resp_rct_3 <- data.frame(value = as.vector(t(df3)))
syn_resp_rct_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(syn_resp_ob, k=4)
syn_resp_ob_1 <- data.frame(value = as.vector(t(df1)))
syn_resp_ob_2 <- data.frame(value = as.vector(t(df2)))
syn_resp_ob_3 <- data.frame(value = as.vector(t(df3)))
syn_resp_ob_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(syn_resp_ex, k=4)
syn_resp_ex_1 <- data.frame(value = as.vector(t(df1)))
syn_resp_ex_2 <- data.frame(value = as.vector(t(df2)))
syn_resp_ex_3 <- data.frame(value = as.vector(t(df3)))
syn_resp_ex_4 <- data.frame(value = as.vector(t(df4)))



mean(sim_resp_rct_1$value)
mean(sim_resp_rct_2$value)
mean(sim_resp_rct_3$value)
mean(sim_resp_rct_4$value)
sd(sim_resp_rct_1$value)
sd(sim_resp_rct_2$value)
sd(sim_resp_rct_3$value)
sd(sim_resp_rct_4$value)
range(sim_resp_rct_1$value)
range(sim_resp_rct_2$value)
range(sim_resp_rct_3$value)
range(sim_resp_rct_4$value)

mean(sim_resp_ob_1$value)
mean(sim_resp_ob_2$value)
mean(sim_resp_ob_3$value)
mean(sim_resp_ob_4$value)
sd(sim_resp_ob_1$value)
sd(sim_resp_ob_2$value)
sd(sim_resp_ob_3$value)
sd(sim_resp_ob_4$value)
range(sim_resp_ob_1$value)
range(sim_resp_ob_2$value)
range(sim_resp_ob_3$value)
range(sim_resp_ob_4$value)

mean(sim_resp_ex_1$value)
mean(sim_resp_ex_2$value)
mean(sim_resp_ex_3$value)
mean(sim_resp_ex_4$value)
sd(sim_resp_ex_1$value)
sd(sim_resp_ex_2$value)
sd(sim_resp_ex_3$value)
sd(sim_resp_ex_4$value)
range(sim_resp_ex_1$value)
range(sim_resp_ex_2$value)
range(sim_resp_ex_3$value)
range(sim_resp_ex_4$value)




mean(syn_resp_rct_1$value)
mean(syn_resp_rct_2$value)
mean(syn_resp_rct_3$value)
mean(syn_resp_rct_4$value)
sd(syn_resp_rct_1$value)
sd(syn_resp_rct_2$value)
sd(syn_resp_rct_3$value)
sd(syn_resp_rct_4$value)
range(syn_resp_rct_1$value)
range(syn_resp_rct_2$value)
range(syn_resp_rct_3$value)
range(syn_resp_rct_4$value)

mean(syn_resp_ob_1$value)
mean(syn_resp_ob_2$value)
mean(syn_resp_ob_3$value)
mean(syn_resp_ob_4$value)
sd(syn_resp_ob_1$value)
sd(syn_resp_ob_2$value)
sd(syn_resp_ob_3$value)
sd(syn_resp_ob_4$value)
range(syn_resp_ob_1$value)
range(syn_resp_ob_2$value)
range(syn_resp_ob_3$value)
range(syn_resp_ob_4$value)

mean(syn_resp_ex_1$value)
mean(syn_resp_ex_2$value)
mean(syn_resp_ex_3$value)
mean(syn_resp_ex_4$value)
sd(syn_resp_ex_1$value)
sd(syn_resp_ex_2$value)
sd(syn_resp_ex_3$value)
sd(syn_resp_ex_4$value)
range(syn_resp_ex_1$value)
range(syn_resp_ex_2$value)
range(syn_resp_ex_3$value)
range(syn_resp_ex_4$value)


test_resp_rct_1 <- tsum.test(
  mean.x = mean(sim_resp_rct_1$value), s.x = sd(sim_resp_rct_1$value), n.x = length(sim_resp_rct_1$value),
  mean.y = mean(syn_resp_rct_1$value), s.y = sd(syn_resp_rct_1$value), n.y = length(syn_resp_rct_1$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_rct_2 <- tsum.test(
  mean.x = mean(sim_resp_rct_2$value), s.x = sd(sim_resp_rct_2$value), n.x = length(sim_resp_rct_2$value),
  mean.y = mean(syn_resp_rct_2$value), s.y = sd(syn_resp_rct_2$value), n.y = length(syn_resp_rct_2$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_rct_3 <- tsum.test(
  mean.x = mean(sim_resp_rct_3$value), s.x = sd(sim_resp_rct_3$value), n.x = length(sim_resp_rct_3$value),
  mean.y = mean(syn_resp_rct_3$value), s.y = sd(syn_resp_rct_3$value), n.y = length(syn_resp_rct_3$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_rct_4 <- tsum.test(
  mean.x = mean(sim_resp_rct_4$value), s.x = sd(sim_resp_rct_4$value), n.x = length(sim_resp_rct_4$value),
  mean.y = mean(syn_resp_rct_4$value), s.y = sd(syn_resp_rct_4$value), n.y = length(syn_resp_rct_4$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)




test_resp_ob_1 <- tsum.test(
  mean.x = mean(sim_resp_ob_1$value), s.x = sd(sim_resp_ob_1$value), n.x = length(sim_resp_ob_1$value),
  mean.y = mean(syn_resp_ob_1$value), s.y = sd(syn_resp_ob_1$value), n.y = length(syn_resp_ob_1$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ob_2 <- tsum.test(
  mean.x = mean(sim_resp_ob_2$value), s.x = sd(sim_resp_ob_2$value), n.x = length(sim_resp_ob_2$value),
  mean.y = mean(syn_resp_ob_2$value), s.y = sd(syn_resp_ob_2$value), n.y = length(syn_resp_ob_2$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ob_3 <- tsum.test(
  mean.x = mean(sim_resp_ob_3$value), s.x = sd(sim_resp_ob_3$value), n.x = length(sim_resp_ob_3$value),
  mean.y = mean(syn_resp_ob_3$value), s.y = sd(syn_resp_ob_3$value), n.y = length(syn_resp_ob_3$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ob_4 <- tsum.test(
  mean.x = mean(sim_resp_ob_4$value), s.x = sd(sim_resp_ob_4$value), n.x = length(sim_resp_ob_4$value),
  mean.y = mean(syn_resp_ob_4$value), s.y = sd(syn_resp_ob_4$value), n.y = length(syn_resp_ob_4$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)




test_resp_ex_1 <- tsum.test(
  mean.x = mean(sim_resp_ex_1$value), s.x = sd(sim_resp_ex_1$value), n.x = length(sim_resp_ex_1$value),
  mean.y = mean(syn_resp_ex_1$value), s.y = sd(syn_resp_ex_1$value), n.y = length(syn_resp_ex_1$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ex_2 <- tsum.test(
  mean.x = mean(sim_resp_ex_2$value), s.x = sd(sim_resp_ex_2$value), n.x = length(sim_resp_ex_2$value),
  mean.y = mean(syn_resp_ex_2$value), s.y = sd(syn_resp_ex_2$value), n.y = length(syn_resp_ex_2$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ex_3 <- tsum.test(
  mean.x = mean(sim_resp_ex_3$value), s.x = sd(sim_resp_ex_3$value), n.x = length(sim_resp_ex_3$value),
  mean.y = mean(syn_resp_ex_3$value), s.y = sd(syn_resp_ex_3$value), n.y = length(syn_resp_ex_3$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ex_4 <- tsum.test(
  mean.x = mean(sim_resp_ex_4$value), s.x = sd(sim_resp_ex_4$value), n.x = length(sim_resp_ex_4$value),
  mean.y = mean(syn_resp_ex_4$value), s.y = sd(syn_resp_ex_4$value), n.y = length(syn_resp_ex_4$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)




# Extract t-test results with p-values, confidence intervals, and effect sizes
ttest_pvalues <- data.frame(
  Test = c("RCT_1", "RCT_2", "RCT_3", "RCT_4", "OB_1", "OB_2", "OB_3", "OB_4", "EX_1", "EX_2", "EX_3", "EX_4"),
  P_Value = c(test_resp_rct_1$p.value, test_resp_rct_2$p.value, test_resp_rct_3$p.value, test_resp_rct_4$p.value,
              test_resp_ob_1$p.value, test_resp_ob_2$p.value, test_resp_ob_3$p.value, test_resp_ob_4$p.value,
              test_resp_ex_1$p.value, test_resp_ex_2$p.value, test_resp_ex_3$p.value, test_resp_ex_4$p.value),
  CI_Lower = c(test_resp_rct_1$conf.int[1], test_resp_rct_2$conf.int[1], test_resp_rct_3$conf.int[1], test_resp_rct_4$conf.int[1],
               test_resp_ob_1$conf.int[1], test_resp_ob_2$conf.int[1], test_resp_ob_3$conf.int[1], test_resp_ob_4$conf.int[1],
               test_resp_ex_1$conf.int[1], test_resp_ex_2$conf.int[1], test_resp_ex_3$conf.int[1], test_resp_ex_4$conf.int[1]),
  CI_Upper = c(test_resp_rct_1$conf.int[2], test_resp_rct_2$conf.int[2], test_resp_rct_3$conf.int[2], test_resp_rct_4$conf.int[2],
               test_resp_ob_1$conf.int[2], test_resp_ob_2$conf.int[2], test_resp_ob_3$conf.int[2], test_resp_ob_4$conf.int[2],
               test_resp_ex_1$conf.int[2], test_resp_ex_2$conf.int[2], test_resp_ex_3$conf.int[2], test_resp_ex_4$conf.int[2]),
  T_Statistic = c(test_resp_rct_1$statistic, test_resp_rct_2$statistic, test_resp_rct_3$statistic, test_resp_rct_4$statistic,
                  test_resp_ob_1$statistic, test_resp_ob_2$statistic, test_resp_ob_3$statistic, test_resp_ob_4$statistic,
                  test_resp_ex_1$statistic, test_resp_ex_2$statistic, test_resp_ex_3$statistic, test_resp_ex_4$statistic)
)
write.csv(ttest_pvalues, file = "CROHNS_ttest_pvalues_ll.csv", row.names = FALSE)




# for smd:
split_to_dfs(rct_age_smd_all, k=4)
rct_age_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_age_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_age_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_age_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(rct_sex_smd_all, k=4)
rct_sex_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_sex_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_sex_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_sex_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(rct_treat_smd_all, k=4)
rct_treat_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_treat_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_treat_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_treat_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(rct_outcome_smd_all, k=4)
rct_outcome_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_outcome_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_outcome_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_outcome_smd_4 <- data.frame(value = as.vector(t(df4)))


split_to_dfs(ob_age_smd_all, k=4)
ob_age_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_age_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_age_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_age_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ob_sex_smd_all, k=4)
ob_sex_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_sex_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_sex_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_sex_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ob_treat_smd_all, k=4)
ob_treat_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_treat_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_treat_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_treat_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ob_outcome_smd_all, k=4)
ob_outcome_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_outcome_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_outcome_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_outcome_smd_4 <- data.frame(value = as.vector(t(df4)))



split_to_dfs(ex_age_smd_all, k=4)
ex_age_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_age_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_age_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_age_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ex_sex_smd_all, k=4)
ex_sex_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_sex_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_sex_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_sex_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ex_treat_smd_all, k=4)
ex_treat_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_treat_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_treat_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_treat_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ex_outcome_smd_all, k=4)
ex_outcome_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_outcome_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_outcome_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_outcome_smd_4 <- data.frame(value = as.vector(t(df4)))



# just median, range, sd

range(rct_age_smd_1$value)
median(rct_age_smd_1$value)
sd(rct_age_smd_1$value)
range(rct_age_smd_2$value)
median(rct_age_smd_2$value)
sd(rct_age_smd_2$value)
range(rct_age_smd_3$value)
median(rct_age_smd_3$value)
sd(rct_age_smd_3$value)
range(rct_age_smd_4$value)
median(rct_age_smd_4$value)
sd(rct_age_smd_4$value)

range(rct_sex_smd_1$value)
median(rct_sex_smd_1$value)
sd(rct_sex_smd_1$value)
range(rct_sex_smd_2$value)
median(rct_sex_smd_2$value)
sd(rct_sex_smd_2$value)
range(rct_sex_smd_3$value)
median(rct_sex_smd_3$value)
sd(rct_sex_smd_3$value)
range(rct_sex_smd_4$value)
median(rct_sex_smd_4$value)
sd(rct_sex_smd_4$value)


range(rct_treat_smd_1$value)
median(rct_treat_smd_1$value)
sd(rct_treat_smd_1$value)
range(rct_treat_smd_2$value)
median(rct_treat_smd_2$value)
sd(rct_treat_smd_2$value)
range(rct_treat_smd_3$value)
median(rct_treat_smd_3$value)
sd(rct_treat_smd_3$value)
range(rct_treat_smd_4$value)
median(rct_treat_smd_4$value)
sd(rct_treat_smd_4$value)



range(rct_outcome_smd_1$value)
median(rct_outcome_smd_1$value)
sd(rct_outcome_smd_1$value)
range(rct_outcome_smd_2$value)
median(rct_outcome_smd_2$value)
sd(rct_outcome_smd_2$value)

range(rct_outcome_smd_3$value)
median(rct_outcome_smd_3$value)
sd(rct_outcome_smd_3$value)
range(rct_outcome_smd_4$value)
median(rct_outcome_smd_4$value)
sd(rct_outcome_smd_4$value)


range(ob_age_smd_1$value)
median(ob_age_smd_1$value)
sd(ob_age_smd_1$value)
range(ob_age_smd_2$value)
median(ob_age_smd_2$value)
sd(ob_age_smd_2$value)
range(ob_age_smd_3$value)
median(ob_age_smd_3$value)
sd(ob_age_smd_3$value)
range(ob_age_smd_4$value)
median(ob_age_smd_4$value)
sd(ob_age_smd_4$value)

range(ob_sex_smd_1$value)
median(ob_sex_smd_1$value)
sd(ob_sex_smd_1$value)
range(ob_sex_smd_2$value)
median(ob_sex_smd_2$value)
sd(ob_sex_smd_2$value)
range(ob_sex_smd_3$value)
median(ob_sex_smd_3$value)
sd(ob_sex_smd_3$value)
range(ob_sex_smd_4$value)
median(ob_sex_smd_4$value)
sd(ob_sex_smd_4$value)


range(ob_treat_smd_1$value)
median(ob_treat_smd_1$value)
sd(ob_treat_smd_1$value)
range(ob_treat_smd_2$value)
median(ob_treat_smd_2$value)
sd(ob_treat_smd_2$value)
range(ob_treat_smd_3$value)
median(ob_treat_smd_3$value)
sd(ob_treat_smd_3$value)
range(ob_treat_smd_4$value)
median(ob_treat_smd_4$value)
sd(ob_treat_smd_4$value)


range(ob_outcome_smd_1$value)
median(ob_outcome_smd_1$value)
sd(ob_outcome_smd_1$value)
range(ob_outcome_smd_2$value)
median(ob_outcome_smd_2$value)
sd(ob_outcome_smd_2$value)
range(ob_outcome_smd_3$value)
median(ob_outcome_smd_3$value)
sd(ob_outcome_smd_3$value)

range(ob_outcome_smd_4$value)
median(ob_outcome_smd_4$value)
sd(ob_outcome_smd_4$value)



range(ex_age_smd_1$value)
median(ex_age_smd_1$value)
sd(ex_age_smd_1$value)
range(ex_age_smd_2$value)
median(ex_age_smd_2$value)
sd(ex_age_smd_2$value)
range(ex_age_smd_3$value)
median(ex_age_smd_3$value)
sd(ex_age_smd_3$value)
range(ex_age_smd_4$value)
median(ex_age_smd_4$value)
sd(ex_age_smd_4$value)

range(ex_sex_smd_1$value)
median(ex_sex_smd_1$value)
sd(ex_sex_smd_1$value)
range(ex_sex_smd_2$value)
median(ex_sex_smd_2$value)
sd(ex_sex_smd_2$value)
range(ex_sex_smd_3$value)
median(ex_sex_smd_3$value)
sd(ex_sex_smd_3$value)
range(ex_sex_smd_4$value)
median(ex_sex_smd_4$value)
sd(ex_sex_smd_4$value)


range(ex_treat_smd_1$value)
median(ex_treat_smd_1$value)
sd(ex_treat_smd_1$value)
range(ex_treat_smd_2$value)
median(ex_treat_smd_2$value)
sd(ex_treat_smd_2$value)
range(ex_treat_smd_3$value)
median(ex_treat_smd_3$value)
sd(ex_treat_smd_3$value)
range(ex_treat_smd_4$value)
median(ex_treat_smd_4$value)
sd(ex_treat_smd_4$value)


range(ex_outcome_smd_1$value)
median(ex_outcome_smd_1$value)
sd(ex_outcome_smd_1$value)
range(ex_outcome_smd_2$value)
median(ex_outcome_smd_2$value)
sd(ex_outcome_smd_2$value)
range(ex_outcome_smd_3$value)
median(ex_outcome_smd_3$value)
sd(ex_outcome_smd_3$value)
range(ex_outcome_smd_4$value)
median(ex_outcome_smd_4$value)
sd(ex_outcome_smd_4$value)

































































# repeat for rs:


# creating scenarios for number of observations
# halve the n due to pre and post being added together later
nsize <- function (scenario, datatype, simsyn) {
  if (scenario == 1) { # scenario one is n = 150 to n = 150
    N <- 75
  }   else if (scenario == 2 ) { # scenario two is n = 500 to n = 500
    N <- 250
  }   else if (scenario == 3 & datatype == "rct") { # scenario three is realistic sample sizes
    N <- 628
  }   else if (scenario == 3 & datatype == "ob"){
    N <- 221
  }   else if (scenario == 3 & datatype == "ex"){
    N <- 322
  }   else if (scenario == 4 & simsyn == "sim"){ # scenario four is n = 150 to n = 500
    N <- 75
  }   else if (scenario == 4 & simsyn == "syn"){
    N <- 250
  }
  return(N)
}



# function for smd calculation:
process_data <- function(df1, df2, var_names, table_list_name, smd_list_prefix, i) {
  # merge df
  df1$Key <- 1:nrow(df1)
  df2$Key <- 1:nrow(df2)
  merged_df <- merge(df1, df2, by = "Key", all = TRUE)
  merged_df <- merged_df[, -which(names(merged_df) == "Key")]

  # make a single df
  combined_df <- data.frame(
    lapply(var_names, function(var) {
      c(merged_df[[paste0(var, ".x")]], merged_df[[paste0(var, ".y")]])
    })
  )
  colnames(combined_df) <- var_names

  # save table
  temp <- get(table_list_name, envir = .GlobalEnv)
  temp[[i]] <- CreateTableOne(vars = var_names, data = combined_df, test = FALSE)
  assign(table_list_name, temp, envir = .GlobalEnv)

  # calculate smd
  for (var in var_names) {
    x <- merged_df[[paste0(var, ".x")]]
    y <- merged_df[[paste0(var, ".y")]]

    pooled_sd <- sqrt((sd(x, na.rm = TRUE)^2 + sd(y, na.rm = TRUE)^2) / 2)
    mean_diff <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
    smd_value <- mean_diff / pooled_sd

    smd_list_name <- paste0(smd_list_prefix, "_", var, "_smd_all")
    smd_list <- get(smd_list_name, envir = .GlobalEnv)
    smd_list[[i]] <- smd_value
    assign(smd_list_name, smd_list, envir = .GlobalEnv)
  }
}


# function for smd calculation across data types
compute_smd_across <- function(df1, df2, vars, prefix, suffix, i) {
  for (var in vars) {
    x1 <- df1[[paste0(var, ".", suffix)]]
    x2 <- df2[[paste0(var, ".", suffix)]]

    sd1 <- sd(x1, na.rm = TRUE)
    sd2 <- sd(x2, na.rm = TRUE)
    pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
    mean_diff <- mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)

    # by summix
    smd_list_name <- if (suffix == "x") {
      paste0(prefix, "_", var, "_smd_sim")
    } else if (suffix == "y") {
      paste0(prefix, "_", var, "_smd_syn")
    } else {
      stop("compute_smd_across: suffix must be 'x' or 'y'")
    }

    smd_list <- get(smd_list_name, envir = .GlobalEnv)
    smd_list[[i]] <- mean_diff / pooled_sd
    assign(smd_list_name, smd_list, envir = .GlobalEnv)
  }
}
# function for making empty variable lists
create_smd_lists <- function(non_combo_prefixes, combo_prefixes, variables, n) {

  for (pref in non_combo_prefixes) {
    for (var in variables) {
      list_name <- paste0(pref, "_", var, "_smd_all")
      assign(list_name, vector("list", length = n), envir = .GlobalEnv)
    }
  }
  # combos for sim/syn
  for (pref in combo_prefixes) {
    for (var in variables) {
      assign(paste0(pref, "_", var, "_smd_sim"), vector("list", length = n), envir = .GlobalEnv)
      assign(paste0(pref, "_", var, "_smd_syn"), vector("list", length = n), envir = .GlobalEnv)
    }
  }
}
# function for making empty table lists
create_table_lists <- function(names, n) {
  for (name in names) {
    assign(name, vector("list", length = n), envir = .GlobalEnv)
  }
}


# function for control response rate

control_response <- function(df, list_name, i) {
  num_responders <- df %>%
    filter(treat == 0, outcome == 1) %>%
    nrow()

  num_control <- df %>%
    filter(treat == 0) %>%
    nrow()

  rate <- ifelse(num_control > 0, num_responders / num_control, NA)

  tmp <- get(list_name, envir = .GlobalEnv)
  tmp[[i]] <- rate
  assign(list_name, tmp, envir = .GlobalEnv)
}


# for loop for sim and analysis

nSim <- 40 # number of simulations # small for testing
set.seed(9999)


#  make empty lists

vars <- c("age", "sex", "treat", "outcome")
non_combo_prefixes <- c("rct", "ob", "ex")                   # these get *_smd_all
combo_prefixes     <- c("rct_ob", "ob_ex", "rct_ex")         # these get *_smd_sim and *_smd_syn


create_smd_lists(non_combo_prefixes, combo_prefixes, vars, nSim)

create_table_lists(c("table_rct_data", "table_ob_data", "table_ex_data"), nSim)

sim_results_rct <- vector("list", length = nSim)
sim_results_ob  <- vector("list", length = nSim)
sim_results_ex  <- vector("list", length = nSim)
syn_results_rct <- vector("list", length = nSim)
syn_results_ob <- vector("list", length = nSim)
syn_results_ex <- vector("list", length = nSim)

sim_resp_rct <- vector("list", length = nSim)
sim_resp_ob  <- vector("list", length = nSim)
sim_resp_ex  <- vector("list", length = nSim)
syn_resp_rct <- vector("list", length = nSim)
syn_resp_ob <- vector("list", length = nSim)
syn_resp_ex <- vector("list", length = nSim)






# list for sum stats

sim_sum_rct_mean <- vector("list", length = nSim)
sim_sum_ob_mean  <- vector("list", length = nSim)
sim_sum_ex_mean  <- vector("list", length = nSim)
syn_sum_rct_mean <- vector("list", length = nSim)
syn_sum_ob_mean <- vector("list", length = nSim)
syn_sum_ex_mean <- vector("list", length = nSim)


sim_sum_rct_med <- vector("list", length = nSim)
sim_sum_ob_med  <- vector("list", length = nSim)
sim_sum_ex_med  <- vector("list", length = nSim)
syn_sum_rct_med <- vector("list", length = nSim)
syn_sum_ob_med <- vector("list", length = nSim)
syn_sum_ex_med <- vector("list", length = nSim)

sim_sum_rct_range <- vector("list", length = nSim)
sim_sum_ob_range  <- vector("list", length = nSim)
sim_sum_ex_range  <- vector("list", length = nSim)
syn_sum_rct_range <- vector("list", length = nSim)
syn_sum_ob_range <- vector("list", length = nSim)
syn_sum_ex_range <- vector("list", length = nSim)

sim_sum_rct_sd <- vector("list", length = nSim)
sim_sum_ob_sd  <- vector("list", length = nSim)
sim_sum_ex_sd  <- vector("list", length = nSim)
syn_sum_rct_sd <- vector("list", length = nSim)
syn_sum_ob_sd <- vector("list", length = nSim)
syn_sum_ex_sd <- vector("list", length = nSim)

# list for dist
sim_kurt_rct <- vector("list", length = nSim)
sim_skew_rct <- vector("list", length = nSim)

sim_kurt_ob <- vector("list", length = nSim)
sim_skew_ob <- vector("list", length = nSim)

sim_kurt_ex <- vector("list", length = nSim)
sim_skew_ex <- vector("list", length = nSim)

syn_kurt_rct <- vector("list", length = nSim)
syn_skew_rct <- vector("list", length = nSim)

syn_kurt_ob <- vector("list", length = nSim)
syn_skew_ob <- vector("list", length = nSim)

syn_kurt_ex <- vector("list", length = nSim)
syn_skew_ex <- vector("list", length = nSim)


scenario <- 1  #  scenario before the loop


# loop for gen simulations


for (i in 1:nSim) {

  if ((i - 1) %% (nSim/4) == 0 & i > 1) {  # change scenario every 10,000 iterations
    scenario <- scenario + 1
  }



  N <- nsize(scenario, "rct", "sim") # selecting scenario


  treat <- sample(c(1, 0), N, replace = TRUE, prob = c(0.63, 0.37))

  # compliance to treatment
  comply <- ifelse(treat == 1, sample(c(1, 0), N, replace = TRUE, prob = c(1, 0)), sample(c(1, 0), N, replace = TRUE, prob = c(0, 1)))

  # other variables
  sex <- rbinom(N, 1, 0.522)
  age <- rtruncnorm(N, a = 18, b = 72, mean = 37.3, sd = 11.8)
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
    mutate(treat_effect = ifelse(post == 0, NA, treat_effect))  # filter treatment effect to apply only to post-treatment group

  # calculate sex_effect
  df_combined <- df_combined %>%
    mutate(sex_effect = ifelse(sex == 1 | sex == 0, 0, NA)) # sex effect not recorded or reported / only baseline / assumed as zero

  # create overall treatment effect including demographic influences
  df_combined <- df_combined %>%
    mutate(overall_treat_effect = ifelse(post == 0, 0, treat_effect + sex_effect))

  # create pre and post outcomes (baseline vs treatment)
  df_combined <- df_combined %>%
    mutate(outcome = ifelse(post == 0, rbinom(nrow(df_combined[df_combined$post == 0, ]), 1, 0.359),
                            sapply(overall_treat_effect, function(prob) rbinom(1, 1, prob))))
  # store dfs in results
  sim_results_rct[[i]] <- df_combined

  # remove df2
  rm(df2)





  # observational
  # treatment

  N <- nsize(scenario, "ob", "sim") # selecting scenario

  treat=sample(c(1,0), N,replace=TRUE, prob=c(.66,.34)) # probability of being exposed to the treatment or not


  # other variables
  sex <- rbinom(N, 1, 0.602)
  age <- rtruncnorm(N, a = 16, b = 65, mean = 38.2, sd = 16.96)
  age <- round(age, 0) # round age to whole numbers
  random_var <- rbinom(N, 1, prob = 0.5)
  post <- rep(0, N)  # pre and post assignment

  # create df combining variables
  ob_df <- data.frame(treat, sex, age, random_var, post)

  # add id
  ob_df <- ob_df %>% dplyr::mutate(id = row_number())

  # create post df
  ob_df2 <- ob_df %>% mutate(post = 1)

  # combine pre and post
  ob_df_combined <- rbind(ob_df, ob_df2)

  # define treatment effect, only for post group
  # define treat
  ob_df_combined <- ob_df_combined %>% mutate(treat_effect=
                                                ifelse(treat==1, .387,
                                                       ifelse(treat==0, .242,
                                                              NA)))%>%
    mutate(treat_effect = ifelse(post == 0, NA, treat_effect))



  # calculate sex_effect and race_effect
  ob_df_combined <- ob_df_combined %>%
    mutate(sex_effect = ifelse(sex == 1 | sex == 0, 0, NA))

  # create overall treatment effect including demographic influences
  ob_df_combined <- ob_df_combined %>%
    mutate(overall_treat_effect = ifelse(post == 0, 0, treat_effect + sex_effect))


  # create pre and post outcomes
  ob_df_combined <- ob_df_combined %>%
    mutate(outcome = ifelse(post == 0, rbinom(sum(post == 0), 1, 0.242),
                            rbinom(sum(post == 1), 1, overall_treat_effect)))
  # store dfs in results
  sim_results_ob[[i]] <- ob_df_combined

  # remove df2
  rm(ob_df2)




  # external
  # treatment
  N <- nsize(scenario, "ex", "sim") # selecting scenario

  outcome=sample(c(1,0), N*2,replace=TRUE, prob=c(.64,.36))

  # other variables
  treat=sample(c(1,0), N*2,replace=TRUE, prob=c(.88,.12))
  sex <- rbinom(N, 1, 0.48)
  age <- rtruncnorm(N*2, a = 15, b = 70, mean = 38, sd = 17.78)
  age <- round(age, 0) # round age to whole numbers
  random_var <- rbinom(N*2, 1, prob = 0.5)



  # create df combining variables
  ex_df <- data.frame(outcome, treat, sex, age, random_var)

  # add id
  ex_df <- ex_df %>% dplyr::mutate(id = row_number())

  # store dfs in results
  sim_results_ex[[i]] <- ex_df











  # begin synthetic data creation:
  # MODEL - rct
  rct_mydata = subset(df_combined, select = c("id", "age", "sex", "treat", "outcome"))

  N <- nsize(scenario, "rct", "syn") # selecting scenario

  # minimum number of observations needed is 10
  rct_mysyn <- syn(rct_mydata,method = "sample", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)

  summary(rct_mysyn)

  # store dfs in results
  syn_results_rct[[i]] <- rct_mysyn$syn




  # MODEL - obs
  ob_mydata = subset(ob_df_combined, select = c("id", "age", "sex", "treat", "outcome"))



  N <- nsize(scenario, "ob", "syn") # selecting scenario

  # minimum number of observations needed is 10
  ob_mysyn <- syn(ob_mydata,method = "sample", k = N, cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)

  summary(ob_mysyn)

  # store dfs in results
  syn_results_ob[[i]] <- ob_mysyn$syn





  # MODEL - ex
  ex_mydata = subset(ex_df, select = c("id", "age", "sex", "treat", "outcome"))


  N <- nsize(scenario, "ex", "syn") # selecting scenario


  print(i) # sim number

  # minimum number of observations needed is 10
  ex_mysyn <- syn(ex_mydata,method = "sample", k = N,  cont.na = NULL, minnumlevels = 10, maxfaclevels = 80)

  summary(ex_mysyn)

  # store dfs in results
  syn_results_ex[[i]] <- ex_mysyn$syn





  # recreate synthetic data in df format
  rct_mysyn <- rct_mysyn$syn
  ob_mysyn <- ob_mysyn$syn
  ex_mysyn <- ex_mysyn$syn

  # dist stats
  sim_kurt_rct[[i]] <- data.frame(age = kurtosis(df_combined$age),
                                  sex = kurtosis(df_combined$sex),
                                  treat = kurtosis(df_combined$treat),
                                  outcome = kurtosis(df_combined$outcome))


  sim_kurt_ob[[i]] <- data.frame(age = kurtosis(ob_df_combined$age),
                                 sex = kurtosis(ob_df_combined$sex),
                                 treat = kurtosis(ob_df_combined$treat),
                                 outcome = kurtosis(ob_df_combined$outcome))



  sim_kurt_ex[[i]] <- data.frame(age = kurtosis(ex_df$age),
                                 sex = kurtosis(ex_df$sex),
                                 treat = kurtosis(ex_df$treat),
                                 outcome = kurtosis(ex_df$outcome))






  sim_skew_rct[[i]] <- data.frame(age = skewness(df_combined$age),
                                  sex = skewness(df_combined$sex),
                                  treat = skewness(df_combined$treat),
                                  outcome = skewness(df_combined$outcome))


  sim_skew_ob[[i]] <- data.frame(age = skewness(ob_df_combined$age),
                                 sex = skewness(ob_df_combined$sex),
                                 treat = skewness(ob_df_combined$treat),
                                 outcome = skewness(ob_df_combined$outcome))



  sim_skew_ex[[i]] <- data.frame(age = skewness(ex_df$age),
                                 sex = skewness(ex_df$sex),
                                 treat = skewness(ex_df$treat),
                                 outcome = skewness(ex_df$outcome))





  syn_kurt_rct[[i]] <- data.frame(age = kurtosis(rct_mysyn$age),
                                  sex = kurtosis(rct_mysyn$sex),
                                  treat = kurtosis(rct_mysyn$treat),
                                  outcome = kurtosis(rct_mysyn$outcome))


  syn_kurt_ob[[i]] <- data.frame(age = kurtosis(ob_mysyn$age),
                                 sex = kurtosis(ob_mysyn$sex),
                                 treat = kurtosis(ob_mysyn$treat),
                                 outcome = kurtosis(ob_mysyn$outcome))



  syn_kurt_ex[[i]] <- data.frame(age = kurtosis(ex_mysyn$age),
                                 sex = kurtosis(ex_mysyn$sex),
                                 treat = kurtosis(ex_mysyn$treat),
                                 outcome = kurtosis(ex_mysyn$outcome))






  syn_skew_rct[[i]] <- data.frame(age = skewness(rct_mysyn$age),
                                  sex = skewness(rct_mysyn$sex),
                                  treat = skewness(rct_mysyn$treat),
                                  outcome = skewness(rct_mysyn$outcome))


  syn_skew_ob[[i]] <- data.frame(age = skewness(ob_mysyn$age),
                                 sex = skewness(ob_mysyn$sex),
                                 treat = skewness(ob_mysyn$treat),
                                 outcome = skewness(ob_mysyn$outcome))



  syn_skew_ex[[i]] <- data.frame(age = skewness(ex_mysyn$age),
                                 sex = skewness(ex_mysyn$sex),
                                 treat = skewness(ex_mysyn$treat),
                                 outcome = skewness(ex_mysyn$outcome))



  # summary stats

  # mean, median, range, sd,

  sim_sum_rct_mean[[i]] <- data.frame(age = mean(df_combined$age),
                                      sex = mean(df_combined$sex),
                                      treat = mean(df_combined$treat),
                                      outcome = mean(df_combined$outcome))


  sim_sum_rct_med[[i]] <- data.frame(age = median(df_combined$age),
                                     sex = median(df_combined$sex),
                                     treat = median(df_combined$treat),
                                     outcome = median(df_combined$outcome))

  sim_sum_rct_range[[i]] <- data.frame(age = range(df_combined$age),
                                       sex = range(df_combined$sex),
                                       treat = range(df_combined$treat),
                                       outcome = range(df_combined$outcome))


  sim_sum_rct_sd[[i]] <- data.frame(age = sd(df_combined$age),
                                    sex = sd(df_combined$sex),
                                    treat = sd(df_combined$treat),
                                    outcome = sd(df_combined$outcome))






  sim_sum_ob_mean[[i]]  <- data.frame(age = mean(ob_df_combined$age),
                                      sex = mean(ob_df_combined$sex),
                                      treat = mean(ob_df_combined$treat),
                                      outcome = mean(ob_df_combined$outcome))


  sim_sum_ob_med[[i]]  <- data.frame(age = median(ob_df_combined$age),
                                     sex = median(ob_df_combined$sex),
                                     treat = median(ob_df_combined$treat),
                                     outcome = median(ob_df_combined$outcome))


  sim_sum_ob_range[[i]]  <- data.frame(age = range(ob_df_combined$age),
                                       sex = range(ob_df_combined$sex),
                                       treat = range(ob_df_combined$treat),
                                       outcome = range(ob_df_combined$outcome))



  sim_sum_ob_sd[[i]]  <- data.frame(age = sd(ob_df_combined$age),
                                    sex = sd(ob_df_combined$sex),
                                    treat = sd(ob_df_combined$treat),
                                    outcome = sd(ob_df_combined$outcome))





  sim_sum_ex_mean[[i]]  <- data.frame(age = mean(ex_df$age),
                                      sex = mean(ex_df$sex),
                                      treat = mean(ex_df$treat),
                                      outcome = mean(ex_df$outcome))


  sim_sum_ex_med[[i]]  <- data.frame(age = median(ex_df$age),
                                     sex = median(ex_df$sex),
                                     treat = median(ex_df$treat),
                                     outcome = median(ex_df$outcome))

  sim_sum_ex_range[[i]]  <- data.frame(age = range(ex_df$age),
                                       sex = range(ex_df$sex),
                                       treat = range(ex_df$treat),
                                       outcome = range(ex_df$outcome))

  sim_sum_ex_sd[[i]]  <- data.frame(age = sd(ex_df$age),
                                    sex = sd(ex_df$sex),
                                    treat = sd(ex_df$treat),
                                    outcome = sd(ex_df$outcome))











  syn_sum_rct_mean[[i]] <- data.frame(age = mean(rct_mysyn$age),
                                      sex = mean(rct_mysyn$sex),
                                      treat = mean(rct_mysyn$treat),
                                      outcome = mean(rct_mysyn$outcome))


  syn_sum_rct_med[[i]] <- data.frame(age = median(rct_mysyn$age),
                                     sex = median(rct_mysyn$sex),
                                     treat = median(rct_mysyn$treat),
                                     outcome = median(rct_mysyn$outcome))

  syn_sum_rct_range[[i]] <- data.frame(age = range(rct_mysyn$age),
                                       sex = range(rct_mysyn$sex),
                                       treat = range(rct_mysyn$treat),
                                       outcome = range(rct_mysyn$outcome))



  syn_sum_rct_sd[[i]] <- data.frame(age = sd(rct_mysyn$age),
                                    sex = sd(rct_mysyn$sex),
                                    treat = sd(rct_mysyn$treat),
                                    outcome = sd(rct_mysyn$outcome))






  syn_sum_ob_mean[[i]] <- data.frame(age = mean(ob_mysyn$age),
                                     sex = mean(ob_mysyn$sex),
                                     treat = mean(ob_mysyn$treat),
                                     outcome = mean(ob_mysyn$outcome))


  syn_sum_ob_med[[i]] <- data.frame(age = median(ob_mysyn$age),
                                    sex = median(ob_mysyn$sex),
                                    treat = median(ob_mysyn$treat),
                                    outcome = median(ob_mysyn$outcome))

  syn_sum_ob_range[[i]] <- data.frame(age = range(ob_mysyn$age),
                                      sex = range(ob_mysyn$sex),
                                      treat = range(ob_mysyn$treat),
                                      outcome = range(ob_mysyn$outcome))



  syn_sum_ob_sd[[i]] <- data.frame(age = sd(ob_mysyn$age),
                                   sex = sd(ob_mysyn$sex),
                                   treat = sd(ob_mysyn$treat),
                                   outcome = sd(ob_mysyn$outcome))








  syn_sum_ex_mean[[i]] <- data.frame(age = mean(ex_mysyn$age),
                                     sex = mean(ex_mysyn$sex),
                                     treat = mean(ex_mysyn$treat),
                                     outcome = mean(ex_mysyn$outcome))


  syn_sum_ex_med[[i]] <- data.frame(age = median(ex_mysyn$age),
                                    sex = median(ex_mysyn$sex),
                                    treat = median(ex_mysyn$treat),
                                    outcome = median(ex_mysyn$outcome))

  syn_sum_ex_range[[i]] <- data.frame(age = range(ex_mysyn$age),
                                      sex = range(ex_mysyn$sex),
                                      treat = range(ex_mysyn$treat),
                                      outcome = range(ex_mysyn$outcome))



  syn_sum_ex_sd[[i]] <- data.frame(age = sd(ex_mysyn$age),
                                   sex = sd(ex_mysyn$sex),
                                   treat = sd(ex_mysyn$treat),
                                   outcome = sd(ex_mysyn$outcome))



  # calculating smd
  process_data(df_combined, rct_mysyn, vars, "table_rct_data", "rct", i)
  process_data(ob_df_combined, ob_mysyn, vars, "table_ob_data", "ob", i)
  process_data(ex_df, ex_mysyn, vars, "table_ex_data", "ex", i)

  # Synthetic vs Synthetic

  rct_ex <- merge(rct_mysyn, ex_mysyn, by = "id", suffixes = c(".x", ".y"))
  rct_ob <- merge(rct_mysyn, ob_mysyn, by = "id", suffixes = c(".x", ".y"))
  ob_ex  <- merge(ob_mysyn, ex_mysyn, by = "id", suffixes = c(".x", ".y"))

  compute_smd_across(rct_ex, rct_ex, vars, "rct_ex", "y", i)
  compute_smd_across(rct_ob, rct_ob, vars, "rct_ob", "y", i)
  compute_smd_across(ob_ex,  ob_ex,  vars, "ob_ex",  "y", i)

  compute_smd_across(rct_ex, rct_ex, vars, "rct_ex", "x", i)
  compute_smd_across(rct_ob, rct_ob, vars, "rct_ob", "x", i)
  compute_smd_across(ob_ex,  ob_ex,  vars, "ob_ex",  "x", i)





  # calculate control response rate:

  control_response(rct_mysyn, "syn_resp_rct", i)
  control_response(ob_mysyn,  "syn_resp_ob",  i)
  control_response(ex_mysyn,  "syn_resp_ex",  i)

  control_response(df_combined,     "sim_resp_rct", i)
  control_response(ob_df_combined,  "sim_resp_ob",  i)
  control_response(ex_df,           "sim_resp_ex",  i)
}









# ANALYSIS AND GRAPHING



# function to split lists
distribute_into_lists <- function(input_list, sizes) {
  if (sum(sizes) != length(input_list)) {
    stop("The sum of sizes must equal the length of the input list")
  }

  result <- vector("list", length(sizes))
  start_index <- 1

  for (i in seq_along(sizes)) {
    end_index <- start_index + sizes[i] - 1
    result[[i]] <- input_list[start_index:end_index]
    start_index <- end_index + 1
  }

  return(result)
}


list_size <- c(nSim/4, nSim/4, nSim/4, nSim/4)
save_dir <- "C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Saved_Rdata/August/CR_RS/"


# save data:

# non-combos
non_combos <- c("rct", "ob", "ex")
variables <- c("age", "sex", "treat", "outcome")
non_combo_names <- unlist(lapply(non_combos, function(pref) {
  paste0(pref, "_", variables, "_smd_all")
}))

# combos
combos <- c("rct_ex", "rct_ob", "ob_ex")
combo_names <- unlist(lapply(combos, function(pref) {
  paste0(rep(pref, each = length(variables)), "_", variables, "_smd_", rep(c("sim", "syn"), each = length(variables) * 1))
}))

# misc
extra_objects <- c(
  "sim_results_rct", "sim_results_ob", "sim_results_ex",
  "syn_results_rct", "syn_results_ob", "syn_results_ex",
  "sim_resp_rct", "sim_resp_ob", "sim_resp_ex",
  "syn_resp_rct", "syn_resp_ob", "syn_resp_ex",
  "table_rct_data", "table_ob_data", "table_ex_data"
)

# final list of all object names
obj_names <- c(extra_objects, non_combo_names, combo_names)

# save all
for (obj_name in obj_names) {

  if (!exists(obj_name)) {
    warning(paste("Object not found:", obj_name))
    next
  }

  obj <- get(obj_name)

  # check for list
  if (!is.list(obj)) {
    warning(paste("Object is not a list:", obj_name))
    next
  }

  obj_split <- distribute_into_lists(obj, list_size)

  for (scenario_index in seq_along(obj_split)) {
    saveRDS(
      obj_split[[scenario_index]],
      file = paste0(save_dir, obj_name, "_scenario", scenario_index, ".Rdata")
    )
  }
}


# dist stats

sim_kurt_rct = do.call(rbind, sim_kurt_rct)
sim_kurt_ob = do.call(rbind,sim_kurt_ob)
sim_kurt_ex = do.call(rbind, sim_kurt_ex)
sim_skew_rct = do.call(rbind, sim_skew_rct)
sim_skew_ob = do.call(rbind,sim_skew_ob)
sim_skew_ex = do.call(rbind, sim_skew_ex)

syn_kurt_rct = do.call(rbind, syn_kurt_rct)
syn_kurt_ob = do.call(rbind,syn_kurt_ob)
syn_kurt_ex = do.call(rbind, syn_kurt_ex)
syn_skew_rct = do.call(rbind, syn_skew_rct)
syn_skew_ob = do.call(rbind,syn_skew_ob)
syn_skew_ex = do.call(rbind, syn_skew_ex)

write.csv(sim_kurt_rct, file = "cr_sim_kurt_rct_rs.csv", row.names = FALSE)
write.csv(sim_kurt_ob, file = "cr_sim_kurt_ob_rs.csv", row.names = FALSE)
write.csv(sim_kurt_ex, file = "cr_sim_kurt_ex_rs.csv", row.names = FALSE)
write.csv(sim_skew_rct, file = "cr_sim_skew_rct_rs.csv", row.names = FALSE)
write.csv(sim_skew_ob, file = "cr_sim_skew_ob_rs.csv", row.names = FALSE)
write.csv(sim_skew_ex, file = "cr_sim_skew_ex_rs.csv", row.names = FALSE)

write.csv(syn_kurt_rct, file = "cr_syn_kurt_rct_rs.csv", row.names = FALSE)
write.csv(syn_kurt_ob, file = "cr_syn_kurt_ob_rs.csv", row.names = FALSE)
write.csv(syn_kurt_ex, file = "cr_syn_kurt_ex_rs.csv", row.names = FALSE)
write.csv(syn_skew_rct, file = "cr_syn_skew_rct_rs.csv", row.names = FALSE)
write.csv(syn_skew_ob, file = "cr_syn_skew_ob_rs.csv", row.names = FALSE)
write.csv(syn_skew_ex, file = "cr_syn_skew_ex_rs.csv", row.names = FALSE)

sim_sum_rct_mean = do.call(rbind, sim_sum_rct_mean)
sim_sum_rct_med = do.call(rbind, sim_sum_rct_med)
sim_sum_rct_range = do.call(rbind, sim_sum_rct_range)
sim_sum_rct_sd = do.call(rbind, sim_sum_rct_sd)

sim_sum_ob_mean = do.call(rbind, sim_sum_ob_mean)
sim_sum_ob_med = do.call(rbind, sim_sum_ob_med)
sim_sum_ob_range = do.call(rbind, sim_sum_ob_range)
sim_sum_ob_sd = do.call(rbind, sim_sum_ob_sd)


sim_sum_ex_mean = do.call(rbind, sim_sum_ex_mean)
sim_sum_ex_med = do.call(rbind, sim_sum_ex_med)
sim_sum_ex_range = do.call(rbind, sim_sum_ex_range)
sim_sum_ex_sd = do.call(rbind, sim_sum_ex_sd)



syn_sum_rct_mean = do.call(rbind, syn_sum_rct_mean)
syn_sum_rct_med = do.call(rbind, syn_sum_rct_med)
syn_sum_rct_range = do.call(rbind, syn_sum_rct_range)
syn_sum_rct_sd = do.call(rbind, syn_sum_rct_sd)

syn_sum_ob_mean = do.call(rbind, syn_sum_ob_mean)
syn_sum_ob_med = do.call(rbind, syn_sum_ob_med)
syn_sum_ob_range = do.call(rbind, syn_sum_ob_range)
syn_sum_ob_sd = do.call(rbind, syn_sum_ob_sd)


syn_sum_ex_mean = do.call(rbind, syn_sum_ex_mean)
syn_sum_ex_med = do.call(rbind, syn_sum_ex_med)
syn_sum_ex_range = do.call(rbind, syn_sum_ex_range)
syn_sum_ex_sd = do.call(rbind, syn_sum_ex_sd)




write.csv(sim_sum_rct_mean, file = "cr_sim_sum_rct_rs_mean.csv", row.names = FALSE)
write.csv(sim_sum_rct_med, file = "cr_sim_sum_rct_rs_med.csv", row.names = FALSE)
write.csv(sim_sum_rct_range, file = "cr_sim_sum_rct_rs_range.csv", row.names = FALSE)
write.csv(sim_sum_rct_sd, file = "cr_sim_sum_rct_rs_sd.csv", row.names = FALSE)

write.csv(sim_sum_ob_mean, file = "cr_sim_sum_ob_rs_mean.csv", row.names = FALSE)
write.csv(sim_sum_ob_med, file = "cr_sim_sum_ob_rs_med.csv", row.names = FALSE)
write.csv(sim_sum_ob_range, file = "cr_sim_sum_ob_rs_range.csv", row.names = FALSE)
write.csv(sim_sum_ob_sd, file = "cr_sim_sum_ob_rs_sd.csv", row.names = FALSE)

write.csv(sim_sum_ex_mean, file = "cr_sim_sum_ex_rs_mean.csv", row.names = FALSE)
write.csv(sim_sum_ex_med, file = "cr_sim_sum_ex_rs_med.csv", row.names = FALSE)
write.csv(sim_sum_ex_range, file = "cr_sim_sum_ex_rs_range.csv", row.names = FALSE)
write.csv(sim_sum_ex_sd, file = "cr_sim_sum_ex_rs_sd.csv", row.names = FALSE)


write.csv(syn_sum_rct_mean, file = "cr_syn_sum_rct_rs_mean.csv", row.names = FALSE)
write.csv(syn_sum_rct_med, file = "cr_syn_sum_rct_rs_med.csv", row.names = FALSE)
write.csv(syn_sum_rct_range, file = "cr_syn_sum_rct_rs_range.csv", row.names = FALSE)
write.csv(syn_sum_rct_sd, file = "cr_syn_sum_rct_rs_sd.csv", row.names = FALSE)

write.csv(syn_sum_ob_mean, file = "cr_syn_sum_ob_rs_mean.csv", row.names = FALSE)
write.csv(syn_sum_ob_med, file = "cr_syn_sum_ob_rs_med.csv", row.names = FALSE)
write.csv(syn_sum_ob_range, file = "cr_syn_sum_ob_rs_range.csv", row.names = FALSE)
write.csv(syn_sum_ob_sd, file = "cr_syn_sum_ob_rs_sd.csv", row.names = FALSE)

write.csv(syn_sum_ex_mean, file = "cr_syn_sum_ex_rs_mean.csv", row.names = FALSE)
write.csv(syn_sum_ex_med, file = "cr_syn_sum_ex_rs_med.csv", row.names = FALSE)
write.csv(syn_sum_ex_range, file = "cr_syn_sum_ex_rs_range.csv", row.names = FALSE)
write.csv(syn_sum_ex_sd, file = "cr_syn_sum_ex_rs_sd.csv", row.names = FALSE)




# plotting smd figures




# plotting smd figures

variables <- c("Age", "Sex", "Treatment", "Outcome")
types <- c("RCT", "Observational", "External")

# build variable names
make_varname <- function(prefix, var, suffix = "") {
  paste0(prefix, "_", tolower(var), "_smd", suffix)
}

# get relevant lists for each suffix and scenario
get_smd_lists <- function(scenario_num, suffix) {
  get_slice <- function(varname) {
    full_list <- get(varname)
    total_len <- length(full_list)
    num_scenarios <- 4
    num_sims_per_scenario <- total_len / num_scenarios

    if (num_sims_per_scenario != floor(num_sims_per_scenario)) {
      stop("Number of simulations per scenario is not an integer. Check data length.")
    }

    start_idx <- (scenario_num - 1) * num_sims_per_scenario + 1
    end_idx <- scenario_num * num_sims_per_scenario
    full_list[start_idx:end_idx]
  }

  if (suffix == "_all") {
    rct_list <- list(
      get_slice(make_varname("rct", "age", suffix)),
      get_slice(make_varname("rct", "sex", suffix)),
      get_slice(make_varname("rct", "treat", suffix)),
      get_slice(make_varname("rct", "outcome", suffix))
    )
    ob_list <- list(
      get_slice(make_varname("ob", "age", suffix)),
      get_slice(make_varname("ob", "sex", suffix)),
      get_slice(make_varname("ob", "treat", suffix)),
      get_slice(make_varname("ob", "outcome", suffix))
    )
    ex_list <- list(
      get_slice(make_varname("ex", "age", suffix)),
      get_slice(make_varname("ex", "sex", suffix)),
      get_slice(make_varname("ex", "treat", suffix)),
      get_slice(make_varname("ex", "outcome", suffix))
    )
  } else if (suffix %in% c("_sim", "_syn")) {
    rct_list <- list(
      get_slice(make_varname("rct_ob", "age", suffix)),
      get_slice(make_varname("rct_ob", "sex", suffix)),
      get_slice(make_varname("rct_ob", "treat", suffix)),
      get_slice(make_varname("rct_ob", "outcome", suffix))
    )
    ob_list <- list(
      get_slice(make_varname("ob_ex", "age", suffix)),
      get_slice(make_varname("ob_ex", "sex", suffix)),
      get_slice(make_varname("ob_ex", "treat", suffix)),
      get_slice(make_varname("ob_ex", "outcome", suffix))
    )
    ex_list <- list(
      get_slice(make_varname("rct_ex", "age", suffix)),
      get_slice(make_varname("rct_ex", "sex", suffix)),
      get_slice(make_varname("rct_ex", "treat", suffix)),
      get_slice(make_varname("rct_ex", "outcome", suffix))
    )
  } else {
    stop("Suffix not recognized.")
  }

  list(RCT = rct_list, Observational = ob_list, External = ex_list)
}

# mean smd plot
create_mean_smd_plot <- function(scenario_num, suffix) {
  all_lists <- get_smd_lists(scenario_num, suffix)
  smd_mean <- do.call(rbind, lapply(names(all_lists), function(type_name) {
    vars_list <- all_lists[[type_name]]
    means <- sapply(vars_list, function(x) {
      v <- unlist(x)
      if (length(v) == 0) return(NA_real_)
      mean(v[is.finite(v)], na.rm = TRUE)
    })
    data.frame(label = variables, type = type_name, SMD = means, stringsAsFactors = FALSE)
  }))

  # drop rows with NA SMD
  smd_mean <- smd_mean[is.finite(smd_mean$SMD), ]
  smd_mean$type <- factor(smd_mean$type, levels = types)

  ggplot(smd_mean, aes(y = label, x = SMD, shape = type)) +
    geom_point(size = 3, position = position_dodgev(height = 0.5)) +
    ggtitle(paste0("SMD - ",
                   ifelse(suffix == "_all", "Original Data",
                          ifelse(suffix == "_sim", "Simulated Combos", "Synthetic Combos")),
                   " - Scenario ", scenario_num)) +
    xlab("Average Standardised Mean Difference (SMD)") +
    scale_shape_manual(values = c(0, 5, 2)) +
    ylab("Variable") +
    xlim(-0.5, 0.5) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"))
}

# color boxplot
create_boxplot_color <- function(scenario_num, suffix) {
  all_lists <- get_smd_lists(scenario_num, suffix)

  rows <- list()
  for (type_name in names(all_lists)) {
    vars_list <- all_lists[[type_name]]
    for (j in seq_along(vars_list)) {
      raw_vals <- unlist(vars_list[[j]])
      raw_vals <- raw_vals[is.finite(raw_vals)]
      if (length(raw_vals) == 0) next
      df_j <- data.frame(
        label = variables[j],
        type = type_name,
        SMD = raw_vals,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- df_j
    }
  }

  if (length(rows) == 0) {
    warning("No finite SMD values found for scenario ", scenario_num, " suffix ", suffix)
    return(ggplot() + ggtitle("No data"))
  }

  boxplot_df <- do.call(rbind, rows)
  boxplot_df$type <- factor(boxplot_df$type, levels = types)
  boxplot_df$label <- factor(boxplot_df$label, levels = variables)

  ggplot(boxplot_df, aes(x = label, y = SMD, fill = type)) +
    geom_boxplot(color = "black", alpha = 0.9, outlier.size = 1, coef = 1.5) +
    coord_cartesian(ylim = c(-1, 1)) +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.background = element_blank()
    ) +
    ggtitle(paste0("Standard Mean Difference by Data Type - ",
                   ifelse(suffix == "_all", "Original Data",
                          ifelse(suffix == "_sim", "Simulated Combos", "Synthetic Combos")),
                   " - Scenario ", scenario_num)) +
    xlab("Variables") +
    guides(fill = guide_legend(title = "Data Type"))
}

# run and save plots
scenarios <- 1:4
suffixes <- c("_all", "_sim", "_syn")
fig_dir <- "C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Aug/CR_RS"

for (suffix in suffixes) {
  mean_plots <- lapply(scenarios, create_mean_smd_plot, suffix = suffix)
  boxplots_color <- lapply(scenarios, create_boxplot_color, suffix = suffix)

  jpeg(file.path(fig_dir, paste0("AllScenarios_mean_smd_plot", suffix, ".jpg")), width = 1200, height = 800)
  grid.arrange(grobs = mean_plots, ncol = 2, nrow = 2)
  dev.off()

  jpeg(file.path(fig_dir, paste0("AllScenarios_smdbox_color", suffix, ".jpg")), width = 1200, height = 800)
  grid.arrange(grobs = boxplots_color, ncol = 2, nrow = 2)
  dev.off()
}












# create df from response list
make_resp_df <- function(resp_list, data_type, method) {
  data.frame(
    rate = unlist(resp_list),
    iteration = seq_along(unlist(resp_list)),
    data_type = data_type,
    method = method
  )
}

# combine all response data into one data frame
resp_df <- dplyr::bind_rows(
  make_resp_df(sim_resp_rct, "RCT", "Sim"),
  make_resp_df(syn_resp_rct, "RCT", "Syn"),
  make_resp_df(sim_resp_ob,  "OB",  "Sim"),
  make_resp_df(syn_resp_ob,  "OB",  "Syn"),
  make_resp_df(sim_resp_ex,  "EX",  "Sim"),
  make_resp_df(syn_resp_ex,  "EX",  "Syn")
)

# assign scenarios to resp df
resp_df <- resp_df %>%
  dplyr::mutate(
    scenario = ceiling(iteration / (nSim / 4))
  )

# summarise mean and se
summary_df <- resp_df %>%
  dplyr::group_by(data_type, method, scenario) %>%
  dplyr::summarise(
    mean_rate = mean(rate, na.rm = TRUE),
    se = sd(rate, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    scenario_label = factor(paste("Scenario", scenario),
                            levels = paste("Scenario", 1:4))
  )

# plot
p <- ggplot(summary_df, aes(x = method, y = mean_rate, color = data_type)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = mean_rate - se, ymax = mean_rate + se),
                position = position_dodge(width = 0.5), width = 0.2) +
  facet_wrap(~scenario_label, ncol = 1) +
  labs(title = "Control Response Rate by Method, Data Type, and Scenario",
       x = "Method",
       y = "Mean Control Response Rate",
       color = "Data Type") +
  scale_color_brewer(palette = "Dark2") +

  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    panel.spacing = unit(1.2, "lines")
  )

ggsave("C:/Users/c3058452/OneDrive - Newcastle University/Work in Progress/Figures_Aug/CR_RS/control_response_rate_means.jpg",
       p, width = 6, height = 12)











# get exact values for reporting:

# function to split response rate values into sep dfs:
split_to_dfs <- function(values, k = 4, prefix = "df") {
  # k is number of groups
  split_values <- split(values, cut(seq_along(values), k, labels = FALSE))

  # loops
  for (i in seq_along(split_values)) {
    df_name <- paste0(prefix, i)
    assign(df_name, data.frame(value = split_values[[i]]), envir = .GlobalEnv)
  }
}

split_to_dfs(sim_resp_rct, k=4)
sim_resp_rct_1 <- data.frame(value = as.vector(t(df1)))
sim_resp_rct_2 <- data.frame(value = as.vector(t(df2)))
sim_resp_rct_3 <- data.frame(value = as.vector(t(df3)))
sim_resp_rct_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(sim_resp_ob, k=4)
sim_resp_ob_1 <- data.frame(value = as.vector(t(df1)))
sim_resp_ob_2 <- data.frame(value = as.vector(t(df2)))
sim_resp_ob_3 <- data.frame(value = as.vector(t(df3)))
sim_resp_ob_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(sim_resp_ex, k=4)
sim_resp_ex_1 <- data.frame(value = as.vector(t(df1)))
sim_resp_ex_2 <- data.frame(value = as.vector(t(df2)))
sim_resp_ex_3 <- data.frame(value = as.vector(t(df3)))
sim_resp_ex_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(syn_resp_rct, k=4)
syn_resp_rct_1 <- data.frame(value = as.vector(t(df1)))
syn_resp_rct_2 <- data.frame(value = as.vector(t(df2)))
syn_resp_rct_3 <- data.frame(value = as.vector(t(df3)))
syn_resp_rct_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(syn_resp_ob, k=4)
syn_resp_ob_1 <- data.frame(value = as.vector(t(df1)))
syn_resp_ob_2 <- data.frame(value = as.vector(t(df2)))
syn_resp_ob_3 <- data.frame(value = as.vector(t(df3)))
syn_resp_ob_4 <- data.frame(value = as.vector(t(df4)))

split_to_dfs(syn_resp_ex, k=4)
syn_resp_ex_1 <- data.frame(value = as.vector(t(df1)))
syn_resp_ex_2 <- data.frame(value = as.vector(t(df2)))
syn_resp_ex_3 <- data.frame(value = as.vector(t(df3)))
syn_resp_ex_4 <- data.frame(value = as.vector(t(df4)))



mean(sim_resp_rct_1$value)
mean(sim_resp_rct_2$value)
mean(sim_resp_rct_3$value)
mean(sim_resp_rct_4$value)
sd(sim_resp_rct_1$value)
sd(sim_resp_rct_2$value)
sd(sim_resp_rct_3$value)
sd(sim_resp_rct_4$value)
range(sim_resp_rct_1$value)
range(sim_resp_rct_2$value)
range(sim_resp_rct_3$value)
range(sim_resp_rct_4$value)

mean(sim_resp_ob_1$value)
mean(sim_resp_ob_2$value)
mean(sim_resp_ob_3$value)
mean(sim_resp_ob_4$value)
sd(sim_resp_ob_1$value)
sd(sim_resp_ob_2$value)
sd(sim_resp_ob_3$value)
sd(sim_resp_ob_4$value)
range(sim_resp_ob_1$value)
range(sim_resp_ob_2$value)
range(sim_resp_ob_3$value)
range(sim_resp_ob_4$value)

mean(sim_resp_ex_1$value)
mean(sim_resp_ex_2$value)
mean(sim_resp_ex_3$value)
mean(sim_resp_ex_4$value)
sd(sim_resp_ex_1$value)
sd(sim_resp_ex_2$value)
sd(sim_resp_ex_3$value)
sd(sim_resp_ex_4$value)
range(sim_resp_ex_1$value)
range(sim_resp_ex_2$value)
range(sim_resp_ex_3$value)
range(sim_resp_ex_4$value)




mean(syn_resp_rct_1$value)
mean(syn_resp_rct_2$value)
mean(syn_resp_rct_3$value)
mean(syn_resp_rct_4$value)
sd(syn_resp_rct_1$value)
sd(syn_resp_rct_2$value)
sd(syn_resp_rct_3$value)
sd(syn_resp_rct_4$value)
range(syn_resp_rct_1$value)
range(syn_resp_rct_2$value)
range(syn_resp_rct_3$value)
range(syn_resp_rct_4$value)

mean(syn_resp_ob_1$value)
mean(syn_resp_ob_2$value)
mean(syn_resp_ob_3$value)
mean(syn_resp_ob_4$value)
sd(syn_resp_ob_1$value)
sd(syn_resp_ob_2$value)
sd(syn_resp_ob_3$value)
sd(syn_resp_ob_4$value)
range(syn_resp_ob_1$value)
range(syn_resp_ob_2$value)
range(syn_resp_ob_3$value)
range(syn_resp_ob_4$value)

mean(syn_resp_ex_1$value)
mean(syn_resp_ex_2$value)
mean(syn_resp_ex_3$value)
mean(syn_resp_ex_4$value)
sd(syn_resp_ex_1$value)
sd(syn_resp_ex_2$value)
sd(syn_resp_ex_3$value)
sd(syn_resp_ex_4$value)
range(syn_resp_ex_1$value)
range(syn_resp_ex_2$value)
range(syn_resp_ex_3$value)
range(syn_resp_ex_4$value)


test_resp_rct_1 <- tsum.test(
  mean.x = mean(sim_resp_rct_1$value), s.x = sd(sim_resp_rct_1$value), n.x = length(sim_resp_rct_1$value),
  mean.y = mean(syn_resp_rct_1$value), s.y = sd(syn_resp_rct_1$value), n.y = length(syn_resp_rct_1$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_rct_2 <- tsum.test(
  mean.x = mean(sim_resp_rct_2$value), s.x = sd(sim_resp_rct_2$value), n.x = length(sim_resp_rct_2$value),
  mean.y = mean(syn_resp_rct_2$value), s.y = sd(syn_resp_rct_2$value), n.y = length(syn_resp_rct_2$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_rct_3 <- tsum.test(
  mean.x = mean(sim_resp_rct_3$value), s.x = sd(sim_resp_rct_3$value), n.x = length(sim_resp_rct_3$value),
  mean.y = mean(syn_resp_rct_3$value), s.y = sd(syn_resp_rct_3$value), n.y = length(syn_resp_rct_3$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_rct_4 <- tsum.test(
  mean.x = mean(sim_resp_rct_4$value), s.x = sd(sim_resp_rct_4$value), n.x = length(sim_resp_rct_4$value),
  mean.y = mean(syn_resp_rct_4$value), s.y = sd(syn_resp_rct_4$value), n.y = length(syn_resp_rct_4$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)




test_resp_ob_1 <- tsum.test(
  mean.x = mean(sim_resp_ob_1$value), s.x = sd(sim_resp_ob_1$value), n.x = length(sim_resp_ob_1$value),
  mean.y = mean(syn_resp_ob_1$value), s.y = sd(syn_resp_ob_1$value), n.y = length(syn_resp_ob_1$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ob_2 <- tsum.test(
  mean.x = mean(sim_resp_ob_2$value), s.x = sd(sim_resp_ob_2$value), n.x = length(sim_resp_ob_2$value),
  mean.y = mean(syn_resp_ob_2$value), s.y = sd(syn_resp_ob_2$value), n.y = length(syn_resp_ob_2$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ob_3 <- tsum.test(
  mean.x = mean(sim_resp_ob_3$value), s.x = sd(sim_resp_ob_3$value), n.x = length(sim_resp_ob_3$value),
  mean.y = mean(syn_resp_ob_3$value), s.y = sd(syn_resp_ob_3$value), n.y = length(syn_resp_ob_3$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ob_4 <- tsum.test(
  mean.x = mean(sim_resp_ob_4$value), s.x = sd(sim_resp_ob_4$value), n.x = length(sim_resp_ob_4$value),
  mean.y = mean(syn_resp_ob_4$value), s.y = sd(syn_resp_ob_4$value), n.y = length(syn_resp_ob_4$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)




test_resp_ex_1 <- tsum.test(
  mean.x = mean(sim_resp_ex_1$value), s.x = sd(sim_resp_ex_1$value), n.x = length(sim_resp_ex_1$value),
  mean.y = mean(syn_resp_ex_1$value), s.y = sd(syn_resp_ex_1$value), n.y = length(syn_resp_ex_1$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ex_2 <- tsum.test(
  mean.x = mean(sim_resp_ex_2$value), s.x = sd(sim_resp_ex_2$value), n.x = length(sim_resp_ex_2$value),
  mean.y = mean(syn_resp_ex_2$value), s.y = sd(syn_resp_ex_2$value), n.y = length(syn_resp_ex_2$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ex_3 <- tsum.test(
  mean.x = mean(sim_resp_ex_3$value), s.x = sd(sim_resp_ex_3$value), n.x = length(sim_resp_ex_3$value),
  mean.y = mean(syn_resp_ex_3$value), s.y = sd(syn_resp_ex_3$value), n.y = length(syn_resp_ex_3$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)

test_resp_ex_4 <- tsum.test(
  mean.x = mean(sim_resp_ex_4$value), s.x = sd(sim_resp_ex_4$value), n.x = length(sim_resp_ex_4$value),
  mean.y = mean(syn_resp_ex_4$value), s.y = sd(syn_resp_ex_4$value), n.y = length(syn_resp_ex_4$value),
  alternative = "two.sided",
  var.equal = FALSE,
  conf.level = 0.95
)


# Extract t-test results with p-values, confidence intervals, and effect sizes
ttest_pvalues <- data.frame(
  Test = c("RCT_1", "RCT_2", "RCT_3", "RCT_4", "OB_1", "OB_2", "OB_3", "OB_4", "EX_1", "EX_2", "EX_3", "EX_4"),
  P_Value = c(test_resp_rct_1$p.value, test_resp_rct_2$p.value, test_resp_rct_3$p.value, test_resp_rct_4$p.value,
              test_resp_ob_1$p.value, test_resp_ob_2$p.value, test_resp_ob_3$p.value, test_resp_ob_4$p.value,
              test_resp_ex_1$p.value, test_resp_ex_2$p.value, test_resp_ex_3$p.value, test_resp_ex_4$p.value),
  CI_Lower = c(test_resp_rct_1$conf.int[1], test_resp_rct_2$conf.int[1], test_resp_rct_3$conf.int[1], test_resp_rct_4$conf.int[1],
               test_resp_ob_1$conf.int[1], test_resp_ob_2$conf.int[1], test_resp_ob_3$conf.int[1], test_resp_ob_4$conf.int[1],
               test_resp_ex_1$conf.int[1], test_resp_ex_2$conf.int[1], test_resp_ex_3$conf.int[1], test_resp_ex_4$conf.int[1]),
  CI_Upper = c(test_resp_rct_1$conf.int[2], test_resp_rct_2$conf.int[2], test_resp_rct_3$conf.int[2], test_resp_rct_4$conf.int[2],
               test_resp_ob_1$conf.int[2], test_resp_ob_2$conf.int[2], test_resp_ob_3$conf.int[2], test_resp_ob_4$conf.int[2],
               test_resp_ex_1$conf.int[2], test_resp_ex_2$conf.int[2], test_resp_ex_3$conf.int[2], test_resp_ex_4$conf.int[2]),
  T_Statistic = c(test_resp_rct_1$statistic, test_resp_rct_2$statistic, test_resp_rct_3$statistic, test_resp_rct_4$statistic,
                  test_resp_ob_1$statistic, test_resp_ob_2$statistic, test_resp_ob_3$statistic, test_resp_ob_4$statistic,
                  test_resp_ex_1$statistic, test_resp_ex_2$statistic, test_resp_ex_3$statistic, test_resp_ex_4$statistic)
)
write.csv(ttest_pvalues, file = "CROHNS_ttest_pvalues_rs.csv", row.names = FALSE)



# for smd:
split_to_dfs(rct_age_smd_all, k=4)
rct_age_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_age_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_age_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_age_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(rct_sex_smd_all, k=4)
rct_sex_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_sex_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_sex_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_sex_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(rct_treat_smd_all, k=4)
rct_treat_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_treat_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_treat_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_treat_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(rct_outcome_smd_all, k=4)
rct_outcome_smd_1 <- data.frame(value = as.vector(t(df1)))
rct_outcome_smd_2 <- data.frame(value = as.vector(t(df2)))
rct_outcome_smd_3 <- data.frame(value = as.vector(t(df3)))
rct_outcome_smd_4 <- data.frame(value = as.vector(t(df4)))


split_to_dfs(ob_age_smd_all, k=4)
ob_age_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_age_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_age_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_age_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ob_sex_smd_all, k=4)
ob_sex_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_sex_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_sex_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_sex_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ob_treat_smd_all, k=4)
ob_treat_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_treat_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_treat_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_treat_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ob_outcome_smd_all, k=4)
ob_outcome_smd_1 <- data.frame(value = as.vector(t(df1)))
ob_outcome_smd_2 <- data.frame(value = as.vector(t(df2)))
ob_outcome_smd_3 <- data.frame(value = as.vector(t(df3)))
ob_outcome_smd_4 <- data.frame(value = as.vector(t(df4)))



split_to_dfs(ex_age_smd_all, k=4)
ex_age_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_age_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_age_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_age_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ex_sex_smd_all, k=4)
ex_sex_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_sex_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_sex_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_sex_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ex_treat_smd_all, k=4)
ex_treat_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_treat_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_treat_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_treat_smd_4 <- data.frame(value = as.vector(t(df4)))
split_to_dfs(ex_outcome_smd_all, k=4)
ex_outcome_smd_1 <- data.frame(value = as.vector(t(df1)))
ex_outcome_smd_2 <- data.frame(value = as.vector(t(df2)))
ex_outcome_smd_3 <- data.frame(value = as.vector(t(df3)))
ex_outcome_smd_4 <- data.frame(value = as.vector(t(df4)))



# just median, range, sd

range(rct_age_smd_1$value)
median(rct_age_smd_1$value)
sd(rct_age_smd_1$value)
range(rct_age_smd_2$value)
median(rct_age_smd_2$value)
sd(rct_age_smd_2$value)
range(rct_age_smd_3$value)
median(rct_age_smd_3$value)
sd(rct_age_smd_3$value)
range(rct_age_smd_4$value)
median(rct_age_smd_4$value)
sd(rct_age_smd_4$value)

range(rct_sex_smd_1$value)
median(rct_sex_smd_1$value)
sd(rct_sex_smd_1$value)
range(rct_sex_smd_2$value)
median(rct_sex_smd_2$value)
sd(rct_sex_smd_2$value)
range(rct_sex_smd_3$value)
median(rct_sex_smd_3$value)
sd(rct_sex_smd_3$value)
range(rct_sex_smd_4$value)
median(rct_sex_smd_4$value)
sd(rct_sex_smd_4$value)


range(rct_treat_smd_1$value)
median(rct_treat_smd_1$value)
sd(rct_treat_smd_1$value)
range(rct_treat_smd_2$value)
median(rct_treat_smd_2$value)
sd(rct_treat_smd_2$value)
range(rct_treat_smd_3$value)
median(rct_treat_smd_3$value)
sd(rct_treat_smd_3$value)
range(rct_treat_smd_4$value)
median(rct_treat_smd_4$value)
sd(rct_treat_smd_4$value)



range(rct_outcome_smd_1$value)
median(rct_outcome_smd_1$value)
sd(rct_outcome_smd_1$value)
range(rct_outcome_smd_2$value)
median(rct_outcome_smd_2$value)
sd(rct_outcome_smd_2$value)

range(rct_outcome_smd_3$value)
median(rct_outcome_smd_3$value)
sd(rct_outcome_smd_3$value)
range(rct_outcome_smd_4$value)
median(rct_outcome_smd_4$value)
sd(rct_outcome_smd_4$value)


range(ob_age_smd_1$value)
median(ob_age_smd_1$value)
sd(ob_age_smd_1$value)
range(ob_age_smd_2$value)
median(ob_age_smd_2$value)
sd(ob_age_smd_2$value)
range(ob_age_smd_3$value)
median(ob_age_smd_3$value)
sd(ob_age_smd_3$value)
range(ob_age_smd_4$value)
median(ob_age_smd_4$value)
sd(ob_age_smd_4$value)

range(ob_sex_smd_1$value)
median(ob_sex_smd_1$value)
sd(ob_sex_smd_1$value)
range(ob_sex_smd_2$value)
median(ob_sex_smd_2$value)
sd(ob_sex_smd_2$value)
range(ob_sex_smd_3$value)
median(ob_sex_smd_3$value)
sd(ob_sex_smd_3$value)
range(ob_sex_smd_4$value)
median(ob_sex_smd_4$value)
sd(ob_sex_smd_4$value)


range(ob_treat_smd_1$value)
median(ob_treat_smd_1$value)
sd(ob_treat_smd_1$value)
range(ob_treat_smd_2$value)
median(ob_treat_smd_2$value)
sd(ob_treat_smd_2$value)
range(ob_treat_smd_3$value)
median(ob_treat_smd_3$value)
sd(ob_treat_smd_3$value)
range(ob_treat_smd_4$value)
median(ob_treat_smd_4$value)
sd(ob_treat_smd_4$value)


range(ob_outcome_smd_1$value)
median(ob_outcome_smd_1$value)
sd(ob_outcome_smd_1$value)
range(ob_outcome_smd_2$value)
median(ob_outcome_smd_2$value)
sd(ob_outcome_smd_2$value)
range(ob_outcome_smd_3$value)
median(ob_outcome_smd_3$value)
sd(ob_outcome_smd_3$value)

range(ob_outcome_smd_4$value)
median(ob_outcome_smd_4$value)
sd(ob_outcome_smd_4$value)



range(ex_age_smd_1$value)
median(ex_age_smd_1$value)
sd(ex_age_smd_1$value)
range(ex_age_smd_2$value)
median(ex_age_smd_2$value)
sd(ex_age_smd_2$value)
range(ex_age_smd_3$value)
median(ex_age_smd_3$value)
sd(ex_age_smd_3$value)
range(ex_age_smd_4$value)
median(ex_age_smd_4$value)
sd(ex_age_smd_4$value)

range(ex_sex_smd_1$value)
median(ex_sex_smd_1$value)
sd(ex_sex_smd_1$value)
range(ex_sex_smd_2$value)
median(ex_sex_smd_2$value)
sd(ex_sex_smd_2$value)
range(ex_sex_smd_3$value)
median(ex_sex_smd_3$value)
sd(ex_sex_smd_3$value)
range(ex_sex_smd_4$value)
median(ex_sex_smd_4$value)
sd(ex_sex_smd_4$value)


range(ex_treat_smd_1$value)
median(ex_treat_smd_1$value)
sd(ex_treat_smd_1$value)
range(ex_treat_smd_2$value)
median(ex_treat_smd_2$value)
sd(ex_treat_smd_2$value)
range(ex_treat_smd_3$value)
median(ex_treat_smd_3$value)
sd(ex_treat_smd_3$value)
range(ex_treat_smd_4$value)
median(ex_treat_smd_4$value)
sd(ex_treat_smd_4$value)


range(ex_outcome_smd_1$value)
median(ex_outcome_smd_1$value)
sd(ex_outcome_smd_1$value)
range(ex_outcome_smd_2$value)
median(ex_outcome_smd_2$value)
sd(ex_outcome_smd_2$value)
range(ex_outcome_smd_3$value)
median(ex_outcome_smd_3$value)
sd(ex_outcome_smd_3$value)
range(ex_outcome_smd_4$value)
median(ex_outcome_smd_4$value)
sd(ex_outcome_smd_4$value)





















































