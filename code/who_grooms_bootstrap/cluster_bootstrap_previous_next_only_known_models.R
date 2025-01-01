## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).
## To get access to the database contact the corresponding author.

## The code below runs the models that look at which males groom. It only includes cases were paternity of the previous and next offspring of the mother are know.

# libraries
packages <- c(
  "doSnow"
  , "foreach"
  , "parallel"
  , "lmerTest"
  , "tidyverse")

## install packages if needed and open libaries
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
}
lapply(packages, require, character.only = TRUE, warn.conflicts = TRUE,
       quietly = FALSE)

## open the file with fuctions
source('./code/3. functions_for_data_prep.R')

setwd('./who_groomed_bootstrap')

latest_version_date <- read_lines("./data/latest_version_date.txt")
who_grooms_bootstrap <- read_csv(
  paste("./data/who_grooms_bootstrap", latest_version_date, ", .csv"))


  temp_data <- who_grooms_bootstrap %>%
    filter(has_consort_data) %>%
    ## basic information
    select(focal, focal_birth, focal_grp,  start, end, nr_days
           , mom, dad, AMales,
           ## ages
           , kid_age, mom_age, AMales_age
           ## early adversity
           , maternal_loss, maternal_rank, maternal_SCI_F
           , density, sibling, drought, cumulative_adversity
           , observer_effort
           ## male details
           , is_dad, mean_ordrank, mean_proprank
           , contains('estrous')
           , proportion_consort_time
           , daily_d_days_rate, offspring_years, nr_pdads
           , previous_kid_with_mom, next_kid_with_mom ## see bootstrap

           ## responses
           , does_groom,  paternal_res_i_adj_zscored)

           ## create random
    mutate(bootstrap_previous =
             if_else(previous_kid_with_mom == .5, previous_random, previous_kid_with_mom),
           bootstrap_next =
             if_else(next_kid_with_mom == .5, next_random, next_kid_with_mom)) %>%
    select(-previous_random, -next_random) %>%
    filter(!is.na(observer_effort)) %>%
    rename(dyad_does_groom = does_groom) %>%
    mutate(consort_pa = estrous_c != 0) %>%
    filter(has_consort_data)

  print(ii)

  temp_data_only_known <- temp_data %>%
    filter(previous_kid_with_mom != 0.5) %>%
    filter(next_kid_with_mom != 0.5)

  ##############################################################################
  ##############################################################################
  ## Full models
  print("who grooms all males only known next and previous")
  # This would be table S4
  ## who grooms consort_prop

  who_grooms_only_known_model_prop_consort_formula <- dyad_does_groom ~ 1 +
    observer_effort +
    kid_age+ mom_age+ AMales_age +
    cumulative_adversity +
    is_dad + mean_ordrank +
    proportion_consort_time +
    daily_d_days_rate + offspring_years + nr_pdads +
    (1|focal) +  (1|AMales)

  who_grooms_only_known_model_prop_consort <- glmer(who_grooms_only_known_model_prop_consort_formula,
                                              data = temp_data
                                              , weights = log(nr_days)
                                              , family = binomial
                                              , control=
                                                glmerControl(optimizer = "bobyqa",
                                                             optCtrl=list(maxfun=2e5))
                                              , na.action = 'na.fail')

  who_grooms_only_known_model_prop_consort _results <- who_grooms_only_known_model_prop_consort %>%
    broom.mixed::tidy() %>%
    left_join(vif.mer(who_grooms_only_known_model_prop_consort), by = 'term') %>%
    mutate(run = paste0("model_", 0))

  write.table(who_grooms_only_known_model_prop_consort_results,
              file = "./data/who_grooms_only_known_model_prop_consort_results_raw_13SEP24.txt"
              ,  append = TRUE, col.names = FALSE)

   ##############################################################################
  ##############################################################################
  #Full model only dad (me)
  ## This would be table 3
  print("full model only dad")
  print("using proportion_consort_time")
  temp_data_only_known_only_dad <- temp_data_only_known %>%
    filter(is_dad == TRUE)

who_grooms_only_known_model_prop_consort_only_dad_formula <- update(who_grooms_only_known_model_prop_consort_formula,
                                              . ~ . - is_dad)

who_grooms_only_known_model_prop_consort_only_dad <- glmer(who_grooms_only_known_model_prop_consort_only_dad _formula,
                                            data = temp_data_only_known_only_dad
                                            , weights = log(nr_days)
                                            , family = binomial
                                            , control=
                                              glmerControl(optimizer = "bobyqa",
                                                           optCtrl=list(maxfun=2e5))
                                            , na.action = 'na.fail')

who_grooms_only_known_model_prop_consort_only_dad_results <- who_grooms_only_known_model_prop_consort_only_dad %>%
    broom.mixed::tidy() %>%
    left_join(vif.mer(who_grooms_only_known_model_prop_consort_only_dad), by = 'term') %>%
    mutate(run = paste0("model_", 0))

  write.table(who_grooms_only_known_model_prop_consort_only_dad_results
              , file = "./data/ who_grooms_only_known_model_prop_consort_only_dad_results_raw_13SEP24.txt"
              ,  append = TRUE, col.names = FALSE)
