## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).
## To get access to the database contact the corresponding author.

## The code below runs the models that look at which males groom. The code is bootstrapped to deal were paternity of the previous and next offspring of the mother are missing.

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

## number of interations
nr_draws = 1000

## Detect how many cores
## If you run on personal computer specify the number of cores
ncores <- detectCores()

setwd('./who_groomed_bootstrap')

latest_version_date <- read_lines("./data/latest_version_date.txt")
who_grooms_bootstrap <- read_csv(
  paste("./data/who_grooms_bootstrap", latest_version_date, ", .csv")
  )

## create the function that will be bootstrapped
## Include all the variables that should be included in the model
run_bootstrap_model <- function() {
   # prep the data
  temp_data <- who_grooms_bootstrap %>%
    ## basic information
    filter(has_consort_data) %>%
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
           , does_groom,  paternal_res_i_adj_zscored

           ## create random variable
           ## This then uses the randomly assigned paternity in
           ## 7. who_groomed bootstrap
           , previous_random := !!sym(paste0("previous_", ii))
           , next_random := !!sym(paste0("next_", ii))) %>%
    ## If patenity is know use that info otherwise you 'random' value
    mutate(bootstrap_previous =
             if_else(previous_kid_with_mom == .5, previous_random, previous_kid_with_mom),
           bootstrap_next =
             if_else(next_kid_with_mom == .5, next_random, next_kid_with_mom)) %>%
    select(-previous_random, -next_random) %>%
    filter(!is.na(observer_effort)) %>%
    rename(dyad_does_groom = does_groom) %>%
    mutate(consort_pa = estrous_c != 0)

  print(ii)

  ##############################################################################
  ##############################################################################
  ## Full models
  print("who grooms full models")

  ## create the function
  who_grooms_full_model_prop_consort_formula <- dyad_does_groom ~ 1 +
    observer_effort +
    kid_age+ mom_age+ AMales_age +
    cumulative_adversity +
    is_dad + mean_ordrank +
    proportion_consort_time +
    daily_d_days_rate + offspring_years + nr_pdads +
    bootstrap_previous + bootstrap_next +
    (1|focal) +  (1|AMales)

  ## run the model
  who_grooms_full_model_prop_consort <- glmer(
    who_grooms_full_model_prop_consort_formula,
    data = temp_data
    , weights = log(nr_days)
    , family = binomial
    , control = glmerControl(optimizer = "bobyqa",
                             optCtrl=list(maxfun=2e5))
    , na.action = 'na.fail')
  # extract data from model
  who_grooms_full_model_prop_consort_results <- who_grooms_full_model_prop_consort %>%
    broom.mixed::tidy() %>%
    left_join(vif.mer(who_grooms_full_model_prop_consort), by = 'term') %>%
    mutate(run = paste0("model_", ii))

  ## write results of i model round to table
  ## round are appeneded
  write.table(who_grooms_full_model_prop_consort_results,
              file = paste("./data/who_grooms_full_model_prop_consort_results_raw_",
                           latest_version_date,".txt")
              ,  append = TRUE, col.names = FALSE)

   ##############################################################################
  ##############################################################################
  ## similiar idea but now only for fathers

  print("full model only dad")
  print("using proportion_consort_time")
  ## update the dataset
  temp_data_only_dad <- temp_data %>%
    filter(is_dad == TRUE)

  ## update the mldep
  who_grooms_full_model_consort_prop_only_dad_formula <-
    update(who_grooms_full_model_prop_consort_formula, . ~ . - is_dad)

  ## run model
who_grooms_full_model_consort_prop_only_dad <-
  glmer(who_grooms_full_model_consort_prop_only_dad_formula,
        data = temp_data_only_dad
        , weights = log(nr_days)
        , family = binomial
        , control = glmerControl(optimizer = "bobyqa",
                                 optCtrl=list(maxfun=2e5))
        , na.action = 'na.fail')

who_grooms_full_model_consort_prop_only_dad_results <- who_grooms_full_model_consort_prop_only_dad %>%
    broom.mixed::tidy() %>%
    left_join(vif.mer(who_grooms_full_model_consort_prop_only_dad), by = 'term') %>%
    mutate(run = paste0("model_", ii))

  write.table(who_grooms_full_model_consort_prop_only_dad_results
              , file = paste0("./data/who_grooms_full_model_consort_prop_",
                              "only_dad_results_raw",
                              latest_version_date,
                              ".txt")
              ,  append = TRUE, col.names = FALSE)
}


## create a cluster adn run the actual models
cl <- makeCluster(ncores,  outfile="")
registerDoSNOW(cl)
clusterExport(cl, list("vif.mer"
                       , "run_bootstrap_model"));

foreach(ii = 1:nr_draws,
        .packages = c('Matrix', 'lmerTest', 'tidyverse')) %dopar% {
          run_bootstrap_model()
        }
stopCluster(cl)

## after running the code more the tables to local comouter

