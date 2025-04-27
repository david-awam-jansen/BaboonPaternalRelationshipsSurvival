## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).

## The code below was written to deal with some of the issue that were raised during review.
## There are some additional analysis, getting AIC values for some of tha analysis and an adaptation to some of the figures.


packages <- c(
    "tidyverse",
    "showtext",
    "cowplot",
    "survival",
    "lmerTest",
    "broom",
    "broom.mixed",
    "purrr",
    "flextable",
    "MuMIn",
    "parallel"
)

## install packages if needed and open libaries
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
}
lapply(packages, require, character.only = TRUE, warn.conflicts = TRUE,
       quietly = FALSE)

source('./code/3. functions_for_data_prep.R')

## Table 1 AIC values
Table1_data <- read_csv("./data/Table1_data.csv")

juvenile_model <- Table1_data %>%
    filter(str_detect(paternal_groom, "paternal")) %>%
    lmer(formula = bond_strength ~ juvenile_age + is_dad + (1|juvenile_id))

juvenile_model_results <- juvenile_model %>%
    broom.mixed::tidy() %>%
    filter(effect == 'fixed')

juvenile_model_results <- juvenile_model_results %>%
    select(term, estimate, std.error, statistic, df, p.value)

names(juvenile_model_results)[2:6] <- c('\U03B2','\U03C3', 't', 'df', 'p')

juv_mod_interpretation = c("",
                           "\U2191 Juvenile age \U2191 bond strength",
                           "\U2191 bond strength = male is father")

juvenile_model_table <- juvenile_model_results %>%
    mutate(term = c(
        "Intercept",
        "Age of juvenile",
        "Is the male the father?" )) %>%
    mutate(interpretation = c("",
                              "\U2191 Juvenile age \U2191 bond strength",
                              "\U2191 bond strength = male is father")) %>%
    mutate(p = format.pval(p, eps = .001, digits = 2)) %>%
    flextable::flextable() %>%
    flextable::colformat_double(j = 2:5,digits = 3, big.mark = "", decimal.mark = ".") %>%
    flextable::colformat_double(j = 6,digits = 3, big.mark = "", decimal.mark = ".") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::align(align = "left", part = "header") %>%
    #flextable::rotate(j = 2:5, rotation="btlr",part="header") %>%
    flextable::align(align = "left", part = "all") %>%
    flextable::width(j = 1, width= 1.7) %>%
    flextable::width(j = 2:4, width= .5) %>%
    flextable::width(j = 5, width= .7) %>%
    flextable::width(j = 7, width= 2)

AIC(juvenile_model
                  , update(juvenile_model, . ~ . - is_dad)
                  , update(juvenile_model, . ~ . - juvenile_age)) %>%
    as_tibble(rownames = "Model") %>%
    mutate(dAIC = round(AIC - AIC(juvenile_model), 2)) %>%
    mutate(AIC = round(AIC, 2)) %>%
    mutate(model = c(
        "Paternal bond strength ~  1 + Age of juvenile + Is male the father",
        "Paternal bond strength ~  1 + Age of juvenile",
        "Paternal bond strength ~  1 + Is male the father")) %>%
    select(model, AIC, dAIC) %>%
    flextable::flextable() %>%
    flextable::colformat_double(j = 2:3,digits = 2,
                                big.mark = "", decimal.mark = ".") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::align(align = "left", part = "header") %>%
    #flextable::rotate(j = 2:5, rotation="btlr",part="header") %>%
    flextable::align(align = "left", part = "all") %>%
    flextable::width(j = 1, width= 6.5)

##############################################################################
## Table 4
Table4_data <- read_csv("./data/Table4_data.csv")

Table4_data_only_known <- Table4_data %>%
    filter(previous_kid_with_mom != 0.5)

dad_overlap_only_known_model_prop_consort_formula <- dad_overlap_years ~ 1 +
    dad_age +
    mom_age +
    male_rank +
    daily_d_days_rate +
    proportion_consort_time +
    nr_pdads +
    previous_kid_with_mom +
    offspring_years +
    cumulative_adversity +
    (1|dad)

dad_overlap_only_known_model_prop_consort <- lmer(dad_overlap_only_known_model_prop_consort_formula
                                                  , data = Table4_data_only_known
                                                  , na.action = 'na.fail')

dad_overlap_only_known_model_prop_consort %>%  tidy()

dad_overlap_only_known_model_prop_consort_AIC <- dad_overlap_only_known_model_prop_consort  %>%  AIC()

dad_overlap_only_known_model_prop_consort_dredge <-
    dredge(dad_overlap_only_known_model_prop_consort)

dad_overlap_best_models <- subset(dad_overlap_only_known_model_prop_consort_dredge, delta <= 2)
full_model_index <- which.max(rowSums(!is.na(dad_overlap_only_known_model_prop_consort_dredge)))

dad_overlap_best_models_tibble <- dad_overlap_best_models %>%
    as_tibble() %>%
    select(AICc, delta, everything()) %>%
    select(-logLik, -df, -weight)

dad_overlap_best_models_tibble %>%
    pivot_longer(names_to = "variable",
             values_to = "value", `(Intercept)`:proportion_consort_time) %>%
    filter(!is.na(value)) %>%
    pivot_wider(names_from = "variable", values_from = "value") %>%
    flextable() %>%
    colformat_double(digits = 3)

dad_overlap_best_models_model.avg <- MuMIn::model.avg(dad_overlap_best_models, revised.var = TRUE)

summary(dad_overlap_best_models_model.avg)$coefmat.subset %>%
    as.data.frame() %>%                # Convert to data frame first
    rownames_to_column(var = "Variable") %>%  # Preserve variable names
    as_tibble() %>%
    flextable() %>%
    colformat_double(digits = 3)
#############################################################################
#############################################################################
## Table S4
## Full models
### This data is not publicaly available.
who_grooms_bootstrap <- read_csv("./data/who_grooms_bootstrap_13SEP24.csv")


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
           , does_groom,  paternal_res_i_adj_zscored) %>%
    filter(!is.na(observer_effort)) %>%
    rename(dyad_does_groom = does_groom) %>%
    mutate(consort_pa = estrous_c != 0)

temp_data_only_known <- temp_data %>%
    filter(previous_kid_with_mom != 0.5) %>%
    filter(next_kid_with_mom != 0.5)

who_grooms_bootstrap %>%
    group_by(focal, is_dad) %>%
     mutate(has_consort_data = !is.na(proportion_consort_time),
           has_paternity = max(previous_kid_with_mom %in% c(0,1) &
               next_kid_with_mom %in% c(0,1))==1) %>%
    select(focal, is_dad, has_consort_data, has_paternity) %>%
    distinct() %>%
    mutate(check =  case_when(
        has_consort_data == TRUE & has_paternity == TRUE ~ "has_both",
        has_consort_data == FALSE & has_paternity == FALSE ~ "missing_both",
        has_consort_data == FALSE & has_paternity == TRUE ~ "missing_consort",
        has_consort_data == TRUE & has_paternity == FALSE ~ "missing_paternity",
        TRUE ~ "check"
    )) %>%
    group_by(is_dad, check) %>%
    summarise(cases = n()) %>%
    mutate(is_dad = if_else(is_dad== TRUE, "dad", "other males")) %>%
    pivot_wider(names_from = is_dad, values_from = cases)

temp_data_only_known_only_dad <- temp_data_only_known %>%
    filter(is_dad == TRUE)

## This will not be in data repository
write_csv(temp_data_only_known, "./data/temp_data_only_known.csv")
write_csv(temp_data_only_known_only_dad, "./data/temp_data_only_known_only_dad")

temp_data_only_known <- read_csv("./data/temp_data_only_known.csv")

who_grooms_only_known_model_prop_consort_formula <- dyad_does_groom ~ 1 +
    observer_effort +
    kid_age+ mom_age+ AMales_age +
    cumulative_adversity +
    is_dad + mean_ordrank +
    proportion_consort_time +
    daily_d_days_rate + offspring_years + nr_pdads +
    previous_kid_with_mom + next_kid_with_mom +
    (1|focal) +  (1|AMales)

who_grooms_only_known_model_prop_consort_only_dad_formula <-
    update(who_grooms_only_known_model_prop_consort_formula,
           . ~ . - is_dad)

who_grooms_only_known_model_prop_consort_only_dad <-
    glmer(who_grooms_only_known_model_prop_consort_only_dad_formula,
          data = temp_data_only_known_only_dad
          , weights = log(nr_days)
          , family = binomial
          , control= glmerControl(optimizer = "bobyqa",
                                  optCtrl=list(maxfun=2e5))
          , na.action = 'na.fail')

who_grooms_only_known_model_prop_consort_only_dad_results <- who_grooms_only_known_model_prop_consort_only_dad %>%
    broom.mixed::tidy() %>%
    filter(effect == "fixed") %>%
    left_join(vif.mer(dad_overlap_only_known_model_prop_consort), by = 'term') %>%
    mutate(run = paste0("model_", 0)) %>%
    mutate(row = NA) %>%
    select(row, everything()) %>%
    mutate(term = map_chr(.x = term, fix_terms)) %>%
    mutate(order = if_else(str_detect(term, "Intercept"), 1, 2)) %>%
    arrange(order, p.value) %>%
    mutate(p.value = format.pval(p.value, eps = .001, digits = 2))

explain_text <- c(
    c(""),
    c("\U2191 maternal age  \U2191 co-residency"),
    c("\U2191 paternal age  \U2191 co-residency"),
    c("\U2191 offspring \U2191 co-residency"),
    c("\U2191 fertile females \U2191 co-residency"),
    c("trend: had offsping \U2191 co-residency"),
    rep(c(""),4)

)

t1 <- dad_overlap_only_known_model_prop_consort_results %>%
    filter(effect == 'fixed') %>%
    select(-effect, -group, -all_of(intersect(names(t1), "model_run"))) %>%
    mutate(term = case_when(str_detect(term , "bootstrap_previous") ~ "previous_kid_with_mom",
                            str_detect(term , "bootstrap_next") ~ "next_kid_with_mom",
                            TRUE ~ term)) %>%
    group_by(term) %>%
    summarise_all(.funs = median) %>%
    mutate(term = map_chr(.x = term, fix_terms)) %>%
    mutate(order = if_else(str_detect(term, "Intercept"), 1, 2)) %>%
    arrange(order, p.value) %>%
    mutate(p.value = format.pval(p.value, eps = .001, digits = 2)) %>%
    select(-order, -VIF) %>%
    mutate(interpretation = explain_text)

t1 <- t1 %>%  select(-df)

names(t1)[2:5] <- c('\U03B2','\U03C3', 't', 'p')

t1 %>%
    flextable() %>%
    colformat_double(digits = 3) %>%
    theme_vanilla() %>%
    fontsize(size = 8, part = "all") %>%
    align(align = "left", part = "header") %>%
    align(j = 2:5, align = "center", part = "body" ) %>%
    width(j = 1, width = 2.2, unit = "in") %>%
    width(j = 2:5, width = .5, unit = "in") %>%
    width(j = 6, width =3, unit = "in") %>%
    height(height = 2,part = "header") %>%
    valign(valign = "bottom", part = "header") %>%
    valign(valign = "top", part = "body")


## some of the models take long
## here we use parallel coding to speed it up
clu <- makeCluster(detectCores() - 2) # Use all cores except 1
clusterExport(clu, c("temp_data_only_known_only_dad",
                    "who_grooms_only_known_model_prop_consort_only_dad"))
clusterEvalQ(clu,library(lme4,logical.return =T))
clusterEvalQ(clu,library(MuMIn,logical.return =T))
who_grooms_only_known_model_prop_consort_only_dad_dredge <-
    dredge(who_grooms_only_known_model_prop_consort_only_dad,
           cluster = clu)
stopCluster(clu)

dad_grooms_best_models <- subset(who_grooms_only_known_model_prop_consort_only_dad_dredge, delta <= 2)

dad_grooms_best_models_tibble <- dad_grooms_best_models %>%
    as_tibble() %>%
    select(AICc, delta, everything()) %>%
    select(-logLik, -df, -weight)

dad_grooms_best_models_tibble %>%
    pivot_longer(names_to = "variable",
                 values_to = "value", `(Intercept)`:proportion_consort_time) %>%
    filter(!is.na(value)) %>%
    mutate(variable = case_when(
        variable == "cumulative_adversity" ~ "cum_adv",
        variable == "proportion_consort_time" ~ "prop_consort",
        variable == "daily_d_days_rate" ~ "d_days",
        variable == "mean"
        TRUE ~ variable)) %>%
    pivot_wider(names_from = "variable", values_from = "value") %>%
    flextable() %>%
    colformat_double(digits = 3)

dad_overlap_best_models_model.avg <- MuMIn::model.avg(dad_overlap_best_models, revised.var = TRUE)

summary(dad_overlap_best_models_model.avg)$coefmat.subset %>%
    as.data.frame() %>%                # Convert to data frame first
    rownames_to_column(var = "Variable") %>%  # Preserve variable names
    as_tibble() %>%
    flextable() %>%
    colformat_double(digits = 3)


dad_grooms_best_models <- subset(who_grooms_only_known_model_prop_consort_only_dad_dredge, delta <= 2)

dad_grooms_best_models_model.avg <- MuMIn::model.avg(dad_grooms_best_models, revised.var = TRUE)

summary(dad_grooms_best_models_model.avg)$coefmat.subset %>%
    as.data.frame() %>%                # Convert to data frame first
    rownames_to_column(var = "Variable") %>%  # Preserve variable names
    as_tibble() %>%
    flextable() %>%
    colformat_double(digits = 3)

dad_overlap_best_models_tibble %>%
    pivot_longer(names_to = "variable",
                 values_to = "value", `(Intercept)`:proportion_consort_time) %>%
    filter(!is.na(value)) %>%
    pivot_wider(names_from = "variable", values_from = "value") %>%
    flextable() %>%
    colformat_double(digits = 3)
#############################################################################
#############################################################################
## table S5

sci <- read_csv("./data/Table_S2_data.csv")


sciF_model <-
    lmer(formula = "SCI_F ~ jDSI_paternal + age_class + mean_rank + (1|focal)",
         na.action = "na.fail", data = sci)

sciM_model <-
    lmer(formula = "SCI_M ~ jDSI_paternal + age_class + mean_rank + (1|focal)",
         na.action = "na.fail", data = sci)

sciF_model_results <- sciF_model %>%
    tidy()   %>%
    mutate(term = case_when(term == "mean_rank" ~ "adult rank",
                            term == "age_class" ~ "adult age",
                            TRUE ~ term)) %>%
    filter(effect == 'fixed') %>%
    mutate(p.value = round(p.value, 3)) %>%
    mutate(p.value = format.pval(p.value, eps = .001, digits = 2)) %>%
    select(term, estimate, std.error, statistic, df, p.value) %>%
    mutate(term = case_when(
        term == "jDSI_paternal" ~ "Average grooming bond strength with father (in the juvenile period)",
        term == "adult_age" ~ "Female age in a given year of adulthood ",
        term == "adult_rank" ~ "Female dominance rank  in a given year of adulthood",
        TRUE ~ term))  %>%
    add_row(term = "Predictors of adult female social connectedness to other adult females",
            .before = 1)

sciM_model_results <- sciM_model %>%
    tidy() %>%
    filter(effect == 'fixed') %>%
    mutate(p.value = round(p.value, 3)) %>%
    mutate(p.value = format.pval(p.value, eps = .001, digits = 2)) %>%
    select(term, estimate, std.error, statistic, df, p.value) %>%
    mutate(term = case_when(
        term == "jDSI_paternal" ~ "Average grooming bond strength with father (in the juvenile period)",
        term == "age_class" ~ "Female age in a given year of adulthood ",
        term == "mean_rank" ~ "Female dominance rank  in a given year of adulthood",
        TRUE ~ term))   %>%
    select(term, estimate, std.error, statistic, df, p.value) %>%
    add_row(term = "Predictors of adult female social connectedness to adult males", .before = 1)

sci_model <- bind_rows(sciF_model_results, sciM_model_results)
names(sci_model)[2:6] <- c('\U03B2','\U03C3', 'r', 'df', 'p')

sci_model %>%
    rename(Term = term) %>%
    flextable() %>%
    merge_at(i = 1, j = NULL, part = "body") %>%
    merge_at(i = 6, j = NULL, part = "body") %>%
    bg(i = 1, bg = 'darkgrey') %>%
    bg(i = 6, bg = 'darkgrey') %>%
    colformat_double(digits = 2) %>%
    fontsize(size = 9, part = "all") %>%
    width(j = 1, width= 3) %>%
    width(j = 2:6, width= .5) %>%
    fontsize(size = 9, part = "all") %>%
    #theme_vanilla() %>%
    align(j = 1, align = "left") %>%
    align(j = 2:6, align = "right") %>%
    align(j = 2:6, align = "center", part = "header") %>%
    align(i = c(1, 6), align = "center") %>%
    autofit()


get_sci_model_details <- function(model, deltaAIC=2) {
    model_dregde <- dredge(model)
    models <- get.models(model_dregde, subset = delta < 2000)

    models %>%
        tibble() %>%
        rename(model = ".") %>%
        mutate(results = map(.x = model, .f = tidy, conf.int = TRUE)) %>%
        mutate(AICc = map(.x = model, .f = AICc)) %>%
        mutate(model_nr = 1:n()) %>%
        unnest(c(results, AICc)) %>%
        filter(effect == "fixed",
               term != "(Intercept)") %>%
        select(-model, -effect, - group) %>%
        mutate_if(is.numeric, round, digits = 3) %>%
        arrange(AICc) %>%
        mutate(dAICc = AICc - min(AICc)) %>%
        mutate(AICc = round(AICc, 2),
               dAICc = round(dAICc, 2)) %>%
        mutate(stat_value =
                   paste0(estimate, "\n (", conf.low," - ", conf.high,")")) %>%
        mutate(term = case_when(
            term == "mean_rank" ~ "Average female dominance rank in a given year of adulthood",
            term == "age_class" ~ "Female age at the start of a given year of adulthood",
            term == "jDSI_paternal" ~ "Mean DSIpaternal \n (in the juvenile period)")) %>%         select(model_nr, term, stat_value, AICc, dAICc) %>%
        pivot_wider(names_from = term, values_from = stat_value) %>%
        select(Model = model_nr,
               "Mean DSIpaternal \n (in the juvenile period)",
               "Female age at the start of a given year of adulthood",
               "Average female dominance rank in a given year of adulthood",
               AICc,
               dAICc
        ) %>%
        filter(dAICc <= deltaAIC) %>%
        flextable() %>%
        width(j = 2:4, 1.5) %>%
        align(j = 2:4, align = "center")
}

get_sci_model_details(sciF_model, deltaAIC = 2)

get_sci_model_details(sciM_model, deltaAIC =10)


#############################################################################
#############################################################################
## table S5
who_grooms_only_known_model_prop_consort <-
    glmer(who_grooms_only_known_model_prop_consort_formula,
          data = temp_data_only_known
          , weights = log(nr_days)
          , family = binomial
          , control= glmerControl(optimizer = "bobyqa",
                                  optCtrl=list(maxfun=2e5))
          , na.action = 'na.fail')

who_grooms_only_known_model_prop_consort_results <- who_grooms_only_known_model_prop_consort %>%
    broom.mixed::tidy() %>%
    filter(effect == "fixed") %>%
    left_join(vif.mer(who_grooms_only_known_model_prop_consort), by = 'term') %>%
    mutate(run = paste0("model_", 0)) %>%
    mutate(row = NA) %>%
    select(row, everything()) %>%
    mutate(coded_term = forcats::fct_relevel(term,
        "(Intercept)", "is_dadTRUE", "kid_age", "AMales_age", "mom_age",
        "mean_ordrank", "daily_d_days_rate", "proportion_consort_time",
        "nr_pdads", "next_kid_with_mom", "previous_kid_with_mom",
        "offspring_years", "cumulative_adversity", "observer_effort")) %>%
    mutate(term = map_chr(.x = term, fix_terms)) %>%
    mutate(order = if_else(str_detect(term, "Intercept"), 1, 2)) %>%
    arrange(order, p.value) %>%
    mutate(p.value = format.pval(p.value, eps = .001, digits = 2))

explain_text <- c(
    c(""),
    c("\U2191 juvenile age  \U2191 grooming"),
    c("\U2191 consort time  \U2191 grooming"),
    c("\U2191 effort  \U2191 grooming"),
    c("\U2191 male age  \U2191 grooming"),
    c("male is father \U2191 grooming"),
    c("\U2191 offspring \U2191 grooming"),
    c("\U2191 high ordinal rank \U2193 grooming"),
    c("\U2191 fertile females rank \U2193 grooming"),
    c("had offsping \U2191 grooming"),
    c("trent: will have offsping \U2193 grooming"),
    rep(c(""),3)

)

t1 <- who_grooms_only_known_model_prop_consort_results %>%
    select(-row, -effect, -group, -run, -order, -VIF) %>%
    mutate(interpretation = explain_text)

names(t1)[2:5] <- c('\U03B2','\U03C3', 't', 'p')

t1 %>%
    arrange(coded_term) %>%
    select(-coded_term) %>%
    flextable() %>%
    colformat_double(digits = 3) %>%
    theme_vanilla() %>%
    fontsize(size = 8, part = "all") %>%
    align(align = "left", part = "header") %>%
    align(j = 2:5, align = "center", part = "body" ) %>%
    width(j = 1, width = 2.2, unit = "in") %>%
    width(j = 2:5, width = .5, unit = "in") %>%
    width(j = 6, width =3, unit = "in") %>%
    height(height = 2,part = "header") %>%
    valign(valign = "bottom", part = "header") %>%
    valign(valign = "top", part = "body")


##  Dredging this model will take a long time.
## I ran it on the ND cluster using 48 cores
## See dredge_tableS5. R fer details.
clu <- makeCluster(4)
clusterExport(clu, c("temp_data_only_known",
                     "who_grooms_only_known_model_prop_consort"))
clusterEvalQ(clu,library(lme4,logical.return =T))
clusterEvalQ(clu,library(MuMIn,logical.return =T))
who_grooms_only_known_model_prop_consort_dredge <-
    dredge(who_grooms_only_known_model_prop_consort,
           cluster = clu)
stopCluster(clu)

## After running download ...
## and open here
S5_table_results <- read_csv("./data/S5_table_results.csv")

who_grooms_best_models <- subset(who_grooms_only_known_model_prop_consort_dredge, delta <= 2)

who_grooms_best_models_tibble <- who_grooms_best_models %>%
    as_tibble() %>%
    select(AICc, delta, everything()) %>%
    select(-logLik, -df, -weight)

who_grooms_best_models_tibble %>%
    pivot_longer(names_to = "variable",
                 values_to = "value", `(Intercept)`:proportion_consort_time) %>%
    filter(!is.na(value)) %>%
    mutate(variable = case_when(
        variable == "cumulative_adversity" ~ "cum_adv",
        variable == "proportion_consort_time" ~ "prop_consort",
        variable == "daily_d_days_rate" ~ "d_days",
        variable == "mean",
        TRUE ~ variable)) %>%
    pivot_wider(names_from = "variable", values_from = "value") %>%
    flextable() %>%
    colformat_double(digits = 3)

who_grooms_best_models_model.avg <- MuMIn::model.avg(who_grooms_best_models, revised.var = TRUE)

summary(who_grooms_best_models_model.avg)$coefmat.subset %>%
    as.data.frame() %>%                # Convert to data frame first
    rownames_to_column(var = "Variable") %>%  # Preserve variable names
    as_tibble() %>%
    flextable() %>%
    colformat_double(digits = 3)
################################################################################
################################################################################
## Fig 2
load("./data/Fig2.RData")
xdata_females_with_social <- read_csv("./data/Table2_data.csv")

text_height = 3.6
leg_x_pos = 0.04
leg_y_pos = .2
legend.spacing.y = 1
legend_spacing = .4
text_size = 12
annodate_fontsize = 6

surv_data_cumpat_long <- surv_data_cumpat %>%
    pivot_longer(cols = -age,
                 names_to = "type",
                 values_to = "predicted.value") %>%
    mutate(type = case_when(
        type == "low_cum_low_pat" ~ "Low DSIpaternal (low ELA)",
        type == "low_cum_high_pat" ~ "High DSIpaternal (low ELA)",
        type == "high_cum_low_pat" ~ "Low DSIpaternal (high ELA)",
        type == "high_cum_high_pat" ~ "High DSIpaternal (high ELA)")) %>%
    mutate(type = forcats::fct_relevel(type,
                                       "High DSIpaternal (low ELA)",
                                       "Low DSIpaternal (low ELA)",
                                       "High DSIpaternal (high ELA)",
                                       "Low DSIpaternal (high ELA)")
    )

surv_cumdpres_data_long <- surv_cumdpres_data %>%
    pivot_longer(cols = -age,
                 names_to = "type",
                 values_to = "predicted.value") %>%
    mutate(type = case_when(
        type == "low_cum_1y" ~ "Short co-residency (low ELA)",
        type == "low_cum_4y" ~ "Long co-residency (low ELA)",
        type == "high_cum_1y" ~ "Short co-residency (high ELA)",
        type == "high_cum_4y" ~ "Long co-residency (high ELA)")) %>%
    mutate(type = forcats::fct_relevel(type,
                                       "Long co-residency (low ELA)",
                                       "Short co-residency (low ELA)",
                                       "Long co-residency (high ELA)",
                                       "Short co-residency (high ELA)")
           )

### Make the actual plot
## This needed a lot of fiddling an I made a function
loadfonts(device = "win")

A <- coxph(data = xdata_females_with_social,
           Surv(statage, adult_survival_status) ~
               cumulative_adversity + jDSI_paternal + dad_overlap_years) %>%
    tidy() %>%
    mutate(term = case_when(term == "cumulative_adversity" ~ "Cumulative adversity",
                            # term == "jDSI_paternal" ~ "Grooming relationship (DSI <sub>paternal</sub>)",
                            term == "jDSI_paternal" ~ "Grooming relationship (DSIpaternal)",
                            term == "dad_overlap_years" ~ "Paternal co-residency")) %>%
    mutate(term = forcats::fct_relevel(term, "Paternal co-residency", after = 0L)) %>%
    mutate(term = forcats::fct_relevel(term, "Cumulative adversity", after = Inf)) %>%
    mutate(low.95 = exp(estimate - 1.96*std.error),
           ,high.95 = exp(estimate + 1.96*std.error)) %>%
    ggplot(aes(y = term, x = exp(estimate))) +
        geom_segment(aes(x = low.95, xend = high.95, yend = term, color = term),
                         size = 5) +
        geom_point(size = 2, color = "black") +
        theme_cowplot(font_size = text_size, font_family = "Times New Roman") +
        labs(x="Hazard ratio") +
        geom_vline(xintercept = 1) +
        scale_color_manual(values = c( "dodgerblue4", "green4", "darkgrey")) +
        annotate("text",x= .83, y=text_height,
                 label="Enhanced survival",
                 family = "Times New Roman",
                 color = "black", size = (text_size-3)/.pt) +
        annotate("text",x= 1.16,y=text_height,
                 label="Reduced survival",
                 family = "Times New Roman",
                 color = "black", size = (text_size-3)/.pt)+
        geom_curve(x = .67, y = text_height - 0.025, xend = .50,
                   yend = text_height - 0.025, curvature = 0,
                   arrow = arrow(length = unit(0.1, "inch")), size = 1,
                color = "black") +
        geom_curve(x = 1.33, y = text_height - 0.025, xend = 1.50,
                   yend = text_height - 0.025,
                    curvature = 0,
                    arrow = arrow(length = unit(0.1, "inch")), size = 1,
                    color = "black")  +
        scale_x_continuous(limits = c(.5, 1.75)) +
        coord_cartesian(ylim=c(1.2,3),clip="off") +
        theme(aspect.ratio = .2,
              legend.position = "none",
              text = element_text(size = text_size)) +

        labs(y="")

A2 <- A +
    #scale_y_discrete(labels = manual_labels) + # Use manual labels
    theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    legend.title = element_text(size =11),
    legend.text = element_text(size =11),
    #axis.text.y = element_markdown(size =11, family = "Times New Roman"),
    plot.margin = unit(c(1, 0, 0.5, 1.5), "cm") # Increased top margin
)

B <- surv_data_cumpat_long %>%
    ggplot() +
        geom_line(aes(x = age,
                      y = predicted.value,
                      colour = type,
                      linetype = type),
                      size = .8) +
        cowplot::theme_cowplot(font_size = text_size, font_family = "Times New Roman") +
        scale_color_manual(
                values = c("green4", "green3", "green4", "green3"),
                labels = c(
                    expression(paste("High DSI"[paternal], " (low ELA)")),
                    expression(paste("Low DSI"[paternal], " (high ELA)")),
                    expression(paste("High DSI"[paternal], " (high ELA)")),
                    expression(paste("Low DSI"[paternal], " (low ELA)"))
                    )) +
        scale_linetype_manual(
                values = c("solid", "solid", "dotted", "dotted"),
                labels = c(
                    expression(paste("High DSI"[paternal], " (low ELA)")),
                    expression(paste("Low DSI"[paternal], " (high ELA)")),
                    expression(paste("High DSI"[paternal], " (high ELA)")),
                    expression(paste("Low DSI"[paternal], " (low ELA)"))
                ))

## cumulative adn dad overlap
C<-  surv_cumdpres_data_long %>%
        ggplot() +
            geom_line(aes(x = age,
                          y = predicted.value,
                          colour = type,
                          linetype = type),
                      size = .8) +
            scale_color_manual(values = c("dodgerblue4", "dodgerblue",
                                          "dodgerblue4", "dodgerblue")) +
            scale_linetype_manual(values = c("solid", "solid", "dotted", "dotted")) +
            cowplot::theme_cowplot(font_size = text_size, font_family = "Times New Roman")

## Most of the custom settings are done in the next part
set_plot_layout <- list(
    scale_x_continuous(expand = c(0, 0), limits = c(0, 27)),
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)),
    cowplot::theme_cowplot(font_size = text_size),
    theme(legend.text = element_text(size = text_size * 0.75),
          legend.position = c(leg_x_pos, leg_y_pos),
          legend.spacing = unit(0.7, "cm"),
          legend.spacing.y = unit(0.7, "cm"),
          legend.key.width = unit(.45, "cm"),
          legend.key.height = unit(legend_spacing, "cm"),
          legend.key.size = unit(1, "cm"),
          legend.margin = margin(t = 5, b = 5, unit = "pt"),
          plot.margin = unit(c(1, 1, .5, 1), "cm"), # Added margin
          text = element_text(family = "Times New Roman", size = text_size),
          axis.text = element_text(size = text_size - 1),
          axis.title.y = element_text(margin = margin(r = 25)),
          plot.title = element_text(hjust = 0.5, size = text_size, face = "bold")
    ),
    labs(x = "Age in years",
         y = "Proportion surviving",
         linetype = "",
         color = "")
    )
B2 <- B  +
    set_plot_layout +
    labs(title = expression(bold("Grooming relationship (DSI"[paternal]*")")))

C2 <- C +
    set_plot_layout +
    labs(title = "Paternal co-residency")

combined <- plot_grid(
    plot_grid(A2 + theme(plot.margin = unit(c(-6, 0, -6, 0), "cm"))
              , NULL
              , labels = c("A.")
              , label_size = text_size
              , label_x = .04
              , label_y = .80
              , hjust = -.5
              , rel_widths = c(1,.1)
              ),
    plot_grid(B2,C2,
              labels = c("B.", "C.")
              , align = "v"
              , label_size = text_size
              , label_x = .1
              , label_y = .95
              , hjust = -.5
              ) +  theme(plot.margin = unit(c(-2, 0, 0, 0), "cm")),
    nrow = 2,
    rel_heights = c(1, .8)) +
    theme(plot.margin = unit(c(-1.5, 0,0, 0), "cm"))

ggsave("./figures/Fig2.jpg", combined, width = 8.5, units = "in", dpi = 600)
