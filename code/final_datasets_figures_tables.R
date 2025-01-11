## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).

## The code below contains an overview of the code to create the final datasets.
## These datasets were used for the final analysises and creating the tables and figures.
## The final datasets are publicly available at https://zenodo.org/records/14590285
## The data in the publicly available has been anonymized.
## The key to get back tp original snames is available up on request.

## libraries
packages <- c(
    "tidyverse"
    , "purrr"
    , "RPostgreSQL"
)

## install packages if needed and open libaries
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
}
lapply(packages, require, character.only = TRUE, warn.conflicts = TRUE,
       quietly = FALSE)

coded_names <- biograph_l %>%
    arrange(birth) %>%
    select(sname, sex) %>%
    filter(!is.na(sname)) %>%
    group_by(sex) %>%
    mutate(ID = paste0(sex, 1:n()))

source('./code/3. functions_for_data_prep.R')
load("./data/data_set_for_dad_analysis_21JUN22.Rdata") ## this can be the older version
load("./data/datasets_for_paper_1MAR24.Rdata")
load("./data/Silk_figure_data_1MAR24.Rdata")

xdata_females_with_social <- xdata_females_with_social %>%
    rename_with(~str_remove(., '.universal')) %>%
    mutate(cumulative_adversity =  if_else(cumulative_adversity > 3, 3, cumulative_adversity))

## function to anonymize the data
replace_names_with_ids <- function(data, coded_names, columns) {
    # Ensure `columns` is a character vector
    if (!is.character(columns)) {
        stop("`columns` must be a character vector of column names.")
    }

    # Create a lookup table for faster matching
    lookup <- setNames(coded_names$ID, coded_names$sname)

    # Replace names with IDs in specified columns
    data %>%
        mutate(across(
            .cols = all_of(columns),
            .fns = ~ lookup[.] %>% replace_na("NA") # Match and replace NA if not found
        ))
}



################################################################################
## Fig 1
## Fig 1A
Fig1A_data <- xdata_females_with_social %>%
    select(focal, dad_overlap_years) %>%
    mutate(focal_order = as.numeric(fct_reorder(focal, dad_overlap_years)))

class_Fig1A_data <- Fig1A_data %>%
    mutate(class = round(dad_overlap_years)) %>%
    group_by(class) %>%
    summarise(cases =n(), .groups = 'drop') %>%
    mutate(class_max = cumsum(cases)) %>%
    mutate(percentage = round(class_max/max(class_max) * 100,0))

# Replace names in the `focal` column
Fig1A_data <- replace_names_with_ids(
    data = Fig1A_data,
    coded_names = coded_names,
    columns = c("focal")
)

write_csv(x = Fig1A_data, file = "./data/data_FigA.csv")

dad_overlap_plot <- ggplot() +
    geom_segment(data = Fig1A_data,
                 aes(x = 0, xend= dad_overlap_years,
                     y=focal_order,
                     yend = focal_order),
                 size = .2) +
    labs(x = "Cumulative duration of co-residency (years)",
         y = "Juvenile female subjects") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_y_continuous(breaks= temp_data_overlap_plot$class_max,
                       labels= paste0(class_Fig1A_data$class_max, " (",
                                      class_Fig1A_data$percentage, '%)'),
                       expand = c(0, 0), limits = c(0, NA)) +
    geom_segment(data = class_Fig1A_data,
                 aes(x=class, xend = class,
                     y = 0, yend = class_max),
                 linetype = "dashed", color = "red", size = 1) +
    geom_segment(data = class_Fig1A_data,
                 aes(x=0, xend = class,
                     y = class_max, yend = class_max),
                 linetype = "dashed", color = "red", size =1) +
    cowplot::theme_cowplot(font_size = 10)

dad_overlap_plot

## Fig 1B
Fig1B_data <- actor_actees_juvAM %>%
    mutate(paternal_grooms =
               if_else(paternal_grooms == "Grooms with dad",
                       "Grooming between juvenile females and their fathers",
                       "Grooming between juvenile females and any adult male")) %>%

    mutate(paternal_grooms =
               fct_relevel(.f = paternal_grooms, "Grooming between juvenile females and their fathers", after = 0L)) %>%
    mutate(dyad_type = case_when(
        dyad_type == "Juvenile groomed AM" ~ "Juv. female grooms any adult male",
        dyad_type == "AM groomed juvenile" ~ "Any adult male grooms juv. female",
        dyad_type == "Juvenile groomed dad" ~ "Daughter grooms father",
        dyad_type == "Dad groomed juvenile" ~ "Father grooms daughter")) %>%
    mutate(dyad_type = forcats::fct_relevel(dyad_type,
                                            "Juv. female grooms any adult male",
                                            "Any adult male grooms juv. female",
                                            "Daughter grooms father",
                                            "Father grooms daughter",
                                            after =0L)) %>%
    filter(!is.na(actor_dad)| !is.na(actee_dad)) %>%
    ungroup() %>%
    mutate(is_focal_juvenile = actor %in% xdata_females_with_social$focal |
               actee %in% xdata_females_with_social$focal) %>%
    mutate(juvenile = if_else(condition = actor_age < 4, actor, actee)) %>%
    select(juvenile, is_focal_juvenile,
           actor, actee, age_class, paternal_grooms,dyad_type) %>%
    group_by(juvenile, age_class, paternal_grooms, dyad_type) %>%
    summarise(cases = n(), .groups = 'drop')

# Replace names in the `focal` column
Fig1B_data <- replace_names_with_ids(
    data = Fig1B_data,
    coded_names = coded_names,
    columns = c("juvenile")
)

write_csv(x = Fig1B_data, file = "./data/data_FigB.csv")

directional_plot <-
    tibble(juvenile = rep(unique(step1$juvenile), each = 16),
           age_class = rep(rep(unique(step1$age_class), each = 4), times = 610),
           paternal_grooms = rep(rep(unique(Fig1B_data$paternal_grooms), each = 2), times = 610*4),
           dyad_type = rep(unique(Fig1B_data$dyad_type), times = 610 * 4)) %>%
    left_join(Fig1B_data) %>%
    mutate(cases = if_else(is.na(cases), 0L, cases)) %>%
    group_by(juvenile, age_class, paternal_grooms) %>%
    mutate(total = sum(cases)) %>%
    filter(total > 0) %>%
    mutate(percentage = cases/sum(cases)) %>%
    group_by(age_class, paternal_grooms, dyad_type) %>%
    summarise(mean.value = mean(percentage, na.rm = TRUE),
              sd.value = sd(percentage, na.rm = TRUE),
              n.value = n(),
              .groups = 'drop') %>%
    mutate(se.value = sd.value/ sqrt(n.value),
           lower.ci.value =
               mean.value - qt(1 - (0.05 / 2), n.value - 1) * se.value,
           upper.ci.value =
               mean.value + qt(1 - (0.05 / 2), n.value - 1) * se.value) %>%
    ggplot(aes(x=age_class,
               y=mean.value,
               group = dyad_type,
               fill = dyad_type,
               color = dyad_type)) +
    geom_point(position=position_dodge(.5), size = 5) +
    geom_errorbar(aes(ymin=lower.ci.value, ymax=upper.ci.value),
                  position = position_dodge(.5))	+
    geom_line(linetype = 'dashed', size = .5, position = position_dodge(.5)) +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(. ~ paternal_grooms, ncol = 1) +
    scale_fill_manual(values = c("#2c7bb6", "#abd9e9",
                                 "#d7191c", "#fdae61"),
                      breaks = c("Father grooms daughter",
                                 "Daughter grooms father",
                                 "Any adult male grooms juv. female",
                                 "Juv. female grooms any adult male")) +
    scale_color_manual(values = c("#2c7bb6", "#abd9e9",
                                  "#d7191c", "#fdae61"),
                       breaks = c("Father grooms daughter",
                                  "Daughter grooms father",
                                  "Any adult male grooms juv. female",
                                  "Juv. female grooms any adult male")) +
    cowplot::theme_cowplot(font_size = 10) +
    guides(fill = 'none',
           color = guide_legend(ncol = 2, size = 5)) +
    labs(x= "Juvenile female age (years)",
         y = "Proportions of grooming initiated") +
    theme(strip.text.x = element_text(size = 8, angle = 0),
          legend.text=element_text(size=rel(0.6)),
          legend.position = "bottom",
          legend.title = element_blank(),
          plot.margin = margin(t = 5, r = -20, b = 5, l = 5))

directional_plot

## To create publication ready figures font sizes may have to be adapted
fig1 <- cowplot::plot_grid(dad_overlap_plot, directional_plot, align = "v", labels = "AUTO", rel_widths = c(.3, .3))
fig1

################################################################################
## Table 1
juvenile_values <- read_csv("./data/juvenile_values.csv") %>%
    mutate(age_class = floor(juvenile_age))

Table1_data <- juvenile_values %>%
    inner_join(biograph_l %>%
                   filter(sex == "F") %>%
                   select(juvenile_id = sname)) %>%
    filter(focal %in% xdata_females_with_social$focal &
               focal == juvenile_id) %>%
    filter(str_detect(paternal_groom, "paternal"))

Table1_data <- replace_names_with_ids(
    data = Table1_data,
    coded_names = coded_names,
    columns = c("focal", "juvenile_id", "sname")
)

write_csv(Table1_data, "./data/Table1_data.csv")

juvenile_model <- Table1_data %>%
    filter(str_detect(paternal_groom, "paternal")) %>%
    lmer(formula = bond_strength ~ juvenile_age + is_dad + (1|juvenile_id)) %>%
    tidy() %>%
    filter(effect == 'fixed')


juvenile_model <- juvenile_model %>%
    select(term, estimate, std.error, statistic, df, p.value)

names(juvenile_model)[2:6] <- c('\U03B2','\U03C3', 't', 'df', 'p')

juv_mod_interpretation = c("",
                           "\U2191 Juvenile age \U2191 bond strength",
                           "\U2191 bond strength = male is father")

juvenile_model %>%
    mutate(term = c(
        "Intercept",
        "Age of juvenile",
        "Is the male the father?" )) %>%
    mutate(interpretation = c("",
                              "\U2191 Juvenile age \U2191 bond strength",
                              "\U2191 bond strength = male is father")) %>%
    mutate(p = format.pval(p, eps = .001, digits = 2)) %>%
    flextable() %>%
    colformat_double(j = 2:5,digits = 3, big.mark = "", decimal.mark = ".") %>%
    colformat_double(j = 6,digits = 3, big.mark = "", decimal.mark = ".") %>%
    fontsize(size = 9, part = "all") %>%
    align(align = "left", part = "header") %>%
    #flextable::rotate(j = 2:5, rotation="btlr",part="header") %>%
    align(align = "left", part = "all") %>%
    fontsize(size = 8, part = "all") %>%
    width(j = 1, width= 1.7) %>%
    width(j = 2:4, width= .5) %>%
    width(j = 5, width= .7) %>%
    width(j = 7, width= 2)

################################################################################
## Table 2
Table2_data <- xdata_females_with_social %>%
    select(focal, statage, adult_survival_status, cumulative_adversity, dad_overlap_years,
           jDSI_paternal, jDSI_M, jDSI_Mde)

Table2_data <- replace_names_with_ids(
    data = Table2_data,
    coded_names = coded_names,
    columns = c("focal")
)
write_csv(Table2_data, "./data/Table2_data.csv")

set_overall <- tibble(formula =
                          c(c("cumulative_adversity"),
                            c("jDSI_paternal"),
                            c("cumulative_adversity + jDSI_paternal"),
                            c("cumulative_adversity + jDSI_Mtop"),
                            c("cumulative_adversity + jDSI_Mde_top"),
                            c("cumulative_adversity + jDSI_paternal + jDSI_Mde_top"),
                            c("cumulative_adversity + dad_overlap_years"),
                            c("cumulative_adversity + dad_overlap_years + jDSI_paternal")
                          )) %>%
    mutate(model = map(.x = formula, .f = get_coxph_model)) %>%
    mutate(confint = map(.x = model, .f =get_confint)) %>%
    mutate(results = map(.x = model, .f = tidy)) %>%
    mutate(model_details = map(.x = model, .f = glance)) %>%
    mutate(AICc = map_dbl(.x = model, .f = MuMIn::AICc)) %>%
    mutate(model_check = map_dbl(.x = model, .f = get_coxzph)) %>%
    unnest(cols = c(results, confint))

## Note for this to work dataset has to be called xdata_females_with_social
## The reduced dataset can be found in Table2_data

set_overall <- tibble(formula =
										c(c("cumulative_adversity"), ##E
											c("jDSI_paternal"),
											c("cumulative_adversity + jDSI_paternal"),  ##B
											c("cumulative_adversity + jDSI_Mtop"),
                    	c("cumulative_adversity + jDSI_Mde_top"), ## F
                      c("cumulative_adversity + jDSI_paternal + jDSI_Mde_top"), ## D
											c("cumulative_adversity + dad_overlap_years"), ##C
											c("cumulative_adversity + dad_overlap_years + jDSI_paternal")  ##A
											)) %>%
  mutate(model = map(.x = formula, .f = get_coxph_model)) %>%
  mutate(confint = map(.x = model, .f =get_confint)) %>%
  mutate(results = map(.x = model, .f = tidy)) %>%
	mutate(model_details = map(.x = model, .f = glance)) %>%
  mutate(AICc = map_dbl(.x = model, .f = MuMIn::AICc)) %>%
	mutate(model_check = map_dbl(.x = model, .f = get_coxzph)) %>%
  unnest(cols = c(results, confint))

## The actual table in the paper was mannually edited.
## The data needed for the table can be extracted usinh the following function.
make_flextable(set_overall)

################################################################################
#Fig 2.

## There is quiet a bit of prep needed to get to generate Fig 2.
# First run two survival models
coxph_cumpat_fit <- coxph(Surv(statage, adult_survival_status) ~ cumulative_adversity + jDSI_paternal,  data = Table2_data)
coxph_cumdpres_fit <- coxph(Surv(statage, adult_survival_status) ~ cumulative_adversity + dad_overlap_years,  data = xdata_females_with_social)

## get the  low, mean and high values of jDSI_paternal
## These are the mean = 1SD, the mean, the mean + 1SD for jDSO_paternal
## For overlap I set it at 1 and 4 years
## I am looking at cummulative adversity of 1 and 3

Low_Pat = mean(xdata_females_with_social$jDSI_paternal, na.rm = TRUE) -
    sd(xdata_females_with_social$jDSI_paternal, na.rm =  TRUE)
Mean_Pat = mean(xdata_females_with_social$jDSI_paternal, na.rm =  TRUE)
High_Pat = mean(xdata_females_with_social$jDSI_paternal, na.rm =  TRUE) +
    sd(xdata_females_with_social$jDSI_paternal, na.rm = TRUE)

## use these values to create survival curves.

fit_cumpat_surv <- survfit(coxph_cumpat_fit
                           , newdata=data.frame(cumulative_adversity = c(1, 1, 3, 3),
                                                jDSI_paternal = c(Low_Pat, High_Pat, Low_Pat, High_Pat)))

fit_cumdpres_surv <- survfit(coxph_cumdpres_fit
                             , newdata=data.frame(cumulative_adversity = c(0, 0, 3, 3),
                                                  dad_overlap_years = c(1,4,1,4)))
## get some summary data
## the order is the order in the new_data (so check if needded)

## I extracted the 0.5 survival point for the low and high jDSI pater
medium_low_Pat_high_cum = as.numeric(summary(fit_cumpat_surv)$table[,'median'][3])
medium_low_Pat_low_cum = as.numeric(summary(fit_cumpat_surv)$table[,'median'][1])
medium_high_Pat_high_cum = as.numeric(summary(fit_cumpat_surv)$table[,'median'][4])
medium_high_Pat_low_cum = as.numeric(summary(fit_cumpat_surv)$table[,'median'][2])

medium_low_pres_high_cum = as.numeric(summary(fit_cumdpres_surv)$table[,'median'][3])
medium_low_pres_low_cum = as.numeric(summary(fit_cumdpres_surv)$table[,'median'][4])
medium_high_pres_high_cum = as.numeric(summary(fit_cumdpres_surv)$table[,'median'][1])
medium_high_pres_low_cum = as.numeric(summary(fit_cumdpres_surv)$table[,'median'][2])

cumpat_survival_difference_lowadversity <- medium_high_Pat_low_cum - medium_low_Pat_low_cum
cumpat_survival_difference_highadversity <- medium_high_Pat_high_cum - medium_low_Pat_high_cum
cumpat_survival_difference <- medium_high_Pat_low_cum - medium_low_Pat_high_cum

## put it all together
surv_data_cumpat <- data.frame(age = fit_cumpat_surv$time,
                               low_cum_low_pat = fit_cumpat_surv$surv[,1],
                               low_cum_high_pat = fit_cumpat_surv$surv[,2],
                               high_cum_low_pat = fit_cumpat_surv$surv[,3],
                               high_cum_high_pat = fit_cumpat_surv$surv[,4])

surv_cumdpres_data <- data.frame(
    age = fit_cumdpres_surv$time,
    low_cum_1y = fit_cumdpres_surv$surv[,1],
    low_cum_4y = fit_cumdpres_surv$surv[,2],
    high_cum_1y = fit_cumdpres_surv$surv[,3],
    high_cum_4y = fit_cumdpres_surv$surv[,4])

## Then make long format for ggplot

surv_data_cumpat_long <- surv_data_cumpat %>%
    pivot_longer(cols = -age,
                 names_to = "type",
                 values_to = "predicted.value") %>%
    mutate(type = case_when(
        type == "low_cum_low_pat" ~ "Low cum. adversity with weak paternal bond",
        type == "low_cum_high_pat" ~ "Low cum. adversity with strong paternal bond",
        type == "high_cum_low_pat" ~ "High cum. adversity with weak paternal bond",
        type == "high_cum_high_pat" ~ "High cum. adversity with strong paternal bond"))

surv_cumdpres_data_long <- surv_cumdpres_data %>%
    pivot_longer(cols = -age,
                 names_to = "type",
                 values_to = "predicted.value") %>%
    mutate(type = case_when(
        type == "low_cum_1y" ~ "Low cum. adversity \n with 1 year",
        type == "low_cum_4y" ~ "Low cum. adversity \n with 4 year",
        type == "high_cum_1y" ~ "High cum. adversity \n 1 year",
        type == "high_cum_4y" ~ "High cum. adversity \n 4 year"))


## I extracted the 0.5 survival point for the low and high jDSI pater
medium_low_Res = as.numeric(summary(fit_cumres_surv)$table[,'median'][3])
medium_high_Res = as.numeric(summary(fit_cumres_surv)$table[,'median'][2])

survival_difference_Res = medium_high_Res - medium_low_Res
cumres_survival_difference_Res = medium_high_Res - medium_low_Res

cumres_survival_difference_lowadversity <- summary(fit_cumres_surv)$table[,'median'][2] -  summary(fit_cumres_surv)$table[,'median'][1]

cumres_survival_difference_highadversity <- summary(fit_cumres_surv)$table[,'median'][4] -  summary(fit_cumres_surv)$table[,'median'][3]

text_height = 3.6
A <- coxph(data = Fig2_data, Surv(statage, adult_survival_status) ~ cumulative_adversity + jDSI_paternal + dad_overlap_years) %>%
    tidy() %>%
    mutate(term = case_when(term == "cumulative_adversity" ~ "cumulative adversity",
                            term == "jDSI_paternal" ~ "paternal dyadic bond strength",
                            term == "dad_overlap_years" ~ "co-residency with father")) %>%
    mutate(term = forcats::fct_relevel(term, "cumulative adversity", after = Inf)) %>%
    mutate(low.95 = exp(estimate - 1.96*std.error),
           high.95 = exp(estimate + 1.96*std.error),
           low.99 = exp(estimate - 2.575*std.error),
           high.99 = exp(estimate + 2.575*std.error)) %>%
    ggplot(aes(y = term, x = exp(estimate))) +
    geom_segment(aes(x = low.99, xend = high.99, yend = term),
                 linewidth =  5, color = "dodgerblue", alpha = .7) +
    geom_segment(aes(x = low.95, xend = high.95, yend = term),
                 linewidth = 5, color = "dodgerblue4") +
    geom_point(size = 2,, color = "black") +
    theme_cowplot() +
    # scale_x_continuous(trans='log10',
    # 									 breaks = c(.5, .7, 1, 1,3, 1.5, 2),
    # 									 limits = c(.45, 2)
    labs(x="Hazard ratio") +
    geom_vline(xintercept = 1) +
    annotate("text",x= .8,y=text_height,label="Enhanced survival", color = "gray")+
    annotate("text",x= 1.2,y=text_height,label="Reduces survival", color = "gray")+
    geom_curve(x = .63, y = text_height - 0.025, xend = .55, yend = text_height - 0.025, curvature = 0,
               arrow = arrow(length = unit(0.08, "inch")), size = 1,
               color = "gray") +
    geom_curve(x = 1.36, y = text_height - 0.025, xend = 1.44, yend = text_height - 0.025, curvature = 0,
               arrow = arrow(length = unit(0.08, "inch")), size = 1,
               color = "gray")  +
    #scale_y_discrete(expand = expand_scale(mult = c(0.5, .5))) +
    coord_cartesian(ylim=c(1.2,3),clip="off") +
    theme(aspect.ratio = .2) +
    labs(y="")

## survival plotea and paternal
B <- ggplot() +
    geom_line(data = surv_data_cumpat_long, aes(x = age,
                                                y = predicted.value,
                                                colour = type,
                                                linetype = type),
              size = 1.2) +
    geom_segment(aes(x = 0, xend = medium_high_Pat_low_cum,
                     y  = 0.5, yend = 0.5),
                 color = 'black', linetype = 'dashed') +


    geom_segment(aes(x = medium_high_Pat_low_cum, xend = medium_high_Pat_low_cum,
                     y  = 0, yend = 0.5),
                 color = 'black', linetype = 'dashed') +
    geom_segment(aes(x = medium_low_Pat_low_cum, xend = medium_low_Pat_low_cum,
                     y  = 0, yend = 0.5),
                 color = 'black', linetype = 'dashed') +
    geom_segment(aes(x = medium_high_Pat_low_cum - 0.3, xend = medium_low_Pat_low_cum + 0.3,
                     y  = 0.07, yend = 0.07),
                 size = .2, arrow = arrow(length = unit(0.05, "inches"))) +
    geom_segment(aes(x = medium_low_Pat_low_cum + 0.3, xend = medium_high_Pat_low_cum - 0.3,
                     y  = 0.07, yend = 0.07),
                 size = .2, arrow = arrow(length = unit(0.05, "inches"))) +
    annotate("text",x = medium_low_Pat_low_cum + (medium_high_Pat_low_cum - medium_low_Pat_low_cum)/2, y = 0.03, size = 2,
             label = paste0(round(medium_high_Pat_low_cum - medium_low_Pat_low_cum,2), " y")) +



    geom_segment(aes(x = 0, xend = medium_high_Pat_high_cum,
                     y  = 0.5, yend = 0.5),
                 color = 'black', linetype = 'dotdash') +
    geom_segment(aes(x = medium_low_Pat_high_cum, xend = medium_low_Pat_high_cum,
                     y  = 0, yend = 0.5),
                 color = 'black', linetype = 'dotted') +
    geom_segment(aes(x = medium_high_Pat_high_cum, xend = medium_high_Pat_high_cum,
                     y  = 0, yend = 0.5),
                 color = 'black', linetype = 'dotted') +
    geom_segment(aes(x = medium_low_Pat_high_cum + 0.3, xend = medium_high_Pat_high_cum -0.3,
                     y  = 0.06, yend = 0.06),
                 size = .2, arrow = arrow(length = unit(0.05, "inches"))) +
    geom_segment(aes(x = medium_high_Pat_high_cum - 0.3, xend = medium_low_Pat_high_cum + 0.3,
                     y  = 0.06, yend = 0.06),
                 size = .2, arrow = arrow(length = unit(0.05, "inches"))) +
    annotate("text",x = medium_low_Pat_high_cum + (medium_high_Pat_high_cum - medium_low_Pat_high_cum)/2
             , y = 0.03, size = 2, size = 5,
             label = paste0(round(medium_low_Pat_low_cum - medium_low_Pat_high_cum,2), " y")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 31)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(legend.position = "none",
          legend.text=element_text(size=7),
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    scale_color_manual(values = c("firebrick", "firebrick", "dodgerblue4", "dodgerblue4")) +
    scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted")) +
    labs(x = "Age",
         y = "Survival change",
         linetype = "",
         color = "")  +
    guides(colour = guide_legend(nrow = 2))


C <- ggplot() +
    geom_line(data = surv_cumdpres_data_long, aes(x = age,
                                                y = predicted.value,
                                                colour = type,
                                                linetype = type),
              size = 1.2) +
    geom_segment(aes(x = medium_low_pres_low_cum, xend = medium_high_pres_low_cum,
                     y  = 0.5, yend = 0.5),
                 color = 'black', linetype = 'dashed') +
    geom_segment(aes(x = medium_high_pres_high_cum, xend = medium_high_pres_high_cum,
                     y  = 0, yend = 0.5),
                 color = 'black', linetype = 'dashed') +
    geom_segment(aes(x = medium_high_pres_low_cum-.1, xend = medium_high_pres_low_cum-.1,
                     y  = 0, yend = 0.5),
                 color = 'black', linetype = 'dashed') +
    geom_segment(aes(x = 0, xend = medium_low_pres_low_cum,
                     y  = 0.5, yend = 0.5),
                 color = 'black', linetype = 'dotdash') +
    geom_segment(aes(x = medium_low_pres_low_cum, xend = medium_low_pres_low_cum,
                     y  = 0, yend = 0.5),
                 color = 'black', linetype = 'dotted') +
    geom_segment(aes(x = medium_low_pres_high_cum, xend = medium_low_pres_high_cum,
                     y  = 0, yend = 0.5),
                 color = 'black', linetype = 'dotted') +
    geom_segment(aes(x = medium_low_pres_high_cum + 0.3, xend = medium_low_pres_low_cum -0.3,
                     y  = 0.06, yend = 0.06),
                 size = .2, arrow = arrow(length = unit(0.05, "inches"))) +
    geom_segment(aes(x = medium_low_pres_low_cum - 0.3, xend = medium_low_pres_high_cum + 0.3,
                     y  = 0.06, yend = 0.06),
                 size = .2, arrow = arrow(length = unit(0.05, "inches"))) +
    annotate("text",x = medium_low_pres_high_cum + (medium_low_pres_low_cum - medium_low_pres_high_cum)/2,
             y = 0.03, size = 2,
             label = paste0(round(medium_low_pres_low_cum - medium_low_pres_high_cum,2), " y")) +

    geom_segment(aes(x = medium_high_pres_high_cum + 0.3, xend = medium_high_pres_low_cum -0.3,
                     y  = 0.07, yend = 0.07),
                 size = .2, arrow = arrow(length = unit(0.05, "inches"))) +
    geom_segment(aes(x = medium_high_pres_low_cum - 0.3, xend = medium_high_pres_high_cum + 0.3,
                     y  = 0.07, yend = 0.07),
                 size = .2, arrow = arrow(length = unit(0.05, "inches"))) +
    annotate("text",x = medium_high_pres_high_cum + (medium_high_pres_low_cum -
                                                         medium_high_pres_high_cum)/2, y = 0.03, size = 2,
             label = paste0(round(medium_high_pres_low_cum - medium_high_pres_high_cum,2), " y")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 31)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    cowplot::theme_cowplot(font_size = 8) +
    #theme(legend.position = "bottom") +
    theme(#legend.position = c(.6, .92),
        #legend.position = "none",
        legend.text=element_text(size=7),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    scale_color_manual(values = c("firebrick", "firebrick", "dodgerblue4", "dodgerblue4")) +
    scale_linetype_manual(values = c("solid", "dotted", "solid", "dotted")) +
    labs(x = "Age",
         y = "Survival change",
         linetype = "",
         color = "")  +
    guides(colour = guide_legend(nrow = 2))

plot_grid(plot_grid(A, labels = c("A.")
                    , label_size = 8
                    , label_x = 0
                    , label_y = 1
                    , hjust = -.5),
          plot_grid(B, C, rel_heights = c(1, 2),
                    labels = c("B.", "C.")
                    , align = "v"
                    , label_size = 8
                    , label_x = 0
                    , label_y = 1
                    , hjust = -.5
          ), nrow = 2)

################################################################################
## Table 3

who_grooms_bootstrap <- read_csv("./data/who_grooms_bootstrap_13SEP24.csv")
Table3_data <- who_grooms_bootstrap %>%
    filter(is_dad & has_consort_data) %>%
    select(focal
           , kid_age, mom_age, AMales_age
           , mean_ordrank
           , daily_d_days_rate
           , proportion_consort_time
           , nr_pdads
           , next_kid_with_mom
           , previous_kid_with_mom
           , offspring_years
           , cumulative_adversity
           , observer_effort
           , does_groom)

Table3_data <- replace_names_with_ids(
    data = Table3_data,
    coded_names = coded_names,
    columns = c("focal")
)

write_csv(Table3_data, "./data/Table3_data.csv")

who_grooms_full_model_only_dad_interpretation = c(
    c(""),
    c("\U2193 group size \U2191 grooming"),
    c("\U2191 juvenile age  \U2191 grooming"),
    #	c("father = yes  \U2191 grooming"),
    c("\U2191 d-days \U2193 grooming"),
    c("\U2193 male rank \U2191 grooming"),
    c("\U2191 paternal sibling \U2191 grooming"),
    c("check what this would mean"),
    #c("\U2191 males age  \U2191 grooming"),
    c("\U2191 mounts \U2191 grooming"),
    #	c("prior kid = yes \U2191 grooming interactions"),
    c(""),
    c(""),
    c(""),
    c(""),
    c(""),
    c("")
)

make_model_flextable(model = "who_grooms_full_model_only_dad",
                     dataset = who_grooms_bootstrap,
                     explain_text = who_grooms_full_model_only_dad_interpretation)

Table4_data <- who_grooms_bootstrap %>%
    filter(is_dad & has_consort_data) %>%
    inner_join(select(xdata_females_with_social, focal, dad_overlap_years))  %>%
    select(focal
           , mom_age, AMales_age
           , mean_ordrank
           , daily_d_days_rate
           , proportion_consort_time
           , nr_pdads
           , previous_kid_with_mom
           , offspring_years
           , cumulative_adversity
           , dad_overlap_years)

Table4_data <- replace_names_with_ids(
    data = Table4_data,
    coded_names = coded_names,
    columns = c("focal")
)

write_csv(Table4_data, "./data/Table4_data.csv")