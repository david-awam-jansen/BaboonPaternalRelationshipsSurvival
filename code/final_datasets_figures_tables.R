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


## Note some tables and figures were changed after review
## Please see additional_analysis.R for details.

## libraries
packages <- c(
      "broom"
    , "broom.mixed"
    , "cowplot"
    , "extrafont"
    , "lmerTest"
    , "MuMIn"
    , "purrr"
    , "showtext"
    , "survival"
    , "systemfonts"
    , "tidyverse")

## install packages if needed and open libaries
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
}

lapply(packages, require, character.only = TRUE, warn.conflicts = TRUE,
       quietly = FALSE)

## briefly set to main directory
setwd("N:/Other computers/My Laptop (before crash)/DSI_paternal_female_survival")
load("./data/data_set_for_dad_analysis_21JUN22.Rdata") ## this can be the older version
load("./data/datasets_for_paper_1MAR24.Rdata")
load("./data/Silk_figure_data_1MAR24.Rdata")
setwd("N:/My Drive/BaboonPaternalRelationshipsSurvival")

## to import fonts
# This may take a while
#font_import()
#loadfonts()    # Load fonts into the current session

## will be used to rename files
coded_names <- biograph_l %>%
    arrange(birth) %>%
    select(sname, sex) %>%
    filter(!is.na(sname)) %>%
    group_by(sex) %>%
    mutate(ID = paste0(sex, 1:n()))

source('./code/3. functions_for_data_prep.R')
font_add(family = "Times New Roman", regular = "times.ttf")

# Define the formatting and saving function
format_and_save_plot <- function(plot, file_name, width = 6, height = 4, dpi = 300, formats = c("png", "tiff", "eps")) {
    # Load required libraries
    showtext_opts(dpi = dpi)
    showtext_auto(enable = TRUE)


    # Add formatting to the plot while preserving existing theme elements
    formatted_plot <- plot +
        theme(
            text = element_text(family = "Times New Roman"),
            plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 9),
            axis.text = element_text(size = 9),
            plot.tag = element_text(size = 9, family = "Times New Roman", face = "italic"),
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 9)
        )

    # Save the plot in specified formats
    for (format in formats) {
        ggsave(
            filename = paste0(file_name, ".", format),
            plot = formatted_plot,
            width = width,
            height = height,
            dpi = dpi,
            device = format
        )
    }

    return(formatted_plot)
}

xdata_females_with_social <- xdata_females_with_social %>%
    rename_with(~str_remove(., '.universal')) %>%
    mutate(cumulative_adversity =  if_else(cumulative_adversity > 3, 3, cumulative_adversity))

xdata_females_with_social %>%
    select(sname = focal, name, pid, birth, matgrp, mom, dad) %>%
    write_csv("table_of_focals_Beth_26MAR25.csv")

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

# Replace names in the `focal` column
Fig1A_data <- replace_names_with_ids(
    data = Fig1A_data,
    coded_names = coded_names,
    columns = c("focal")
)

write_csv(x = Fig1A_data, file = "./data/data_FigA.csv")
Fig1A_data <- read_csv("./data/data_Fig1A.csv")

classes_Fig1A_data <- Fig1A_data %>%
    mutate(class = ceiling(dad_overlap_years)) %>%
    group_by(class) %>%
    summarise(cases =n(), .groups = 'drop') %>%
    mutate(class_max = cumsum(cases)) %>%
    mutate(percentage = round(class_max/max(class_max) * 100,0))

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

write_csv(x = Fig1B_data, file = "./data/data_Fig1B.csv")
Fig1B_data <- read_csv(file = "./data/data_Fig1B.csv") %>%
    mutate(paternal_grooms = forcats::fct_relevel(paternal_grooms,
            "Grooming between juvenile females and their fathers",
            "Grooming between juvenile females and any adult male"))

Fig1B_plot_data <-
    tibble(juvenile = rep(unique(Fig1B_data$juvenile), each = 16),
                          age_class = rep(rep(unique(Fig1B_data$age_class),
                                              each = 4), times = 610),
                          paternal_grooms = rep(rep(unique(Fig1B_data$paternal_grooms), each = 2),
                                                times = 610*4),
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
    mutate(dyad_type = case_when(
        dyad_type == "Any adult male grooms juv. female" ~ "Adult male grooms juvenile female",
        dyad_type == "Juv. female grooms any adult male" ~ "Juvenile female grooms adult male",
        dyad_type == "Daughter grooms father" ~ "Juvenile female  grooms father",
        dyad_type == "Father grooms daughter" ~ "Father grooms juvenile female")) %>%
    mutate(dyad_type = forcats::fct_relevel(dyad_type,
                                  "Father grooms juvenile female",
                                  "Juvenile female  grooms father",
                                  "Adult male grooms juvenile female",
                                  "Juvenile female grooms adult male"))

Fig1A <- ggplot() +
    geom_segment(data = Fig1A_data,
                 aes(x = 0, xend= dad_overlap_years,
                     y=focal_order,
                     yend = focal_order),
                 size = .2) +
    labs(x = "Cumulative co-residency with father (years)",
         y = "Juvenile female subjects") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_y_continuous(breaks= classes_Fig1A_data$class_max,
                       labels= paste0(classes_Fig1A_data$class_max, " (",
                                      classes_Fig1A_data$percentage, '%)'),
                       expand = c(0, 0), limits = c(0, NA)) +
    geom_segment(data = classes_Fig1A_data,
                 aes(x=class, xend = class,
                     y = 0, yend = class_max),
                 linetype = "dashed", color = "red", linewidth = 1) +
    geom_segment(data = classes_Fig1A_data,
                 aes(x=0, xend = class,
                     y = class_max, yend = class_max),
                 linetype = "dashed", color = "red", linewidth = 1) +
    cowplot::theme_cowplot(font_family = "Times New Roman") +
    theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        legend.title = element_text(size =11),
        legend.text = element_text(size =11)
    )

ggsave(plot = Fig1A, filename = "./figures/Fig1A.jpg", width = 8.5, dpi = 600)

Fig1B <- Fig1B_plot_data %>%
    ggplot(aes(x=age_class,
               y=mean.value,
               group = dyad_type,
               fill = dyad_type,
               color = dyad_type)) +
    geom_point(position=position_dodge(.5), size = 2) +
    geom_errorbar(aes(ymin=lower.ci.value, ymax=upper.ci.value),
                  position = position_dodge(.5), width = .4)	+
    geom_line(linetype = 'dashed', size = .5, position = position_dodge(.5)) +
    scale_x_discrete(expand = expansion(add = c(0.2, 0.2))) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0.017)) +

    facet_wrap(. ~ paternal_grooms, ncol = 1) +
    scale_fill_manual(values = c("#2c7bb6", "#abd9e9",
                                 "#d7191c", "#fdae61"),
                      breaks = c("Father grooms juvenile female",
                                 "Juvenile female  grooms father",
                                 "Adult male grooms juvenile female",
                                 "Juvenile female grooms adult male")) +
    scale_color_manual(values = c("#2c7bb6", "#abd9e9",
                                  "#d7191c", "#fdae61"),
                       breaks = c("Father grooms juvenile female",
                                  "Juvenile female  grooms father",
                                  "Adult male grooms juvenile female",
                                  "Juvenile female grooms adult male")) +
    cowplot::theme_cowplot(font_family = "Times New Roman") +
    guides(fill = 'none',
           color = guide_legend(ncol = 1, size = 5)) +
    labs(x= "Juvenile female age (years)",
         y = "Proportions of grooming initiated",
         color = "") +
    theme(strip.text.x = element_text(size = 8, angle = 0),
          legend.text=element_text(size=rel(0.6)),
          legend.position = "bottom",
          legend.title = element_blank()) +
    theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        legend.title = element_text(size =11),
        legend.text = element_text(size =11))


ggsave(plot = Fig1B, filename = "./figures/Fig1B.jpg", width = 8.5, dpi = 600)

## To create publication ready figures font sizes may have to be adapted
Fig1 <- cowplot::plot_grid(Fig1A, NULL, Fig1B,
                           align = "hv",
                           labels = c("A", "B"),
                           label_x = 0,
                           label_y = 1,
                           hjust = -.5,
                           ncol = 3,
                           rel_widths = c(1.5, .1, 1.5))

wrapped_Fig1A <- wrap_elements(full = Fig1A)
wrapped_Fig1B <- wrap_elements(full = Fig1B)

# Combine side by side
Fig1 <- wrapped_Fig1A + wrapped_Fig1B +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(family = "Times New Roman"))

# Now save the combined figure
# For some reason the format_and_save+plot does not work for this combined figure
ggsave(plot = Fig1, filename = "./figures/Fig1.jpg", width = 8.5, dpi = 600)

################################################################################
## Table 1
## After review the data for this table changed
## Please see additional  analysis for details.

juvenile_model_values <- read_csv("./data/data_for_juvenile_model.csv") %>%
    inner_join(biograph_l %>%
                   filter(sex == "F") %>%
                   select(juvenile_id = sname)) %>%
    filter(focal %in% xdata_females_with_social$focal &
               focal == juvenile_id) %>%
    filter(str_detect(paternal_groom, "paternal"))

Table1_data <- replace_names_with_ids(
    data = juvenile_model_values,
    coded_names = coded_names,
    columns = c("focal", "juvenile_id", "sname")
)

write_csv(Table1_data, "./data/Table1_data.csv")
Table1_data <- read_csv("./data/Table1_data.csv")

juvenile_model <- Table1_data %>%
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

juvenile_model_results %>%
    mutate(term = c(
        "Intercept",
        "Age of juvenile",
        "Is the male the father?" )) %>%
    mutate(interpretation = c("",
                              "\U2191 Juvenile age \U2191 bond strength",
                              "\U2191 bond strength = male is father")) %>%
    mutate(p = format.pval(p, eps = .001, digits = 1)) %>%
    flextable::flextable() %>%
    flextable::colformat_double(j = 2:5,digits = 3, big.mark = "", decimal.mark = ".") %>%
    flextable::colformat_double(j = 6,digits = 3, big.mark = "", decimal.mark = ".") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::align(align = "left", part = "header") %>%
    #flextable::rotate(j = 2:5, rotation="btlr",part="header") %>%
    flextable::align(align = "left", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::width(j = 1, width= 1.7) %>%
    flextable::width(j = 2:4, width= .5) %>%
    flextable::width(j = 5, width= .7) %>%
    flextable::width(j = 7, width= 2)


## getting some AIC values
AIC(juvenile_model
    , update(juvenile_model, . ~ . - is_dad)
    , update(juvenile_model, . ~ . - juvenile_age)
    ) %>%
    as_tibble(rownames = "Model") %>%
    mutate(dAIC = AIC - AIC(juvenile_model))

################################################################################
## Table 2
## Some of the data in Table2_data will also be used for Fig2
Table2_data <- xdata_females_with_social %>%
    select(focal, statage, adult_survival_status, cumulative_adversity, dad_overlap_years,
           jDSI_paternal, jDSI_M, jDSI_Mde, jDSI_Mtop, jDSI_Mde_top)

Table2_data <- replace_names_with_ids(
    data = Table2_data,
    coded_names = coded_names,
    columns = c("focal")
)
write_csv(Table2_data, "./data/Table2_data.csv")
Table2_data <- read_csv("./data/Table2_data.csv")

## Some of the code depends on the name xdata_females_with_social
xdata_females_with_social <- Table2_data

# In the next block of code survival models are conducted for a range of models
## the Formulat indicates which variables were included
## See ./code/3. functions_for_data_prep.R for the functions
set_overall <- tibble(formula =	c(c("cumulative_adversity"), ##A
                                  c("jDSI_paternal"), ## B
                                  c("cumulative_adversity + jDSI_paternal"),  ##C
                                  c("cumulative_adversity + jDSI_Mtop"), ## D
                                  c("cumulative_adversity + jDSI_Mde_top"), ## E
                                  c("cumulative_adversity + jDSI_paternal + jDSI_Mde_top"), ## F
                                  c("cumulative_adversity + dad_overlap_years"), ##G
								  c("cumulative_adversity + dad_overlap_years + jDSI_paternal") ## H
								  )) %>%
    mutate(model = map(.x = formula, .f = get_coxph_model)) %>% ## run the survival models
    mutate(confint = map(.x = model, .f =get_confint)) %>% ##  get confidence intervals
    mutate(results = map(.x = model, .f = tidy)) %>% ## get tidy model resuls
    mutate(model_details = map(.x = model, .f = glance)) %>% ## get some model variables
    mutate(AICc = map_dbl(.x = model, .f = MuMIn::AICc)) %>% ## get addapted AICc values
	mutate(model_check = map_dbl(.x = model, .f = get_coxzph)) %>% ## model check
    unnest(cols = c(results, confint)) ## extract all results

## During analysis the following formula made the tables.
## Significant variables were highlighted and models were sorted accoding to AIC
make_flextable(set_overall)
## The actual table in the paper was manually edited and the order was changed
################################################################################
#Fig 2.
Table2_data <- read_csv("./data/Table2_data.csv")
xdata_females_with_social <- Table2_data

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

fit_cumpres_surv <- survfit(coxph_cumdpres_fit
                             , newdata=data.frame(cumulative_adversity = c(0, 0, 3, 3),
                                                  dad_overlap_years = c(1,4,1,4)))
## get some summary data
## the order is the order in the new_data (so check if needded)

## I extracted the 0.5 survival point for the low and high jDSI pater
medium_low_Pat_high_cum = as.numeric(summary(fit_cumpat_surv)$table[,'median'][3])
medium_low_Pat_low_cum = as.numeric(summary(fit_cumpat_surv)$table[,'median'][1])
medium_high_Pat_high_cum = as.numeric(summary(fit_cumpat_surv)$table[,'median'][4])
medium_high_Pat_low_cum = as.numeric(summary(fit_cumpat_surv)$table[,'median'][2])

medium_low_pres_high_cum = as.numeric(summary(fit_cumpres_surv)$table[,'median'][3])
medium_low_pres_low_cum = as.numeric(summary(fit_cumpres_surv)$table[,'median'][4])
medium_high_pres_high_cum = as.numeric(summary(fit_cumpres_surv)$table[,'median'][1])
medium_high_pres_low_cum = as.numeric(summary(fit_cumpres_surv)$table[,'median'][2])

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
    age = fit_cumpres_surv$time,
    low_cum_1y = fit_cumpres_surv$surv[,1],
    low_cum_4y = fit_cumpres_surv$surv[,2],
    high_cum_1y = fit_cumpres_surv$surv[,3],
    high_cum_4y = fit_cumpres_surv$surv[,4])

save(surv_data_cumpat, surv_cumdpres_data, file = "./data/Fig2.RData")

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
medium_low_Res = as.numeric(summary(fit_cumpres_surv)$table[,'median'][3])
medium_high_Res = as.numeric(summary(fit_cumpres_surv)$table[,'median'][2])

survival_difference_Res = medium_high_Res - medium_low_Res
cumres_survival_difference_Res = medium_high_Res - medium_low_Res

cumres_survival_difference_lowadversity <- summary(fit_cumpres_surv)$table[,'median'][2] -  summary(fit_cumpres_surv)$table[,'median'][1]

cumres_survival_difference_highadversity <- summary(fit_cumpres_surv)$table[,'median'][4] -  summary(fit_cumpres_surv)$table[,'median'][3]

### Make the actual plot
## This needed a lot of fiddling an I made a function

## Note this figure was changed after review.
## See additional analysis for details
loadfonts(device = "win")
make_plots <- function(text_height, leg_x_pos, leg_y_pos,
                       legend_spacing, legend.spacing.y, text_size,
                       annodate_fontsize, y_axis_margin) {

    A <- coxph(data = xdata_females_with_social,
               Surv(statage, adult_survival_status) ~
                   cumulative_adversity + jDSI_paternal + dad_overlap_years) %>%
        tidy() %>%
        mutate(term = case_when(term == "cumulative_adversity" ~ "cumulative adversity",
                                term == "jDSI_paternal" ~ "paternal dyadic bond strength",
                                term == "dad_overlap_years" ~ "co-residency with father")) %>%
        mutate(term = forcats::fct_relevel(term, "cumulative adversity", after = Inf)) %>%
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
        annotate("text",x= .85, y=text_height,
                 label="Enhanced survival", color = "black",
                 size = (text_size-2)/.pt) +
        annotate("text",x= 1.15,y=text_height,
                 label="Reduced survival",
                 color = "black", size = (text_size-2)/.pt)+
        geom_curve(x = .70, y = text_height - 0.025, xend = .50, yend = text_height - 0.025, curvature = 0,
                   arrow = arrow(length = unit(0.1, "inch")), size = 1,
                   color = "black") +
        geom_curve(x = 1.3, y = text_height - 0.025, xend = 1.50, yend = text_height - 0.025, curvature = 0,
                   arrow = arrow(length = unit(0.1, "inch")), size = 1,
                   color = "black")  +
        scale_x_continuous(limits = c(.5, 1.75)) +
        coord_cartesian(ylim=c(1.2,3),clip="off") +
        theme(aspect.ratio = .2,
              legend.position = "none",
              text = element_text(size = text_size)) +
        labs(y="")

    B <- surv_data_cumpat_long %>%
        mutate(cumulative_level = if_else(str_detect(tolower(type), "low"),
                                          "Low cumulative adversity (1)",
                                          "High cumulative adversity (3)"),
               paternal_level = if_else(str_detect(tolower(type), "weak"),
                                        "Weak paternal bond (bottom 25%)",
                                        "Strong paternal bond (top 25%)")) %>%
        ggplot() +
        geom_line(aes(x = age,
                      y = predicted.value,
                      colour = paternal_level,
                      linetype = cumulative_level),
                  size = .8) +
        scale_color_manual(values = c("green4", "green2")) +
        scale_linetype_manual(values = c("dotted", "solid")) +
        cowplot::theme_cowplot(font_size = text_size, font_family = "Times New Roman")


    ## cumulative adn dad overlap
    C<-  surv_cumdpres_data_long %>%
        mutate(cumulative_level = if_else(str_detect(tolower(type), "low"),
                                          "Low cumulative adversity (1)",
                                          "High cumulative adversity (3)"),
               paternal_level = if_else(str_detect(tolower(type), "4"),
                                        "Long paternal co-residency (4y)",
                                        "Short paternal co-residency (1y)")) %>%
        ggplot() +
        geom_line(aes(x = age,
                      y = predicted.value,
                      colour = paternal_level,
                      linetype = cumulative_level),
                  size = .8) +
        scale_color_manual(values = c("dodgerblue4", "dodgerblue")) +
        scale_linetype_manual(values = c("dotted", "solid")) +
        cowplot::theme_cowplot(font_size = text_size, font_family = "Times New Roman")

    ## Most of the custom settings are done in the next part
    set_plot_layout <- list(
        scale_x_continuous(expand = c(0, 0), limits = c(0, 27)),
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1)),
        cowplot::theme_cowplot(font_size = text_size),
        theme(legend.text=element_text(size=text_size*0.75),
              legend.position = c(leg_x_pos, leg_y_pos),
              legend.spacing = unit(10, "cm"),  # Adjust spacing between items
              #legend.spacing.y = unit(-.3, "cm"), ## distance between
              legend.spacing.y = unit(1, "cm"), ## distance between
              legend.key.width =  unit(1, "cm"),
              legend.key.height =  unit(legend_spacing, "cm"),
              legend.key.size = unit(0.8, "cm"),  # Consistent key size
              legend.margin = margin(t = -5, b = -5, unit = "pt"),
              plot.margin = unit(c(0, 0, 0, 0,0), "cm"),
              text = element_text(family = "Times New Roman", size = text_size),
              axis.text = element_text(size = text_size-1),
              plot.title = element_text(hjust = 0.5, size = text_size,  face = "bold"),
              axis.title.y = element_text(margin = y_axis_margin)),
        labs(x = "Age in years",
             y = "Proportion surviving",
             linetype = "",
             color = "")#,
        #guides(colour = guide_legend(nrow = 2),
         #      linetype = guide_legend(override.aes = list(size = 1, n.dots = 2)))
    )

        A2 <- A +  theme(
        text = element_text(family = "Times New Roman"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        legend.title = element_text(size =11),
        legend.text = element_text(size =11),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
        )

    B2 <- B  + set_plot_layout +
        labs(title = expression("Grooming relationship (DSI"[paternal]*")"))

    C2 <- C + set_plot_layout +
        labs(title = "Paternal co-residency")

    combined <- plot_grid(plot_grid(A2, NULL, labels = c("A.")
                                    , label_size = text_size
                                    , label_x = 0
                                    , label_y = .8
                                    , hjust = -.5
                                    , rel_widths = c(1,.1)),
                          plot_grid(B2, C2,
                                    labels = c("B.", "C.")
                                    , align = "v"
                                    , label_size = text_size
                                    , label_x = 0
                                    , label_y = 1.05
                                    , hjust = -.5
                          ), nrow = 2, rel_heights = c(1, 3))

    ggsave("./figures/Fig2.jpg", combined, width = 8.5, height = 6, units = "in", dpi = 1200)



}

text_height = 3.6
leg_x_pos = 0.01
leg_y_pos = .21
legend.spacing.y = 1
legend_spacing = .4
text_size = 12
annodate_fontsize = 6
y_axis_margin <- margin(0, 0,0, 0)

make_plots(text_height, leg_x_pos, leg_y_pos, legend_spacing, legend.spacing.y, text_size, annodate_fontsize,  y_axis_margin)
## Plot looks ok, but needs some layout fiddling. It might need patchwork instead of cowplot.

################################################################################
## Table 3

who_grooms_bootstrap <- read_csv("./data/who_grooms_bootstrap_13SEP24.csv")

Table3_data <- who_grooms_bootstrap %>%
    filter(is_dad & has_consort_data) %>%
    select(focal
           , AMales
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
           , does_groom
           , nr_days
           )

Table3_data <- replace_names_with_ids(
    data = Table3_data,
    coded_names = coded_names,
    columns = c("focal", "AMales")
)

write_csv(Table3_data, "./data/Table3_data.csv")
Table3_data <- read_csv("./data/Table3_data.csv")

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

## The table for this model is made by a function.
## The function works with model results from the bootstrapping
## It does model averaging ect. You have to give it the name of model,
## the dataset and an vector with interpretations.


make_model_flextable(model = "who_grooms_full_model_only_dad",
                     dataset = who_grooms_bootstrap,
                     explain_text = who_grooms_full_model_only_dad_interpretation)

################################################################################
## Table 4
dad_overlap_only_know_data <- read_csv("./data/dad_overlap_only_know_data.csv")

Table4_data <- dad_overlap_only_know_data %>%
    select(focal
           , mom_age, dad_age,
           , male_rank
           , daily_d_days_rate
           , proportion_consort_time
           , nr_pdads
           , previous_kid_with_mom
           , offspring_years
           , cumulative_adversity
           , dad_overlap_years
           , dad)


Table4_data <- replace_names_with_ids(
    data = Table4_data,
    coded_names = coded_names,
    columns = c("focal", "dad")
)

write_csv(Table4_data, "./data/Table4_data.csv")
Table4_data <- read_csv("./data/Table4_data.csv")

dad_overlap_only_known_model_prop_consort_formula <- dad_overlap_years ~ 1 +
    mom_age+ dad_age +
    male_rank +
    daily_d_days_rate +
    proportion_consort_time +
    nr_pdads +
    previous_kid_with_mom +
    offspring_years +
    cumulative_adversity +
    (1|dad)

dad_overlap_only_known_model_prop_consort <- lmer(dad_overlap_only_known_model_prop_consort_formula
                                                  , data = Table4_data
                                                  , na.action = 'na.fail')

dad_overlap_only_known_model_prop_consort_results <- dad_overlap_only_known_model_prop_consort %>%
    broom.mixed::tidy() %>%
    left_join(vif.mer(dad_overlap_only_known_model_prop_consort), by = 'term') %>%
    mutate(run = paste0("model_", 0)) %>%
    mutate(row = NA) %>%
    select(row, everything()) %>%
    mutate(term = if_else(str_detect(term, "rank"), "dad_rank", term))

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