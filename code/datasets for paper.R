## old fiel version
## library
library(lubridate)
library(tidyverse)

source('./code/3. functions_for_data_prep.R')

## this the last time the babase data was downloaded

## from babase
#biograph_l <- read_csv("./data/biograph_l.csv")
#ranks_l <- read_csv("./data/ranks_l.csv")
#parents_l <- read_csv("./data/parents_l.csv")
#actor_actees_l <- read_rds("./data/actor_actee_l.RDS")
## the next file contains data from babase this can be the older
## THis mainly consist of babase data
load("./data/data_set_for_dad_analysis_21JUN22.Rdata") ## this can be the older version
members_l <- read_csv("./data/members_l.csv") %>%
	arrange(sname, date) %>%
	select(-sys_period)

members_AM_l <- members_l %>%  ## Only ranked males
	#select(-sys_period) %>%
	inner_join(rankdates_l) %>%
	filter(date > ranked) %>%
	select(sname, grp, date, ranked)

members_juveniles_l <- members_l %>%
	inner_join(select(biograph_l, sname, birth), by = 'sname') %>%
	mutate(age = as.numeric(date- birth)/365.25) %>%
	filter(age < 4) %>%
	inner_join(select(parents_l, sname = kid, dad))

## preparing the juveinle early adversity dataset
xdata_ea_raw <- read_csv(file = paste("./data/", "ea_dataset_less_restricted_",
																			"21JUN22", ".csv", sep="")) %>%
	mutate(end_juvenile = birth + years(4) - days(1)) %>%
	mutate(dad_overlap = pmap(.l = list(sname, dad, birth, end_juvenile),
														.f = get_overlap)) %>%
	unnest(cols = c(dad_overlap)) %>%
	rename_with(stringr::str_replace,
							pattern = "partner", replacement = "dad")

xdata_ea <- xdata_ea_raw %>%
	inner_join(select(parents_l, sname = kid, mom, dad),
						 by = c("sname", "mom", "dad")) %>%
	mutate(statage = as.numeric((statdate - birth) / 365.25)) %>%
	mutate(included_cases = complete.cases(
		sname, birth, matgrp, bstatus, statage,
		mom, dad,
		maternal_loss, maternal_SCI_F,
		maternal_rank, competing_sibling, density, drought)) %>%
	mutate(survival_status_at_age_4 = case_when(
		ymd(birth) + years(4) < ymd(statdate) ~ 0,
		ymd(birth) + years(4) > ymd(statdate) &
			status %in% c(0, 1) ~ 1,
		ymd(birth) + years(4) > ymd(statdate) &
			status %in% c(2, 3) ~ 0),
		survival_status_at_age_4_desc = case_when(
			ymd(birth) + years(4) < ymd(statdate) ~ 0,
			ymd(birth) + years(4) > ymd(statdate) &
				status %in% c(0, 1) ~ 1,
			ymd(birth) + years(4) > ymd(statdate) &
				status %in% c(2, 3) ~ 999),
		survival_status_at_age_4_desc = factor(
			survival_status_at_age_4_desc,
			labels = c("Alive", "Died", "Censored")
		),
		age_at_age_4 = ifelse(statage <= 4, statage, 4),
		adult_survival_status = ifelse(status == 1, 1, 0),
		## create binary variables using Tung eta al cutoffs
		maternal_rank_binary = maternal_rank >= 12,
		maternal_SCI_binary = maternal_SCI_F <= as.numeric(
			quantile(xdata_ea_raw$maternal_SCI_F, na.rm = TRUE, probs = 0.25)),
		density_binary = density > 35,
		sex = factor(sex, labels = c("Females", "Males"))) %>%
	group_by(sname) %>%
	mutate(cumulative_adversity = sum(maternal_loss + maternal_rank_binary +
																			maternal_SCI_binary + competing_sibling +  density_binary + drought)) %>%
	ungroup()

xdata_ea_named <- xdata_ea %>%
	mutate(sex= factor(sex, labels = c("Females", "Males")),
				 maternal_loss = factor(maternal_loss,
				 											 levels = c(FALSE, TRUE),
				 											 labels = c("Mom alive", "Mom dies")),
				 maternal_SCI_F_binary = factor(maternal_SCI_F_binary,
				 															 levels = c(FALSE, TRUE),
				 															 labels = c("Socially connected", "Socially isolated")),
				 maternal_rank_binary = factor(maternal_rank_binary,
				 															levels = c(FALSE, TRUE),
				 															labels = c("Normal rank", "Low rank")),
				 competing_sibling = factor(competing_sibling,
				 													 levels = c(FALSE, TRUE),
				 													 labels = c("Absent", "Com sib. present")),
				 drought = factor(drought,
				 								 levels = c(FALSE, TRUE),
				 								 labels = c("Wet", "Dry")),
				 density_binary = factor(density_binary,
				 												levels = c(FALSE, TRUE),
				 												labels = c("Normal group", "Large group")),
				 first_born = factor(first_born,
				 										levels = c(FALSE, TRUE),
				 										labels = c("Experienced mother", "First born")))

## social indexes
dsi_data_raw <- read_rds(paste0("./data/", "dsi_data_",
													 "1MAR24", ".rds")) %>%
	select(focal = sname, focal_grp = grp, everything())
# This file is double nested. DSI values are inside zscored_values
# and specific dyad values for the focal are inside focal data

## same as above butnow universal zacored values
dsi_universal <- read_rds(paste("./data/", "dsi_universal_",
															 "1MAR24", ".rds", sep="")) %>%
	select(focal, focal_grp = grp, everything()) %>%
	select(-sname)
## Like above this file has a nested data. However this file also has the dyad data for all the other possible dyads.

dsi_data <- dsi_data_raw%>%
	unnest(cols = c(zscored_values))

## same as above butnow universal zacored values
dsi_universal <- read_rds(paste("./data/", "dsi_universal_",
																"1MAR24", ".rds", sep="")) %>%
	select(focal, focal_grp = grp, everything()) %>%
	select(-sname)

universal_zscored_values <- read_rds("./data/universal_zscored_values_1MAR24.RDS")

paternal_data <- universal_zscored_values%>%
	select(universal_zscored_values) %>%
	unnest(cols = c(universal_zscored_values)) %>%
	select(focal, row_nr = my_row, paternal_position, AF_partners,
				 AM_partners, nr_AM, nr_AM_partners, mom, dad, mom_present, mom_groomed,
				 dad_present, dad_groomed, dad_groom_class = groom_type)

## I forgot to add Mtop and Mde_top to the updated DSI calculation part on
## cluster. I manually make them here
## Mtop is the top grooming partner strength
## Mde_top is the top male grooming partner strength that is not the dad

DSI_M_top_values <- dsi_universal %>%
	select(focal, focal_grp, age_class, row_nr, universal_zscored_resvalues) %>%
	filter(age_class < 4) %>%
	unnest(cols = c(universal_zscored_resvalues)) %>%
	filter(focal_check == TRUE) %>%
	group_by(focal, focal_grp, age_class, row_nr) %>%
	arrange(-universal_zscore, paternal_groom) %>%
	filter(str_detect(paternal_groom, "paternal")) %>%
	slice(1) %>%
	ungroup() %>%
	select(row_nr, DSI_Mtop = universal_res_value)

DSI_Mde_top_values <- dsi_universal %>%
	select(focal, focal_grp, age_class,  row_nr, universal_zscored_resvalues) %>%
	filter(age_class < 4) %>%
	unnest(cols = c(universal_zscored_resvalues)) %>%
	filter(focal_check == TRUE) %>%
	group_by(focal, focal_grp, age_class, row_nr) %>%
	arrange(-universal_zscore, paternal_groom) %>%
	filter(str_detect(paternal_groom, "non_paternal")) %>%
	slice(1) %>%
	ungroup() %>%
	select(row_nr, DSI_Mde_top = universal_res_value)

dsi_data_combined <- dsi_data %>%
	left_join(select(dsi_universal, row_nr, contains("DSI"))) %>%
	left_join(select(paternal_data, everything())) %>%
	left_join(DSI_M_top_values) %>%
	left_join(DSI_Mde_top_values)

write_csv(dsi_data_combined, "./data/dsi_data_combined.csv")
write_csv(paternal_data, "./data/paternal_data.csv")


social_indexes_males <- dsi_data_combined %>%
	mutate(age = as.numeric((start - birth)/365.25)) %>%
	filter(age < 4) %>%
	select(focal, focal_grp, start, end, age, contains("DSI"), everything())

write_csv(social_indexes_males, "./data/social_indexes_males.csv")

male_values <- dsi_data_combined %>%
	select(focal, focal_grp, start, end, focal_data) %>%
	unnest(cols = c(focal_data)) %>%
	filter(focal_check == TRUE & str_detect(paternal_groom, "paternal")) %>%
	rename(AMales = sname)

juvenile_values <- dsi_universal %>%
	unnest(cols = universal_zscored_resvalues) %>%
	select(row_nr, focal, focal_start = start, juvenile_id = partner, sname,
				 paternal_groom, birth,  bond_strength = universal_zscore, i_adj) %>%
	mutate(juvenile_age = as.numeric((focal_start - birth)/365.25),
				 is_dad = if_else(paternal_groom == "paternal",
				 								 1 ,0))

juvenile_social_indexes_males = social_indexes_males  %>%
	pivot_longer(names_to = "index", values_to = "value", cols = starts_with("DSI")) %>%
	group_by(age_class, index) %>%
	mutate(value.scaled = get_zscore(value)) %>%
	group_by(focal, index) %>%
	summarise(jmean = mean(value.scaled, na.rm = TRUE),
						nr_years = sum(!is.na(value.scaled)),
						.groups = 'drop') %>%
	filter(nr_years >= 2) %>%
	select(-nr_years) %>%
	ungroup() %>%
	pivot_wider(names_from = index, values_from = jmean, names_prefix = "j")

# dad_top_groomer <- dsi_data_combined %>%
# 	inner_join(select(xdata_females_with_social	, focal)) %>%
# 	#inner_join(select(biograph_l, focal, birth)) %>%
# 	mutate(age = round(as.numeric((start - birth)/365.25), 0)) %>%
# 	select(row_nr, focal, focal_grp, start, end, age, age_class, focal_data) %>%
# 	mutate(top_dyads = map2(.x = focal_data, .y="paternal_res_i_adj", .f= get_top_grooms),
# 				 has_real_partners = map_lgl(.x = top_dyads, .f = has_real_partners),
# 				 nr_dyads = map_dbl(.x = focal_data, .f = check_nr_dyads)) %>%
# 	left_join(select(dsi_universal, row_nr, dad_present)) %>%
# 	mutate(is_dad_top = map2_lgl(.x = top_dyads, .y="paternal_res_i_adj_zscored",
# 															 .f = get_dad_top)) %>%
# 	select(-contains("DSI")) %>%
# 	filter(focal %in% xdata_females_with_social$focal) %>%
# 	mutate(dad_status =
# 				 	case_when(dad_present == FALSE ~ "Dad not present",
# 				 						has_real_partners == FALSE ~ "No adult male grooming partners",
# 				 						is_dad_top == TRUE ~  "Dad is the top grooming partner",
# 				 						is_dad_top == FALSE ~  "Dad is the not top partner")) %>%
# 	mutate(dad_status =
# 				 	forcats::fct_relevel(dad_status,
# 				 											 "Dad is the top grooming partner", after =0L))
#
# write_csv(dad_top_groomer, "./data/dad_top_groomer.csv")

juvenile_social_indexes  <- dsi_data_combined %>%
	filter(age_class < 4) %>%
	pivot_longer(names_to = "index", values_to = "value", contains('DSI')) %>%
	group_by(age_class, index) %>%
	mutate(value.scaled = get_zscore(value)) %>%
	group_by(focal, index) %>%
	summarise(jmean = mean(value.scaled, na.rm = TRUE),
						nr_years = sum(!is.na(value.scaled)),
						.groups = 'drop') %>%
	group_by(focal) %>%
	mutate(nr_years = max(nr_years)) %>%
	filter(nr_years >= 2) %>%
	select(-nr_years) %>%
	ungroup() %>%
	pivot_wider(names_from = index, values_from = jmean, names_prefix = "j")

grooms_by_dads <- dsi_universal %>%
	select(focal, focal_dad = dad, focal_grp,
				 focal_start = start, birth, dad_groomed,
				 , universal_zscored_resvalues) %>%
	unnest(universal_zscored_resvalues) %>%
	filter(sname == focal_dad |
				 	focal_check == TRUE & paternal_groom == "paternal") %>%
	filter(dad_groomed) %>%
	select(focal, focal_dad, focal_start, birth, juvenile_id = partner, paternal_groom,
				 bond_strength = universal_zscore, dad_groomed)

write_csv(grooms_by_dads, file = "./data/grooms_by_dads.csv")

## preparing xdata with social
xdata_females_with_social <- xdata_ea %>%
	rename(focal = sname) %>%
	filter(sex == 'Females')  %>%
	filter(age_at_age_4 >= 4) %>%
	filter(!is.na(dad) & !is.na(mom)) %>%
	inner_join(juvenile_social_indexes, by = "focal")

xdata_females_with_social_named <- xdata_ea_named %>%
	rename(focal = sname) %>%
	filter(sex == 'Females')  %>%
	filter(age_at_age_4 >= 4) %>%
	filter(!is.na(dad) & !is.na(mom)) %>%
	inner_join(juvenile_social_indexes, by = "focal")

dad_data <- xdata_females_with_social %>%
	select(focal, birth, dad, dad_overlap_years) %>%
	mutate(dad_data = pmap(.l = list(focal, dad, birth), .f = get_some_dad_data)) %>%
	unnest(cols = c(dad_data), keep_empty = TRUE)  %>%
	mutate(father_rankdata = map2(.x = dad, .y = rnkdate, .f = get_next_dad_rank)) %>%
	unnest(cols = c(father_rankdata), keep_empty = TRUE) %>%
	mutate(end_juv = birth + lubridate::years(4) - lubridate::days(1)) %>%
	mutate(age_of_juvenile = as.numeric(last_date - birth)/365.25,
				 days_to_next_grp = as.numeric(next_rnkdate - last_date),
				 days_to_statdate = as.numeric(dad_statdate - last_date)) %>%
	mutate(reason = case_when(
		end_juv == last_date ~ "reached end of juvenile period",
		!is.na(next_ranked_grp) ~ "moved to other grp",
		is.na(next_ranked_grp) & dad_status == 1 ~ "father died",
		is.na(next_ranked_grp) & dad_status != 0 ~ paste0("unknow/censored"),
		TRUE ~ "gone before birth"))


dad_groom_status <- dsi_data_combined %>%
	filter(focal %in% xdata_females_with_social$focal) %>%
	select(row_nr, age_class, , focal_data, dad_present, dad_groomed) %>%  unnest(cols = focal_data) %>%
	filter(focal_check) %>%
	filter(str_detect(paternal_groom, "paternal")) %>%
	group_by(row_nr) %>%
	#arrange(-{{ variable }}) %>%
	arrange(row_nr, -res_i_adj) %>%
	select(row_nr, age_class, , sname, paternal_groom, i_total, res_i_adj, dad_present, dad_groomed) %>%
	filter(sname != 'XXX') %>% ## double check if this is needed
	mutate(position = cume_dist(res_i_adj)) %>%
	group_by(row_nr, age_class, dad_present, dad_groomed) %>%
	summarise(has_grooming_partners = sum(i_total > 0.00274, na.rm = TRUE) > 0,
						is_dad_top=sum(paternal_groom == "paternal" & position == 1) != 0) %>%
	ungroup() %>%
	mutate(dad_status =
			 	case_when(dad_present == FALSE ~ "Dad not present",
			 						has_grooming_partners == FALSE ~ "No adult male grooming partners",
			 						is_dad_top == FALSE ~  "Dad is the not top partner",
			 						has_grooming_partners == TRUE &is_dad_top == TRUE ~  "Dad is the top grooming partner",
			 						TRUE ~ "check")) %>%
	mutate(dad_status =
				 	forcats::fct_relevel(dad_status,
				 											 "Dad is the top grooming partner", after =0L))

dad_groom_status_step1 <- juvenile_values %>%
	filter(focal %in% xdata_females_with_social$focal) %>%
	filter(juvenile_id == focal) %>%
	filter(str_detect(paternal_groom, "paternal")) %>%
	distinct() %>%
	inner_join(dsi_data %>%
						 	filter(focal %in% xdata_females_with_social$focal) %>%
						 	select(row_nr, focal, focal_grp, focal_data) %>%
						 	unnest(cols = c(focal_data)) %>%
						 	filter(focal_check & str_detect(paternal_groom, "paternal")) %>%
						 	select(row_nr, focal, focal_grp, sname, i_total, paternal_groom)) %>%
	arrange(row_nr, -bond_strength, paternal_groom) %>%
	group_by(row_nr, focal, focal_grp) %>%
	mutate(position = cume_dist(bond_strength))

dad_groom_status_step2 <- dad_groom_status_step1 %>%
	select(row_nr, focal, focal_start, juvenile_age, sname, paternal_groom, i_total, position) %>%
	filter(paternal_groom == "paternal") %>%
	mutate(real_groomer = i_total >= 1 & !is.na(i_total))

dad_groom_status_step3 <- dad_groom_status_step1 %>%
	inner_join(dad_groom_status_step2) %>%
  filter(paternal_groom == "paternal") %>%
	mutate(real_groomer = i_total >= 1 & !is.na(i_total),
				 dad_present = sname != 'XXX') %>%
	ungroup() %>%
	mutate(is_dad_top = position == 1) %>%
	select(row_nr, focal, focal_start, juvenile_age, position, real_groomer, dad_present, is_dad_top)

dad_groom_status <- dad_groom_status_step3 %>%
	inner_join(dad_groom_status_step3 %>%
						 	group_by(row_nr) %>%
						 	summarise(has_grooming_partners = sum(real_groomer))) %>%
	mutate(dad_status =
				 	case_when(dad_present == FALSE ~ "Dad not present",
				 						has_grooming_partners == 0 ~ "No adult male grooming partners",
				 						is_dad_top == FALSE ~  "Dad is the not top partner",
				 						is_dad_top == TRUE ~  "Dad is the top grooming partner",
				 						TRUE ~ "check")) %>%
	mutate(dad_status =
				 	forcats::fct_relevel(dad_status,
				 											 "Dad is the top grooming partner", after =0L)) %>%
	mutate(age = round(juvenile_age))

dad_data <- xdata_females_with_social %>%
    select(focal, birth, dad, dad_overlap_years) %>%
    mutate(dad_data = pmap(.l = list(focal, dad, birth), .f = get_some_dad_data)) %>%
    unnest(cols = c(dad_data), keep_empty = TRUE)  %>%
    mutate(father_rankdata = map2(.x = dad, .y = rnkdate, .f = get_next_dad_rank)) %>%
    unnest(cols = c(father_rankdata), keep_empty = TRUE) %>%
    mutate(end_juv = birth + lubridate::years(4) - lubridate::days(1)) %>%
    mutate(age_of_juvenile = as.numeric(last_date - birth)/365.25,
           days_to_next_grp = as.numeric(next_rnkdate - last_date),
           days_to_statdate = as.numeric(dad_statdate - last_date)) %>%
    mutate(reason = case_when(
        end_juv == last_date ~ "reached end of juvenile period",
        !is.na(next_ranked_grp) ~ "moved to other grp",
        is.na(next_ranked_grp) & dad_status == 1 ~ "father died",
        is.na(next_ranked_grp) & dad_status != 0 ~ paste0("unknow/censored"),
        TRUE ~ "gone before birth"))

write_rds(dad_data, "./data/dad_data.rds")

temp3 <- xdata_females_with_social %>%
    arrange(dad_overlap_years) %>%
    mutate(ordered_sname = seq(1:nrow(xdata_females_with_social))) %>%
    inner_join(select(dad_data, focal, reason))

silk_figure_data <- Silk_fig4_step1  %>%
    filter(focal %in% xdata_females_with_social$focal) %>%  ## do we want to restrict?
    inner_join(Silk_grooming_step2) %>%
    select(focal, focal_grp, age_class, nr_days, contains('mean'), paternal, non_paternal) %>%
    pivot_longer(names_to = 'index', values_to = 'value', mean_nr_adult_males:non_paternal) %>%
    group_by(age_class, index) %>%
    summarise(mean.value = mean(value, na.rm = TRUE),
              sd.value = sd(value, na.rm = TRUE),
              n.value = n(),
              .groups = 'drop') %>%
    mutate(se.value = sd.value/ sqrt(n.value),
           lower.ci.value =
               mean.value - qt(1 - (0.05 / 2), n.value - 1) * se.value,
           upper.ci.value =
               mean.value + qt(1 - (0.05 / 2), n.value - 1) * se.value)  %>%
    mutate(type = if_else(condition = str_detect(index, "paternal"),
                          true = "Number of grooming dyads",
                          false = "Number of males present")) %>%
    mutate(index = case_when(
        index == "mean_dad_presence" ~ "Probability of father being around",
        index == "mean_nr_adult_males" ~ "Number of adult males present in group",
        index == "mean_nr_adult_males_not_father" ~
            "Number of adult males present in group (excluding father)",
        index == "paternal" ~ "Probability of a grooming dyad with father",
        index == "non_paternal" ~
            "Number of grooming dyads with other adult males")) %>%
    mutate(type = forcats::fct_relevel(type, "Number of males present", after = 0L)) %>%
    mutate(age_class = case_when(age_class == 1 ~ "0-1",
                                 age_class == 2 ~ "1-2",
                                 age_class == 3 ~ "2-3",
                                 age_class == 4 ~ "3-4"
    ))






save(social_indexes_males, male_values, dad_groom_status,
		 juvenile_values, juvenile_social_indexes_males,
		 juvenile_social_indexes, grooms_by_dads,
		 xdata_females_with_social, xdata_females_with_social_named,
		 dad_data,
		 file = "./data/datasets_for_paper_1MAR24.Rdata")

## preparting some data based onm Silk et al papers
Silk_grooming_step1 <- dsi_data %>%
	mutate(grooming = map(.x = focal_data, .f = get_grooming_partners))

Silk_grooming_step2 <- Silk_grooming_step1 %>%
	unnest(cols = c(grooming), keep_empty = TRUE) %>%
	select(focal, age_class, focal_grp, start, end, paternal, non_paternal) %>%
	mutate_at(vars(paternal:non_paternal), ~replace(., is.na(.), 0))

Silk_fig4_step1 <- dsi_data %>%
	inner_join(select(parents_l, focal = kid, dad)) %>%
	select(focal, focal_grp, start, end, dad, ranked) %>%
	mutate(AM_count = pmap(.l = list(focal, focal_grp, start, end, dad), .f = get_AM_count)) %>%
	unnest(cols = c(AM_count))

Silk_fig4_step2 <- Silk_fig4_step1 %>%
	inner_join(Silk_grooming_step2) %>%
	select(focal, focal_grp, age_class, nr_days, contains('mean'), paternal, non_paternal) %>%
	pivot_longer(names_to = 'index', values_to = 'value',
							 mean_nr_adult_males:non_paternal) %>%
	group_by(focal, age_class, index) %>%
	summarise(value = weighted.mean(x = value, w = nr_days, na.rm = TRUE),
						.groups='drop')

Silk_fig4_step2_wide <- Silk_fig4_step2 %>%
	pivot_wider(names_from = index, values_from = value) %>%
	mutate(AM_grooming = paternal + non_paternal)

Silk_fig4_step3 <- Silk_fig4_step2 %>%
	group_by(age_class, index) %>%
	summarise(value = mean(value, na.rm = FALSE)
						, .groups='drop') %>%
	pivot_wider(names_from = 'index', values_from = 'value')

save(Silk_fig4_step1,
		 #Silk_fig4_step1_wide,
		 Silk_fig4_step2,
		 Silk_fig4_step2_wide,
		 Silk_fig4_step3,
		 #Silk_fig4_step3b,
		 #Silk_fig4_step4,
		 Silk_grooming_step1,
		 Silk_grooming_step2,
		 file = paste0("./data/Silk_figure_data_", "1MAR24", ".Rdata"))

actor_actees_l <- actor_actees_l %>%
	filter(act == 'G ') %>%
	inner_join(select(parents_l, actor = kid, actor_dad = dad)) %>%
	inner_join(select(parents_l, actee = kid, actee_dad = dad)) %>%
	inner_join(select(biograph_l,
										actor = sname, actor_birth = birth, actor_sex = sex)) %>%
	inner_join(select(biograph_l,
										actee = sname, actee_birth = birth, actee_sex = sex)) %>%
	left_join(select(rankdates_l, actor = sname, actor_ranked = ranked)) %>%
	left_join(select(rankdates_l, actee = sname, actee_ranked = ranked)) %>%
	mutate(actor_age = (date - actor_birth)/365.25,
         actee_age = (date - actee_birth)/365.25) %>%
	mutate(actor_is_ranked = date >= actor_ranked,
				 actee_is_ranked = date >= actee_ranked)

actor_actees_juvAM <- actor_actees_l %>%
	filter((actor_is_ranked & actee_age < 4 ) |
				 	(actee_is_ranked & actor_age < 4 )) %>%
	mutate(dyad_type = case_when(actor == actee_dad ~ "Dad groomed juvenile",
															 actee == actor_dad ~ "Juvenile groomed dad",
															 actor_age < 4 & actee_is_ranked ~ "Juvenile groomed AM",
															 actee_age < 4 & actor_is_ranked ~ "AM groomed juvenile",
															 TRUE ~ "check")) %>%
	mutate(paternal_grooms = case_when(
		str_detect(dyad_type, regex('dad', ignore_case = T)) ~ 'Grooms with dad',
		str_detect(dyad_type, "AM") ~ 'Grooms with adult males')) %>%
	group_by(iid) %>%
	mutate(juvenile_age = min(actor_age, actee_age)) %>%
	mutate(age_class = case_when(
		juvenile_age <= 1 ~ "0-1",
		juvenile_age <= 2 ~ "1-2",
		juvenile_age <= 3 ~ "2-3",
		juvenile_age <= 4 ~ "3-4"))

write_csv(actor_actees_juvAM, './data/actor_actees_juvAM.csv')

observer_effort_values <-dsi_data %>%
	ungroup() %>%
	mutate(observer_effort = map_dbl(.x = focal_data,
																		.f = get_mean_AF_log2OE)) %>%
	#mutate(observer_effort = map_dbl(.x = focal_data,
	#.f = get_mean_observer_effort)) %>% ## very similiar as above
	select(focal, focal_grp, start, end, observer_effort)

## get bootstrap results
results_files_path = './data/bootstrap_results/raw_files/'
result_files <- list.files(path = results_files_path, pattern = "best_model")

for(ii in 1:length(result_files)) {
	xdata <- read_delim(file = paste0(results_files_path, result_files[ii]),
											delim = " ", col_names = FALSE, show_col_types = FALSE)

	if(ncol(xdata) == 10) names(xdata) <- c("row", "effect", "group",  "term", "estimate" ,"std.error", "statistic", "p.value", "VIF", "model_run")
	if(ncol(xdata) == 11) names(xdata) <- c("row", "effect", "group",  "term", "estimate" ,"std.error", "statistic", "df", "p.value", "VIF", "model_run")

	xdata <- xdata %>%
		filter(effect == 'fixed') %>%
		select(-effect, -group, -row, -model_run) %>%
		mutate(term = case_when(str_detect(term , "bootstrap_previous") ~ "previous_kid_with_mom",
														str_detect(term , "bootstrap_next") ~ "next_kid_with_mom",
														TRUE ~ term))
	new_filename <- str_replace(string = result_files[ii], pattern = "_raw.text",
															replacement = ".csv")
	new_summary <- str_replace(string = result_files[ii], pattern = "_raw.text",
														 replacement = "_summary.csv")
	write_csv(xdata, paste0("./data/bootstrap_results/clean_files/", new_filename))
	print(paste0(new_filename, " has been saved"))

	xdata %>%
		group_by(term) %>%
		summarise_all(.funs = mean) %>%
		mutate(sig = p.value < 0.05) %>%
		write_csv(paste0("./data/bootstrap_results/", new_summary))
}






#
#

#
# dad_top_groomer <- dsi_data_combined %>%
# 	#inner_join(select(biograph_l, sname, birth)) %>%
# 	#mutate(age = as.numeric((start - birth)/365.25)) %>%
# 	#rename(focal = sname) %>%
# 	#unnest(cols = c(zscored_values)) %>%
# 	mutate(top_dyads = map2(.x = focal_data, .y="res_i_adj_zscored", .f= get_top_grooms),
# 				 has_real_partners = map_lgl(.x = top_dyads, .f = has_real_partners),
# 				 nr_dyads = map_dbl(.x = top_dyads, .f = check_nr_dyads)) %>%
# 	#left_join(select(dsi_universal, row_nr, dad_present)) %>%
# 	mutate(is_dad_top = map2_lgl(.x = top_dyads, .y="res_i_adj_zscored", .f = get_dad_top))

#

# xdata_females_with_social <- read_csv(paste0('./data/xdata_females_with_social_', latest_version_date, '.csv'))
#
# dsi_data <- read_rds(paste0("./data/dsi_data_", latest_version_date,".rdsxxx"))
# dsi_universal_data <- read_rds(paste0("./data/dsi_universal_", latest_version_date,".rds"))

# dsi_data_combined <- dsi_data %>%
# unnest(cols = c(zscored_values)) %>%
# 	left_join(dsi_universal_data) %>%
# 	select(-sname) %>%
# 	select(focal, focal_grp = grp, everything())
#
# social_indexes_males <- dsi_data_combined %>%
# 		mutate(age = as.numeric((start - birth)/365.25)) %>%
# 		filter(age < 4) %>%
#   	select(focal, focal_grp, start, end, age, contains("DSI"), everything())
#
# male_values <- dsi_data_combined %>%
# 	select(focal, focal_grp, start, focal_data) %>%
# 	unnest(cols = c(focal_data)) %>%
# 	filter(focal_check == TRUE & str_detect(paternal_groom, "paternal")) %>%
# 	rename(AMales = sname)
#
# juvenile_social_indexes_males = social_indexes_males  %>%
#   pivot_longer(names_to = "index", values_to = "value", cols = starts_with("DSI")) %>%
# 	group_by(age_class, index) %>%
# 	mutate(value.scaled = zscore(value)) %>%
# 	group_by(focal, index) %>%
# 	summarise(jmean = mean(value.scaled, na.rm = TRUE),
# 						nr_years = sum(!is.na(value.scaled)),
# 						.groups = 'drop') %>%
# 	filter(nr_years >= 2) %>%
# 	select(-nr_years) %>%
# 	ungroup() %>%
# 	pivot_wider(names_from = index, values_from = jmean, names_prefix = "j")

# observer_effort_values <-dsi_data %>%
# 	unnest(cols = c(zscored_values)) %>%
# 	ungroup() %>%
# 	mutate(observer_effort = map_dbl(.x = focal_data, .f = get_mean_observer_effort)) %>%
# 	select(focal = sname, focal_grp = grp, start, end, observer_effort)
#


# temp_files <- my_subset_files %>%
# 	mutate(file_path = paste0("./data/my_subsets/", file)) %>%
# 	separate(file, sep = "_", into = c(NA, NA, "row_nr", NA)) %>%
# 	mutate(row_nr = as.numeric(row_nr)) %>%
# 	mutate(data = map(.x = file_path, .f = read.table, header = TRUE)) %>%
# 	left_join(select(iyol_sub, row_nr, focal = sname, grp, start, end),
# 						by = "row_nr") %>%
# 	group_by(row_nr, focal, grp, start, end) %>%
# 	mutate(zscored = map(.x = data, .f = make_zscored_values))


# temp_files %>%
# 	ungroup() %>%
# 	slice(1) %>%
# 	unnest(cols = c(zscored))
# 	mutate(top_dyads = map(.x = zscored, .f= get_top_grooms),
# 				 has_real_partners = map_lgl(.x = top_dyads, .f = has_real_partners),
# 				 nr_dyads = map_dbl(.x = top_dyads, .f = check_nr_dyads)) %>%



# DSI_zscored_values = my_subset_files %>%
# 	select(-data, file_path) %>%
# 	mutate(DSI_values = map(.x = zscored, .f = get_zscored_DSI_values)) %>%
# 	select(-zscored) %>%
# 	unnest(cols = c(DSI_values)) %>%
# 	rename_at(vars(starts_with("DSI")), function(x) paste0(x, "_zscored"))
#
# write_csv(DSI_zscored_values, paste0("./data/", "DSI_zscored_values_",
#  																		latest_version_date, ".csv", sep=""))
#
# DSI_zscored_values <- read_csv(paste0("./data/", "DSI_zscored_values_",
# 																												 latest_version_date, ".csv", sep="")) %>%
# 	inner_join(select(biograph_l, focal = sname, focal_birth = birth)) %>%
# 	mutate(age = as.numeric((start - focal_birth)/365.25))

## These values are based on universal zscored values
# dsi_paternal <- read_rds(paste("./data/", "dsi_paternal_",
# 													 latest_version_date, ".rds", sep=""))

# model_results <- str_remove_all(string = list.files(path = './data/bootstrap_results/', pattern = ".csv"), pattern = "_results.csv")


