## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).
## To get access to the database contact the corresponding author.

library(tidyverse)

source('./code/3. functions_for_data_prep.R')

## open datasets
latest_version_date <- read_lines("./data/latest_version_date.txt")
load(paste0("./data/data_set_for_dad_analysis_", latest_version_date,".Rdata"))
load(paste0("./data/datasets_for_paper_", latest_version_date, ".Rdata"))
load(paste0("./data/Silk_figure_data_", latest_version_date, ".Rdata"))

## ***Important note*** The DSI values in this analysis are based on the universal zscored values. ## In the data these are DSI_variable.universal.
## To make naming in the paper easier I renamed all of these to DSI_variable (i.e. I removed the .universal)

xdata_females_with_social <- xdata_females_with_social %>%
	rename_with(~str_remove(., '.universal')) %>%
	mutate(cumulative_adversity =  if_else(cumulative_adversity > 3, 3, cumulative_adversity))


#members_raw <- load('./data/members_raw.Rdata')
members_l <- read_csv('./data/members_l.csv')
members_AM_l <- members_l %>%  ## Only ranked males
	select(-sys_period) %>%
	inner_join(rankdates_l) %>%
	filter(date > ranked) %>%
	select(sname, grp, date, ranked)

members_juveniles_l <- members_l %>%
	select(-sys_period)
	inner_join(select(biograph_l, sname, birth), by = 'sname') %>%
	mutate(age = as.numeric(date- birth)/365.25) %>%
	filter(age < 4) %>%
	inner_join(select(parents_l, sname = kid, dad))


## the starting point here are the focal females (per year/group)
who_grooms_step1 <- social_indexes_males %>%
	select(focal, focal_grp, start, end, contains('DSI')) %>%
	inner_join(select(xdata_females_with_social,
	                  focal, mom, dad, maternal_loss, maternal_rank, maternal_SCI_F,
	                  density, sibling, drought, cumulative_adversity), by = "focal") %>%
	inner_join(select(biograph_l, focal = sname, focal_birth = birth), by = "focal") %>%
	inner_join(select(biograph_l, mom = sname, mom_birth = birth), by = "mom") %>%
	mutate(kid_age = as.numeric(start - focal_birth)/365.25) %>%
	mutate(mom_age = as.numeric(start - mom_birth)/365.25) %>%  # mom age
    ## How many days is the focal present
	mutate(nr_focal_days = pmap_dbl(
		.l = list(focal, focal_grp, start, end),
		.f = get_focal_days)) %>%
	ungroup()

who_grooms_step2 <- who_grooms_step1 %>%
	inner_join(observer_effort_values) %>%
	ungroup()

who_grooms_step3a <- who_grooms_step2 %>%
	## in the time slot we are looking at which which males were present
    mutate(all_males = pmap(.l = list(focal_grp, start, end), .f = get_all_males)) %>%
	unnest(cols = c(all_males)) %>%
    ## now we have a list that has a rows per focal female and every male present
	rename(AMales = sname)

who_grooms_step3b <- who_grooms_step3a  %>%
    ## here we add the dydic values for the focal-male dyad
	left_join(male_values %>%
	              mutate(does_groom = i_total > 1/365.25) %>%
	              select(focal_check, focal, focal_grp, start, AMales,
	                     paternal_groom, paternal_res_i_adj_zscored, does_groom) %>%
	              filter(focal_check == TRUE),
	          by = c("focal", "focal_grp", "start", "AMales")) %>%
	filter(!is.na(paternal_res_i_adj_zscored))

who_grooms_step3c <- who_grooms_step3b %>%
		group_by(focal, focal_grp, start) %>%
	mutate(nr_males = n_distinct(AMales),
				 nr_grooms = sum(!is.na(paternal_data))) %>%
	group_by(AMales) %>%
	mutate(is_dad = AMales == dad) %>%
	inner_join(select(biograph_l, AMales = sname, AMales_birth = birth)) %>%
	mutate(AMales_age = as.numeric((start - AMales_birth)/365.25)) %>%
	ungroup()

## what was the rank for the male
who_grooms_step4 <- who_grooms_step3c %>%
	mutate(AMales_rank = pmap(.l = list(AMales, focal_grp, start, end),
														.f = get_mean_rank)) %>%
	unnest(cols = c(AMales_rank), keep_empty = TRUE)

who_grooms_step5 <- who_grooms_step4 %>%
    ## was the male present in period focal female was concived and was maiting behaviour observed
	left_join(select(potential_dads_l, focal = kid, AMales = pdad,
									 pdad_status = status,
									 estrous_presence, estrous_me, estrous_c),
						by = c('AMales', "focal")) %>%
	mutate(estrous_presence =if_else(is.na(estrous_presence), 0, estrous_presence),
				 estrous_c =if_else(is.na(estrous_c), 0, estrous_c),
				 estrous_me =if_else(is.na(estrous_me), 0, estrous_me)) %>%
    ## how many fertile females are around
	mutate(daily_d_days_rate = pmap_dbl(.l = list(AMales, focal_grp, start, end),
	                                    .f = get_d_days_rate)) %>%
    ## how many potential dads are around
	mutate(nr_pdads = map_dbl(.x = focal, .f = get_nr_pdads)) %>%
    ## how many offspring does male have in the group (in offspring years)
	mutate(offspring_years = pmap_dbl(
		.l = list(focal, AMales, focal_grp, start, end),
		.f = get_offspring_years))

who_grooms_step6 <- who_grooms_step5 %>%
	ungroup() %>%
	mutate(kids = pmap(.l = list(focal, focal_birth, mom, AMales),
										 .f = get_kids)) %>%
	unnest(cols = c(kids), keep_empty = TRUE)

who_grooms_step6 %>%
	#select(-zscore) %>%
	filter(if_any(everything(), ~is.na(.))) %>%
	select_if(function(x) any(is.na(x)))

## We no longer look at hybrid but can be included ny uncommenting.
# who_grooms_step7 <- who_grooms_step6 %>%
# 	left_join(select(hybridgene_scores_l, AMales = sname, hybrid_score = score)) %>%
# 	left_join(select(anubis_estimates, AMales = sname, anubis_admix, notes) %>%
# 							filter(is.na(notes) & !is.na(anubis_admix)))

who_grooms_step8 <- who_grooms_step6 %>%
	select(focal, focal_birth, focal_grp, mom, dad, AMales, start, end, nr_days
				 ## ages
				 , kid_age, mom_age, AMales_age,
				 , observer_effort
				 ## early adversity
				 , maternal_loss, maternal_rank, maternal_SCI_F, density,
				 sibling, drought, cumulative_adversity
				 ## male details
				 , is_dad, mean_ordrank, mean_proprank,
				 , contains('estrous'), , daily_d_days_rate, offspring_years, nr_pdads
				 , previous_kid_with_mom, next_kid_with_mom ## see bootstrap
				 ## for subset
				 ##, hybrid_score, anubis_admix
				 ## responses
				 , does_groom,  paternal_res_i_adj_zscored) %>%
	filter(!is.na(mean_proprank)) %>%  ## There are a few issues with male ranks
	filter(nr_days >= 30) ## Males have to be in group for at least 30 days

## Not all previous or next siblings have know dads.
## We initially assigned a .5 value to the binary variable if a malke was the dad
## We then decided to use bootstrapping to use 'randonmly' assigned paterity.
## Only males that were potential dads can be assigned paternity

nr_draws = 1000

who_grooms_step_previous <- who_grooms_step8 %>%
	filter(previous_kid_with_mom == .5) %>%
	group_by(focal, focal_birth, focal_grp, start, end, mom) %>%
	nest() %>%
	ungroup() %>%
	mutate(period = "previous") %>%
	mutate(get_previous = pmap(.l = list(focal, focal_birth, focal_grp, mom, period),
														 .f = get_maternal_sibling)) %>%
	unnest(cols= c(get_previous), keep_empty = TRUE) %>%
	mutate(male_counts = pmap(.l = list(offspring_sname, focal_grp, start, end),
														.f = get_male_present_count)) %>%
	unnest(cols= c(male_counts)) %>%
	#mutate(pdads_sname = map2(.x = data, .y = period, .f = get_pdads)) %>%
	mutate(pdads_sname = pmap(.l = list(offspring_sname, focal_grp, start, end),
														.f = get_pdad_sname)) %>%
	mutate(previous_random_draws = pmap(
	    .l = list(nr_potential_dads, pdad_present, pdads_sname,
	              draws = nr_draws, period),
		.f = get_random_draws)) %>%
	select(focal, focal_grp, start, end, previous_random_draws)

who_grooms_step_next <- who_grooms_step8 %>%
	filter(next_kid_with_mom == .5) %>%
	group_by(focal, focal_birth, focal_grp, start, end, mom) %>%
	nest() %>%
	ungroup() %>%
	mutate(period = "next") %>%
	mutate(get_next = pmap(.l = list(focal, focal_birth, focal_grp, mom, period),
												 .f = get_maternal_sibling)) %>%
	unnest(cols= c(get_next), keep_empty = TRUE) %>%
	mutate(male_counts = pmap(.l = list(offspring_sname, focal_grp, start, end),
														.f = get_male_present_count)) %>%
	unnest(cols= c(male_counts)) %>%
	#mutate(pdads_sname = map2(.x = data, .y = period, .f = get_pdads)) %>%
	mutate(pdads_sname = pmap(.l = list(offspring_sname, focal_grp, start, end),
														.f = get_pdad_sname)) %>%
	mutate(next_random_draws = pmap(
	    .l = list(nr_potential_dads, pdad_present,  pdads_sname,
	              draws = nr_draws, period),
	    .f = get_random_draws)) %>%
	select(focal, focal_grp, start, end, next_random_draws)

who_grooms_step_combined <- who_grooms_step8 %>%
	group_by(focal, focal_birth, focal_grp, start, end, mom) %>%
	nest() %>%
	left_join(who_grooms_step_previous) %>%
	left_join(who_grooms_step_next)

who_grooms_bootstrap <- who_grooms_step_combined %>%
	mutate(data = pmap(.l = list(data, previous_random_draws, next_random_draws),
										 .f = assign_random_draws)) %>%
	unnest(cols = c(data)) %>%
	select(-previous_random_draws, -next_random_draws) %>%
	filter(!is.na(mean_proprank)) %>%
	filter(nr_days >= 30)

who_grooms_bootstrap <- who_grooms_bootstrap %>%
	inner_join(xdata_females_with_social %>%
						 	select(focal, competing_sibling, contains("binary")))

write_csv(who_grooms_bootstrap,
          paste0("./data/who_grooms_bootstrap_", latest_version_date, ".csv"))


## adding the new consort prop data (calculated by Beth)
consort_prop <- read_csv("./data/ConsortTime_final_1Sep2024.csv")


selected_consort_focals <- unique(consort_prop$focal)

who_grooms_bootstrap <- who_grooms_bootstrap %>%
	left_join(select(consort_prop, focal, mom, AMales = actor, proportion_consort_time)) %>%
	mutate(has_consort_data = focal %in% selected_consort_focals) %>%
	mutate(proportion_consort_time = case_when(
	    focal %in% selected_consort_focals & is.na(proportion_consort_time) ~ 0,
	    !(focal %in% selected_consort_focals) ~ NA,
	     TRUE ~ proportion_consort_time)) %>%
	write_csv(paste0("./data/who_grooms_bootstrap_", latest_version_date, ".csv"))


# move who_grooms_bootstrap.csv to cluster and run cluster_bootstrap_previous_next_code.R
# and get the results back here

result_files <- list.files(path = './data/bootstrap_results/raw_files/', pattern = "full")

for(ii in 1:length(result_files)) {
	xdata <- read_delim(file = paste0('./data/bootstrap_results/raw_files/',
																		result_files[ii]),
											delim = " ", col_names = FALSE)


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
	summary_filename <- str_replace(string = result_files[ii], pattern = "_raw.text",
															replacement = "_summary.csv")

	write_csv(xdata, paste0("./data/bootstrap_results/clean_files/", new_filename))

	xdata %>%
		group_by(term) %>%
		summarise_all(.funs = mean) %>%
		mutate(is_sig = p.value < 0.05) %>%
		write_csv(paste0("./data/bootstrap_results/summary_files/", summary_filename ))

	print(paste0(new_filename, " has been saved and summarized"))
}
