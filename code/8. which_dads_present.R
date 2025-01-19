## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).

## The code below contains was used to generate the data for table 4.
## In this table we invesigate what variables prediced how long baboon fathers overlap with daughters/

## The final datasets are publicly available at https://zenodo.org/records/14590285
## The data that is publicly available has been anonymized.
## The key to get back to original snames is available up on request.


packages =
	c("broom"
		, "broom.mixed"
		, "cowplot"
		, "flextable"
		, "forcats"
		, "ftExtra"
		, "ggfortify"
		, "ggplot2"
		, "grid"
		, "kableExtra"
		, "lmerTest"
		, "lubridate"
		, "magick"
		, "MuMIn"
		, "patchwork"
		, "officer"
		, "purrr"
		, "rptR"
		, "scales"
		, "survival"
		, "survminer"
		, "survMisc"
		, "tidyverse"
	)

# install packages if needed and open libaries
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
	install.packages(setdiff(packages, rownames(installed.packages()))
									 , dependencies = TRUE,
									 repos = "http://cran.us.r-project.org")
}
lapply(packages, require, character.only = TRUE, warn.conflicts = TRUE, quietly = TRUE)

## This code uses a lots of functions inside tibbles using the purrr package
## The functions can be found here.
source('./code/3. functions_for_data_prep.R')
## This mainly consist of babase data
load("./data/data_set_for_dad_analysis_21JUN22.Rdata") ## this can be the older version
load("./data/datasets_for_paper_1MAR24.Rdata")
load("./data/Silk_figure_data_1MAR24.Rdata")

## ***Important note*** The DSI values in this analysis are based on the universal zscored values.
## In the data these are DSI_variable.universal.
## To make naming in the paper easier I renamed all of these to DSI_variable (i.e. I removed the .universal)

xdata_females_with_social <- xdata_females_with_social %>%
	rename_with(~str_remove(., '.universal')) %>%
	mutate(cumulative_adversity =  if_else(cumulative_adversity > 3, 3, cumulative_adversity))


## load data needed for this analysis
biograph_l <- read_csv("./data/biograph_l.csv")
parents_l <- read_csv("./data/parents_l.csv")
rankdates_l <- read_csv("./data/data/rankdates_l.csv")
load("./data/data_set_for_dad_analysis_21JUN22.Rdata") ## this can be the older version
consort_time <- read_csv("./data/ConsortTime_final_1Sep2024.csv")

members_l <- read_csv("./data/members_l.csv")
paternal_data <- read_csv("./data/paternal_data.csv")
rank_l <- read_csv("./data/ranks_l.csv")

## for this analysis we looked variables related to the dad..
## We focused on the year before the male left.

step1 <- xdata_females_with_social %>% ## 8 & 9
	select(focal, birth, mom, dad, dad_overlap_years, cumulative_adversity) %>%
	inner_join(select(parents_l, focal = kid, zdate)) %>%
	mutate(overlap_data = pmap( .l = list(focal, dad, birth, birth + years(4) - days(1))
															, .f = get_overlap)) %>% ## how long did they overlap
	unnest(cols = c(overlap_data), keep_empty = TRUE) %>%
	inner_join(select(biograph_l, mom = sname, mom_birth = birth)) %>%
	inner_join(select(biograph_l, dad = sname, dad_birth = birth)) %>%
	mutate(mom_age = as.numeric(first_breakup - mom_birth)/365.25,
				 dad_age = as.numeric(first_breakup - dad_birth)/365.25) %>%
	mutate(male_rank = pmap(.l = list(dad, last_overlap_grp, first_breakup),
													.f = get_male_last_rank)) %>% ## get the male o
	unnest(cols = c(male_rank), keep_empty = TRUE) %>%
	mutate(daily_d_days_rate = pmap_dbl(.l = list(dad, NA, first_breakup %m-% months(1) %m+% days(1),
																								first_breakup, per_group = TRUE),
																			.f = get_d_days_rate)) %>% ## measure of access to fertile females
	left_join(select(consort_time, focal, dad = actor, proportion_consort_time)) %>% ## 4
	## filter out cases were we don't have consort data time daugher was conceived
	filter(focal %in% unique(consort_time$focal)) %>%
	## if a male did not have any consort he would jave an NA which in this case would be a zero
	## We already filtered out the juvenile females for whom no consort data is available
	mutate(proportion_consort_time = if_else(is.na(proportion_consort_time), 0, proportion_consort_time)) %>%
	mutate(nr_pdads = map_dbl(.x = focal, .f = get_nr_pdads)) %>% ##  how many potential dads were around at conceiving
	mutate(kids = pmap(.l = list(focal, birth, mom, dad),
										 .f = get_kids)) %>% ## ## get the previous and next offspring of the mother
	unnest(cols = c(kids), keep_empty = TRUE) %>%
	mutate(offspring_years = pmap_dbl(
		.l = list(focal, dad, focal_grp = NA, first_breakup - years(1) + days(1),
							first_breakup, per_group = FALSE),
				.f = get_offspring_years)) ## how many paternal offspring does the make have


## checking for missing data
step1 %>%
	select(mom_age, dad_age, male_rank, daily_d_days_rate, proportion_consort_time, nr_pdads, offspring_years,
	 			 previous_kid_with_mom, next_kid_with_mom, offspring_years) %>%
	filter(if_any(everything(), ~is.na(.))) %>%
	select_if(function(x) any(is.na(x)))

# We have some missing data related to ranks

## First we will work on the cases were the next offspring of the mother is know
dad_overlap_only_know_data <- step1 %>%
	filter(previous_kid_with_mom != 0.5) %>%
	filter(!is.na(male_rank)) %>%
	select(focal,
				 dad_overlap_years,
				 mom_age, dad_age,
				 male_rank,
				 daily_d_days_rate,
				 proportion_consort_time,
				 nr_pdads,
				 previous_kid_with_mom,
				 offspring_years,
				 cumulative_adversity,
				 dad)

## save the data
dad_overlap_only_know_data %>%
	write_csv(file = "./data/dad_overlap_only_know_data.csv")

## this is the model
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

## run the model
dad_overlap_only_known_model_prop_consort <- lmer(dad_overlap_only_known_model_prop_consort_formula,
																									data = dad_overlap_only_know_data
 																									, na.action = 'na.fail')

## extract the model data
dad_overlap_only_known_model_prop_consort_results <- dad_overlap_only_known_model_prop_consort %>%
	broom.mixed::tidy() %>%
	left_join(vif.mer(dad_overlap_only_known_model_prop_consort), by = 'term') %>%
	mutate(run = paste0("model_", 0)) %>%
	mutate(row = NA) %>%
	select(row, everything()) %>%
	mutate(term = if_else(str_detect(term, "rank"), "dad_rank", term))

dad_overlap_only_known_model_prop_consort_results



## We can also try and bootstrap the cases were we don't know the previous offspring
## see who bootstrapped for more details

nr_draws = 1000

step1_previous <- step1 %>%
	filter(previous_kid_with_mom == .5) %>%
	group_by(focal, dad, birth, start, end, mom) %>%
	nest() %>%
	ungroup() %>%
	mutate(period = "previous") %>%
	mutate(get_previous = pmap(.l = list(focal, birth, focal_grp=NA, mom, period),
														 .f = get_maternal_sibling)) %>%
	unnest(cols= c(get_previous), keep_empty = TRUE) %>%
	mutate(male_counts = pmap(.l = list(offspring_sname, dad, start, end),
														.f = get_male_present_count_overall)) %>%
	unnest(cols= c(male_counts)) %>%
	mutate(pdads_sname = pmap(.l = list(offspring_sname, dad, start, end),
														.f = get_pdad_sname_overall)) %>%
	mutate(previous_random_draws = pmap(.l = list(nr_potential_dads, nr_males, pdad_present, pdads_sname, draws = nr_draws, period),
																			.f = get_random_draws)) %>%
	select(focal, start, end, previous_random_draws)

step1_bootstrap <- step1_previous %>%
	mutate(data = pmap(.l = list(data, previous_random_draws, next_random_draws),
										 .f = assign_random_draws_overall)) %>%
	unnest(cols = c(data)) %>%
	select(-previous_random_draws, -next_random_draws)

## save this data and run on parrall computer
## see who bootstrapped for details
