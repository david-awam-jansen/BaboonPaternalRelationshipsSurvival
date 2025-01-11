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

members_l <- read_csv("./data/members_l.csv")

#load() ## TBD see datasets for papers
paternal_data <- read_csv("./data/paternal_data.csv")

#members_raw <- load('./data/members_raw.Rdata')
biograph_l <- read_csv("./data/biograph_l.csv")
#members_l <- read_csv('./data/members_l.csv')
parents_l <- read_csv("./data/parents_l.csv")
rankdates_l <- read_csv("./data/data/rankdates_l.csv")
load("./data/data_set_for_dad_analysis_21JUN22.Rdata") ## this can be the older version
consort_time <- read_csv("./data/ConsortTime_final_1Sep2024.csv")
rank_l <- read_csv("./data/ranks_l.csv")

make_overlap_plot <- function(focal, mom, dad, start, focal_grp) {
	start = ymd(start)
	pre = start - years(1) - months(6)
	end = start + years(4) - days(1)

	zdate_focal = parents_l %>%
		filter(kid == focal) %>%
		select(zdate) %>%
		pull()


	members_l %>%
		mutate(grp = as.factor(grp)) %>%
		filter(sname %in% c(focal, mom, dad)) %>%
		filter(date > pre) %>%
		filter(date < end) %>%
		mutate(sname = case_when(sname == focal ~ paste(sname, "(focal)"),
														 sname == mom ~ paste(sname, "(mom)"),
														 TRUE ~ paste(sname, "(dad)"
														 ))) %>%
		mutate(is_focal_grp = grp == focal_grp) %>%
		ggplot(aes(x = date, y = sname, color = grp, size = is_focal_grp)) +
		geom_point() +
		geom_vline(xintercept = start, color = "black") +
		geom_vline(xintercept = zdate_focal, linetype = "dashed") +
		cowplot::theme_cowplot()

}

step1 <- xdata_females_with_social %>% ## 8 & 9
	select(focal, birth, mom, dad, dad_overlap_years, cumulative_adversity) %>%
	inner_join(select(parents_l, focal = kid, zdate)) %>%
	mutate(overlap_data = pmap( .l = list(focal, dad, birth, birth + years(4) - days(1))
															, .f = get_overlap)) %>%
	unnest(cols = c(overlap_data), keep_empty = TRUE) %>%
	inner_join(select(biograph_l, mom = sname, mom_birth = birth)) %>% ## 1
	inner_join(select(biograph_l, dad = sname, dad_birth = birth)) %>% ## 1
	mutate(mom_age = as.numeric(first_breakup - mom_birth)/365.25,
				 dad_age = as.numeric(first_breakup - dad_birth)/365.25) %>%
	mutate(male_rank = pmap(.l = list(dad, last_overlap_grp, first_breakup),
													.f = get_male_last_rank)) %>% ## 2
	unnest(cols = c(male_rank), keep_empty = TRUE) %>%
	mutate(daily_d_days_rate = pmap_dbl(.l = list(dad, NA, first_breakup %m-% months(1) %m+% days(1),
																								first_breakup, per_group = TRUE),
																			.f = get_d_days_rate)) %>% ## 3
	filter(focal %in% unique(consort_time$focal)) %>%
	left_join(select(consort_time, focal, dad = actor, proportion_consort_time)) %>% ## 4
	mutate(proportion_consort_time = if_else(is.na(proportion_consort_time), 0, proportion_consort_time)) %>%
	mutate(nr_pdads = map_dbl(.x = focal, .f = get_nr_pdads)) %>% ## 5
	mutate(kids = pmap(.l = list(focal, birth, mom, dad),
										 .f = get_kids)) %>% ## 6
	unnest(cols = c(kids), keep_empty = TRUE) %>%
	mutate(offspring_years = pmap_dbl(
		.l = list(focal, dad, focal_grp = NA, first_breakup - years(1) + days(1),
							first_breakup, per_group = FALSE),
				.f = get_offspring_years)) ## 7

## dealing with some missing ranks
step1 %>%
	select(mom_age, dad_age, male_rank, daily_d_days_rate, proportion_consort_time, nr_pdads, offspring_years,
	 			 previous_kid_with_mom, next_kid_with_mom, offspring_years) %>%
	filter(if_any(everything(), ~is.na(.))) %>%
	select_if(function(x) any(is.na(x)))




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
	#mutate(pdads_sname = map2(.x = data, .y = period, .f = get_pdads)) %>%
	mutate(pdads_sname = pmap(.l = list(offspring_sname, dad, start, end),
														.f = get_pdad_sname_overall)) %>%
	mutate(previous_random_draws = pmap(.l = list(nr_potential_dads, nr_males, pdad_present, pdads_sname, draws = nr_draws, period),
																			.f = get_random_draws)) %>%
	select(focal, start, end, previous_random_draws)

step1_next <- step1 %>%
	mutate(start = first_breakup - years(1) + days(1),
				 end = first_breakup) %>%
	filter(next_kid_with_mom == .5) %>%
	group_by(focal, dad, birth, start, end, mom) %>%
	nest() %>%
	ungroup() %>%
	mutate(period = "next") %>%
	mutate(get_previous = pmap(.l = list(focal, birth, focal_grp=NA, mom, period),
														 .f = get_maternal_sibling)) %>%
	unnest(cols= c(get_previous), keep_empty = TRUE) %>%
	mutate(male_counts = pmap(.l = list(offspring_sname, dad, start, end),
														.f = get_male_present_count_overall)) %>%
	unnest(cols= c(male_counts)) %>%
	#mutate(pdads_sname = map2(.x = data, .y = period, .f = get_pdads)) %>%
	mutate(pdads_sname = pmap(.l = list(offspring_sname, dad, start, end),
														.f = get_pdad_sname_overall)) %>%
	mutate(next_random_draws = pmap(.l = list(nr_potential_dads, nr_males, pdad_present, pdads_sname, draws = nr_draws, period),
																			.f = get_random_draws)) %>%
	select(focal, start, end, next_random_draws)

step1_combined <- step1 %>%
	mutate(start = first_breakup - years(1) + days(1),
				 end = first_breakup) %>%
	group_by(focal, birth, start, end, mom) %>%
	nest() %>%
	left_join(step1_previous) %>%
	left_join(step1_next)

step1_bootstrap <- step1_combined %>%
	mutate(data = pmap(.l = list(data, previous_random_draws, next_random_draws),
										 .f = assign_random_draws_overall)) %>%
	unnest(cols = c(data)) %>%
	select(-previous_random_draws, -next_random_draws)

step1_bootstrap %>%
	write_csv("./data/step1_28SEP24.csv")

step1 <- read_csv("./data/step1_28SEP24.csv")

dad_overlap_only_known_model_prop_consort_formula <- dad_overlap_years ~ 1 +
	mom_age+ dad_age +
	cumulative_adversity +
	male_rank +
	proportion_consort_time +
	daily_d_days_rate +
	offspring_years +
	nr_pdads +
	previous_kid_with_mom +
	#next_kid_with_mom +
	(1|dad)

# step1 %>%
# 	filter(previous_kid_with_mom != 0.5 #& next_kid_with_mom != 0.5
# 				 ) %>%
# 	filter(!is.na(rank)) %>%
# 	select(mom_age,  dad_age, cumulative_adversity, rank, proportion_consort_time,
# 				 daily_d_days_rate, offspring_years, nr_pdads, previous_kid_with_mom, next_kid_with_mom) %>%
# 	View()


185 - step1 %>%
	filter(previous_kid_with_mom != 0.5) %>% ##
	filter(!is.na(male_rank)) %>%  ##
	nrow()

dad_overlap_only_know_data <- step1 %>%
	filter(previous_kid_with_mom != 0.5) %>%
	filter(!is.na(male_rank))

step1 %>%
	filter(is.na(male_rank)) %>%
	select(focal, dad, last_overlap_grp, birth, zdate, first_breakup)


dad_overlap_only_know_data %>%
	select(dad) %>%
	n_distinct()


dad_overlap_only_known_model_prop_consort <- lmer(dad_overlap_only_known_model_prop_consort_formula,
																									data = dad_overlap_only_know_data
#																									, family = poisson
#																									, control=
#																										glmerControl(optimizer = "bobyqa",
#																																 optCtrl=list(maxfun=2e5))
																									, na.action = 'na.fail')

dad_overlap_only_known_model_prop_consort_results <- dad_overlap_only_known_model_prop_consort %>%
	broom.mixed::tidy() %>%
	left_join(vif.mer(dad_overlap_only_known_model_prop_consort), by = 'term') %>%
	mutate(run = paste0("model_", 0)) %>%
	mutate(row = NA) %>%
	select(row, everything()) %>%
	mutate(term = if_else(str_detect(term, "rank"), "dad_rank", term))

write.table(x = dad_overlap_only_known_model_prop_consort_results,
					"./data/bootstrap_results/dad_overlap_only_known_model_prop_consort_results.txt")







make_model_flextable("dad_overlap_only_known_model_prop_consort", dataset = step1, explain_text = explain_text)
model <- "dad_overlap_only_known_model_prop_consort"

dad_overlap_only_known_model_prop_consort_results %>%
	arrange(p.value) %>%
	flextable() %>%
	colformat_double(digits = 3)



t1 <- read.table(paste0('./data/bootstrap_results/', model, "_results.txt")) %>%
	as_tibble() %>%
	mutate(term = if_else(str_detect(term, "rank"), "dad_rank", term))

explain_text <- c(
	c(""),
	c("\U2191 mothers' age  \U2191 longer co-residency"),
	c("\U2191 fathers' age  \U2191 longer co-residency"),
	c("\U2191 paternal offspring \U2191 longer co-residency"),
	c("\U2191 increased maiting opportunities \U2191 longer co-residency"),
	rep(c(""),5)

)

t1 <- t1 %>%
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

if(!'df' %in% names(t1)) t1 <- t1 %>% add_column(df = NA)

t1 <- t1 %>%  select(-df)

names(t1)[2:5] <- c('\U03B2','\U03C3', 't', 'p')


t1 %>%
	select(-run) %>%
	mutate(term = if_else(term == "rank", "Rank of father", term)) %>%
	flextable() %>%
	colformat_double(digits = 3) %>%
	theme_vanilla() %>%
	# add_footer_row(values = paste0('There are ',
	# 															 nrow(dataset), " data points in this analyis with ",
	# 															 n_distinct(dataset$AMales), " adult males."),
	# 							 colwidths = ncol(t1)) %>%
	fontsize(size = 8, part = "all") %>%
	align(align = "left", part = "header") %>%
	align(j = 2:5, align = "center", part = "body" ) %>%
	#flextable::rotate(j = 2:5, rotation="btlr",part="header") %>%
	width(j = 1, width = 2.2, unit = "in") %>%
	width(j = 2:5, width = .5, unit = "in") %>%
	width(j = 6, width =3, unit = "in") %>%
	height(height = 2,part = "header") %>%
	valign(valign = "bottom", part = "header") %>%
	valign(valign = "top", part = "body")






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
