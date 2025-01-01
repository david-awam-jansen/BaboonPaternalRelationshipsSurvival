## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).

## The code below creates the eealy adversity data.
## For details on the early adversity see Tung, J., Archie, E.A., Altmann, J. &amp; Alberts, S.C. 2016 Cumulative early adversity predicts longevity in wild baboons. Nature Communications 7, 11181.
## To get access to the database contact the corresponding author.

## list of packages needed
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

latest_version_date <- read_lines("./data/latest_version_date.txt")
load(paste0("./data/data_set_for_dad_analysis_", latest_version_date, ".RData"))

xdata_l <- biograph_l %>%
  filter(## exclude animals born before the onset of monitoring
         birth > '1971-07-01'
         ## exclude animals born before rainfall estimates
         & birth > '1976-08-15'
         & matgrp < 3             ## exclude lodge group
         ## exclude proton's group because the data are so thin
         & matgrp != 1.300
         & bstatus == 0           ## only use animals with good age estimates
         & statdate > birth       ## exclude animals who lived for < 1 day
         ## exclude individuals that died before they could be sexed
         & sex != 'U'

  ) %>%
  filter(birth + years(4) < dmy(latest_version_date)) %>%
  inner_join(select(parents_l, sname = kid, mom, dad)) %>%
  mutate(birth_rnkdate = floor_date(birth, 'month'))

## i. Density
## Number of adults in group on date of birth
## Only use if group was a was a study group
## Note this is stricter the Tung et al 2016 analysis
## Let for instance look at an example of filtering out members data based on
## groups history and behave gaps

members_l <- members_l %>%
  inner_join(select(biograph_l, sname, sex), by = c("sname")) %>%
  left_join(select(maturedates_l, sname, matured), by = c("sname"))

adult_group_sizes_l <- members_l %>%
  filter(date > matured) %>%
  select(grp, date, sname) %>%
  inner_join(xdata_l %>%
               select(grp = matgrp, date = birth) %>%
               distinct(),
             by = c("grp", "date")) %>%
  group_by(grp, date) %>%
  summarise(nr_adults = n()) %>%
  arrange(desc(date)) %>%
	ungroup()

## create a temp tibble
xdata_l <- xdata_l %>%
  left_join(adult_group_sizes_l %>%
  						select(matgrp = grp, birth = date, density = nr_adults),
  					by = c("birth", "matgrp")) %>%
	filter(!is.na(density)) ## exclude individuals for who mo density

## ii. maternal rank
maternal_ranks_l <- ranks_l %>%
  filter(rnktype == 'ADF') %>%
  select(sname, grp, rnkdate, rank) %>%
	 inner_join(xdata_l %>%
  					 	select(sname = mom, kid = sname, grp = matgrp,
  					 				 rnkdate = birth_rnkdate),
  					 by = c("sname", "grp", "rnkdate")) %>%
  select(mom = sname, sname = kid, maternal_rank = rank)

xdata_l <- xdata_l %>%
	left_join(select(maternal_ranks_l, mom, sname, maternal_rank)) %>%
	filter(!is.na(maternal_rank))

## iii. Rainfall and drought
get_rainfall <- function(focal_birth) {
  rainfall_l %>%
    filter(date >= focal_birth,
           date <= focal_birth + years(1)) %>%
    select(rain) %>%
    pull() %>%
    sum(na.rm = TRUE)
}

xdata_l <- xdata_l %>%
  mutate(rainfall = map_dbl(.x = birth, .f = get_rainfall)) %>%
  mutate(drought = rainfall <= 200)

## iv. sibling
get_siblings <- function(focal_birth, focal_mom) {
  sibblings_temp <- pregdata_l %>%
    filter(mom == focal_mom) %>%
    filter(birth > ymd(focal_birth),
           birth <= ymd(focal_birth) + months(18)) %>%
    select(sibling = kid, sibling_birth = birth, sibling_dad = dad)

    if(nrow(sibblings_temp) == 0) {
      tibble(sibling = NA, sibling_birth = NA
             , sibling_dad = NA
             , competing_sibling = 0)
    } else(
        sibblings_temp %>%
      mutate(competing_sibling = 1)
    )
}

xdata_l <- xdata_l %>%
	 mutate(sibling = map2(.x = birth, .y = mom, .f = get_siblings)) %>%
	 unnest(cols = c(sibling), keep_empty = TRUE)

xdata_l <- xdata_l %>%
  mutate(age_at_birth_sibling = (as.numeric((ymd(sibling_birth) - ymd(birth))/365.25)))

## v. maternal loss
xdata_l <- xdata_l %>%
  inner_join(maternal_loss_l)

## vi. maternal SCI to be added
## check https://github.com/amboseli/ramboseli/blob/master/documentation/sociality-indices.md fow these are calculated
maternal_sci_values <- read_csv(paste0("./data/maternal_sci_values_31MAY22.csv"))

xdata_l <- xdata_l %>%
	left_join(select(maternal_sci_values, mom = sname, sname = kid, years, contains("SCI"))) %>%
	rename_with(.fn = ~ gsub("SCI", "maternal_SCI", .x, fixed = TRUE),
							.col = starts_with("SCI")) %>%
	filter(!is.na(maternal_SCI_F))

## For the early adversity analysis we often use a binary variable
## For example low versus high density
xdata_l <- xdata_l %>%
    mutate(
        maternal_rank_binary = maternal_rank >= quantile(maternal_rank, prob=c(0.75)),
        density_binary = density >= quantile(density, prob=c(0.75)),
        maternal_SCI_F_binary =  maternal_SCI_F <= quantile(maternal_SCI_F, prob=c(0.25)))

## parity
xdata_l <- xdata_l %>%
    ## Every offspring of a female is given a pid.
    ## This numeric part of this code would be the number of the offspring
	mutate(parity = readr::parse_number(pid),
				 first_born = parity == 1)

## 2 of the groups lived in a habitat that had a low avaiability of food.
## After the moved they had a better habitat.
## All other groups were always in this higher quality habitat
xdata_l <- xdata_l %>%
    mutate(habitat_quality = case_when(
        matgrp == 1 & year(ymd(birth)) <= 1987 ~  "low",
		matgrp == 2 & year(ymd(birth)) <= 1991 ~ "low",
		TRUE ~ "high"),
		habitat_quality = factor(habitat_quality))

## The cumulative adversity is the sum of the binary adversities
xdata_l <- xdata_l %>%
	mutate(cumulative_adversity =
				 	sum(density_binary + drought + competing_sibling +
				 			maternal_loss + maternal_rank_binary +
				 			    maternal_SCI_F_binary
				 				))

## save thea early adversity data set
write.csv(x = xdata_l,
					file = paste("./data/", "ea_dataset_less_restricted_", latest_version_date, ".csv", sep=""))

## Quick overview of where there is missing data
desc_missing <- NULL

numbers <- xdata_l %>%
	dplyr::group_by(sex) %>%
	summarise(number = n())

p1 <- xdata_l %>%
	select(sex,dad, density, maternal_rank, maternal_loss, maternal_SCI_F, competing_sibling, drought, cumulative_adversity) %>%
	dplyr::group_by(sex) %>%
	naniar::miss_var_summary() %>%
	#mutate((n_miss/pct_miss) * 100)
	dplyr::mutate(variable = factor(variable,
	                                levels = (sort(unique(variable),
	                                               decreasing = TRUE))))  %>%
	# forcats::fct_relevel(rev(desc_missing), after = Inf) %>%
	# forcats::fct_relevel(rev(custom_order), after = Inf)) %>%
	left_join(numbers, by = c("sex")) %>%
	mutate(number = paste(n_miss, " out of ", number, sep = "\n")) %>%
	ggplot(aes(sex, variable, fill = pct_miss)) +
	geom_tile(color = "black", size = 2)  +
	geom_text(aes(label = number), color = "black", size = 3)  +
	#facet_wrap(~sex) +
	scale_fill_gradient2(low = "white", high = "firebrick",
											 limits = c(0, 100), name = "% Miss") +
	theme_minimal() + theme(legend.position = "bottom") +
	ggtitle("Missing data overview for included juveniles")

p1

ggsave(plot = p1, filename = "./results/overview_of_missing_early_adversity.jpeg",
			 width =5,height =10)
