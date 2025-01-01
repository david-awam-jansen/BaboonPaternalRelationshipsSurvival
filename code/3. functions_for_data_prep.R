## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).

## The code below contains an overview of custom functions.
## These functions were used for the waggling of data, analysis, creating of tables and figures
## Some of these were written before changes were main to tidyverse packages.
## They need to be alter to work with newer packages.

## They are also not in the most logical order, but I attempted to group them roughly per analysis

## get a zscore
get_zscore <-  function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)

## calculate the VIF factor for a linear mixed model
### Based on code by A.  Franks
## https://hlplab.wordpress.com/2011/02/24/diagnosing-collinearity-in-lme4/
vif.mer <- function (fit) {
    ## adapted from rms::vif

    v <- vcov(fit)
    nam <- names(fixef(fit))

    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
        v <- v[-(1:ns), -(1:ns), drop = FALSE]
        nam <- nam[-(1:ns)]
    }

    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam

    v %>% as_tibble() %>%
        mutate(term = nam) %>%
        rename(VIF = value)
}

get_focal_days  <- function(focal, focal_grp, start, end) {
    start_date <- ymd(start)
    end_date <- ymd(end)

    members_partner_overlap_l %>%
        filter(sname == focal
               & date >= start_date
               & date <= end_date
               & (grp == focal_grp |  grpofresidency== focal_grp )) %>%
        group_by(sname) %>%
        summarise(nr_days = n(), .groups = 'drop') %>%
        pull()
}

get_mean_observer_effort <-function(df) {

    if(nrow(df) > 1) {
        xx <- df %>%
            ungroup() %>%
            select(mean_AF_log2OE) %>%
            filter(!is.na(mean_AF_log2OE))

        mean(xx$mean_AF_log2OE, na.rm = TRUE)
    } else {
        NA
    }
}

## Note.

## For cases where e.g. the father was not present a 'fake' dad was created that had very low grooming effort (e.g. 1 groom per year). These fake males were given the sname XXX so I could easliy keep track of them. Something similar was also  done for the calculations of the SCI values. As the male was not really present some values such as observer effort could be calculated. The functions below help to deal with some of these issue and/or filter out the fake males.

get_mean_AF_log2OE <- function(df) {
		df %>%
			filter(dyad_type == "AF-JUV") %>%
			group_by(partner, partner_grp) %>%
			summarise(mean_AF_log2OE = log2(mean(OE)), .groups = "drop") %>%
			pull(mean_AF_log2OE)
}
## Are any of the groomig partners real?
has_real_partners <- function(df) {
	df %>%
		summarise(nr_partners = sum(grooming_partner != 'XXX') > 0 ) %>%
		pull()
}

## extract the real dyads
check_nr_dyads <- function(df) {
    df %>%
        filter(focal_check) %>%
        filter(str_detect(paternal_groom, "paternal")) %>%
        filter(i_total > 1/365 & sname != 'XXX') %>%
        nrow()
}

## is the dad the top dyadic partner
get_dad_top <- function(df, variable) {

	if(nrow(df) == 0) {

		FALSE
	} else {
		df %>%
			arrange({{variable}}, paternal_groom) %>%
			slice(1) %>%
			mutate(dad_top = paternal_groom == 'paternal') %>%
			select(dad_top) %>%
			pull()
	}
}

## turns the regression values calculated as part of the dyadic strengths values into zscores (per dryad).
make_zscored_values <- function(df) {
    df %>%
        filter(res_i_adj > -9999) %>%
        group_by(dyad_type) %>%
        mutate_at(vars(starts_with('res_i')),
                  .funs = list(zscore_i_adj = ~zscore(.))) %>%
        ungroup() %>%
        select(focal_check, tmp, sname, sname_sex, sname_sex_class, sname_grp,
               partner, partner_sex, partner_sex_class,  partner_grp,
               dyad, OE, log2OE, partner_mom, partner_dad, dyad_type,
               res_i_adj, zscore_i_adj)
}


## how long did the focal overlap with the partner
## In most cases how long did an juvenile female overlap with a male
## This can be done per group or overall
## It also deals with the fact that males can leave the group/population for a while and then return. Here we decided that if a male left for more then 30 days he had left and that overlap had ended.

get_overlap <- function(focal, focal_partner, start, end, per_group = FALSE) {
	start = ymd(start)
	end = ymd(end)

	df <- members_l %>%
		filter(sname == focal) %>%
		arrange(date) %>%
		filter(between(date, start, end)) %>%
		select(date, focal = sname, focal_grp = grp) %>%
		full_join(members_l %>%
								filter(sname == focal_partner) %>%
								filter(between(date, start, end)) %>%
								select(date, partner = sname, partner_grp = grp),
							by = "date")  %>%
		mutate(same_grp = focal_grp == partner_grp) %>%
		filter(same_grp)

	if(nrow(df) != 0) {
	    first_breakup <- df %>%
	        filter(same_grp == TRUE) %>%
		    mutate(previous_date = lag(date, default = first(date)),
		           nr_days_since_last = as.numeric(date - previous_date)) %>%
		    filter(nr_days_since_last > 30) %>%
		    slice(1) %>%
		    select(previous_date) %>%
		    pull()

		last_overlap_date <- df %>%
		    filter(same_grp) %>%
			tail(1) %>%
			select(date) %>%
			pull()

		if(length(first_breakup) == 0) first_breakup <- last_overlap_date
		if(is.na(first_breakup)) first_breakup <- last_overlap_date

		last_overlap_grp <- df %>%
		    filter(date == first_breakup) %>%
			select(focal_grp) %>%
			pull()

		total_number_overlap_days <- df %>%
		    filter(same_grp) %>%
			nrow()

		first_number_overlap_days <- df %>%
		    filter(same_grp & date < first_breakup) %>%
			nrow()

		df %>%
		    summarise(
		        partner_overlap_years = sum(same_grp, na.rm = TRUE)/365.25,
		        focal_days_known_grp = sum(!is.na(focal_grp))/365.25,
		        partner_overlap_proportional =
		            partner_overlap_years/partner_overlap_years,
		        last_overlap_date = last_overlap_date,
				first_breakup = first_breakup,
				last_overlap_grp = last_overlap_grp,
				total_number_overlap_days = total_number_overlap_days,
				first_number_overlap_days = first_number_overlap_days,
				days_since_birth = as.numeric(last_overlap_date - start))
		} else {
		    zdate_grp = parents_l %>%
		        filter(focal == kid) %>%
				select(dadgrp) %>%
				pull()

		    df <- members_l %>%
		        filter(sname == focal_partner &
		               grp == zdate_grp &
					   date < start) %>%
		        mutate(focal_grp = grp) %>%
				mutate(same_grp = TRUE) ## not really, but needed for later steps

			last_overlap_date <- df %>%
				filter(same_grp) %>%
				tail(1) %>%
				select(date) %>%
				pull()

			last_overlap_grp <- df %>%
				filter(same_grp) %>%
				tail(1) %>%
				select(focal_grp) %>%
				pull()

			number_overlap_days <- df %>%
				filter(same_grp) %>%
				nrow()

			df %>%
				summarise(
				    partner_overlap_years = 0,
					focal_days_known_grp = NA,
					partner_overlap_proportional = 0,
					last_overlap_date = last_overlap_date,
					first_breakup = last_overlap_date,
					last_overlap_grp = last_overlap_grp)
		}
}

## This function was used in some initial analysis and data exploration.
## The code above deals better with males leaving

get_overlap_per_group <- function(focal, focal_partner, start, end) {
	start = ymd(start)
	end = ymd(end)

	df <- members_l %>%
		filter(sname == focal) %>%
		arrange(date) %>%
		filter(between(date, start, end)) %>%
		select(date, focal = sname, focal_grp = grp) %>%
		full_join(members_l %>%
								filter(sname == focal_partner) %>%
								filter(between(date, start, end)) %>%
								select(date, partner = sname, partner_grp = grp),
							by = "date")  %>%
		mutate(same_grp = focal_grp == partner_grp)

	df %>%
		group_by(focal_grp) %>%
		summarise(
		    partner_overlap_years = sum(same_grp, na.rm = TRUE)/365.25,
		    focal_years_known_grp = sum(!is.na(focal_grp))/365.25,
		    partner_overlap_proportional = partner_overlap_years/focal_years_known_grp,
							last_overlap_day = max(date))
}

## This plot helped to sort some issues
check_overlap_plot <- function(focal, focal_partner, start, end) {
	start = ymd(start)
	end = ymd(end)

	df <- members_partner_overlap_l %>%
		filter(sname == focal) %>%
		arrange(date) %>%
		filter(between(date, start, end)) %>%
		select(date, focal = sname, focal_grpofresidency = grpofresidency) %>%
		full_join(members_partner_overlap_l %>%
								filter(sname == focal_partner) %>%
								filter(between(date, start, end)) %>%
								select(date, partner = sname,
								       partner_grp = grp,
								       partner_grpofresidency = grpofresidency),
								by = "date")  %>%
	    mutate(same_grp = focal_grpofresidency == partner_grpofresidency |
					 	focal_grpofresidency == partner_grp)

	df %>%
	    select(date, same_grp, focal_grpofresidency, dad_grpofresidency) %>%
	    pivot_longer(names_to = "sname", values_to = "grp",
	                 focal_grpofresidency:dad_grpofresidency) %>%
		inner_join(select(groups_l, grp = gid, social_group = name), by = "grp") %>%
		mutate(sname = if_else(str_detect(sname, 'focal'), 'Juvenile', 'Dad')) %>%
		mutate(age = as.numeric(date - focal_birth)/365.25) %>%
		ggplot(aes(x = age, y = sname, color = social_group)) +
		geom_point() +
		geom_vline(xintercept = df1$first_time_not_same_grp, color = 'red', size = 2) +
		geom_vline(xintercept = df1$last_time_not_same_grp, color = 'green', size = 2) +
		cowplot::theme_cowplot() +
		labs(x = "Age of juvenile in years",
				 y ="") +
		theme(legend.position = "bottom")
}


## The next functions where used to get the data for the analysis related to who groomed and the grooming strengths.

## get list of all males and how many days they were present
get_all_males <- function(focal_grp, start, end) {
    start_date <- ymd(start)
    end_date <- ymd(end)

    members_AM_l %>%
        filter(date >= start_date
               & date <= end_date
               & grp == focal_grp) %>%
        group_by(sname) %>%
        summarise(nr_days = n(), .groups = 'drop')

}

## Get mean rank of the male
get_mean_rank <- function(AMales, focal_grp, start, end) {
    start_date <- ymd(start)
    end_date <- ymd(start)

    proportional_ranks_l %>%
        filter(sname == AMales) %>%
        filter(rnktype == 'ADM') %>%
        filter(rnkdate >= floor_date(ymd(start), "month") &
                   rnkdate <= floor_date(ymd(end), "month")) %>%
        select(ordrank, proprank) %>%
        summarise(mean_ordrank = mean(ordrank),
                  mean_proprank = round(mean(proprank), 2),
                  .groups = 'drop')
}


#Get number of fertile feamles (expressed im days)
get_d_days_rate <- function(AMales, focal_grp, start, end, per_group = TRUE) {
    start <- ymd(start)
    end <- ymd(end)

    if(per_group==TRUE) {

        members_AM_l %>%
            filter(sname == AMales
                   #			 & 	grp == focal_grp
                   & date >=start
                   & date <= end) %>%
            left_join(d_days_l, by = c("grp", "date")) %>%
            summarise(daily_d_days_rate = sum(nr_d_days, na.rm = TRUE)/n()) %>%
            pull()
    }else{
        members_AM_l %>%
            filter(sname == AMales
                   & date >=start
                   & date <= end) %>%
            left_join(d_days_l, by = c("grp", "date")) %>%
            summarise(daily_d_days_rate = sum(nr_d_days, na.rm = TRUE)/n()) %>%
            pull()

    }
}

## Get the number of kids a male has in a group
get_kids <- function(focal, focal_birth, focal_mom, AMales) {
    df <- parents_l %>%
        arrange(zdate) %>%
        select(kid, mom, dad) %>%
        inner_join(select(biograph_l, kid = sname, kid_birth = birth), by = "kid") %>%
        filter(mom == focal_mom) %>%
        filter(kid != focal) %>%
        mutate(age_diff = as.numeric(kid_birth - ymd(focal_birth)))

    previous_kid_with_mom_df <- df %>%
        filter(age_diff < 0)

    if(nrow(previous_kid_with_mom_df) >= 1) {
        previous_kid_with_mom_df <- previous_kid_with_mom_df %>%
            filter(age_diff == max(age_diff))

        previous_pdads <- potential_dads_l %>%
            filter(kid %in% previous_kid_with_mom_df$kid) %>%
            select(pdad) %>%
            pull()

        previous_kid_with_mom <-
            tibble(previous_kid_with_mom = case_when(
                nrow(previous_kid_with_mom_df) == 0 ~ 0,
                !is.na(previous_kid_with_mom_df$dad) & previous_kid_with_mom_df$dad == AMales ~ 1,
                !is.na(previous_kid_with_mom_df$dad) & previous_kid_with_mom_df$dad != AMales ~ 0,
                is.na(previous_kid_with_mom_df$dad) & AMales %in% c(previous_pdads) ~ 0.5,
                is.na(previous_kid_with_mom_df$dad) & !(AMales %in% c(previous_pdads)) ~ 0,
                TRUE ~ 99))

        previous_kin = tibble(previous_sibling = previous_kid_with_mom_df$kid,
                              previous_kid_with_mom)
        } else {
            previous_kin  =tibble(previous_sibling = NA, previous_kid_with_mom = 0)
    }

    next_kid_with_mom_df <- df %>%
        filter(age_diff > 0)

    if(nrow(next_kid_with_mom_df) >= 1) {
        next_kid_with_mom_df <- next_kid_with_mom_df %>%
            filter(age_diff == min(age_diff))

        next_pdads <- potential_dads_l %>%
            filter(kid %in% next_kid_with_mom_df$kid) %>%
            select(pdad) %>%
            pull()

        next_kid_with_mom <-
            tibble(next_kid_with_mom = case_when(
                !is.na(next_kid_with_mom_df$dad) & next_kid_with_mom_df$dad == AMales ~ 1,
                !is.na(next_kid_with_mom_df$dad) & next_kid_with_mom_df$dad != AMales ~ 0,
                is.na(next_kid_with_mom_df$dad) & AMales %in% c(next_pdads) ~ 0.5,
                is.na(next_kid_with_mom_df$dad) & !(AMales %in% c(next_pdads)) ~ 0,
                nrow(next_kid_with_mom_df) == 0 ~ 0,
                TRUE ~ 99))

        next_kin = tibble(next_sibling = next_kid_with_mom_df$kid, next_kid_with_mom)

        } else {

        next_kin  =tibble(next_sibling = NA, next_kid_with_mom = 0)

        }

    bind_cols(previous_kin,next_kin)
}

## Check for maternal sibbling in group
get_maternal_sibling <- function(focal_sname, focal_birth, focal_grp, focal_mom, period) {

	df <- parents_l %>%
		inner_join(select(biograph_l, kid = sname, kid_birth = birth), by = "kid") %>%
		filter(mom == focal_mom) %>%
		arrange(kid_birth)

	if(period == "previous") {
		df2 <- df %>%
			filter(kid_birth < focal_birth) %>%
			#filter(kid != focal_sname) %>%
			filter(kid_birth == max(kid_birth))
	} else {
		df2 <- df %>%
			filter(kid_birth >= focal_birth) %>%
			filter(kid != focal_sname) %>%
			filter(kid_birth == min(kid_birth))
	}

	df2 %>%
		inner_join(potential_dads_l %>%
							 	group_by(kid) %>%
							 	summarise(nr_potential_dads = n())
							 ##,.groups = 'drop'
							 , by = "kid") %>%
		select(offspring_sname = kid, nr_potential_dads)
}


## How many paternal offspring are present in the period of overlap.
## It is expressed as offspring years to deal with offspring maturing, leaving or dying.
get_offspring_years <- function(focal_sname, AMales, focal_grp, start, end, per_group = TRUE) {
	start <- lubridate::ymd(start)
	end <- lubridate::ymd(end)

	if(per_group == TRUE) {
		members_juveniles_l  %>%
			filter(dad == AMales &
					 	grp == focal_grp) %>%
			filter(date >= start & date <= end) %>%
			select(sname, date, birth) %>%
			filter(sname != focal_sname) %>%
			arrange(date) %>%
			nrow()/365.25
	}else{
		members_juveniles_l  %>%
			filter(dad == AMales) %>%
			filter(date >= start & date <= end) %>%
			select(sname, date, birth) %>%
			filter(sname != focal_sname) %>%
			arrange(date) %>%
			nrow()/365.25
	}
}


## How many potential dads are there
get_nr_potential_dads <- function(focal_sibling) {
    potential_dads_l %>%
        filter(kid == focal_sibling) %>%
        n_distinct()
}

## Get the number of potential dads
get_nr_pdads <- function(focal) {
    potential_dads_l %>%
        filter(kid == focal) %>%
        group_by(kid) %>%
        summarise(nr_potential_dads = n()
                  ,.groups = 'drop'
                  , by = "kid") %>%
        select(nr_potential_dads) %>%
        pull()
}

## how many males for who groomed analysis
get_male_present_count <-  function(focal_kid, focal_grp, start, end) {

	start_date <- ymd(start)
	end_date <- ymd(end)

	members_AM_l %>%
		filter(date >= start_date
					 & date <= end_date
					 & grp == focal_grp) %>%
		distinct(sname) %>%
		left_join(select(potential_dads_l, kid, sname = pdad) %>%
								filter(kid == focal_kid),
							by = 'sname') %>%
		summarise(nr_males = n(),
							pdad_present = sum(!is.na(kid)==TRUE))
}

## get snames of potential dads
get_pdad_sname <-  function(focal_kid, focal_grp, start, end) {

	start_date <- ymd(start)
	end_date <- ymd(end)

	members_AM_l %>%
		filter(date >= start_date
					 & date <= end_date
					 & grp == focal_grp) %>%
		distinct(sname) %>%
		inner_join(select(potential_dads_l, kid, sname = pdad) %>%
								filter(kid == focal_kid),
							by = 'sname') %>%
		pull(sname)
}

## how many adult males are present
get_AM_count <- function(focal, focal_grp, start, end, focal_dad) {
    start <- lubridate::ymd(start)
    end <- lubridate::ymd(end)
    members_l %>%
        inner_join(members_l %>%
                       filter(sname == focal
                              & grp == focal_grp
                              & date >= start
                              & date <= end
                       ) %>%
                       select(grp, date),
                   by = c("grp", "date")) %>%
        inner_join(select(rankdates_l, sname, ranked)) %>%
        filter(date > ranked) %>%
        group_by(date) %>%
        summarise(nr_adult_males = n()
                  , nr_adult_males_not_dad = sum(sname != focal_dad)
                  , dad_presence = sum(sname == focal_dad)
                  , .groups = 'drop') %>%
        summarise(nr_days = n()
                  , max_nr_adult_males = max(nr_adult_males)
                  , max_nr_adult_males_not_dad = max(nr_adult_males_not_dad)
                  , mean_nr_adult_males = mean(nr_adult_males)
                  , mean_nr_adult_males_not_dad = mean(nr_adult_males_not_dad)
                  , mean_dad_presence = mean(dad_presence)
                  , .groups = 'drop')
}

## how many grooing partners does a focal have
get_grooming_partners <- function(df) {
    df %>%
        filter(focal_check
               & sname != 'XXX'
               & str_detect(paternal_groom, "paternal")
        ) %>%
        filter(i_total > 0.00274) %>%
        group_by(paternal_groom) %>%
        summarize(N = n(),
                  .groups = 'drop') %>%
        pivot_wider(names_from = 'paternal_groom', values_from = 'N')
}

## gets some dad regarding fathers
get_some_dad_data <- function(focal, focal_dad, focal_birth) {
    members_l %>%
        filter(date >= focal_birth &
                   date < focal_birth + lubridate::years(4) &
                   sname == focal) %>%
        select(grp, date, sname) %>%
        inner_join(members_l %>%
                       filter(date >= focal_birth &
                                  date <= focal_birth + lubridate::years(4) &
                                  sname == focal_dad) %>%
                       select(grp, date, focal_dad = sname),
                   by = c("grp", "date")) %>%
        arrange(date) %>%
        filter(!is.na(focal_dad)) %>%
        tail(1) %>%
        inner_join(biograph_l %>%
                       filter(sname == focal_dad) %>%
                       select(focal_dad = sname, dad_statdate = statdate, dad_status = status),
                   by = "focal_dad") %>%
        mutate(time_to_statdate = as.numeric(dad_statdate - date)/365.25) %>%
        mutate(rnkdate = lubridate::floor_date(date, 'month')) %>%
        select(last_date = date, dad_statdate, dad_status, time_to_statdate, rnkdate)
}

## What is the rank of the dad in next group
get_next_dad_rank <- function(focal_father, last_rnkdate) {
    ranks_l %>%
        filter(rnktype == 'ADM') %>%
        filter(sname == focal_father) %>%
        filter(rnkdate > last_rnkdate) %>%
        arrange(rnkdate) %>%
        slice(1) %>%
        select(next_ranked_grp = grp, next_rnkdate = rnkdate, rank)

}

## what was the rank of the male when he left
get_male_last_rank <- function(AMales, last_overlap_grp, first_breakup){

    ranks_l %>%
        filter(rnktype == "ADM") %>%
        filter(sname == AMales) %>%
        filter(grp == last_overlap_grp) %>%
        filter(floor_date(ymd(first_breakup), "month") == rnkdate) %>%
        select(rank) %>%
        pull()
}

## These next functions are used for the bootstrapping to handle missing cases in regards to the paternity of the previous or next offspring.
get_random_draws <- function(nr_potential_dads, nr_males, pdad_present, pdads_sname, draws = nr_draws, period)  {


	paternity <- c(1, rep(0, times = nr_potential_dads -1))

	xx <-replicate(draws,sample(x = paternity, size = pdad_present,replace = FALSE))

	print("2")
	if(pdad_present > 1) {
		xx <- xx %>%
			as_tibble() %>%
			mutate(AMales = c(pdads_sname))
	} else {
		xx <- t(xx) %>%
			as_tibble() %>%
			mutate(AMales = c(pdads_sname))
	}

	names(xx)[1:draws] <- paste0(period, "_", seq(1:draws))
	return(xx)
} ## end random draw


assign_random_draws <- function(df, previous_set, next_set) {
	if(!is.null(previous_set)) {
		df <- df %>%
			left_join(previous_set, by = 'AMales')
	}

	if(!is.null(next_set)) {
		df <- df %>%
			left_join(next_set, by = 'AMales')
	}
	return(df)
}

assign_random_draws_overall <- function(df, previous_set, next_set) {
	if(!is.null(previous_set)) {
		df <- df %>%
			left_join(select(previous_set, dad = AMales, everything()), by = 'dad')
	}

	if(!is.null(next_set)) {
		df <- df %>%
			left_join(select(next_set, dad = AMales), by = 'dad')
	}
	return(df)
}


## The next set of functions were used to run models, created tables and figures.
## The made it easier to automate some of the analysis.
## There are often used in a purrr pipeline

## This function helps with making the names in tables etc more human readable.
fix_terms <-function(term) {
	term %>%
		str_replace("kid_age", "Age of juvenile") %>%
		str_replace("dad_age", "Age of dad") %>%
		str_replace("mom_age", "Maternal age") %>%
		str_replace("maternal SCI_M", "Maternal SCI_M") %>%
		str_replace("cumulative_adversity", "Cumulative early life adversity") %>%
		str_replace("maternal_SCI_F_binaryTRUE", "Socially isolated mom") %>%
		str_replace("maternal_SCI_F_binary", "Socially isolated mom") %>%
		str_replace("maternal_SCI_F", "Socially isolated mom") %>%
		str_replace("maternal_lossTRUE", "Maternal loss") %>%
		str_replace("maternal_loss", "Maternal loss") %>%
		str_replace("maternal_rank_binaryTRUE", "Low ranking mom") %>%
		str_replace("maternal_rank_binary", "Low ranking mom") %>%
		str_replace("maternal_rank", "Low ranking mom") %>%
		str_replace("density_binaryTRUE", "Born in large group") %>%
		str_replace("density_binary", "Born in large group") %>%
		str_replace("density", "Born in large group") %>%
		str_replace("siblingTRUE", "Has competing maternal sibling") %>%
		str_replace("droughtTRUE", "Born in drought") %>%
		str_replace("drought", "Born in drought") %>%
		str_replace("maternal_proprank", "Maternal rank (at birth)") %>%
		str_replace("paternal_sibling_years", "Number of paternal siblings") %>%
		str_replace("offspring_years", "Number of kids male has in group") %>%
		str_replace("full_sibling_years", "Number of full siblings") %>%
		str_replace("dad_rank", "Ordinal rank of dad") %>%
		str_replace("male_rank", "Ordinal rank of male") %>%
		str_replace("hybrid_score", "Hybrid score of of the male") %>%
		str_replace("anubis_admix", "Anubis admix score of the male") %>%
		str_replace("nr_potential_dads", "Number of potential dads") %>%

		str_replace("nr_d_days", "Daily rate of d days") %>%
		str_replace("previous_sibling", "Overlap with previous maternal sibling") %>%
		str_replace("next_sibling", "Overlap with next maternal sibling") %>%
		str_replace("focal_age", "Age class of the juvenile") %>%
		str_replace("is_dadTRUE", "Male is the dad") %>%
		str_replace("AMales_age", "Age of male") %>%
		str_replace("mean_ordrank", "Rank of the male") %>%
		str_replace("mean_proprank", "Proportional rank of the male") %>%
		str_replace("previous_kid_with_mom", "Mother had a kid with male") %>%
		str_replace("next_kid_with_mom", "Mother will have a kid with male") %>%
		str_replace("estrous_presence", "Potential dad") %>%
		str_replace("proportion_consort_time", "Proportion of consort time") %>%
		str_replace("daily_d_days_rate", "Daily rate of d days") %>%
		str_replace("nr_pdads", "Number of potential dads at conception") %>%str_replace("jDSI_paternal", "Bond strength with father") %>%
		str_replace("jDSI_maternal", "Bond strength with mother") %>%
		str_replace("jDSI_Mtop", "Top bond strength with any adult male (including father)") %>%
		str_replace("jDSI_Mde_top", "Top bond strength with any adult male (excluding father)") %>%
		str_replace("dad_overlap_l", "Years of co-residency with father")  %>%
		str_replace("dad_overlap_years", "Years of co-residency with father")  %>%
		str_replace("cumulative_adversity:jDSI_paternal","Cumulative early life adversity X bond strength with father") %>%
	  str_replace("cumulative_adversity:dad_overlap_years","Cumulative early life adversity X years of co-residency with father") %>%
		str_replace("observer_effort", "Observer effort") %>%
		str_replace("poly\\(Age of male, 2\\)1", "First polynomial of male age") %>%
		str_replace("poly\\(Age of male, 2\\)2", "Second polynomial of male age") %>%
		str_replace('estrous_c', "Male had consort with mother during fertile period") %>%
		str_replace('competing_sibling', "Competing sibling") %>%
		# str_replace('estrous_me', "The number of mounts and ejaculation interactions observerd between the mother and the male during the fertile period ")
	str_replace('estrous_me', "The number of mounts")
}


## This function does a survival
get_coxph_model <- function(response) {
	coxph(data=xdata_females_with_social, as.formula(paste0("Surv(statage, adult_survival_status) ~ ",																													response)))
}

## Get the zph value of a survival model
get_coxzph <- function(x) {
	tail(cox.zph(x)$table, 1)[,3]
}

## get the confidence intervals of a model

get_confint <- function(model) {
	model %>%
		confint() %>%
		as_tibble() %>%
		setNames(c("conf.low", "conf.high"))
}


## I crated several different functions to create tables.
## These were used for differt models
## Many of these were not used for the final tables in the paper, but rather for the Rmarkdown version of the paper.

make_flextable <-function(df) {

	t1 <- df %>%
		mutate(HR = round(exp(estimate), 3),
					 conf.low = round(exp(conf.low), 3),
					 conf.high = round(exp(conf.high), 3),
					 sig_stars = case_when(p.value < 0.001  ~ "***",
					 											p.value < 0.01  ~ "***",
					 											p.value < 0.05  ~ "*",
					 											p.value < 0.01  ~ ".",
					 											TRUE ~ ""),
					 p.value = round(p.value, 3),
					 #sig_stars = paste("^", "**", sig_stars,"**", "^", sep =""),
					 HR = paste(HR, sig_stars, "\n(", conf.low, "-", conf.high,")", sep = "")) %>%
		select(formula, term, HR, AICc, model_check) %>%
		#arrange(AICc) %>%
		mutate(min_AICc = min(AICc)) %>%
		mutate(dAICc = AICc - min_AICc) %>%
		select(formula, term, HR, AICc, dAICc, model_check) %>%
		#mutate(term = str_remove(term, "TRUE")) %>%
		mutate(term = map_chr(.x = term, .f = fix_terms)) %>%
		mutate(formula =
					 	str_replace(formula, "maternal_loss", "maternal loss"),
					 formula =
					 	str_replace(formula, "jDSI_M","xxx"),
					 formula =
					 	str_replace(formula, "xxxtop",
					 							"strength of strongest juvenile dyadic strength with top male
                    (including dad; top male jDSI)"),
					 formula =
					 	str_replace(formula, "xxxde_top",
					 							"strength of strongest juvenile dyadic strength with top male
                    (excluding dad; top non-dad jDSI)"),
					 formula =
					 	str_replace(formula, "xxxde","juvenile dyadic strength with top3 males (excluding dad; jDSI_M)"),
					 formula =
					 	str_replace(formula, "xxx","juvenile dyadic strength with top3 males (jDSI_M)"),

					 formula =
					 	str_replace(formula, "cumulative_adversity","cumulative adversity"),
					 formula =
					 	str_replace(formula, "jDSI_paternal","juvenile dyadic strength with dad (jDSI with dad)"),
					 formula =
					 	str_replace(formula, "jSCI_M","xxx"),
					 formula =
					 	str_replace(formula, "xxxde","juvenile social connectedness to males excluding dad (jSCI_Mde)"),
					 formula =
					 	str_replace(formula, "xxx","juvenile social connectedness to males (jSCI_M)"),
					 formula =
					 	str_replace(formula, "dad_overlap_l","proportion of juvenile period shared with dad in same grp"),
					 formula =
					 	str_replace(formula, "dad_overlap","proportions of juvenile period shared with dad in same grp"),
					 formula = paste0('Adult female survival - ', formula)) %>%
		group_by(term) %>%
		pivot_wider(names_from = term, values_from = HR, values_fill = NA_character_) %>%
		select(!(AICc:model_check), AICc, dAICc) %>%
		mutate(formula = paste0("", formula)) %>%
		rename(model = formula) %>%
		select(model, contains(c("Cum", "jDSI with dad")), everything())

	col.index <- matrix(ncol =2, nrow = 0)
	#names(col.index) <- c("row", "col")
	for(ii in 1:nrow(t1)) {
		cols <- which((str_detect(string = t1[ii, ], pattern = "\\*", negate = FALSE)))
		temp <- as.matrix(tibble(row = rep(ii, times = length(cols)),
														 col = cols))
		col.index <- rbind(col.index, temp)
	}

	#marking = LETTERS[dim(t1)[1] : 1]
	marking = LETTERS[1:dim(t1)[1]]
	models<- t1$model

	ft <- t1 %>%
		mutate(model = marking) %>%
		flextable() %>%
		flextable::compose(part = "header", j = "dAICc",
											 value = as_paragraph("\U0394", "AICc"))

	for(i in 1:nrow(col.index)) {
		ft <- ft %>%
			bg(i = col.index[i, 1], j = col.index[i, 2], bg = 'yellow')
	}

	ft <- ft %>%
		add_header_row(values = c("", "Hazard ratio (95% CI)", "Model parameters"),
									 colwidths = c(1, ncol(t1) - 3, 2),
									 top = TRUE) %>%
		colformat_double(j = (ncol(t1)- 3):ncol(t1), digits = 3) %>%
		#align(align = "right", part = "footer") %>%
		theme_vanilla()
	#
	ft <- ft %>%
		colformat_double(digits = 2) %>%
		align(align = "left", part = "all") %>%
		align(align = "center", part = "header") %>%
		align(j = nrow(t1), align = "center") %>%
		#align(j = nrow(t1) + 1, align = "center") %>%
		fontsize(size = 7, part = "all") #%>%

	ft %>%
		fontsize(size = 7, part = "all")

	# 	width(j = 1, width = .1) %>%
	#   width(j = (ncol(t1)- 2):ncol(t1), width = .4) %>%
	# 	width(j = (2:(ncol(t1)- 3)), width = 1)
}


## To make tables of the bootstrap models
make_model_flextable <- function(model, dataset, explain_text = NA_character_) {
	t1 <- read.table(paste0('./data/bootstrap_results/', model, "_results.txt")) %>%
		as_tibble()

	if(ncol(t1) == 9) names(t1) <- c("row", "effect", "group",  "term", "estimate" ,"std.error", "statistic", "p.value", "VIF", "model_run")
	if(ncol(t1) == 10) names(t1) <- c("row", "effect", "group",  "term", "estimate" ,"std.error", "statistic", "p.value", "VIF", "model_run")

	if(ncol(t1) == 11) names(t1) <- c("row", "effect", "group",  "term", "estimate" ,"std.error", "statistic", "df", "p.value", "VIF", "model_run")


		t1 <- t1 %>%
		filter(effect == 'fixed') %>%
		select(-row, -effect, -group, -all_of(intersect(names(t1), "model_run"))) %>%
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
		width(j = 1, width = 2, unit = "in") %>%
		width(j = 2:5, width = .5, unit = "in") %>%
		width(j = 6, width =3, unit = "in") %>%
		height(height = 2,part = "header") %>%
		valign(valign = "bottom", part = "header") %>%
		valign(valign = "top", part = "body")
}



