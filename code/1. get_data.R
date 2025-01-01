## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).

## The code below was written to download the data from babase.
## To get access to the database contact the corresponding author.

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

## Get current date to keep track of when the data was downloaded
current_date <- toupper(format(today(), '%d%b%y'))

# Create a connection to the database
babase <- DBI::dbConnect(
	RPostgreSQL::PostgreSQL(),
	host = "localhost",
	port = 22222,
	user = "Username",  ## This is the babase username
	dbname = "babase",
	password = "Password") ## This is the babase password

## Collect the following datasets.
## For more details on these datasets see https://papio.biology.duke.edu/pocket_reference.pdf

## The underscore l indicates that the data was actually downloaded locally.

actor_actees_l <- tbl(babase, 'actor_actees') %>% collect()
behave_gaps_l <- tbl(babase, 'behave_gaps') %>% collect()
biograph_l <- tbl(babase, 'biograph') %>% collect()
groups_l<- tbl(babase, 'groups') %>% collect()
## In the final analysis we did not use hybrid_scores
#hybridgene_scores_l <- tbl(babase, 'hybridgene_scores') %>% collect()
maturedates_l <- tbl(babase, 'maturedates') %>%  collect()
parents_l <- tbl(babase, 'parents') %>% collect()
potential_dads_l <- tbl(babase, 'potential_dads') %>% collect()
rankdates_l <- tbl(babase, 'rankdates') %>% collect()
ranks_l <- tbl(babase, 'ranks') %>% collect()

## To get a measure of the daily rate of fertile females we calculated a vairable called d_days.
d_days_l <- tbl(babase, 'members') %>%
	inner_join(tbl(babase, 'biograph') %>%
						 	filter(sex == 'F') %>%
						 	select(sname), by = "sname") %>%
	inner_join(tbl(babase, 'repstats') %>%
						 	select(sname, date, repstate = state)
						 , by = c("sname", "date")) %>%
	left_join(tbl(babase, 'cycstats'), by = c("sname", "date")) %>%
	left_join(tbl(babase, 'cycpoints_cycles'),
						by = c("sname", "date")) %>%
	filter(code == "D") %>%
	group_by(grp, date) %>%
	summarise(nr_d_days = n()) %>%
	collect()

## Due to a various reasons there are periods were quality of data collection was insufficient. The code below will deal with this by using behave_gaps and groups history.

groups_history <- tbl(babase, "groups_history") %>%
    select(grp = gid, permanent, impermanent, cease_to_exist, last_reg_census) %>%
    group_by(grp) %>%
    mutate(last_date = pmin(impermanent,
                            cease_to_exist,
                            last_reg_census,
                            "2020-07-01")) %>%
    ungroup()

groups_history_l <- groups_history %>%  collect()

## Not all behave gaps effected the data we needed to exclude and David and Beth selected which periods should be excluded.

selected_behave_gaps <- read_csv("./data/selected_behave_gaps.csv")

excluded_members <- selected_behave_gaps %>%
	filter(exclude_group_size == 'YES') %>%
	select(bgid) %>%
	pull()

excluded_ranks <- selected_behave_gaps %>%
    filter(exclude_rank == 'YES') %>%
    select(bgid) %>%
    pull()

## Get members data
members <- tbl(babase, 'members') %>%
	filter(grp < 3) %>%
	select(membid, grp, sname, date, grpofresidency)

members <- members %>%
	left_join(groups_history) %>%
	filter(date >= permanent & date <= last_date)

members_bg <- members %>%
	dplyr::left_join(tbl(babase, "behave_gaps"), by = "grp") %>%
	dplyr::filter(date >= gap_start & date <= gap_end) %>%
	select(membid, bgid) %>%
	filter(bgid %in% excluded_members)

members <- members %>%
	anti_join(members_bg, by = 'membid')

members_l <- members %>%
	collect()

## and create some subsets of the members datasets

members_AM_l <- members_l %>%  ## Only ranked males
	inner_join(rankdates_l) %>%
	filter(date > ranked) %>%
	select(sname, grp, date, ranked)

members_juveniles_l <- members_l %>%
	inner_join(select(biograph_l, sname, birth), by = 'sname') %>%
	mutate(age = as.numeric(date- birth)/365.25) %>%
	filter(age < 4) %>%
	inner_join(select(parents_l, sname = kid, dad))

members_partner_overlap_l <- tbl(babase, 'members') %>%
	collect() %>%
	arrange(sname, date)

## get rank data
ranks <- tbl(babase, 'ranks') %>%
	filter(grp < 3) %>%
	select(rnkid, grp, sname, rnkdate, rnktype, rank)

ranks <- ranks %>%
	inner_join(select(tbl(babase, 'biograph'), sname, sex)) %>%
	left_join(select(tbl(babase, 'maturedates'), sname, matured))

ranks_l <- ranks %>%
	inner_join(groups_history) %>%
	filter(rnkdate >= permanent & rnkdate <= last_date)

ranks_bg <- ranks_l %>%
	inner_join(behave_gaps_l) %>%
	select(rnkid, bgid, rnkdate, gap_start, gap_end) %>%
 	filter(rnkdate > gap_start & rnkdate < gap_end) %>%
	select(rnkid, bgid) %>%
	filter(bgid %in% excluded_ranks)

ranks_l <- ranks_l %>%
	anti_join(ranks_bg, by = 'rnkid')

proportional_ranks_l <- tbl(babase, "proportional_ranks") %>%
	collect()

proportional_ranks_l <-proportional_ranks_l %>%
	inner_join(select(ranks_l, rnkid))


## For the next 3 datasets I used excisting SQL queries.

rainfall_l <- tbl(babase,
                  sql("SELECT hydroyear(wrdaytime::date),
                       date_part('year', wrdaytime::date) AS year
                       , date_part('month', wrdaytime::date) AS month
                       , wrdaytime::date AS date, rain
                       FROM wreadings w
                       LEFT JOIN raingauges r ON r.wrid = w.wrid")) %>%
    collect()

pregdata_l <- tbl(babase,
                  sql("SELECT b.pid, kid AS sname, kid, mom, dad
                      , zdate AS conceive_date, b.birth, bstatus, sex
                      , statdate, status, b.dcause, d.nature, d.agent
                      , pregs.parity, cpt.date AS resume_date
                      FROM parents p
                      JOIN biograph b ON p.kid = b.sname
                      LEFT JOIN pregs ON b.pid = pregs.pid
                      LEFT JOIN cycpoints cpz ON cpz.cpid = pregs.conceive
                      LEFT JOIN cycpoints cpt ON cpt.cpid = pregs.resume
                      LEFT JOIN dcauses d ON d.dcause = b.dcause")) %>%
    collect() %>%
	mutate(birth = ymd(birth))

maternal_loss_l <- tbl(babase,
                       sql("SELECT b.sname, p.mom
                       , (b2.statdate - b.birth)/365.25 AS maternal_loss_age
                       , (b2.statdate - b.birth)/365.25 <= 4 AS maternal_loss
                       FROM biograph b
                       JOIN parents p ON p.kid = b.sname
                       JOIN biograph b2 ON b2.sname = p.mom")) %>%
    collect()

current_date <- toupper(format(today(), '%d%b%y'))


## save only the data sets that was locally downloaded.
save(list = ls(pattern = "_l"),
		 file = paste0("./data/data_set_for_dad_analysis_",
		               current_date, ".RData"))

write_lines(current_date, "./data/latest_version_date.txt")