## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).

## The code below was written to calculate the actual DSI values
## It sources the file that has that has the actual code for the calculations.
## This code was based on social index code in the Ramboseli package.
## It was addapted a lot to enable the calculation of social index values for juveniles.

## To get access to the database contact the corresponding author.

## This code needs a lot of computer power and should be run on a cluster.
## It is possible to run a subset on personal device.

list.of.packages <- list("foreach", "doSNOW", "parallel", "tidyverse",
                         "lubridate", "dbplyr", "purrrlyr",# "RPostgreSQL",
                         "zoo")

## install packages if needed and open libaries
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
}
lapply(packages, require, character.only = TRUE, warn.conflicts = TRUE,
       quietly = FALSE)

latest_version_date <- read_lines("./data/latest_version_date.txt")
load(paste0("./data/data_set_for_social_index_analysis_", latest_version_date,".RData"))

## some data prep
adult_members <- members_l %>%
	select(grp, date, sname, sex) %>%
	filter(grp < 3) %>%
	left_join(select(maturedates_l, sname, matured), by = "sname") %>%
	left_join(select(rankdates_l, sname, ranked), by = "sname") %>%
	mutate(is_adult = case_when(sex == "F" & matured <= date ~ TRUE,
															sex == "M" & ranked <= date ~ TRUE)) %>%
	filter(is_adult == TRUE)

ea_dataset_less <- read_csv(file = paste("./data/", "ea_dataset_less_restricted_", latest_version_date, ".csv", sep=""))

iyol_sub <- iyol_sub %>%
	inner_join(select(ea_dataset_less, sname)) %>%
	filter(days_present > 59) %>%
	filter(age_start_yrs < 4) %>%
	arrange(sname, start, grp) %>%
	mutate(row_nr = row_number()) %>%
	filter(sex_class == "JUV")

## uncomment for testing
#iyol_sub <- iyol_sub %>% sample_n(96)

## This could help if there are errors
dyadic_index_messages <- tibble(row_nr = iyol_sub$row_nr, message = NA)

## Source the actual code for the calculations
source('./code/social_indexes/sociality_indices.R')

## Detect how many cores
## If you run on personal computer specify the number of cores
ncores <- detectCores()

dsi_data <- dyadic_index(my_iyol = iyol_sub, biograph_l = biograph_l,
		members_l = members_l, focals_l = focals_l, females_l = females_l,
		interactions_l = grooming_l,
		min_cores_days = 30,  ## individuals have to live together for at least 30 days
		parallel = TRUE, ncores = ncores,
		directional = FALSE)

dsi_data <- dsi_data %>%
	select(row_nr, everything())

write_rds(dsi_data, paste0("./data/dsi_data_", latest_version_date, ".rds"))

## This file can be useful to see why there are NA values
write.table(dyadic_index_messages, paste0("./data/dyadic_index_messages_",
                                          latest_version_date, ".txt"))

## During the initial writing of the dyadic code we discovered that there are cases were the observer effort correction didn't work as expected. We found some bias related to the fact that observer effort that are a by-product of either especially small groups, or groups who were observed few times in a given period. We can control for this using a “universal” slope. See methods in Rosenbaum et al 2020 for details (https://www.pnas.org/doi/full/10.1073/pnas.2004524117#sec-4).

## Get universal values
universal_values <-  dsi_data %>%
	rename(focal = sname) %>%
	unnest(cols = c(zscored_values)) %>%
	unnest(cols = c(focal_data)) %>%
	mutate(group = dyad_type) %>%
	filter(sname != 'XXX') %>%
	group_by(group) %>%
	nest()  %>%
	mutate(regression = purrr::map(data,get_sub_regression)) %>%
	mutate(model_results = regression %>%  map(.,.f = broom::tidy)) %>%
	unnest(cols = c(model_results)) %>%
	select(dyad_type = group,  term, estimate) %>%
	pivot_wider(names_from = term, values_from = estimate)

## deal with missing data
no_juvenile_focals <- read.table(paste0("./data/dyadic_index_messages_", latest_version_date, ".txt"),
col.names=c("row_nr", "message")) %>%
	filter(message != "Focal data should be available") %>%
	pull(row_nr)

my_subset_files <- tibble(file = list.files(path = "./data/my_subsets/", pattern = "my_subset_\\d+"))  %>%
	separate(file, sep = "_", into = c(NA, NA, "row_nr", NA), remove = FALSE) %>%
	filter(!(row_nr %in% no_juvenile_focals)) %>%
	mutate(row_nr = as.numeric(row_nr))

print("Do I get here")

## Link the subset (individaul year files to orgional dataset)
iyol_sub <-iyol_sub %>%
	left_join(my_subset_files)

## calcuilate the actual values
universal_zscored_values <- dyadic_index.universal(
    my_iyol = iyol_sub
    , biograph_l = biograph_l
    , members_l = members_l
    , focals_l = focals_l
    , females_l = females_l
    , interactions_l = grooming_l
    , min_cores_days = 30
    , parallel = TRUE, ncores = ncores)

print("can I get to this step")
write_rds(universal_zscored_values,
          paste0("./data/universal_zscored_values_", latest_version_date, ".rds"))

dsi_universal_data <- universal_zscored_values %>%
	select(row_nr, universal_zscored_values)

print("I am get close")
dsi_universal_files <- read.table(
    file = paste0("./data/dyadic_index_messages_", latest_version_date, ".txt"),
    col.names=c("row_nr", "message")) %>%
    left_join(select(iyol_sub, -file)) %>%
	left_join(tibble(file = list.files(path = "./data/my_subsets/",
	                                   pattern = "my_subset_\\d+"))  %>%
	              separate(file, sep = "_", into = c(NA, NA, "row_nr", NA),
	                       remove = FALSE) %>%
	              mutate(row_nr = as.numeric(row_nr)))

dsi_universal <- dsi_universal_files %>%
	left_join(dsi_universal_data) %>%
	unnest(cols = c(universal_zscored_values))

write_rds(x = dsi_universal,
          file = paste0("./data/dsi_universal_", latest_version_date, ".rds"))

print("I did it") ## it feels really good if you get here :-)