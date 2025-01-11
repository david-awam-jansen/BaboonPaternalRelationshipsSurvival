Sys.setenv(TZ = 'UTC')
list.of.packages <- list("foreach", "doSNOW", "parallel", "tidyverse",
												 "lubridate", "dbplyr", "purrrlyr",# "RPostgreSQL",
												 "zoo")
new.packages <- list.of.packages[!(list.of.packages %in%
																	 	installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(unlist(new.packages))
lapply(list.of.packages,function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

latest_version_date <- read_lines("./data/latest_version_date.txt")
load(paste0("./data/data_set_for_social_index_analysis_", latest_version_date,".RData"))
parents_l <- read_csv("./data/parents_l.csv")
members_l <- read_csv("./data/members_l.csv")
maturedates_l <- read_csv("./data/maturedates_l.csv")
rankdates_l <- read_csv("./data/rankdates_l.csv")
dsi_data <- read_rds('./data/DSI_zscored_values_21JUN22.rds')
iyol_sub <- read_csv("./data/iyol_sub.csv")
ea_dataset_less <- read_csv("./data/ea_dataset_less.csv")

source('./code/social_indexes/sociality-indices.R')

iyol_sub <- iyol_sub %>%
	inner_join(select(ea_dataset_less, sname)) %>%
	filter(days_present > 59) %>%
	filter(age_start_yrs < 4) %>%
	arrange(sname, start, grp) %>%
	mutate(row_nr = row_number())

adult_members <- members_l %>%
	select(grp, date, sname, sex) %>%
	filter(grp < 3) %>%
	left_join(select(maturedates_l, sname, matured), by = "sname") %>%
	left_join(select(rankdates_l, sname, ranked), by = "sname") %>%
	mutate(is_adult = case_when(sex == "F" & matured <= date ~ TRUE,
															sex == "M" & ranked <= date ~ TRUE)) %>%
	filter(is_adult == TRUE)

dsi_data <- read_rds(paste0("./data/dsi_data_", latest_version_date,".RDS"))

universal_paternal_values <-  dsi_data %>%
	rename(focal = sname) %>%
	unnest(cols = c(subset)) %>%
	unnest(cols = c(focal_data)) %>%
	mutate(group = dyad_type) %>%
	filter(sname != 'XXX') %>%
	group_by(group) %>%
	nest()  %>%
	mutate(regression = purrr::map(data,get_sub_regression)) %>%
	mutate(model_results = regression %>%  map(.,.f = broom::tidy)) %>%
	unnest(model_results) %>%
	select(dyad_type = group,  term, estimate) %>%
	pivot_wider(names_from = term, values_from = estimate)

  no_juvenile_focals <- read.table('./data/dyadic_index_messages_21JUN22.txt',
																 col.names=c("row_nr", "message")) %>%
  	filter(message == "No juvenile dyads present") %>%
  	pull(row_nr)

  my_subset_files <- tibble(file = list.files(path = "./data/my_subsets/", pattern = "my_subset_\\d+"))  %>%
  	separate(file, sep = "_", into = c(NA, NA, "row_nr", NA), remove = FALSE) %>%
  	filter(!(row_nr %in% no_juvenile_focals))

    my_subset_file <-  my_subset_file[1:30, ]

  cl <- makeCluster(4, outfile = "")
    registerDoSNOW(cl)
    pb <- txtProgressBar(min = 0, max = nrow(my_subset_files), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    clusterExport(cl, list("add_variables"
  											 , "get_functions"
  											 , "universal_paternal_values"
  											 , "iyol_sub"));
    opts <- list(progress = progress)
    subset <- foreach(i = 1:nrow(my_subset_files), .options.snow = opts,
    									.packages = c('tidyverse')) %dopar% {
    										add_variables(my_subset_files[i, ]$file)
    									}
    close(pb)
    stopCluster(cl)

dsi_paternal <- add_column(my_subset_files, subset) %>%
	unnest(cols = c(subset))

write_rds(dsi_paternal, paste0("./data/dsi_paternal_", latest_version_date, ".rds"))

