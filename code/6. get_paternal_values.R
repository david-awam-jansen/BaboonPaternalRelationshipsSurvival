## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).

## The code below was written to calculate the actual paternal DSI values
## It sources the file that has that has the actual code for the calculations.
## The calculation of the paternal DSI values was developed for this paper.

## To get access to the database contact the corresponding author.

## This code needs a lot of computer power and should be run on a cluster.
## It is possible to run a subset on personal device.


list.of.packages <- list("foreach"
                         , "doSNOW"
                         , "parallel"
                         , "tidyverse"
                         , "lubridate"
                         , "dbplyr"
                         , "purrrlyr"
                         ,  "zoo")

## install packages if needed and open libaries
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
}
lapply(packages, require, character.only = TRUE, warn.conflicts = TRUE,
       quietly = FALSE)

latest_version_date <- read_lines("./data/latest_version_date.txt")
load(paste0("./data/data_set_for_social_index_analysis_", latest_version_date,".RData"))

dsi_data <- read_rds('./data/DSI_zscored_values.rds')
iyol_sub <- read_csv("./data/iyol_sub.csv")
ea_dataset_less <- read_csv("./data/ea_dataset_less.csv")

## These should have been loaded with data_set_for_social_index_analysis
# parents_l <- read_csv("./data/parents_l.csv")
# members_l <- read_csv("./data/members_l.csv")
# maturedates_l <- read_csv("./data/maturedates_l.csv")
# rankdates_l <- read_csv("./data/rankdates_l.csv")


## Source the actual code for the calculations
## For more details on the actual steps of the calculation see the code
source('./code/social_indexes/sociality_indices.R')

iyol_sub <- iyol_sub %>%
	inner_join(select(ea_dataset_less, sname)) %>%
	filter(days_present > 59) %>%
	filter(age_start_yrs < 4) %>%
	arrange(sname, start, grp) %>%
	mutate(row_nr = row_number())

dsi_data <- read_rds(paste0("./data/dsi_data_", latest_version_date,".RDS"))

## During the initial writing of the dyadic code we discovered that there are cases were the observer effort correction didn't work as expected. We found some bias related to the fact that observer effort that are a by-product of either especially small groups, or groups who were observed few times in a given period. We can control for this using a “universal” slope. See methods in Rosenbaum et al 2020 for details (https://www.pnas.org/doi/full/10.1073/pnas.2004524117#sec-4).

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

## Detect how many cores
## If you run on personal computer specify the number of cores
ncores <- detectCores()

cl <- makeCluster(ncores, outfile = "out.log")
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