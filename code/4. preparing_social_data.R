## This code is part of the analysis for the Early-life paternal relationships predict adult female survival in wild baboons paper.

## The code was created by David Jansen
## Archie Lab; University of Notre Dame
## david.awam.jansen@gmail.com

## The corresponding author of the paper is Elizabeth Archie (earchie@nd.edu).

## The code below was written to prep the data for the calculations of social indexes.
## To get access to the database contact the corresponding author.

list.of.packages <- list("foreach"
                         , "doSNOW"
                         , "parallel"
                         , "tidyverse"
                         , "lubridate"
                         , "dbplyr"
                         , "purrrlyr"
                         , "RPostgreSQL"
                         , "zoo"
                         , "ramboseli")


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(unlist(new.packages))
lapply(list.of.packages, require, character.only = T)

## preparing social dataset

babase <- DBI::dbConnect(
	RPostgreSQL::PostgreSQL(),
	host = "localhost",
	port = 22222,
	user = "jansen",
	dbname = "babase",
	password = "Bab00n3455")

source("./code/social_indexes/biographical-data.R")

## although some of the data was already downloaded earlier.
## See 1. get_data.R
## some of it is redone here.
## This is do make the code from the ramboseli package work

# Create a connection to the database
babase <- DBI::dbConnect(
    RPostgreSQL::PostgreSQL(),
    host = "localhost",
    port = 22222,
    user = "Username",  ## This is the babase username
    dbname = "babase",
    password = "Password") ## This is the babase password




# Get local copy of biograph table
biograph_l <- collect(tbl(babase, "biograph"))

# Make a members subset that excludes behavioral observation gaps
members_l <- subset_members(babase, .adults_only = FALSE)

# Subset other data sets used for sociality indices
focals_l <- subset_focals(babase, members_l)
females_l <- subset_females(members_l)

# Grooming
grooming_l <- subset_interactions(babase, members_l, my_acts = c("G"), .adults_only = FALSE)

# Make an individual-year-of-life data set for adults
iyol <- make_iyol(babase, members_l, focals_l, grooming_l, .adults_only = FALSE)

## Restrict to groups where the animal was present for at least 60 days
iyol_sub <- iyol %>%
	filter(days_present >= 60)

current_date <- toupper(format(today(), '%d%b%y'))
local_files <- ls(pattern = "_l")

save(list = c(local_files, "iyol", "iyol_sub"),
		 file = paste0("./data/data_set_for_social_index_analysis_", current_date, ".RData"))

## This data will be moved to a cluster to calculate the actual social indexes.




