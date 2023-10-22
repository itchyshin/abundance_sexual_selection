## Script is to get some initial summaries
## of relative abundance from eBird
## In our other project we have shown that relative abunance from eBird
## is pretty robust measured different ways (i.e., predicted from a model compared with a simple mean of abundances)
## this script interacts with a BigQuery Database

## packages
library(readr)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(tidyr)


# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa')

# testing first
relative_abund_summary.test <- ebird %>%
  dplyr::select(COUNTRY, COMMON_NAME, SCIENTIFIC_NAME, 
                OBSERVATION_COUNT, CATEGORY) %>%
  dplyr::filter(CATEGORY %in% c("species", "issf", "domestic")) %>%
  dplyr::filter(COUNTRY == "Portugal") %>%
  dplyr::filter(OBSERVATION_COUNT != "X") %>%
  mutate(OBSERVATION_COUNT=as.numeric(as.character(OBSERVATION_COUNT))) %>%
  group_by(COUNTRY, COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarise(mean_abund=mean(OBSERVATION_COUNT)) %>%
  collect(n=Inf)

total_sampling_events.test <- ebird %>%
  dplyr::select(COUNTRY, OBSERVATION_COUNT, CATEGORY, SAMPLING_EVENT_IDENTIFIER) %>%
  dplyr::filter(CATEGORY %in% c("species", "issf", "domestic")) %>%
  dplyr::filter(COUNTRY == "Portugal") %>%
  dplyr::filter(OBSERVATION_COUNT != "X") %>%
  mutate(OBSERVATION_COUNT=as.numeric(as.character(OBSERVATION_COUNT))) %>%
  distinct() %>%
  group_by(COUNTRY) %>%
  summarise(N=n()) %>%
  collect(n=Inf)

# calculating the information
relative_abund_summary.total <- ebird %>%
  dplyr::select(STATE_CODE, COUNTRY, COMMON_NAME, SCIENTIFIC_NAME, 
                OBSERVATION_COUNT, CATEGORY) %>%
  dplyr::filter(CATEGORY %in% c("species", "issf", "domestic")) %>%
  dplyr::filter(OBSERVATION_COUNT != "X") %>%
  mutate(OBSERVATION_COUNT=as.numeric(as.character(OBSERVATION_COUNT))) %>%
  group_by(STATE_CODE, COUNTRY, COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarise(mean_abund=mean(OBSERVATION_COUNT),
            Number_obs=n()) %>%
  collect(n=Inf)


saveRDS(relative_abund_summary.total, "Data from eBird/ebird_relative_abundance_summary.RDS")


total_sampling_events.total <- ebird %>%
  dplyr::select(COUNTRY, STATE_CODE, OBSERVATION_COUNT, 
                CATEGORY, SAMPLING_EVENT_IDENTIFIER) %>%
  dplyr::filter(CATEGORY %in% c("species", "issf", "domestic")) %>%
  dplyr::filter(OBSERVATION_COUNT != "X") %>%
  mutate(OBSERVATION_COUNT=as.numeric(as.character(OBSERVATION_COUNT))) %>%
  distinct() %>%
  group_by(STATE_CODE) %>%
  summarise(N=n()) %>%
  collect(n=Inf)

saveRDS(total_sampling_events.total, "Data from eBird/ebird_state_code_total_checklists.RDS")







