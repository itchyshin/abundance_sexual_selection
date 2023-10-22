# get gridded summary relative abundance
# data from bigquery

## packages
library(readr)
library(DBI)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)


# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa')

# create locality_id and grid_id table
grids <- tbl(con, 'grid_5_joined_with_sites')


  
# calculate the relative abundance for each species for each grid
# which is the sum of all abundances/the total time spent surveying in a grid
grid_minutes <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID, COMMON_NAME, 
                OBSERVATION_COUNT, OBSERVATION_DATE, SCIENTIFIC_NAME, DURATION_MINUTES) %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  group_by(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID) %>%
  summarize(checklist_minutes=sum(DURATION_MINUTES)) %>%
  left_join(., grids, by="LOCALITY_ID") %>%
  group_by(grid_id) %>%
  summarize(total_minutes=sum(checklist_minutes)) %>%
  collect(n=Inf)
  
grid_species_abund <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID, COMMON_NAME, 
                OBSERVATION_COUNT, OBSERVATION_DATE, SCIENTIFIC_NAME, DURATION_MINUTES) %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  mutate(OBSERVATION_COUNT = as.numeric(as.character(OBSERVATION_COUNT))) %>%
  group_by(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, LOCALITY_ID) %>%
  summarize(total_abund=sum(OBSERVATION_COUNT)) %>%
  left_join(., grids, by="LOCALITY_ID") %>%
  group_by(grid_id, COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(total_abund=sum(total_abund)) %>%
  collect(n=Inf)

grid_rel_abund_summary <- grid_species_abund %>%
  left_join(., grid_minutes, by="grid_id") %>%
  mutate(rel_abund=total_abund/total_minutes)

# now calculate the relative frequency for each species
# the number of occurrences/total checklists
grid_checklists <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID, OBSERVATION_DATE) %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  group_by(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID) %>%
  summarize(number_checklists=n()) %>%
  left_join(., grids, by="LOCALITY_ID") %>%
  group_by(grid_id) %>%
  summarize(total_checklists=sum(number_checklists)) %>%
  collect(n=Inf)

grid_species_occurences <- ebird %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID, COMMON_NAME, 
                OBSERVATION_COUNT, OBSERVATION_DATE, SCIENTIFIC_NAME, DURATION_MINUTES) %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  group_by(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, LOCALITY_ID) %>%
  summarize(total_occurrences=n()) %>%
  left_join(., grids, by="LOCALITY_ID") %>%
  group_by(grid_id, COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(total_occurrences=sum(total_occurrences)) %>%
  collect(n=Inf)

grid_rel_freq_summary <- grid_species_occurences %>%
  left_join(., grid_checklists, by="grid_id") %>%
  mutate(rel_frequency=(total_occurrences/total_checklists)*100) %>%
  dplyr::select(grid_id, rel_frequency, COMMON_NAME)

grid_rel_summary <- grid_rel_abund_summary %>%
  left_join(., grid_rel_freq_summary)

saveRDS(grid_rel_summary, "Data/relative_grids_summary.RDS")



