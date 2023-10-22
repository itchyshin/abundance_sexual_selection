# a script to combine monthly summaries into one RDS for analyses


library(dplyr)
library(tidyr)
library(readr)

#i read in the month files manually (probably and easier way to do it that i don't know)
#creating and saving a data frame of the joined months
#this is the data frame that is too large to upload (i just realised)
january <- readRDS("Data/ebird_monthly_grid_summaries/january.RDS") %>%
  dplyr::filter(mean_abund>0)
february <- readRDS("Data/ebird_monthly_grid_summaries/february.RDS") %>%
  dplyr::filter(mean_abund>0)
march <- readRDS("Data/ebird_monthly_grid_summaries/march.RDS") %>%
  dplyr::filter(mean_abund>0)
april <- readRDS("Data/ebird_monthly_grid_summaries/april.RDS") %>%
  dplyr::filter(mean_abund>0)
may <- readRDS("Data/ebird_monthly_grid_summaries/may.RDS") %>%
  dplyr::filter(mean_abund>0)
june <- readRDS("Data/ebird_monthly_grid_summaries/june.RDS") %>%
  dplyr::filter(mean_abund>0)
july <- readRDS("Data/ebird_monthly_grid_summaries/july.RDS") %>%
  dplyr::filter(mean_abund>0)
august <- readRDS("Data/ebird_monthly_grid_summaries/august.RDS") %>%
  dplyr::filter(mean_abund>0)
september <- readRDS("Data/ebird_monthly_grid_summaries/september.RDS") %>%
  dplyr::filter(mean_abund>0)
october <- readRDS("Data/ebird_monthly_grid_summaries/october.RDS") %>%
  dplyr::filter(mean_abund>0)
november <- readRDS("Data/ebird_monthly_grid_summaries/november.RDS") %>%
  dplyr::filter(mean_abund>0)
december <- readRDS("Data/ebird_monthly_grid_summaries/december.RDS") %>%
  dplyr::filter(mean_abund>0)


month.joined <- january %>%
  bind_rows(., february) %>%
  bind_rows(., march) %>%
  bind_rows(., april) %>%
  bind_rows(., may) %>%
  bind_rows(., june) %>%
  bind_rows(., july) %>%
  bind_rows(., august) %>%
  bind_rows(., september) %>%
  bind_rows(., october) %>%
  bind_rows(., november) %>%
  bind_rows(., december)


McQueen <- read.csv(file="Data/Seasonal_Plumage.csv",header=TRUE)

clements <- read_csv("Data/clements_clean.csv")


month.data <- month.joined %>%
  ungroup() %>%
  left_join(., clements) %>%
  dplyr::filter(order=="Passeriformes") %>%
  left_join(McQueen, by = "TipLabel") %>%
  filter(complete.cases(sexual_dichromatism)) %>%
  mutate(., sexual_dichromatism = log10(sexual_dichromatism))

month.data$migration <- str_replace_all(month.data$migration, "0", "resident")
month.data$migration <- str_replace_all(month.data$migration, "1", "migrant")
month.data$migration <- str_replace_all(month.data$migration, "2", "migrant")

#final result

saveRDS(month.data, "Data/grid_ebird_by_month_joined_with_McQueen.RDS")
