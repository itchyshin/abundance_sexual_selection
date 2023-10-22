# preparation monthly

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

#i read in the month files manually (probably and easier way to do it that i don't know)
#creating and saving a data frame of the joined months
#this is the data frame that is too large to upload (i just realised)
january <- readRDS("Data/grid_densities_per_month/january_densities.RDS") 
january$MONTH <- "january"
february <- readRDS("Data/grid_densities_per_month/february_densities.RDS") 
february$MONTH <- "february"
march <- readRDS("Data/grid_densities_per_month/march_densities.RDS")
march$MONTH <- "march"
april <- readRDS("Data/grid_densities_per_month/april_densities.RDS") 
april$MONTH <- "april"
may <- readRDS("Data/grid_densities_per_month/may_densities.RDS") 
may$MONTH <- "may"
june <- readRDS("Data/grid_densities_per_month/june_densities.RDS")
june$MONTH <- "june"
july <- readRDS("Data/grid_densities_per_month/july_densities.RDS") 
july$MONTH <- "july"
august <- readRDS("Data/grid_densities_per_month/august_densities.RDS") 
august$MONTH <- "august"
september <- readRDS("Data/grid_densities_per_month/september_densities.RDS") 
september$MONTH <- "september"
october <- readRDS("Data/grid_densities_per_month/october_densities.RDS") 
october$MONTH <- "october"
november <- readRDS("Data/grid_densities_per_month/november_densities.RDS") 
november$MONTH <- "november"
december <- readRDS("Data/grid_densities_per_month/december_densities.RDS") 
december$MONTH <- "december"


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
  group_by(ebird_COMMON_NAME, MONTH, grid_id) %>%
  slice(1) %>%
  filter(complete.cases(sexual_dichromatism)) %>%
  mutate(., sexual_dichromatism = log10(sexual_dichromatism))

month.data$migration <- str_replace_all(month.data$migration, "0", "resident")
month.data$migration <- str_replace_all(month.data$migration, "1", "migrant")
month.data$migration <- str_replace_all(month.data$migration, "2", "migrant")


# get other bits of info
dat1 <- readRDS("Data/grid_ebird_by_month_joined_with_mcqueen.RDS")
dat <- distinct(dat1[,c(1, 2, 7, 5)])

# combining area_data and density data
month.data %>% left_join(., dat, by = c("ebird_COMMON_NAME", "MONTH", "grid_id")) -> dat_final

# final data
saveRDS(dat_final, "Data/grid_ebird_by_month_joined_with_McQueen_imputed.RDS")


