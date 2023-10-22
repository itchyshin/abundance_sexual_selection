# packages
library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(tidyr)
library(stringr)
library(lme4)
library(tibble)


# read in data for relative abundance in grids
rel_abund_dat <- readRDS("Data/relative_grids_summary.RDS") %>%
  rename(TipLabel=SCIENTIFIC_NAME) %>%
  mutate(TipLabel=gsub(" ", "_", TipLabel))

# read in color dat
dale <- read_csv(file="Data/Dale.csv")
mcqueen <- read_csv(file="Data/Seasonal_Plumage.csv")


# join ebird with mcqueen dichromatism
ebird_mcqueen <- rel_abund_dat %>%
  inner_join(., mcqueen, by="TipLabel") %>%
  mutate(log.rel_abund=log(rel_abund),
         log.rel_freq=log(rel_frequency))


model <- lmer(log.rel_abund ~ log_mass*scale(sexual_dichromatism) + (1+scale(sexual_dichromatism)|grid_id), data = ebird_mcqueen)
summary(model)
hist(ranef(model)$grid_id[,2])

model2 <- lmer(log.rel_freq ~ log_mass*scale(sexual_dichromatism) + (1+scale(sexual_dichromatism)|grid_id), data = ebird_mcqueen)
summary(model2)
hist(ranef(model2)$grid_id[,2])


ranef_df1 <- as.data.frame(ranef(model)$grid_id) %>%
  rownames_to_column(var="grid_id")

ranef_df2 <- as.data.frame(ranef(model2)$grid_id) %>%
  rownames_to_column(var="grid_id")


# read in and download spatial data
# download some natural earth data
countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()

graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()

# read in underlying map of grids for 5 degree grid cells
grid_5_degree <- st_read("Data/grid_5_degree.geojson")

grid_response_summary <- grid_5_degree %>%
  rename(grid_id=ID) %>%
  mutate(grid_id=as.character(as.integer(grid_id))) %>%
  left_join(., ranef_df1, by="grid_id") %>%
  dplyr::select(-`(Intercept)`) %>%
  rename(rel_abund_response=`scale(sexual_dichromatism)`) %>%
  left_join(., ranef_df2, by="grid_id") %>%
  dplyr::select(-`(Intercept)`) %>%
  rename(rel_freq_response=`scale(sexual_dichromatism)`)

ggplot()+
  geom_sf(data = bb, col = "grey20", fill = "transparent")+
  geom_sf(data= grid_response_summary, aes(fill=rel_abund_response))+
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4)+
  theme_minimal()+
  scale_fill_gradient2()+
  theme(axis.text = element_blank())+
  theme(legend.position="bottom")

ggplot()+
  geom_sf(data = bb, col = "grey20", fill = "transparent")+
  geom_sf(data= grid_response_summary, aes(fill=rel_freq_response))+
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4)+
  theme_minimal()+
  scale_fill_gradient2()+
  theme(axis.text = element_blank())+
  theme(legend.position="bottom")

ebird_mcqueen <- ebird_mcqueen %>%
  mutate(migration=as.character(as.numeric(migration)))

mig_model <- lmer(log.rel_abund ~ log_mass*scale(sexual_dichromatism) + (1+scale(sexual_dichromatism)|grid_id) + (1+scale(sexual_dichromatism)|migration), data = ebird_mcqueen)
summary(mig_model)
hist(ranef(model)$grid_id[,2])
