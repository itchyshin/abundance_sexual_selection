library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(sf)
library(rnaturalearth)
data.1 <- readRDS("Data/grid_ebird_by_month_joined_with_mcqueen.RDS")
grids <- st_read("Data/grid_5_degree_with_clim.geojson")

##function

migration_month_function <- function(migration_value){
  
  dat.migration <- data.1 %>%
    dplyr::filter(migration == migration_value)
  
  month_function <- function(MONTH_value){
    
    dat.month <- dat.migration %>%
      dplyr::filter(MONTH == MONTH_value)
    
    model.month <- lmer(log10(mean_abund) ~ log_mass*scale(sexual_dichromatism) + (1 + scale(sexual_dichromatism)|grid_id), data = dat.month)
    
    ranef_df1 <- as.data.frame(ranef(model.month)$grid_id) %>% 
      rownames_to_column(var="grid_id") %>%
      mutate(intercept = `(Intercept)` + fixef(model.month)[1], slope = `scale(sexual_dichromatism)` + fixef(model.month)[3]) %>% 
      mutate(Month = MONTH_value) %>%
      mutate(Migration = migration_value)
  }
  
  modelling_results_month <- bind_rows(lapply(unique(dat.migration$MONTH), month_function))
  
}

modelling_results_migration_month <- bind_rows(lapply(unique(data.1$migration), migration_month_function))

##join with grids

data.model <- grids %>%
  rename(grid_id=ID) %>%
  mutate(grid_id=as.character(as.integer(grid_id))) %>%
  rename(temperature = mat) %>%
  rename(longitude = centroid_lat) %>%
  rename(latitude = centroid_lng) %>%
  rename(precipitation = prec_ann) %>%
  left_join(modelling_results_migration_month, by = 'grid_id')
  
#model temperature

model.temp <- lm(slope ~ temperature*Migration, data = data.model)
hist(model.temp$residuals)
plot(model.temp, which = 2)
summary(model.temp)

ggplot(data.model, aes(temperature, slope, colour = Migration, fill = Migration)) +
  geom_point() +
  geom_smooth(method= "lm" , se=FALSE, col="red") +
  facet_wrap(~ Month, ncol = 4)
  
#model latitude

model.latitude <- lm(slope ~ latitude*Migration, data = data.model)
hist(model.latitude$residuals)
plot(model.latitude, which = 2)
summary(model.latitude)

ggplot(data.model, aes(latitude, slope, colour = Migration, fill = Migration)) +
  geom_point() +
  geom_smooth(method= "lm" , se=FALSE, col="red") +
  facet_wrap(~ Month, ncol = 4)
#model precepitation

model.precipitation <- lm(slope ~ precipitation*Migration, data = data.model)
hist(model.precipitation$residuals)
plot(model.precipitation, which = 2)
summary(model.precipitation)

ggplot(data.model, aes(precipitation, slope, colour = Migration, fill = Migration)) +
  geom_point() +
  geom_smooth(method= "lm" , se=FALSE, col="red") +
  facet_wrap(~ Month, ncol = 4)
















