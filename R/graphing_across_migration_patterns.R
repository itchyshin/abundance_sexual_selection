## This is a script to read in data ready to be modelled
## and then model it at different levels

# packages
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)


# read in data
data.1 <- readRDS("Data/grid_ebird_joined_with_mcqueen.RDS")


# modelling function
migration_function <- function(migration_value){
  
  dat <- data.1 %>%
    dplyr::filter(migration == migration_value)
  
  model.migration <- lmer(log10(mean_abund) ~ log_mass*scale(sexual_dichromatism) + (1 + scale(sexual_dichromatism)|grid_id), data = dat)
  
  # I recommend chaning the name not to use "(Intercept)" "scale(sexual_dichromatism)"
  ranef_df1 <- as.data.frame(ranef(model.migration)$grid_id) %>% 
    rownames_to_column(var="grid_id") %>%
    # fixef(model.migration)[1] == population intercept and poluation slope == fixef(model.migration)[3]
    mutate(intercept = `(Intercept)` + fixef(model.migration)[1], slope = `scale(sexual_dichromatism)` + fixef(model.migration)[3]) %>% 
    mutate(Migration=migration_value)
  
}


# now "apply" this function over each level of
# migration values
# and then save the results in one dataframe
modelling_results_migration <- bind_rows(lapply(unique(data.1$migration), migration_function))

# plotting graphs

library(sf)
library(rnaturalearth)

grids <- st_read("Data/grid_5_degree_with_clim.geojson")

countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()

graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()

modelling_results_migration_grid <- grids %>%
  rename(grid_id = ID) %>%
  mutate(grid_id = as.character(as.integer(grid_id))) %>%
  left_join(modelling_results_migration) %>%
  mutate(Migration = as.character(as.integer(Migration))) %>%
  dplyr::filter(complete.cases(Migration))

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= modelling_results_migration_grid, aes(fill = slope)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom") +
  facet_wrap(~ Migration, ncol=1)





