library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

data.1 <- readRDS("Data/grid_ebird_by_month_joined_with_mcqueen.RDS")

month_function <- function(MONTH_value){
  
  dat <- data.1 %>%
    dplyr::filter(MONTH == MONTH_value)
  
  model.month <- lmer(log10(mean_abund) ~ log_mass*scale(sexual_dichromatism) + (1 + scale(sexual_dichromatism)|grid_id), data = dat)
  
  ranef_df1 <- as.data.frame(ranef(model.month)$grid_id) %>% 
    rownames_to_column(var="grid_id") %>%
    mutate(intercept = `(Intercept)` + fixef(model.month)[1], slope = `scale(sexual_dichromatism)` + fixef(model.month)[3]) %>% 
    mutate(Month = MONTH_value)
}

modelling_results_month <- bind_rows(lapply(unique(data.1$MONTH), month_function))



#graphing results

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

modelling_results_month_grid <- grids %>%
  rename(grid_id = ID) %>%
  mutate(grid_id = as.character(as.integer(grid_id))) %>%
  left_join(modelling_results_month) %>%
  mutate(MONTH = as.character(as.integer(Month))) %>%
  dplyr::filter(complete.cases(Month))

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= modelling_results_month_grid, aes(fill = slope)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom") +
  facet_wrap(~ Month, ncol=4)
