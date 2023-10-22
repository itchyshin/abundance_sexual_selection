library(sf)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(scales)

data.1 <- readRDS("Data/grid_ebird_by_month_joined_with_mcqueen.RDS")
data.1$MONTH <- factor(data.1$MONTH, levels=c("january", "february",
                                              "march", "april", "may",
                                              "june", "july", "august",
                                              "september", "october", 
                                              "november", "december"),
                       labels = c("january", "february",
                                  "march", "april", "may",
                                  "june", "july", "august",
                                  "september", "october", 
                                  "november", "december"))
data.1$migration <- factor(data.1$migration, levels = c("resident", "migrant"),
                           labels = c("resident", "migrant"))


grids <- st_read("Data/grid_5_degree_with_clim.geojson")
countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()
graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()
bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()

###overall

species.richness.overall <- data.1 %>%
  group_by(grid_id) %>%
  summarise(n_distinct(TipLabel)) %>%
  rename("Species" = "n_distinct(TipLabel)")

ggplot(species.richness.overall, aes(grid_id, Species)) + geom_point()

overall.summary <- as.data.frame(t(unclass(summary(species.richness.overall$Species))))

#map

species.richness.overall.map <- grids %>%
  rename(grid_id=ID) %>%
  left_join(species.richness.overall) 

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= species.richness.overall.map, aes(fill = Species)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradientn(limits = c(0,800), colors = c("white", "blue")) +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom")


###months/migration

month_migration_function <- function(MONTH_value){
  
  dat.month <- data.1 %>%
    dplyr::filter(MONTH == MONTH_value)
  
  migration_function <- function(migration_value){
    
    dat.migration <- dat.month %>%
      dplyr::filter(migration == migration_value) %>%
      group_by(grid_id) %>%
      summarise(n_distinct(TipLabel)) %>%
      rename("Species" = "n_distinct(TipLabel)") %>%
      mutate(Migration = migration_value) %>%
      mutate(Month = MONTH_value)
    
  }
  
  modelling_results_migration <- bind_rows(lapply(unique(dat.month$migration), migration_function))
  
}

species.richness.migration.months <- bind_rows(lapply(unique(data.1$MONTH), month_migration_function))

ggplot(species.richness.migration.months, aes(grid_id, Species)) + 
  geom_point() +
  facet_grid(rows = vars(Migration), cols = vars(Month))

#summary table

migration_month_function.2 <- function(migration_value){
  
  data.migration <- species.richness.migration.months %>%
    dplyr::filter(Migration == migration_value)
  
  month_function.2 <- function(Month_value){
    
    data.month <- data.migration %>%
      dplyr::filter(Month == Month_value) 
    
    month.summary <- as.data.frame(t(unclass(summary(data.month$Species)))) %>%
      mutate(Migration = migration_value) %>%
      mutate(Month = Month_value)
  
  }
  
  summary_results_month <- bind_rows(lapply(unique(data.migration$Month), month_function.2))  
  
}

migration.month.summary <- bind_rows(lapply(unique(species.richness.migration.months$Migration), migration_month_function.2))

#map

species.richness.migration.months.map <- grids %>%
  rename(grid_id=ID) %>%
  left_join(species.richness.migration.months) %>%
  drop_na(Migration) %>%
  drop_na(Month)

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= species.richness.migration.months.map, aes(fill = Species)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradientn(limits = c(0,800), colors = c("white", "blue")) +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom") +
  facet_grid(rows = vars(Migration), cols = vars(Month))

#separating so its easier

graph.resident <- filter(species.richness.migration.months.map, Migration == "resident")

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= graph.resident, aes(fill = Species)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradientn(limits = c(0,800), colors = c("white", "blue")) +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom") +
  facet_wrap(~Month, ncol = 4)

graph.migrant <- filter(species.richness.migration.months.map, Migration == "migrant")

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= graph.migrant, aes(fill = Species)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradientn(limits = c(0,800), colors = c("white", "blue")) +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom") +
  facet_wrap(~Month, ncol = 4)

