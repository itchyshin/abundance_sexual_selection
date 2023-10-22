## read in packages and data

library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

data.1 <- readRDS("Data/grid_ebird_by_month_joined_with_mcqueen.RDS")
hist(data.1$sexual_dichromatism)

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

## functioning

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
str(modelling_results_migration_month)

## grouping months

modelling_results_migration_month <- modelling_results_migration_month %>%
  mutate(Grouped_Month = Month)

modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "december", "Dec-Feb")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "january", "Dec-Feb")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "february", "Dec-Feb")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "march", "Mar-May")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "april", "Mar-May")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "may", "Mar-May")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "june", "Jun-Aug")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "july", "Jun-Aug")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "august", "Jun-Aug")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "september", "Sept-Nov")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "october", "Sept-Nov")
modelling_results_migration_month$Grouped_Month <- str_replace_all(modelling_results_migration_month$Grouped_Month, "november", "Sept-Nov")


## graphing

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

## overall

migration_month_grids <- grids %>%
  rename(grid_id = ID) %>%
  mutate(grid_id = as.character(as.integer(grid_id))) %>%
  left_join(modelling_results_migration_month) %>%
  dplyr::filter(complete.cases(Month))

str(migration_month_grids)

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data = migration_month_grids, aes(fill = slope)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom") +
  facet_grid(rows = vars(Month), cols = vars(Migration))

## serparating migration patterns

migration_month_grids_resident <- grids %>%
  rename(grid_id = ID) %>%
  mutate(grid_id = as.character(as.integer(grid_id))) %>%
  left_join(modelling_results_migration_month) %>%
  dplyr::filter(complete.cases(Month)) %>%
  dplyr::filter(Migration == "resident")


ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data = migration_month_grids_resident, aes(fill = slope)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom") +
  facet_wrap(~ Month, ncol = 4)




migration_month_grids_migrant <- grids %>%
  rename(grid_id = ID) %>%
  mutate(grid_id = as.character(as.integer(grid_id))) %>%
  left_join(modelling_results_migration_month) %>%
  dplyr::filter(complete.cases(Month)) %>%
  dplyr::filter(Migration == "migrant")

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data = migration_month_grids_migrant, aes(fill = slope)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom") +
  facet_wrap(~ Month, ncol = 4)

# only get the grids with 12 months of data

grids_to_select <- modelling_results_migration_month %>%
  group_by(grid_id, Migration) %>%
  summarize(N=n()) %>%
  dplyr::filter(N==12) %>%
  mutate(include="Yes")

## temporal plot

modelling_results_migration_month %>%
  left_join(grids_to_select) %>%
  dplyr::filter(complete.cases(include)) %>%
  ggplot(., aes(x=factor(Month, levels=c("january", "february",
                                         "march", "april", "may",
                                         "june", "july", "august",
                                         "september", "october", 
                                         "november", "december")), 
                y=slope, group=grid_id))+
  geom_point(color="gray30", size=0.002)+
  geom_line(color="gray10", size=0.002)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(panel.grid=element_blank())+
  facet_wrap(~Migration, ncol=1, scales="free")+
  xlab("")+
  ylab("Slope")+
  geom_smooth(aes(color=factor(Migration), group=Migration))+
  labs(color="Migration")





