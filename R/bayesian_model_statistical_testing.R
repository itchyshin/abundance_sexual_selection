library(sf)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(MCMCglmm)
library(lme4)
library(metafor)

data.full <- readRDS("Data/bayesian_results_months_migration.RDS")
grids <- st_read("Data/grid_5_degree_with_clim.geojson")

countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()

graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()

#temopral plot (for reference)

grids_to_select <- data.full %>%
  group_by(grid, Migration) %>%
  summarize(N=n()) %>%
  dplyr::filter(N==12) %>%
  mutate(include="Yes")

data.full %>%
  left_join(grids_to_select) %>%
  dplyr::filter(complete.cases(include)) %>%
  ggplot(., aes(x=factor(Month, levels=c("January", "February",
                                         "March", "April", "May",
                                         "June", "July", "August",
                                         "September", "October", 
                                         "November", "December")), 
                y=estimate_slope, group=grid))+
  geom_point(color="gray30", size=0.002)+
  geom_line(color="gray10", size=0.002)+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1.3)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(panel.grid=element_blank())+
  facet_wrap(~Migration, ncol=1, scales="free")+
  xlab("")+
  ylab("Slope")+
  geom_smooth(aes(color=factor(Migration), group=Migration))+
  labs(color="Migration")

#function for comparing month/migration different to zero

migration_month_function <- function(Migration_value){
  
  data.migration <- data.full %>%
    dplyr::filter(Migration == Migration_value)
  
    month_function <- function(Month_value){
      
      data.month <- data.migration %>%
        dplyr::filter(Month == Month_value)
      
      model <- rma(yi = estimate_slope, sei = se_slope, data = data.month)
      
      test_df <- as.data.frame(coef(summary(model))) %>% 
        rownames_to_column(var="Month") %>%
        mutate(Month = Month_value) %>%
        mutate(Migration = Migration_value) %>%
        select(Migration, Month, pval, estimate, se) %>%
        mutate(., Significance = if_else(pval < 0.05, "Significant", "Insignificant"))
    
    }

test_results_month <- bind_rows(lapply(unique(data.migration$Month), month_function))  

}

test_results_final <- bind_rows(lapply(unique(data.full$Migration), migration_month_function))

# jan - M vs resident

data.full %>% 
  filter(Month == "January") %>%  
  group_by(grid) %>% 
  summarise(slope_diff = estimate_slope[Migration == "Resident"] - estimate_slope[Migration == "Migrant"], 
            slope_diff_se = sqrt(se_slope[Migration == "Resident"]^2 + se_slope[Migration == "Migrant"]^2)) -> test_data

model <- rma(yi = slope_diff, sei = slope_diff_se, data = test_data)
summary(model)

#function for comparing residents - migrants

difference_migration_function <- function(Month_value){
  
  data.month.by.month <- data.full %>%
    dplyr::filter(Month == Month_value) %>%
    group_by(grid) %>% 
    summarise(slope_diff = estimate_slope[Migration == "Resident"] - estimate_slope[Migration == "Migrant"], 
              slope_diff_se = sqrt(se_slope[Migration == "Resident"]^2 + se_slope[Migration == "Migrant"]^2))
  
  model <- rma(yi = slope_diff, sei = slope_diff_se, data = data.month.by.month)
  
  test_df <- as.data.frame(coef(summary(model))) %>% 
    rownames_to_column(var="Month") %>%
    mutate(Month = Month_value) %>%
    select(Month, pval, estimate, se) %>%
    mutate(., Significance = if_else(pval < 0.05, "Significant", "Insignificant"))
    
}

test_results_difference <- bind_rows(lapply(unique(data.full$Month), difference_migration_function))


#two-factor ANOVA then TUKEY

aov.test <- aov(estimate_slope ~ Migration*Month, data = data.full)
summary(aov.test)
TukeyHSD(aov.test) 

#one-factor for migrants

data.migrants <- filter(data.full, Migration == "Migrant")
aov.migrant.test <- aov(estimate_slope ~ Month, data = data.migrants)
TukeyHSD(aov.migrant.test)










