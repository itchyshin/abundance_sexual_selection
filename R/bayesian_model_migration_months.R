library(sf)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(MCMCglmm)
library(lme4)

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

##########################
# let's make this Bayesian
##########################

# creating z transformed sexual dichromatism
data.1$z_ln_sexual_dichromatism <- scale(data.1$sexual_dichromatism)
data.1$z_log_mass <- scale(data.1$log_mass)
data.1$z_log_mean_abund <- scale(log(data.1$mean_abund))
# taking NA values and the colum called "family" - MCMC
data.2 <- data.1[ -which(is.na(data.1$log_mass)== T), -10]

#  We are not using this prior - it is a bit informative
#prior0 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = diag(2), nu = 2)))
# We use something else - nonish prior - I do not think whicher we use will not matter too much
prior <-list(R =list(V = 1, nu = 0.002),
             G = list(G1 =list(V =diag(2), nu = 2,
                               alpha.mu =rep(0, 2),
                               alpha.V=diag(25^2, 2, 2))))
# everything is normalised (log or ln) and z-transformed (scaled)
# we need to run this again

month_migration_function <- function(MONTH_value){
  
  dat.month <- data.2 %>%
    dplyr::filter(MONTH == MONTH_value)

  migration_function <- function(migration_value){
  
      dat.migration <- dat.month %>%
        dplyr::filter(migration == migration_value)
  
      model.1b <- MCMCglmm(z_log_mean_abund~ 
                           z_log_mass*z_ln_sexual_dichromatism, 
                         random = ~ us(1+z_ln_sexual_dichromatism):grid_id, 
                         data = dat.migration, 
                         prior = prior, 
                         pr=TRUE, 
                         pl=TRUE,
                         nitt=13000*5, thin=10*5, burnin=3000*5
      )

       # building the table
                         
        n_grid <- length(unique(dat.migration$grid_id))
                             
        stats_table <- summary(model.1b$Sol)$statistics
                             
        mean_intercept <- as.numeric(stats_table[1,1]) 
        sd_intercept <- as.numeric(stats_table[1,2]) 
        mean_slope <- as.numeric(stats_table[3,1])
        sd_slope <- as.numeric(stats_table[3,2]) 
        deviations <- stats_table[5:(n_grid*2+4),1:2]
        deviations <- as.data.frame(deviations) 
        deviations$grid <- row.names(deviations) # making row name a column
        rownames(deviations) <- NULL
        deviations_df <- as_tibble(deviations)
        deviations_df %>% transmute(grid = str_replace_all(grid, "z_ln_sexual_dichromatism.grid_id.", ""),
                                    grid = str_replace_all(grid, paste0("\\(","Intercept", "\\)",".grid_id."), ""),
                                    type = rep(c("intercept","slope"), each = n_grid),
                                    estimate = if_else(type == "intercept", Mean + mean_intercept, Mean + mean_slope), 
                                    se = if_else(type == "intercept", sqrt(SD^2 + sd_intercept^2), sqrt(SD^2 + sd_slope^2))) -> type_df
                             
         # reshaping data
                             
        type_df %>%
          select(-se) %>%
          tidyr::spread(type, estimate) %>%
          rename(estimate_intercept = intercept) %>%
          rename(estimate_slope = slope)-> estimate_df
                             
        type_df %>%
          select(-estimate) %>%
          spread(type, se) %>%
          rename(se_intercept = intercept) %>%
          rename(se_slope = slope) %>%
          left_join(estimate_df, by = "grid") %>%
          mutate(Migration = migration_value) %>%
          mutate(Month = MONTH_value) -> final_df
          
 
  }

  modelling_results_migration <- bind_rows(lapply(unique(dat.month$migration), migration_function))
  
}

modelling_results_month_migration <- bind_rows(lapply(unique(data.2$MONTH), month_migration_function))

#saveRDS(modelling_results_month_migration, "Data/bayesian_results_months_migration.RDS")
modelling_results_month_migration <- readRDS("Data/bayesian_results_months_migration.RDS")

# graphing

grids <- st_read("Data/grid_5_degree_with_clim.geojson")

countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()

graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()

ebird.grid.model.2b <- grids %>%
  rename(grid=ID) %>%
  mutate(grid=as.character(as.integer(grid))) %>%
  left_join(modelling_results_month_migration) %>%
  rename(Latitude = centroid_lng) %>%
  rename(Longitude = centroid_lat) %>%
  drop_na(Migration) %>%
  drop_na(Month)

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= ebird.grid.model.2b, aes(fill = estimate_slope)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradientn(limits = c(-0.04941918, 0.04665269), colors = c("red", "white", "blue")) +
  theme(axis.text = element_blank()) +
  theme(legend.position = "bottom") +
  facet_grid(rows = vars(Migration), cols = vars(Month))

# separating for ease of reading

filter.months <- c("January", "May")

graph.resident <- ebird.grid.model.2b %>%
  filter(., Migration == "Resident")
  
ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= graph.resident, aes(fill = estimate_slope)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradientn(limits = c(-0.05, 0.05), colors = c("red", "white", "blue")) +
  theme(axis.text = element_blank()) +
  theme(text = element_text(size=13)) +
  theme(legend.position = "none") +
  labs(fill = "Regression Value of Sexual Dichromatism \nand Relative Abundance") +
  ggtitle("a) Resident Passerines") +
  guides(fill = guide_colorbar(barwidth = 10)) +
  theme(plot.title = element_text(hjust=0)) +
  facet_wrap(~Month, ncol = 4)

ggsave(file='Final_Full_Resident_Graph.png', width = 10, height=10, dpi=600)

graph.migrant <- ebird.grid.model.2b %>%
  filter(., Migration == "Migrant") 

ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= graph.migrant, aes(fill = estimate_slope)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradientn(limits = c(-0.05, 0.05), colors = c("red", "white", "blue")) +
  theme(axis.text = element_blank()) +
  theme(text = element_text(size=13)) +
  theme(legend.position = "bottom", legend.spacing.x = unit(1, 'cm')) +
  labs(fill = "Slope Values of Relative Abundance \nRegressed on Sexual Dichromatism") +
  ggtitle("b) Migrant Passerines") +
  theme(plot.title = element_text(hjust=0)) +
  guides(fill = guide_colorbar(barwidth = 13)) +
  facet_wrap(~Month, ncol = 4)

ggsave(file='Final_Full_Migrant_Graph.png', width = 10, height=10, dpi=600)

# testing across longitude

lat.investigate <- ebird.grid.model.2b %>%
  group_by(Latitude, Migration, Month) %>%
  summarise(Slope = mean(estimate_slope))

ggplot(lat.investigate, aes(Latitude, Slope)) + 
  geom_point() +
  scale_fill_gradientn(limits = c(-0.04941918, 0.04665269)) +
  facet_grid(rows = vars(Migration), cols = vars(Month))

# only get the grids with 12 months of data

grids_to_select <- modelling_results_month_migration %>%
  group_by(grid, Migration) %>%
  summarize(N=n()) %>%
  dplyr::filter(N==12) %>%
  mutate(include="Yes")

## temporal plot

modelling_results_month_migration %>%
  left_join(grids_to_select) %>%
  dplyr::filter(complete.cases(include)) %>%
  ggplot(., aes(x=factor(Month, levels=c("January", "February",
                                         "March", "April", "May",
                                         "June", "July", "August",
                                         "September", "October", 
                                         "November", "December")), 
                y=estimate_slope, group=grid))+
  #geom_point(color="gray30", size=0.002)+
  geom_line(color="gray10", size=0.002)+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1.3)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(panel.grid=element_blank())+
  facet_wrap(~Migration, ncol=1, scales="fixed")+
  xlab("")+
  ylab("Regression Value of Sexual Dichromatism and Relative Abundance")+
  geom_smooth(aes(color=factor(Migration), group=Migration))+
  labs(color="Migration") +
  theme(legend.position="none")



##########################################
## checking if the slopes are signficantly different from 0
# an example
modelling_results_month_migration %>% filter(Migration == "migrant", Month == "january") -> test_dat
modelling_results_month_migration %>% filter(Migration == "migrant", Month == "august") -> test_dat2
modelling_results_month_migration %>% filter(Migration == "migrant", Month == "may") -> test_dat3
library(metafor)

# jan
model <- rma(yi = estimate_slope, sei = se_slope, data = test_dat)
test.df<-coef(summary(model))

summary(model)
# aug
model2 <- rma(yi = estimate_slope, sei = se_slope, data = test_dat2)

summary(model2)

# may
model3 <- rma(yi = estimate_slope, sei = se_slope, data = test_dat3)

summary(model3)


# jan - M vs resident

modelling_results_month_migration %>% filter(Month == "january") %>%  group_by(grid) %>% 
  summarise(slope_diff = estimate_slope[Migration == "resident"] - estimate_slope[Migration == "migrant"], 
            slope_diff_se = sqrt(se_slope[Migration == "resident"]^2 + se_slope[Migration == "migrant"]^2)) -> test_dat4

# may
model4 <- rma(yi = slope_diff, sei = slope_diff_se, data = test_dat4)

summary(model4)

#######################
# Shinichi's addition
########################
dat.month <- data.2 %>%
  dplyr::filter(MONTH == "april")

dat.migration <- dat.month %>%
  dplyr::filter(migration == "migrant")


dat.month2 <- data.2 %>%
  dplyr::filter(MONTH == "may")

dat.migration2 <- dat.month2 %>%
  dplyr::filter(migration == "migrant")


dat.month3 <- data.2 %>%
  dplyr::filter(MONTH == "june")

dat.migration3 <- dat.month3 %>%
  dplyr::filter(migration == "migrant")



model1 <- MCMCglmm(z_log_mean_abund~ 
              z_log_mass*z_ln_sexual_dichromatism, 
            random = ~ us(1+z_ln_sexual_dichromatism):grid_id, 
            data = dat.migration, 
            prior = prior, 
            pr=TRUE, 
            pl=TRUE,
            nitt=13000, thin=10, burnin=3000)

summary(model1)
summmary(model1$Sol)

slope_april <- model1$Sol[,3]
plot(slope_april)


model2 <- MCMCglmm(z_log_mean_abund~ 
                     z_log_mass*z_ln_sexual_dichromatism, 
                   random = ~ us(1+z_ln_sexual_dichromatism):grid_id, 
                   data = dat.migration2, 
                   prior = prior, 
                   pr=TRUE, 
                   pl=TRUE,
                   nitt=13000, thin=10, burnin=3000)
summary(model2)
summary(model2$Sol)
slope_may <- model2$Sol[,3]

# comparing april vs may

summary(slope_may - slope_april)
plot(slope_may - slope_april)

##########
model3 <- MCMCglmm(z_log_mean_abund~ 
                     z_log_mass*z_ln_sexual_dichromatism, 
                   random = ~ us(1+z_ln_sexual_dichromatism):grid_id, 
                   data = dat.migration2, 
                   prior = prior, 
                   pr=TRUE, 
                   pl=TRUE,
                   nitt=13000, thin=10, burnin=3000)
summary(model3)
summary(model3$Sol)
