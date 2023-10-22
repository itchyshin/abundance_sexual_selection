library(sf)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(MCMCglmm)
library(lme4)
library(stringr)

data <- readRDS("Data/grid_ebird_by_month_joined_with_mcqueen.RDS")

#TODO - Joshn -- I do not think this is correct way of doing group_by
data.1 <- data %>%
  group_by(grid_id, ebird_COMMON_NAME, sexual_dichromatism, log_mass) %>%
  summarise(mean_abund = mean(mean_abund)) %>% ungroup()

# Inded, we get slightly different results for this
test <- data %>%
  group_by(grid_id, ebird_COMMON_NAME) %>%
  summarise(mean_abund = mean(mean_abund), sexual_dichromatism = mean(sexual_dichromatism), log_mass = mean(sexual_dichromatism))  %>% ungroup()

##########################
# let's make this Bayesian
##########################

# creating z transformed sexual dichromatism
data.1$z_ln_sexual_dichromatism <- scale(data.1$sexual_dichromatism)
data.1$z_log_mass <- scale(data.1$log_mass)
data.1$z_log_mean_abund <- scale(log(data.1$mean_abund))
# taking NA values and the colum called "family" - MCMC
data.2 <- data.1[ -which(is.na(data.1$log_mass)== T),]

####### we can just load in the model and skip the next few steps #######

# takes a long time to run the model

#  We are not using this prior - it is a bit informative
#prior0 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = diag(2), nu = 2)))
# We use something else - nonish prior - I do not think whicher we use will not matter too much
prior <-list(R =list(V = 1, nu = 0.002),
             G = list(G1 =list(V =diag(2), nu = 2,
                               alpha.mu =rep(0, 2),
                               alpha.V=diag(25^2, 2, 2))))
# everything is normalised (log or ln) and z-transformed (scaled)
# we need to run this again
model.1b <- MCMCglmm(z_log_mean_abund~ 
                       z_log_mass*z_ln_sexual_dichromatism, 
                     random = ~ us(1+z_ln_sexual_dichromatism):grid_id, 
                     data = data.2, 
                     prior = prior, 
                     pr=TRUE, 
                     pl=TRUE,
                     nitt=13000*5, thin=10*5, burnin=3000*5
)
summary(model.1b)
summary(model.1b$Sol)

##model might need to be resaved
##saveRDS(model.1b, file = "./Objects/model.1b.RDS") 

model.1b <- readRDS("./Objects/model.1b.RDS")

# building the table

n_grid <- length(unique(data.2$grid_id))

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
  left_join(estimate_df, by = "grid") -> final_df

# graphing

grids <- st_read("Data/grid_5_degree_with_clim.geojson")

ebird.grid.model.1b <- grids %>%
  rename(grid=ID) %>%
  mutate(grid=as.character(as.integer(grid))) %>%
  left_join(final_df)


countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()

graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()


ggplot() +
  geom_sf(data = bb, col = "grey20", fill = "transparent") +
  geom_sf(data= ebird.grid.model.1b, aes(fill = estimate_slope)) +
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4) +
  theme_minimal() +
  scale_fill_gradientn(limits = c(-0.05, 0.05), colors = c("red", "white", "blue")) +
  theme(axis.text = element_blank()) +
  theme(text = element_text(size=15)) +
  theme(legend.position = "bottom", legend.spacing.x = unit(1, 'cm')) +
  guides(fill = guide_colorbar(barwidth = 10)) +
  labs(fill = "Slope Values of Relative Abundance \nRegressed on Sexual Dichromatism")

ggsave(file='Overall_Graph_Final.png', width = 10, height=8, dpi=600)




