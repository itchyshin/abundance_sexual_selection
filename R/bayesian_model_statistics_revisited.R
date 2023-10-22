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


data.1$z_ln_sexual_dichromatism <- scale(data.1$sexual_dichromatism)
data.1$z_log_mass <- scale(data.1$log_mass)
data.1$z_log_mean_abund <- scale(log(data.1$mean_abund))

data.2 <- data.1[ -which(is.na(data.1$log_mass)== T), -10]

prior <-list(R =list(V = 1, nu = 0.002),
             G = list(G1 =list(V =diag(2), nu = 2,
                               alpha.mu =rep(0, 2),
                               alpha.V=diag(25^2, 2, 2))))

##test

data.april.migrant <- data.2 %>%
  dplyr::filter(MONTH == "april") %>%
  dplyr::filter(migration == "migrant")

model1 <- MCMCglmm(z_log_mean_abund~ 
                     z_log_mass*z_ln_sexual_dichromatism, 
                   random = ~ us(1+z_ln_sexual_dichromatism):grid_id, 
                   data = data.april.migrant, 
                   prior = prior, 
                   pr=TRUE, 
                   pl=TRUE,
                   nitt=13000, thin=10, burnin=3000)

summary(model1)
slope.april.migrant <- model1$Sol
stats_table <- as.data.frame(t(summary(model1$Sol[,3])$quantiles))

##function?

month_migration_function <- function(MONTH_value){
  
  dat.month <- data.2 %>%
    dplyr::filter(MONTH == MONTH_value)
  
  migration_function <- function(migration_value){
    
    dat.migration <- dat.month %>%
      dplyr::filter(migration == migration_value)
    
    model <- MCMCglmm(z_log_mean_abund~ 
                         z_log_mass*z_ln_sexual_dichromatism, 
                       random = ~ us(1+z_ln_sexual_dichromatism):grid_id, 
                       data = dat.migration, 
                       prior = prior, 
                       pr=TRUE, 
                       pl=TRUE,
                       nitt=13000, thin=10, burnin=3000)
    
    stats_table <- as.data.frame(t(summary(model$Sol[,3])$quantiles)) 
  }
  
  results_migration <- bind_rows(lapply(unique(dat.month$migration), migration_function))
  
}

results_migration_month <- bind_rows(lapply(unique(data.2$MONTH), month_migration_function))
    
    
    
    
    
    
    
