#########
# test
########

# TODO - getting it work with one tree

library(tidyverse)
#library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(MCMCglmm)
library(lme4)
library(stringr)
library(here)
library(ape)
library(metafor)
library(magrittr)
library(purrr)
library(dplyr)
library(parallel)
library(tictoc)
library(furrr)
#library(future)
#future::plan(multisession, workers = 18)
#library(sf)
#library(geosphere)
library(orchaRd)


# bird abundance data
dat <- readRDS(here("Data","grid_ebird_by_month_joined_with_McQueen_imputed.RDS"))

# This has 1000 trees - we will use the first 100
tree <- read.tree(here("Data", "phylo", "phy.tre"))

#tree <- tree[1:10]
tree <- tree[11:20]
#tree <- tree[21:30] 
#tree <- tree[31:40] 
#tree <- tree[41:50] 

# distance matrix between grids
distance_km <- readRDS(here("Data", "distance_km.RDS"))

# Indeed, we get slightly different results for this
# average of 12 month - different grids
dat_ave <- dat %>%
  group_by(grid_id, ebird_COMMON_NAME) %>%
  summarise(density = mean(density), dens_v = mean(dens_se^2), 
            sexual_dichromatism = dplyr::first(sexual_dichromatism), 
            log_mass = dplyr::first(log_mass), migration = dplyr::first(migration), 
            Phylogeny = dplyr::first(TipLabel))  %>% tidyr::nest()

# setting out data frame

num_tree <- 10

res_list_slope <- vector(mode = "list", length = num_tree)
res_list_overall <- vector(mode = "list", length = num_tree)

tic()

for(i in 1:num_tree){
  
  # getting ride of trees and models every time
  res_list_slope[[i]] <- df2
  
  # function to prepare correlation matrix
  cor_tree_fun <- function(df, tree){
    trimed_tree <- drop.tip(tree, which((tree$tip.label %in% df$Phylogeny) == FALSE))
    cor_tree <- vcv(trimed_tree, corr=T)
    cor_tree
  }
  
  # function to catch warning and put NA
  cor_tree_poss <- possibly(.f = cor_tree_fun, otherwise = NA)
  
  # getting cor_tree for all grids
  # at least 50 trees 
  # TODO change it to mulit-core
  # system.time(dat_ave %<>% mutate(cor_tree = purrr::map(data, 
  #                                                      cor_tree_poss, tree = tree[[1]])))
  # tic()
  #dat_ave$cor_tree <- mclapply(dat_ave$data, cor_tree_poss, tree = tree[[i]], mc.cores = 18)
  #toc()
  
  dat_ave$cor_tree <- mclapply(dat_ave$data, cor_tree_poss, tree = tree[[i]], mc.cores = 18)
  
  # the function to fit: phylogenetic meta-analysis 
  # try “control=list(optimizer=“optim”, optmethod=“BFGS”)” if this “control=list(optimizer=“optim”, optmethod=“Nelder-Mead”)” fails
  mod_fun1 <- function(df, cor_tree){
    rma.mv(yi = density, V = dens_v, 
           mod = ~ sexual_dichromatism + log_mass,
           random = list(~1|Phylogeny, ~1|ebird_COMMON_NAME),
           R= list(Phylogeny = cor_tree),
           test = "t",
           sparse = TRUE,
           control = list(optimizer = "optim", optmethod="Nelder-Mead"),
           data = df)
  }
  
  # migration (for migrants)
  mod_fun2 <- function(df, cor_tree){
    rma.mv(yi = density, V = dens_v, 
           mod = ~ sexual_dichromatism*migration + log_mass,
           random = list(~1|Phylogeny, ~1|ebird_COMMON_NAME),
           R= list(Phylogeny = cor_tree),
           test = "t",
           sparse = TRUE,
           control = list(optimizer = "optim", optmethod="Nelder-Mead"),
           data = df)
  }
  # migration (for residents)
  mod_fun3 <- function(df, cor_tree){
    rma.mv(yi = density, V = dens_v, 
           mod = ~ sexual_dichromatism*relevel(factor(migration), ref = "resident") + log_mass,
           random = list(~1|Phylogeny, ~1|ebird_COMMON_NAME),
           R= list(Phylogeny = cor_tree),
           test = "t",
           sparse = TRUE,
           control = list(optimizer = "optim", optmethod="Nelder-Mead"),
           data = df)
  }
  
  # keeps going even with Error etc
  poss_mod_fun1 <- possibly(.f = mod_fun1, otherwise = NULL)
  poss_mod_fun2 <- possibly(.f = mod_fun2, otherwise = NULL)
  poss_mod_fun3 <- possibly(.f = mod_fun3, otherwise = NULL)
  
  # running the model for all the grids
  # TODO change it to mulit-core
  # running with 18 cores
  
  # #tic()
  # list1 <- mcmapply(poss_mod_fun1,
  #          dat_ave$data, dat_ave$cor_tree, 
  #          mc.preschedule = FALSE, mc.cores = 18)
  # 
  # list2 <- mcmapply(poss_mod_fun2,
  #          dat_ave$data, dat_ave$cor_tree, 
  #          mc.preschedule = FALSE, mc.cores = 18)
  # 
  # list3 <-mcmapply(poss_mod_fun3,
  #         dat_ave$data, dat_ave$cor_tree, 
  #         mc.preschedule = FALSE, mc.cores = 18)
  # 
  # 
  # #length(list1)
  # #length(list2)
  # #length(list3)
  # 
  # dat_ave$model1 <- list1
  # dat_ave$model2 <- list2
  # dat_ave$model3 <- list3
  # 
  #toc()
  
  #slower
  system.time(dat_ave %<>% mutate(model1 = purrr::map2(data, cor_tree, poss_mod_fun1),
                                  model2 = purrr::map2(data, cor_tree, poss_mod_fun2),
                                  model3 = purrr::map2(data, cor_tree, poss_mod_fun3)))
  
  # replacing NULL with NA
  null2 <- map_lgl(dat_ave$model2, ~is.null(.x) == TRUE)
  null3 <- map_lgl(dat_ave$model3, ~is.null(.x) == TRUE)
  dat_ave$model2[null2] <- NA
  dat_ave$model3[null3] <- NA
  
  # compact() - getting ride of NULL from model 1
  leave1 <- map_lgl(dat_ave$model1, ~is.null(.x) == FALSE)
  
  #nrow(dat_ave) - (sum(leave))
  dat_ave <- dat_ave[leave1, ]
  
  # function to get slopes + SE for sexual dichromatism
  # if it does not work, put NA
  s_est <- function(mod){mod$b[[2]]}
  s_error <- function(mod){mod$se[[2]]}
  
  b_fun_poss1 <- possibly(.f = s_est, otherwise = NA_real_)
  
  
  se_fun_poss1 <- possibly(.f = s_error, otherwise = NA_real_)
  
  #getting slopes and SE
  dat_ave %>% transmute(grid_id = grid_id,
                        b_overall = purrr::map_dbl(model1, b_fun_poss1),
                        b_migrant = purrr::map_dbl(model2, b_fun_poss1), 
                        b_resident = purrr::map_dbl(model3,b_fun_poss1), 
                        se_overall = purrr::map_dbl(model1,se_fun_poss1), 
                        se_migrant = purrr::map_dbl(model2,se_fun_poss1), 
                        se_resdident = purrr::map_dbl(model3,se_fun_poss1)
  ) -> df2
  
  # saving df2 - we get 50 df2
  # res_list_slope[[i]] <- 
  
  # spatial model 
  df2$const <- 1 # not quite sure what it means for now
  distance <- distance_km[df2$grid_id,df2$grid_id]
  
  
  mod_overall <- rma.mv(yi = b_overall, V = se_overall^2, 
                        random = list(~grid_id|const),
                        struct = "SPGAU",
                        dist = list(as.matrix(distance)), # we need a distance matrix
                        test = "t",
                        sparse = TRUE,
                        control = list(optimizer = "optim", optmethod="Nelder-Mead"),
                        data = df2)
  
  # migrant
  mod_migrant <- rma.mv(yi = b_migrant, V = se_migrant^2, 
                        random = list(~grid_id|const),
                        struct = "SPGAU",
                        dist = list(as.matrix(distance)), # we need a distance matrix
                        test = "t",
                        sparse = TRUE,
                        control = list(optimizer = "optim", optmethod="Nelder-Mead"),
                        data = df2)
  
  mod_resident <- rma.mv(yi = b_resident, V = se_resdident^2, 
                         random = list(~grid_id|const),
                         struct = "SPGAU",
                         dist = list(as.matrix(distance)), # we need a distance matrix
                         test = "t",
                         sparse = TRUE,
                         control = list(optimizer = "optim", optmethod="Nelder-Mead"),
                         data = df2)
  
  # putting 3 models together 
  # res_list_overall[[i]] <- list(mod_overall, mod_migrant, mod_resident)
  res_list_overall[[i]] <- list(mod_overall, mod_migrant, mod_resident)
  
}

toc()

saveRDS(res_list_slope, here("Data", "res_list_slope2.RDS"))
saveRDS(res_list_overall, here("Data", "res_list_overall2.RDS"))

