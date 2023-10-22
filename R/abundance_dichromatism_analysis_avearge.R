###################
# Overall analysis
###################

# loading packages
#library(sf)
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
#library(sf)
#library(geosphere)
library(orchaRd)
#library(furrr) # this seems very slow somehow
#future::plan(multiprocess)

# bird abundance data
data <- readRDS(here("Data","grid_ebird_by_month_joined_with_McQueen_imputed.RDS"))

# This has 1000 trees - we will use the first 100
tree <- read.tree(here("Data", "phylo", "phy.tre"))
#tree <- tree[1:50]

# distance matrix between grids
distance_km <- readRDS(here("Data", "distance_km.RDS"))

# Indeed, we get slightly different results for this
# average of 12 month - different grids
# TODO - check correalitons between dichro and dimorph

# > cor.test(test$sexual_dichromatism, test$sexual_dimorphism)
# 
# Pearson's product-moment correlation
# 
# data:  test$sexual_dichromatism and test$sexual_dimorphism
# t = 7.2763, df = 45534, p-value = 3.485e-13
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.02490236 0.04325076
# sample estimates:
#        cor 
# 0.03407943 

dat_ave <- data %>%
  # this data set has some missing values 
  # a lot of no sex specific mass - only 240 are missing - if they do not have sex specific - can we do it without
  # TODO - we can do ones with 
  left_join(., read_csv("Data/cleaned_body_size_data.csv", na = c(-999, "NA")) %>%
              mutate(sexual_dimorphism = male_body_mass_g/female_body_mass_g) %>%
              dplyr::select(COMMON_NAME, sexual_dimorphism) %>%
              rename(ebird_COMMON_NAME=COMMON_NAME)) %>%
  group_by(grid_id, ebird_COMMON_NAME) %>%
  summarise(density = mean(density), 
            dens_v = mean(dens_se^2), 
            sexual_dichromatism = dplyr::first(sexual_dichromatism),
            sexual_dimorphism = dplyr::first(sexual_dimorphism),
            log_mass = dplyr::first(log_mass), 
            migration = dplyr::first(migration), 
            Phylogeny = dplyr::first(TipLabel)) #-> test #%>% tidyr::nest()

#cor.test(test$sexual_dichromatism, test$sexual_dimorphism)


# setting out data frame

num_tree <- 50

res_list <- vector(mode = "list", length = num_tree)

for(i in 1:num_tree){

# function to prepare correlation matrix
cor_tree_fun <- function(df, tree){
  trimed_tree <- drop.tip(tree, which((tree$tip.label %in% df$Phylogeny) == FALSE))
  cor_tree <- vcv(trimed_tree,corr=T)
  cor_tree
}

# function to catch warning and put NA
cor_tree_poss <- possibly(.f = cor_tree_fun, otherwise = NA)

# getting cor_tree for all grids
# at least 50 trees 
# TODO - this is were we can go tree[[1]] -- tree[[2]]
dat_ave %<>% mutate(cor_tree = purrr::map(data, cor_tree_poss, tree = tree[[i]]))

# the function to fit: phylogenetic meta-analysis 
# try “control=list(optimizer=“optim”, optmethod=“BFGS”)” if this “control=list(optimizer=“optim”, optmethod=“Nelder-Mead”)” fails
mod_fun1 <- function(df, cor_tree){
  rma.mv(yi = density, V = dens_v, 
         mod = ~ sexual_dichromatism + log_mass,
         random = list(~1|Phylogeny, ~1|ebird_COMMON_NAME),
         R= list(Phylogeny = cor_tree),
         test = "t",
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
         control = list(optimizer = "optim", optmethod="Nelder-Mead"),
         data = df)
}

# keeps going even with Error etc
poss_mod_fun1 <- possibly(.f = mod_fun1, otherwise = NULL)
poss_mod_fun2 <- possibly(.f = mod_fun2, otherwise = NULL)
poss_mod_fun3 <- possibly(.f = mod_fun3, otherwise = NULL)

# running the model for all the grids
dat_ave %<>% mutate(model1 = purrr::map2(data, cor_tree, poss_mod_fun1), 
                    model2 = purrr::map2(data, cor_tree, poss_mod_fun2), 
                    model3 = purrr::map2(data, cor_tree, poss_mod_fun3))

# compact() - getting ride of NULL
leave <- map_lgl(dat_ave$model, ~is.null(.x) == FALSE)
#nrow(dat_ave) - (sum(leave))
dat_ave <- dat_ave[leave, ]

# function to get slopes + SE for sexual dichromatism
# if it does not work, put NA
b_fun_poss1 <- possibly(.f = function(mod){
  mod$b[[2]]}, 
  otherwise = NA)

se_fun_poss1 <- possibly(.f = function(mod){
  mod$se[[2]]},
  otherwise = NA)

#getting slopes and SEs
dat_ave %>% transmute(grid_id,
                      b_overall = purrr::map_dbl(model1, b_fun_poss1),
                      b_migrant = purrr::map_dbl(model2, b_fun_poss1), 
                      b_resident = purrr::map_dbl(model3,b_fun_poss1), 
                      se_overall = purrr::map_dbl(model1,se_fun_poss1), 
                      se_migrant = purrr::map_dbl(model2,se_fun_poss1), 
                      se_resdident = purrr::map_dbl(model3,se_fun_poss1)
) -> df2

# spatial model 
df2$const <- 1 # not quite sure what it means for now
distance <- distance_km[df2$grid_id,df2$grid_id]


mod_overall <- rma.mv(yi = b_overall, V = se_overall^2, 
                        random = list(~grid_id|const),
                        struct = "SPGAU",
                        dist = list(as.matrix(distance)), # we need a distance matrix
                        test = "t",
                        control = list(optimizer = "optim", optmethod="Nelder-Mead"),
                        data = df2)

# migrant
mod_migrant <- rma.mv(yi = b_migrant, V = se_migrant^2, 
                       random = list(~grid_id|const),
                       struct = "SPGAU",
                       dist = list(as.matrix(distance)), # we need a distance matrix
                       test = "t",
                       control = list(optimizer = "optim", optmethod="Nelder-Mead"),
                       data = df2)

mod_resident <- rma.mv(yi = b_resident, V = se_resdident^2, 
                        random = list(~grid_id|const),
                        struct = "SPGAU",
                        dist = list(as.matrix(distance)), # we need a distance matrix
                        test = "t",
                        control = list(optimizer = "optim", optmethod="Nelder-Mead"),
                        data = df2)

# putting 3 models together 
res_list[[i]] <-list(mod_overall, mod_migrant, mod_resident) 

}


