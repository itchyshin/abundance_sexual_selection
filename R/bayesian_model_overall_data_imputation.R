# TODO separate the function to get slopes from the figure
# TODO 2 big functions - the one for all estimates and the other for each of 12 months (2 separeate R files?)

#########
library(sf)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(MCMCglmm)
library(lme4)
library(stringr)
library(phylolm)
library(here)
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(metafor)
library(magrittr)
library(purrr)
#library(furrr) # this seems very slow somehow
#future::plan(multiprocess)

data <- readRDS("Data/grid_ebird_by_month_joined_with_McQueen_imputed.RDS")

# just one tree
# phy.tre - got it from Corey - has 1000 trees
tree <- read.tree(here("Data", "phylo", "phy.tre"))

dim(data)
length(unique(factor(data$grid_id)))


#TODO - Josh -- I do not think this is correct way of doing group_by
# anyway it is all good - I will fix it 

# data.1 <- data %>%
#   group_by(grid_id, ebird_COMMON_NAME, sexual_dichromatism, log_mass) %>%
#   summarise(mean_abund = mean(mean_abund)) %>% ungroup()

#TODO 12 month
#TODO 50 trees 
#TODO spatially explicit model - need a distance matrix of different grids

# Indeed, we get slightly different results for this
# average of 12 month - different grids
dat_ave <- data %>%
  group_by(grid_id, ebird_COMMON_NAME) %>%
  summarise(density = mean(density), dens_v = mean(dens_se^2), sexual_dichromatism = dplyr::first(sexual_dichromatism), log_mass = dplyr::first(log_mass), migration = dplyr::first(migration), Phylogeny = dplyr::first(TipLabel))  %>% nest()

# the number of data for each grid
#purrr::map_int(dat_ave$data, nrow) -> n_grid
#n_grid <- tibble(n_grid = n_grid)
#dat_ave <- cbind(dat_ave, n_grid)

# NOW -- took out this and left all the models which can coverge
# some grid has very small sample size - we should have at least 20 species
#sum(n_grid > 19)
#dat_ave %<>% filter(n_grid > 19)

# function to prepare correlation matrix

cor_tree_fun <- function(df, tree){
  trimed_tree <- drop.tip(tree, which((tree$tip.label %in% df$Phylogeny) == FALSE))
  cor_tree <- vcv(trimed_tree,corr=T)
  cor_tree
}

cor_tree_poss <- possibly(.f = cor_tree_fun, otherwise = NA)

# getting cor_tree for all grids
# at least 50 trees but probably just 100 - more we use it more efficient for RRs
# TODO - this is were we can go tree[[1]] -- tree[[2]]
dat_ave %<>% mutate(cor_tree = purrr::map(data, cor_tree_poss, tree = tree[[100]]))

# function to fit metafor models
# try “control=list(optimizer=“optim”, optmethod=“BFGS”)” if this “control=list(optimizer=“optim”, optmethod=“Nelder-Mead”)” fails
# TODO -- we can model migration to see
# TODO -- put it back - we should go back to a normal model
mod_fun1 <- function(df, cor_tree){
  rma.mv(yi = density, V = dens_v, 
         mod = ~ sexual_dichromatism + log_mass,
         random = list(~1|Phylogeny, ~1|ebird_COMMON_NAME),
         R= list(Phylogeny = cor_tree),
         test = "t",
         control = list(optimizer = "optim", optmethod="Nelder-Mead"),
         data = df)
}

# migration is in it
#
mod_fun2 <- function(df, cor_tree){
  rma.mv(yi = density, V = dens_v, 
         mod = ~ sexual_dichromatism*migration + log_mass,
         random = list(~1|Phylogeny, ~1|ebird_COMMON_NAME),
         R= list(Phylogeny = cor_tree),
         test = "t",
         control = list(optimizer = "optim", optmethod="Nelder-Mead"),
         data = df)
}

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
poss_mod_fun1 <- possibly(.f = mod_fun2, otherwise = NULL)
poss_mod_fun2 <- possibly(.f = mod_fun3, otherwise = NULL)

# running the model for all the grids
#test<-dat_ave[100:102,]
dat_ave %<>% mutate(model = purrr::map2(data, cor_tree, poss_mod_fun1), model2 = purrr::map2(data, cor_tree, poss_mod_fun2))

# compact() - getting ride of NULL
# TODO - we probably should count how many NULL
leave <- map_lgl(dat_ave$model, ~is.null(.x) == FALSE)
#nrow(dat_ave) - (sum(leave))
dat_ave <- dat_ave[leave, ]

# function to get slopes for sexual dichromatism

# if it does not work, put NA
b_fun_poss1 <- possibly(.f = function(mod){
  mod$b[[2]]}, 
  otherwise = NA)

# b_fun_poss2 <- possibly(.f = function(mod){
#   mod$b[[2]] + mod$b[[5]]}, 
#   otherwise = NA)

# if warning, put NA
se_fun_poss1 <- possibly(.f = function(mod){
  mod$se[[2]]},
  otherwise = NA)

# I think this is worng and overestiamte for SE
# TODO - fix this 
# se_fun_poss2 <- possibly(.f = function(mod){
#   sqrt(mod$se[[2]]^2 + mod$se[[5]]^2)},
#   otherwise = NA)
 
dat_ave %>% transmute(grid_id, 
                      b_migrant = purrr::map_dbl(model, b_fun_poss1), 
                      b_resident = purrr::map_dbl(model2,b_fun_poss1), 
                      se_migrant = purrr::map_dbl(model,se_fun_poss1), 
                      se_resdident = purrr::map_dbl(model2,se_fun_poss1)
                      ) -> df2

# TODO - need to count which one failed???
# TODO - we need to get this for 50 trees
# TODO - we need to get this for 12 months + all_combined

# overall negative
mod_migrant <- rma(yi = b_migrant, sei = se_migrant, 
            test = "t",
            data = df2)

summary(mod_migrant)

distance <- res[df2$grid_id,df2$grid_id]
# km
distance2 <- res[df2$grid_id,df2$grid_id]/1000
# spatial model 
df2$const <- 1 # not quite sure what it means for now
mod_migrant2 <- rma.mv(yi = b_migrant, V = se_migrant^2, 
              random = list(~grid_id|const),
              struct = "SPGAU",
              dist = list(as.matrix(distance2)), # we need a distance matrix
              test = "t",
              control = list(optimizer = "optim", optmethod="Nelder-Mead"),
              data = df2)

summary(mod_migrant2)


# overall positive
mod_resident <- rma(yi = b_resident, sei = se_resdident,
            test = "t",
            data = df2)

summary(mod_resident)

mod_resident2 <- rma.mv(yi = b_resident, V = se_resdident^2, 
                       random = list(~grid_id|const),
                       struct = "SPGAU",
                       dist = list(as.matrix(distance2)), # we need a distance matrix
                       test = "t",
                       control = list(optimizer = "optim", optmethod="Nelder-Mead"),
                       data = df2)

summary(mod_resident2)


################################
# simpler model ################
################################


mod_fun1 <- function(df, cor_tree){
  rma.mv(yi = density, V = dens_v, 
         mod = ~ sexual_dichromatism + log_mass,
         random = list(~1|Phylogeny, ~1|ebird_COMMON_NAME),
         R= list(Phylogeny = cor_tree),
         test = "t",
         control = list(optimizer = "optim", optmethod="Nelder-Mead"),
         data = df)
}

# keeps going even with Error etc
poss_mod_fun1 <- possibly(.f = mod_fun1, otherwise = NULL)

# running the model for all the grids
dat_ave %<>% mutate(model = purrr::map2_dbl(data, cor_tree, poss_mod_fun1, .progress = TRUE))

# the slope for average density
dat_ave %>% transmute(grid_id, estimate_slope = purrr::map_dbl(model,~.x$b[[2]]), se_slope = purrr::map_dbl(model,~.x$se[[2]])) -> final_df

#TODO - we can do this for all 12 months
#TODO - we can do this for 50 trees

#TODO - we should average out everything 
mod0 <- rma(yi = estimate_slope, sei = se_slope, 
              test = "t",
              data = final_df)

summary(mod0)

# TODO

# getting gird_id and res alligned
distance <- res[final_df$grid_id,final_df$grid_id]
# km
distance2 <- res[final_df$grid_id,final_df$grid_id]/1000
# spatial model 
final_df$const <- 1 # not quite sure what it means for now
mod <- rma.mv(yi = estimate_slope, V = se_slope^2, 
              random = list(~grid_id|const),
              struct = "SPGAU",
              dist = list(as.matrix(distance2)), # we need a distance matrix
              test = "t",
              control = list(optimizer = "optim", optmethod="Nelder-Mead"),
              data = final_df)

summary(mod)

##########################
# let's make this Bayesian
##########################

# # creating z transformed sexual dichromatism
# data.1$z_ln_sexual_dichromatism <- scale(data.1$sexual_dichromatism)
# data.1$z_log_mass <- scale(data.1$log_mass)
# data.1$z_log_mean_abund <- scale(log(data.1$mean_abund))
# # taking NA values and the colum called "family" - MCMC
# data.2 <- data.1[ -which(is.na(data.1$log_mass)== T),]
# 
# ####### we can just load in the model and skip the next few steps #######
# 
# # takes a long time to run the model
# 
# #  We are not using this prior - it is a bit informative
# #prior0 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = diag(2), nu = 2)))
# # We use something else - nonish prior - I do not think whicher we use will not matter too much
# prior <-list(R =list(V = 1, nu = 0.002),
#              G = list(G1 =list(V =diag(2), nu = 2,
#                                alpha.mu =rep(0, 2),
#                                alpha.V=diag(25^2, 2, 2))))
# # everything is normalised (log or ln) and z-transformed (scaled)
# # we need to run this again
# model.1b <- MCMCglmm(z_log_mean_abund~ 
#                        z_log_mass*z_ln_sexual_dichromatism, 
#                      random = ~ us(1+z_ln_sexual_dichromatism):grid_id, 
#                      data = data.2, 
#                      prior = prior, 
#                      pr=TRUE, 
#                      pl=TRUE,
#                      nitt=13000*5, thin=10*5, burnin=3000*5
# )
# summary(model.1b)
# summary(model.1b$Sol)
# 
# ##model might need to be resaved
# ##saveRDS(model.1b, file = "./Objects/model.1b.RDS") 
# 
# model.1b <- readRDS("./Objects/model.1b.RDS")
# 
# # building the table
# 
# n_grid <- length(unique(data.2$grid_id))
# 
# stats_table <- summary(model.1b$Sol)$statistics
# 
# mean_intercept <- as.numeric(stats_table[1,1]) 
# sd_intercept <- as.numeric(stats_table[1,2]) 
# mean_slope <- as.numeric(stats_table[3,1])
# sd_slope <- as.numeric(stats_table[3,2]) 
# deviations <- stats_table[5:(n_grid*2+4),1:2]
# deviations <- as.data.frame(deviations) 
# deviations$grid <- row.names(deviations) # making row name a column
# rownames(deviations) <- NULL
# deviations_df <- as_tibble(deviations)
# deviations_df %>% transmute(grid = str_replace_all(grid, "z_ln_sexual_dichromatism.grid_id.", ""),
#                             grid = str_replace_all(grid, paste0("\\(","Intercept", "\\)",".grid_id."), ""),
#                             type = rep(c("intercept","slope"), each = n_grid),
#                             estimate = if_else(type == "intercept", Mean + mean_intercept, Mean + mean_slope), 
#                             se = if_else(type == "intercept", sqrt(SD^2 + sd_intercept^2), sqrt(SD^2 + sd_slope^2))) -> type_df
# 
# # reshaping data
# 
# type_df %>%
#   select(-se) %>%
#   tidyr::spread(type, estimate) %>%
#   rename(estimate_intercept = intercept) %>%
#   rename(estimate_slope = slope)-> estimate_df
# 
# type_df %>%
#   select(-estimate) %>%
#   spread(type, se) %>%
#   rename(se_intercept = intercept) %>%
#   rename(se_slope = slope) %>%
#   left_join(estimate_df, by = "grid") -> final_df

# graphing

grids <- st_read("Data/grid_5_degree_with_clim.geojson")

#library(geosphere)
#distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
library(dplyr)
library(sf)
library(geosphere)

dat <- st_read("Data/grid_5_degree.geojson") %>%
  st_set_geometry(NULL) %>%
  rename(centroid_lng=centroid_lat,
         centroid_lat=centroid_lng)
# getting the distrance matrix
res <- distm(dat[c("centroid_lng","centroid_lat")],dat[c("centroid_lng","centroid_lat")])
rownames(res) <- dat$ID
colnames(res) <- dat$ID
res <- res/1000

saveRDS(res,here("Data", "distance_km.RDS"))

ebird.grid.model.1b <- grids %>%
  rename(grid_id=ID) %>%
  #mutate(grid_id=as.character(as.integer(grid_id))) %>%
  left_join(final_df)

options("download.file.method" = "libcurl")
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

# some tests ------------

# testing with one grid
dat_1424 <- data %>% filter(grid_id == 1534)
dat_1424$tip.label <- dat_1424$TipLabel
dat_1424$Phylogeny <- dat_1424$TipLabel
dat_1424$Species <- dat_1424$ebird_COMMON_NAME
dat_1424$obs <- factor(1:nrow(dat_1424))
tree <- read.tree(here("Data", "phylo", "AllBirdsEricson1_summary.tre"))

#onetree <- tree[[1]]

match(dat_1424$tip.label, tree$tip.label)

trimed_tree <- drop.tip(tree, which((tree$tip.label %in% dat_1424$tip.label) == FALSE))
cor_tree <- vcv(trimed_tree,corr=T)
# mod <- phylolm(density ~ sexual_dichromatism + log_mass,
#                        data = dat_1424,
#                        phy = trimed_tree,
#                        model = "BM")

#mod <- gls(density ~ scale(sexual_dichromatism)*scale(log_mass), correlation = corBrownian(phy = trimed_tree),
#                 data = dat_1424)
#summary(mod)

mod <- rma.mv(yi = density, V = dens_se^2, 
              mod = ~ scale(sexual_dichromatism) + scale(log_mass),
              random = list(~1|MONTH, ~1|Phylogeny, ~1|Species, ~1| obs),
              R= list(Phylogeny = cor_tree),
              test = "t",
              data = dat_1424)

summary(mod)

res <-lm(density ~ scale(sexual_dichromatism) + scale(log_mass), data = dat_1424)
summary(res)

# making 572 grids to work on 
dat <- data %>% group_by(grid_id) %>% nest()

library(phylolm)



