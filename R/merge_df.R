library(pacman)

p_load(tidyverse,
       here,
       miWQS)

# putting 50 result

df1 <- readRDS(here("Data", "res_list_slope1.RDS"))
df2 <- readRDS(here("Data", "res_list_slope2.RDS"))
df3 <- readRDS(here("Data", "res_list_slope3.RDS"))
df4 <- readRDS(here("Data", "res_list_slope4.RDS"))
df5 <- readRDS(here("Data", "res_list_slope5.RDS"))

dat <- list_rbind(c(df1, df2, df3, df4, df5))

dat %>% ungroup %>%  select("grid_id",
                            "b_overall", "se_overall",
                            "b_migrant", "se_migrant",
                            "b_resident","se_resdident") %>% 
  nest_by(grid_id) -> dat_nest

dat_nest$data[[1]] %>% as.matrix %>% t %>% 
  as.vector -> test


####

my_array <- array(test, dim = c(3, 2, 50))

# combining all the results from 50 models using Rubin's rules
pool.mi(my_array) 


# example 

# combining the results
est_50 <- map_dbl(ma_50, ~ .x$b[[1]])
se_50 <-  map_dbl(ma_50, ~ .x$se)
df_50 <- c(rbind(est_50, se_50))

# combine the data frames into an array
my_array <- array(df_50, dim = c(3, 2, 50))

# combining all the results from 50 models using Rubin's rules
pool.mi(my_array) 