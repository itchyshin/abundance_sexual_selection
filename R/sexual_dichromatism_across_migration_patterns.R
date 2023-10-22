
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(beeswarm)
library(phylolm)
library(ape)

data.1 <- readRDS("Data/grid_ebird_by_month_joined_with_mcqueen.RDS")
data.1$migration <- str_replace_all(data.1$migration, "migrant", "Migrant")
data.1$migration <- str_replace_all(data.1$migration, "resident", "Resident")

data.1$migration <- factor(data.1$migration, levels = c("Resident", "Migrant"),
                           labels = c("Resident", "Migrant"))
str(data.1)
hist(data.1$sexual_dichromatism)

## sexual dichromatism

data.1 %>%
  group_by(migration) %>%
  summarise(n_distinct(TipLabel))

data.unique <- data.1 %>%
  dplyr::select(ebird_COMMON_NAME, migration, sexual_dichromatism) %>%
  distinct()

data.unique %>%
  ggplot(., aes(x = migration, y = sexual_dichromatism, fill = migration)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  theme(legend.position = "none") +
  labs(y = "Sexual Dichromatism Score", x ="")

ggsave(file='Final_Dichromatism_Graph.png', width = 6, height=6, dpi=600)
  
beeswarm(sexual_dichromatism ~ migration, data = data.unique, col = c('red', 'blue'), method = "hex", cex = 0.43)

t.test(sexual_dichromatism ~ migration, data = data.unique, var.equal = TRUE)

## phylogeny and test
#download.file("https://cloudstor.aarnet.edu.au/plus/s/pysK99S3w6Ft4ae/download","Data/phylo/birdtree.tre")

data.test <- data.1 %>%
  dplyr::select(TipLabel, migration, sexual_dichromatism) %>%
  distinct() %>%
  rename(., tip.label = TipLabel) %>%
  remove_rownames() %>%
  column_to_rownames(var = "tip.label")

tree <- read.tree("Data/phylo/birdtree.tre")

run_phylo_model <- function(i){
  onetree <- tree[[i]]
  #run phylo model here
  model.phylo <- phylolm(sexual_dichromatism ~ migration,
                         data = data.test,
                         phy = onetree,
                         model = "BM")
  
  #extract coefficient and SE
  return(coef(model.phylo))
}

phylo_model_out <- map_df(1:1000,run_phylo_model)

quantile(phylo_model_out$migrationresident, c(0.025, 0.975))
hist(phylo_model_out$migrationresident)

## slope

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

modelling_results_migration_month %>%
  group_by(Migration) %>%
  summarise(mean(slope))

modelling_results_migration_month %>%
  ggplot(., aes(x = Migration, y = slope, fill = Migration)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  theme(legend.position = "none")

t.test(slope ~ Migration, data = modelling_results_migration_month, var.equal = TRUE)



