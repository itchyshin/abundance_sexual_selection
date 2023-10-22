library(sf)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(dplyr)
library(MCMCglmm)
library(lme4)

modelling_results_month_migration <- readRDS("Data/bayesian_results_months_migration.RDS")
grids_to_select <- modelling_results_month_migration %>%
  group_by(grid, Migration) %>%
  summarize(N=n()) %>%
  dplyr::filter(N==12) %>%
  mutate(include="Yes")

modelling_results_month_migration %>%
  left_join(grids_to_select) %>%
  dplyr::filter(complete.cases(include)) %>%
  #dplyr::filter(Migration=="Resident") %>%
  ggplot(., aes(x=factor(Month, levels=c("January", "February",
                                         "March", "April", "May",
                                         "June", "July", "August",
                                         "September", "October", 
                                         "November", "December")), 
                y=estimate_slope, group=grid))+
  #geom_point(color="gray30", size=0.002)+
  geom_polygon(color="gray40", size=0.25, fill = NA) +
  geom_segment(aes(x="January", xend="December", y=0,yend=0), color = "red", linetype = "dashed", size = 1.5)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_line(colour="gray80", size = (0.8)))+
  facet_wrap(~Migration, ncol=2, scales="fixed")+
  xlab("")+
  ylab("Slope Values of Relative Abundance \nRegressed on Sexual Dichromatism")+
  geom_smooth(aes(color=factor(Migration), group=Migration))+
  labs(color="Migration") +
  theme(legend.position="none") +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_text(face = "bold")) +
  coord_polar()

ggsave(file='Final_Temporal_Graph.png', width = 15, height=8, dpi=600)
