model <- lmer(log10(mean_abund) ~ log_mass*scale(sexual_dichromatism) + (1+scale(sexual_dichromatism)|grid_id), data = data_joined2)
summary(model)
hist(ranef(model)$grid_id[,2])

library(sf)
library(tidyverse)
library(rnaturalearth)

grids <- st_read("Data/grid_5_degree.geojson")

ranef_df1 <- as.data.frame(ranef(model)$grid_id) %>%
  rownames_to_column(var="grid_id")


idk <- grids %>%
  rename(grid_id=ID) %>%
  mutate(grid_id=as.character(as.integer(grid_id))) %>%
  left_join(ranef_df1)


# read in and download spatial data
# download some natural earth data
countries <- ne_countries(returnclass = "sf") %>%
  st_geometry()

graticules <- ne_download(type = "graticules_15", category = "physical",
                          returnclass = "sf") %>%
  st_geometry()

bb <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf") %>%
  st_geometry()


ggplot()+
  geom_sf(data = bb, col = "grey20", fill = "transparent")+
  geom_sf(data= idk, aes(fill=`scale(sexual_dichromatism)`))+
  geom_sf(data = countries, fill = "grey80", col = "black", lwd = 0.3, alpha=0.4)+
  theme_minimal()+
  scale_fill_gradient2()+
  theme(axis.text = element_blank())+
  theme(legend.position="bottom")


ggplot(idk, aes(x=centroid_lat, y=`scale(sexual_dichromatism)`))+
  geom_point()+
  geom_smooth(method="lm")

