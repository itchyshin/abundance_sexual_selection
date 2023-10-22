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

grids <- st_read("Data/grid_5_degree_with_clim.geojson")
 
data.clean <- grids %>%
  rename(grid_id=ID) %>%
  rename(Latitude = centroid_lng) %>%
  rename(Longitude = centroid_lat) %>%
  left_join(data.1) %>%
  group_by(ebird_COMMON_NAME, migration) %>%
  summarise(Lat = mean(Latitude))%>%
  filter(., migration == "migrant")

sum(data.clean$Lat > 0)
sum(data.clean$Lat < 0)
