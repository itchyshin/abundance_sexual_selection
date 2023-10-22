library(sf)
library(tidyverse)
library(rnaturalearth)
library(raster)

grids <- st_read("Data/grid_5_degree.geojson")
tmin<-getData('worldclim',var="tmin",res=2.5)
tmax<-getData('worldclim',var="tmax",res=2.5)
prec<-getData('worldclim',var="prec",res=2.5)

mat <- mean((tmax+tmin)/2)/10
plot(mat)

prec_ann <- sum(prec)
plot(prec_ann)

grids$mat<-raster::extract(mat,cbind(grids$centroid_lat,grids$centroid_lng))
grids$prec_ann<-raster::extract(prec_ann,cbind(grids$centroid_lat,grids$centroid_lng))

st_write(grids,"Data/grid_5_degree_with_clim.geojson")
