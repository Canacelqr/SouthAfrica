library(tidyverse)
library(rgdal)
library(sf)
library(foreach)
library(doParallel)
registerDoParallel(2)
library(fasterize)
library(sp)
library(raster)
library(vegan)
library(rasterVis)
library(rgdal)
library(rgeos)
library(picante)
library(betapart)
library(CommEcol)
library(leaflet)
##Import data
datadir="C:/Users/qliao/Desktop/southafrica/"
files=data.frame(
  path=list.files(datadir, recursive=T, pattern="shp$"),stringsAsFactors = F)%>%
  mutate(file=basename(path),
         species=gsub(".shp","",file),
         family=gsub("[/].*$","",path))
all_species = foreach(i=1:nrow(files),.combine=rbind,.packages = c("dplyr","sf")) %dopar% {
  sp=read_sf(file.path(datadir,files$path[i]))%>%
    select(-1)%>%
    mutate(family=files$family[i],species=files$species[i]) %>% st_set_crs(4326)
  #return(st_crs(sp))
  return(sp)
}
# Resolution for the species richness map
cfr_bbox=st_bbox(all_species)
r <- raster(
  xmn=cfr_bbox$xmin, 
  xmx=cfr_bbox$xmax, 
  ymn=cfr_bbox$ymin, 
  ymx=cfr_bbox$ymax,
  res=0.1)

# Resolution for the beta diversity map
cfr_bbox=st_bbox(all_species)
r2 <- raster(
  xmn=cfr_bbox$xmin, 
  xmx=cfr_bbox$xmax, 
  ymn=cfr_bbox$ymin, 
  ymx=cfr_bbox$ymax,
  res=0.1)

# raster_alphadiversity<- fasterize(all_species, r, background = 0, by = "species")
# all_alphadiversity=as.data.frame(raster_alphadiversity)
# raster_alphadiversity1<- fasterize(all_species, r, background = 0, by = NULL,fun="sum")

##Create Function of beta diversity
betagrid2<-function(data, radius, phylotree, phylobeta=F, index="sorensen"){
  mean_turnover<-numeric(length(data[,1]))
  mean_nestedness<-numeric(length(data[,1]))
  mean_beta<-numeric(length(data[,1]))
  for(i in 1:length(data[,1])){
    adj<-select.window(xf=data[i,1], yf=data[i,2], radius, xydata=data)[,-c(1,2)]
    res<-beta.pair(as.matrix(adj), index.family=index)
    mean_turnover[i]<-mean(as.matrix(res[[1]])[2:length(as.matrix(res[[1]])[,1]),1],na.rm=TRUE)
    mean_nestedness[i]<-mean(as.matrix(res[[2]])[2:length(as.matrix(res[[2]])[,1]),1],na.rm=TRUE)
    mean_beta[i]<-mean(as.matrix(res[[3]])[2:length(as.matrix(res[[3]])[,1]),1],na.rm=TRUE)
  }
  return(data.frame(cell=row.names(data), mean_turnover, mean_nestedness, mean_beta))
}

##Calculate beta diversity
#process data
raster_alphadiversity3<- fasterize(all_species, r2, background = 0, by = "species")
all_alphadiversity2=as.data.frame(raster_alphadiversity3)
raster_alphadiversity4<- fasterize(all_species, r2, background = 0, by = NULL,fun="sum")

temp <- data.matrix(all_alphadiversity2)
shape<- cbind.data.frame(coordinates(raster_alphadiversity3),temp)
shape <- data.matrix(shape)
row_name <- rep(1:length(shape[,1]),1)
row.names(shape) <- row_name

# #**Result**
# ##**Map of Sum Number Species**
# raster_spatial_richness<- fasterize(all_species, r, background = 0, fun = "sum", by = NULL)
# plot(raster_spatial_richness)

# ##**Map of Spatial Richness**
# spatialrichness<- all_alphadiversity %>%specnumber()
# alpha_diversity=r
# values(alpha_diversity)<-spatialrichness
# mask=calc(raster_alphadiversity1, function(x)ifelse(x>0,1,NA))
# alpha_diversity2=mask(alpha_diversity,mask)
# gplot(alpha_diversity2)+
#   geom_tile(aes(fill=value))+scale_fill_gradientn(colours=c("brown","red","yellow","green","dark green"))

##**Map of Spatial Pattern of Beta Diversity**
results <- betagrid2(data=shape, radius=0.05)
library(rgeos)
shape$betadiv <- results[,4]
beta_diversity=r2
crs(beta_diversity)<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 
values(beta_diversity)<-shape$betadiv
my.colors = colorRampPalette(c("pink","red","gold","yellow","darkgreen","green","blue","darkblue","violet","black"))
plot(beta_diversity, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)

##**Remove the Boundary**
library(sp)
buffer_beta_diversity<- rasterToPolygons(beta_diversity>-Inf,dissolve = T, fun = function(x) !is.na(x), n=4, na.rm=TRUE)
pc <- spTransform( buffer_beta_diversity, CRS( "+init=epsg:4326" ) )
crs(buffer_beta_diversity)<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
buffer_beta_diversity<-gBuffer(pc, byid=FALSE, id=NULL, width=-0.1, quadsegs=5, capStyle="ROUND")
crop_betadiversity<- mask(beta_diversity, buffer_beta_diversity)
my.colors = colorRampPalette(c("brown","red","yellow","green","dark green"))
plot(crop_betadiversity,col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)

##**Map of Species Turnover**
shape$turnover <- results[,2]
turn_over=r2
values(turn_over)<-shape$turnover
my.colors = colorRampPalette(c("brown","red","yellow","green","dark green"))
plot(turn_over, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)

##**Map of Species Nestedness**
shape$nestedness <- results[,3]
nested_ness=r2
values(nested_ness)<-shape$nestedness
my.colors = colorRampPalette(c("brown","red","yellow","green","dark green"))
plot(nested_ness, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)