load("C:/Users/Anis/Documents/Projet d'étude/Documents d'étude/tacsatEflalo_2.RData")
tacsatEflalo_2= dplyr::filter(tacsatEflalo_2, SI_STATE == 1)


install.packages("rnaturalearth")
install.packages("sf")
install.packages("mapdata")
install.packages("ggplot2")
install.packages("ggspatial")
install.packages("rgeos")
library(rnaturalearth)
library(sf)
library(mapdata)
library(ggplot2)
library(ggspatial)
library(rgeos)

# setwd("C:/R_projects/phd_zfh_baptiste_alglave/examples/map/scientific_data")
setwd("C:/Users/Anis/Documents/Projet d'étude/scientific_data")
# C:\Users\Anis\Documents\Projet d'étude\scientific_data.zip

load("res/ptsurvey_wgs84_full.RData")
load("map/strata_ZN.RData")
load("map/strata_ZCC.RData")
load("map/strata_ZCL.RData")
load("map/strata_ZS.RData")

mapBase <- ne_countries(scale = "medium", returnclass = "sf")
# mapBase <- map("worldHires", fill = T, plot = F)
# mapBase <- st_as_sf(mapBase) # %>% filter(ID %in% c("France","Spain"))

# for (j in 23:50)
# {
plot <- ggplot() +
  geom_sf(data = strata_ZN,
          alpha = 0.35,
          size = 0.75,
          color="#7FC97F",
          fill ="#7FC97F")+
  geom_sf(data = strata_ZS,
          alpha = 0.35,
          size = 0.75,
          color="#FFED6F",
          fill="#FFED6F")+
  geom_sf(data = strata_ZCC,
          alpha = 0.35,
          size = 0.75,
          color="white",
          fill="white")+
  geom_sf(data = strata_ZCL,
          alpha = 0.35,
          size = 0.75,
          color="#EF8A62",
          fill="#EF8A62")+
  geom_sf(data = mapBase)+
  
  geom_point(data = tacsatEflalo_2,
             aes(x=SI_LONG,y=SI_LATI,col = tacsatEflalo_2[,35]))+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  theme(aspect.ratio = 1,legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "skyblue"))+
  scale_color_gradient(low = "white", high = "red")+
  xlab("")+ylab("")+
  annotation_scale(location = "bl", line_width = .5) +
  annotation_north_arrow(location = "tl", height = unit(0.7, "cm"), width = unit(0.7, "cm"))

plot
# png(file="C:/Users/Anis/Documents/Projet d'étude/saving_plot[j].png")
# plot[j]
# dev.off()
# }

# Total

# Total<-sum(tacsatEflalo_2[,j],j=23:50,na.rm = FALSE)
# View(Total)

# Total<-colSums(tacsatEflalo_2[,c(23:50)],na.rm = FALSE)
# View(Total)

colSums (tacsatEflalo_2[,c(23:50)], na.rm = TRUE, dims = 1)
Total<-rowSums (tacsatEflalo_2[,c(23:50)], na.rm = TRUE, dims = 1)
View(Total)

mapBase <- ne_countries(scale = "medium", returnclass = "sf")
# mapBase <- map("worldHires", fill = T, plot = F)
# mapBase <- st_as_sf(mapBase) # %>% filter(ID %in% c("France","Spain"))

# for (j in 23:50)
# {
plot <- ggplot() +
  # geom_sf(data = strata_ZN,
  #         alpha = 0.35,
  #         size = 0.75,
  #         color="#7FC97F",
  #         fill ="#7FC97F")+
  # geom_sf(data = strata_ZS,
  #         alpha = 0.35,
  #         size = 0.75,
  #         color="#FFED6F",
  #         fill="#FFED6F")+
  # geom_sf(data = strata_ZCC,
#         alpha = 0.35,
#         size = 0.75,
#         color="white",
#         fill="white")+
# geom_sf(data = strata_ZCL,
#         alpha = 0.35,
#         size = 0.75,
#         color="#EF8A62",
#         fill="#EF8A62")+
# geom_sf(data = mapBase)+

geom_point(data = tacsatEflalo_2,
           aes(x=SI_LONG,y=SI_LATI,col = Total))+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  theme(aspect.ratio = 1,legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "skyblue"))+
  scale_color_gradient(low = "white", high = "red")+
  xlab("")+ylab("")+
  annotation_scale(location = "bl", line_width = .5) +
  annotation_north_arrow(location = "tl", height = unit(0.7, "cm"), width = unit(0.7, "cm"))

plot


Mean<-rowMeans (tacsatEflalo_2[,c(23:50)], na.rm = TRUE, dims = 1)
View(Mean)
str(Mean)

mapBase <- ne_countries(scale = "medium", returnclass = "sf")
# mapBase <- map("worldHires", fill = T, plot = F)
# mapBase <- st_as_sf(mapBase) # %>% filter(ID %in% c("France","Spain"))

# for (j in 23:50)
# {

load("C:/Users/Anis/Documents/Projet d'étude/Documents d'étude/ICES_rect.RData")
ICES_rect
plot <- ggplot() +
  # geom_sf(data = strata_ZN,
  #         alpha = 0.35,
  #         size = 0.75,
  #         color="#7FC97F",
  #         fill ="#7FC97F")+
  # geom_sf(data = strata_ZS,
  #         alpha = 0.35,
  #         size = 0.75,
  #         color="#FFED6F",
  #         fill="#FFED6F")+
  # geom_sf(data = strata_ZCC,
#         alpha = 0.35,
#         size = 0.75,
#         color="white",
#         fill="white")+
# geom_sf(data = strata_ZCL,
#         alpha = 0.35,
#         size = 0.75,
#         color="#EF8A62",
#         fill="#EF8A62")+

geom_sf(data = ICES_rect, fill= "skyblue")+
  geom_sf(data = mapBase)+
  geom_point(data = tacsatEflalo_2,
             aes(x=SI_LONG,y=SI_LATI,col = Mean))+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  theme(aspect.ratio = 1,legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "skyblue"))+
  scale_color_gradient(low = "white", high = "red")+
  xlab("")+ylab("")+
  annotation_scale(location = "bl", line_width = .5) +
  annotation_north_arrow(location = "tl", height = unit(0.7, "cm"), width = unit(0.7, "cm"))

plot

summary(tacsatEflalo_2[,c(23:50)])

# load("C:/Users/Anis/Downloads/ICES_rect.RData")
# load("C:/Users/Anis/Documents/Projet d'étude/Documents d'étude/ICES_rect.RData")
# ICES_rect


# Nombre de pings par séquences de pêche : 
View(table(tacsatEflalo_2$FT_REF))

View(table(tacsatEflalo_2$SI_DATE))


# Moyenne et écart-type du nombre de pings par séquence de pêche
mean(table(tacsatEflalo_2$FT_REF))
sd(table(tacsatEflalo_2$FT_REF))


# barplot(tacsatEflalo_2$SI_DATE)

View(table(tacsatEflalo_2$FT_REF,tacsatEflalo_2$SI_DATE))

head(table(tacsatEflalo_2$FT_REF,tacsatEflalo_2$SI_DATE))
 
tacsatEflalo_2.1=data.frame(table(tacsatEflalo_2$FT_REF,tacsatEflalo_2$SI_DATE))

# Création du sous tableau extrait constitué uniquement des données dont la date est au 01/11/2018
tacsatEflalo_2.2=data.frame(tacsatEflalo_2$SI_DATE == "01/11/2018", tacsatEflalo_2$SI_LATI,tacsatEflalo_2$SI_LONG)

# Représentation des positions de tous les bateaux sur la journée du 01/11/2018 :
plot_3<-ggplot() +
  geom_sf(data = ICES_rect, fill= "skyblue")+
  geom_sf(data = mapBase)+
  geom_point(data = tacsatEflalo_2.2,
             aes(x = tacsatEflalo_2$SI_LONG, y = tacsatEflalo_2$SI_LATI, col = tacsatEflalo_2$VE_REF))+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)

plot_3  

# Utiliser la fonction over du package sp pour croiser les données gps avec les carrés statistique, et alors on saura dans quel carré stat, 
#  tombe chaque ping.

tacsatEflalo_2$SI_STATE

install.packages("raster")
library(raster)
grid_projection <- "+proj=longlat +datum=WGS84"
ptvms_wgs84 <- SpatialPointsDataFrame(coords=cbind(tacsatEflalo_2$SI_LONG,tacsatEflalo_2$SI_LATI),
                                      data=tacsatEflalo_2,
                                      proj4string=crs(grid_projection))

# intersect of grid and VMS data
ICES_rect_2=as(ICES_rect, "Spatial")

ptvms_wgs84_1 <- over(ptvms_wgs84, ICES_rect_2)

ptvms_wgs84_1
ptvms_wgs84_2 <- cbind(ptvms_wgs84@data, ptvms_wgs84_1)
ptvms_wgs84_2
###############################################################################################################################


ICES_rect

plot_3<-ggplot() +
   geom_sf(data = ICES_rect, fill= "skyblue")+
   geom_sf(data = mapBase)+
   geom_point(data = ptvms_wgs84_2,
                   +                aes(x = ptvms_wgs84_2$SI_LONG, y = ptvms_wgs84_2$SI_LATI, col = ptvms_wgs84_2$VE_REF))+
   coord_sf(xlim = c(-3,-2), ylim = c(46,46.5), expand = FALSE)
 
plot_3 


# Création du sous tableau extrait constitué uniquement des données dont la date est au 01/11/2018 et de VE_REF = 16428
# tacsatEflalo_2.3=data.frame(tacsatEflalo_2$SI_DATE == "01/11/2018" , tacsatEflalo_2$SI_LATI,tacsatEflalo_2$SI_LONG)
tacsatEflalo_2.3 = subset(tacsatEflalo_2, tacsatEflalo_2$VE_REF== "16428")
View(tacsatEflalo_2.3)

plot_5<-ggplot() +
  geom_sf(data = ICES_rect, fill= "skyblue")+
  geom_sf(data = mapBase)+
  geom_point(data = tacsatEflalo_2.3,
             +                aes(x = SI_LONG, y = SI_LATI))+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)

plot_5 

View(tacsatEflalo_2.3)
names(tacsatEflalo_2.3)
View(tacsatEflalo_2.2)

#  Avec ptvms_wgs84_2 : ça marche
plot_3<-ggplot() +
  geom_sf(data = ICES_rect, fill= "skyblue")+
  geom_sf(data = mapBase)+
  geom_point(data = ptvms_wgs84_2,
             aes(x = SI_LONG, y = SI_LATI, col = VE_REF))+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)

plot_3 

# 
ptvms_wgs84_2.2 = subset(ptvms_wgs84_2, ptvms_wgs84_2$VE_REF== "16428")

plot_4<-ggplot() +
  geom_sf(data = ICES_rect, fill= "skyblue")+
  geom_sf(data = mapBase)+
  geom_point(data = ptvms_wgs84_2.2,
             aes(x = SI_LONG, y = SI_LATI, col = VE_REF))+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)

plot_4

# Selection des 4 carrés stats où il évolue

plot_4<-ggplot() +
  geom_sf(data = ICES_rect, fill= "skyblue")+
  geom_sf(data = mapBase)+
  geom_point(data = ptvms_wgs84_2.2,
             aes(x = SI_LONG, y = SI_LATI, col = as.factor(LE_KG_Solea_solea)))+
  facet_wrap(.~SI_DATE)+
  coord_sf(xlim = c(-3,-1), ylim = c(45.5,46.5), expand = FALSE)

plot_4

# Borne gauche : [min si_long] et Borne droite : [max si_long] + 1
# Evolution d'un bateau sur un mois
# Over : répartition des pings sur un carré stat

?over

ptvms_wgs84 <- SpatialPointsDataFrame(coords=cbind(tacsatEflalo_2$SI_LONG,tacsatEflalo_2$SI_LATI),
                                      data=tacsatEflalo_2,
                                      proj4string=crs(grid_projection))

# intersect of grid and VMS data
ICES_rect_2=as(ICES_rect, "Spatial")
ptvms_wgs84_1 <- over(ptvms_wgs84, ICES_rect_2)


ptvms_wgs84_1
ptvms_wgs84_2 <- cbind(ptvms_wgs84@data, ptvms_wgs84_1)
library(dplyr)
ptvms_wgs84_2 %>% group_by(ICESNAME_2,SI_DATE,VE_REF) %>% dplyr::summarize(n = n())


# filter(ptvms_wgs84_2, ptvms_wgs84_2$ICESNAME == "23E7" )
# Représenter et résumer l'info :
boxplot() etc...

