library(sf)
library(maps)
library(mapdata)
library(ggspatial)

plot(st_geometry(points_sf_4))

# world <- ne_countries(scale = "medium", returnclass = "sf")
mapBase <- map("worldHires", fill = T, plot = F)
mapBase <- st_as_sf(mapBase) %>% filter(ID %in% c("France","Spain"))

VMS_data <- points_sf %>%
  filter(month %in% c(10:12)) %>%
  filter(LE_MET_level6 %in% c("OTB_DEF_>=70_0","OTB_CEP_>=70_0","OTT_DEF_>=70_0")) %>%
  mutate(MET = str_sub(LE_MET_level6,1,3))


plot_figure <- ggplot()+
  theme_bw()+
  geom_sf(data=VMS_data[which(VMS_data$MET == "OTT"),],
          col="grey",alpha=0.05,shape=16,size=1)+
  geom_sf(data=VMS_data[which(VMS_data$MET == "OTB"),],
          col="bisque",alpha=0.1,shape=16,size=1)+
  geom_sf(data=ObsMer_sf[which(ObsMer_sf$foCatEu6 == "OTT_DEF_>=70_0"),],
          aes(size = CPUE),col = "darkgrey")+
  geom_sf(data=ObsMer_sf[which(ObsMer_sf$foCatEu6 != "OTT_DEF_>=70_0"),],
          aes(size = CPUE),col = "orange")+
  geom_sf(data = st_as_sf(survey_layer_3),col="lightblue1",alpha=0,size=2)+
  geom_point(data=ptsurvey_wgs84_sf_4,
             aes(x = x, y = y,size = CatchWgt_spp),
             col = "skyblue")+
  theme(legend.position = "none")+
  geom_sf(data = ICES_rect_2,alpha=0)+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-5.5,0.25), ylim = c(43,48), expand = FALSE)+
  annotation_scale(location = "bl", line_width = .5) +
  annotation_north_arrow(location = "tl", height = unit(0.7, "cm"), width = unit(0.7, "cm"))+
  xlab("")+ylab("")

# x11()
# plot_figure

ggsave("C:/R_projects/phd_manuscript/Introduction/figures/several-datasources.png",
       width = 1.25*5,height = 7.5)

library(cowplot)
library(ggpubr)

toto <- VMS_data %>%
  group_by(ICESNAME_2) %>%
  dplyr::summarise(land =sum(LE_KG_Solea_solea))
  
toto_2 <- st_join(ICES_rect_2,toto,left=T) %>%
  filter(!is.na(land)) %>%
  mutate(land = land)

vms_plot <- ggplot()+
  theme_bw()+
  geom_sf(data=VMS_data[which(VMS_data$MET == "OTB"),],
          col="darkgrey",shape=16,size=0.1,alpha=0.25)+
  scale_fill_distiller(palette="Spectral")+
  geom_sf(data = ICES_rect_2,alpha=0)+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-5.5,0.25), ylim = c(43,48), expand = FALSE)+
  theme(legend.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5,size = 10),
        plot.title = element_text(face = "bold",hjust = 0.5),
        legend.position = "none")+
  ggtitle("VMS data",subtitle = "Estimated fishing locations")

logbooks_plot <- ggplot()+
  theme_bw()+
  geom_sf(data = ICES_rect_2,alpha=0)+
  geom_sf(data=toto_2,aes(fill=land))+
  scale_fill_distiller(palette = "Spectral")+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-5.5,0.25), ylim = c(43,48), expand = FALSE)+
  theme(legend.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5,size = 10),
        plot.title = element_text(face = "bold",hjust = 0.5),
        legend.position = "none")+
  ggtitle("Logbook data",subtitle = "Resolution: 'statistical rectangle x day x fishing trip x gear'")

legend_logbooks_plot <- as_ggplot(cowplot::get_legend(logbooks_plot+theme(legend.position="bottom",legend.title = element_blank())))

vms_logbook_plot <- ggplot()+
  theme_bw()+
  geom_sf(data=VMS_data[which(VMS_data$MET == "OTB"),],
          aes(col=log(LE_KG_Solea_solea+1),size=0.5),
          alpha=0.25,shape=16,size=0.1)+
  scale_color_distiller(palette = "Spectral")+
  geom_sf(data = ICES_rect_2,alpha=0)+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-5.5,0.25), ylim = c(43,48), expand = FALSE)+
  theme(legend.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5,size = 10),
        plot.title = element_text(face = "bold",hjust = 0.5),
        legend.position = "none")+
  ggtitle("Combined 'VMS x logbook' data")

legend_vms_logbook_plot <- as_ggplot(cowplot::get_legend(vms_logbook_plot+theme(legend.position="bottom",legend.title = element_blank())))

test <- cowplot::plot_grid(vms_plot,logbooks_plot,vms_logbook_plot,
                   NULL,legend_logbooks_plot,legend_vms_logbook_plot,
                   nrow = 2,align = "h",rel_heights = c(1,0.2))

ggsave("C:/R_projects/phd_manuscript/Introduction/figures/vms-logbook-combi.png",
       width = 12,height = 6)
