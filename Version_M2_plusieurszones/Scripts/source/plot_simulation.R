###################
## Simulation plots
###################

## Simulation boxplots for comparison of models
#----------------------------------------------

load(paste0(results_file,"draft/paper_3/res/Results_full_multi_square.RData"))

# Filter simulations that have converged
Results_plot <- Results_2 %>%
  filter(converge == 0) %>% 
  filter(simu_type == "Unsampled Rectangles") %>%  # and the one that will be plotted
  filter(Model %in% c("Integrated model","Scientific model"))

# And those that will be plotted
Results_plot$obs <- NA
Results_plot$obs[which(Results_plot$aggreg_obs == T)] <- "Joint COS model"
Results_plot$obs[which(Results_plot$aggreg_obs == F)] <- "Two-step approach"
Results_plot$obs[which(Results_plot$Estimation_model == 2)] <- "Point-level model"
Results_plot$obs <- factor(Results_plot$obs,levels = c("Two-step approach","Joint COS model", "Point-level model"))

# 
Results_2$obs <- NA
Results_2$obs[which(Results_2$aggreg_obs == T)] <- "Joint COS model"
Results_2$obs[which(Results_2$aggreg_obs == F)] <- "Two-step approach"
Results_2$obs[which(Results_2$Estimation_model == 2)] <- "Point-level model"
Results_2$obs <- factor(Results_2$obs,levels = c("Two-step approach","Joint COS model", "Point-level model"))

Results_sci_plot <- Results_plot %>% 
  filter(obs == "Point-level model")

# MSPE
mspe_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = mspe,
                   fill = obs),outliers = F)+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  xlab("")+ylab("MSPE")+
  ylim(0,NA)

# Beta values
beta_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = beta1_est,
                   fill = obs),outliers = F)+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  geom_hline(yintercept = 2,color="red",linetype="dashed")+
  xlab("")+ylab(TeX("$\\beta$"))

# Range values
range_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = Range,
                   fill = obs),outliers = F)+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  geom_hline(yintercept = 0.6,color="red",linetype="dashed")+
  xlab("")+ylim(0,2.5)

# Variance values
Results_2$obs <- factor(Results_2$obs,levels = c("Two-step approach","Joint COS model"))
variance_plot <- ggplot()+
  geom_boxplot(data = Results_2 %>% filter(!is.na(obs)),
               aes(x = obs,
                   y = sigma_com_est,
                   fill = obs),outliers = F)+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  geom_hline(yintercept = 1,color="red",linetype="dashed")+
  xlab("")+
  scale_fill_manual(breaks = c("Two-step approach","Joint COS model"), 
                    values=c("#F8766D", "#00BA38"))+
  ylim(0,2)+
  xlab("")+ylab(TeX("$\\sigma$"))

# zero-inflation parameter
q1_plot <- ggplot()+
  geom_boxplot(data = Results_2 %>% filter(!is.na(obs)),
               aes(x = obs,
                   y = q1_com_est,
                   fill = obs),outliers = F)+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  geom_hline(yintercept = -1,color="red",linetype="dashed")+
  xlab("")+
  scale_fill_manual(breaks = c("Two-step approach","Joint COS model"), 
                    values=c("#F8766D", "#00BA38"))+
  xlab("")+ylab(TeX("$\\xi$"))

legend <- ggpubr::as_ggplot(ggpubr::get_legend(mspe_plot,position="right"))

boxplot_1 <- plot_grid(mspe_plot,NULL,beta_plot,legend,align = "hv",nrow = 1,rel_widths = c(0.3,0.1,0.3,0.3))

ggsave("../../paper_reallocation/images/boxplot_1.png",height = 4*1.25,width = 8*1.25)

legend2 <- ggpubr::as_ggplot(ggpubr::get_legend(variance_plot,position="right"))

boxplot_3 <- plot_grid(variance_plot,NULL,q1_plot,legend2,nrow = 1,rel_widths = c(1,0.1,1,0.6))

ggsave("../../paper_reallocation/images/boxplot_3.png",height = 4,width = 8)


# boxplot_1_add <- plot_grid(mspe_plot,beta_plot,legend,
#                        range_plot,variance_plot,q1_plot,align = "hv",ncol = 3)
# 
# ggsave("../../paper_reallocation/images/boxplot_1_add.png",height = 8,width = 8)


## Simulation prediction maps for comparison of models
#-----------------------------------------------------

## Simulated data
load("draft/paper_3/res/list_input.RData")
loc_x<-list_input$loc_x
S_x<-list_input$S_x
simu_df <- data.frame(loc_x,S_x=S_x)

# rectangles
ICES_rect <- st_read("/media/balglave/Elements/backup_phd/phd_zfh_baptiste_alglave/data/raw/shapefile/ices_rect_stats/ICES_Statistical_Rectangles_Eco.shp")
st_crs(ICES_rect) <- grid_projection
ICES_rect_2 <- ICES_rect %>%
  filter(Ecoregion == "Bay of Biscay and the Iberian Coast")

## Model outputs
load("draft/paper_3/res/sci_df.RData")
load("draft/paper_3/res/int_Yi_df.RData")
load("draft/paper_3/res/int_Dj_df.RData")
load("draft/paper_3/res/com_Yi_df.RData")
load("draft/paper_3/res/list_input.RData")

loc_x<-list_input$loc_x
S_x<-list_input$S_x
c_com_x<-list_input$c_com_x
y_i2<-list_input$y_i2
index_i<-list_input$index_i
y_sci_i<-list_input$y_sci_i
index_sci_i<-list_input$index_sci_i
boats_i<-list_input$boats_i
cov_x_2<-list_input$cov_x_2
ref_level<-list_input$ref_level
spde<-list_input$spde
mesh<-list_input$mesh
n_cells<-list_input$n_cells
sampling<-list_input$sampling

## Simulated data
sci.data_df <- data.frame(cell = index_sci_i,
                          y_sci = y_sci_i) %>%
  inner_join(loc_x[,c("x","y","cell","ICESNAME")])

com.data_df <- data.frame(cell = index_i,
                          y_com = y_i2) %>%
  inner_join(loc_x[,c("x","y","cell","ICESNAME")]) %>%
  group_by(x,y,cell,ICESNAME) %>%
  dplyr::summarise(n=n())

max_val <- log(max(simu_df$S_x,
               sci_df$S_x,
               int_Yi_df$S_x,
               int_Dj_df$S_x))

data_plot <- ggplot()+
  geom_point(data=simu_df,aes(x=x,y=y,col=log(S_x)),shape=15,size=2,alpha=0.25)+
  geom_point(data=simu_df[which(simu_df$ICESNAME %in% com.data_df$ICESNAME),],
             aes(x=x,y=y,col=S_x),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+ # , limits=c(0,max_val+0.01*max_val))+
  geom_point(data=sci.data_df,aes(x=x,y=y),col="red",size=1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank()
  )+
  geom_sf(data=ICES_rect_2,alpha=0)+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

ggsave("../../paper_reallocation/images/data_plot.png",width = 6, height = 6)


simu_plot <- ggplot(simu_df)+
  geom_point(aes(x=x,y=y,col=log(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Simulated latent field")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank())+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

sci_plot <- ggplot(sci_df)+
  geom_point(aes(x=x,y=y,col=log(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Point-level model")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank())+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

int_Yi_plot <- ggplot(int_Yi_df)+
  geom_point(aes(x=x,y=y,col=log(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Two-step approach")+ # integrated model
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank())+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

int_Dj_plot <- ggplot(int_Dj_df)+
  geom_point(aes(x=x,y=y,col=log(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Joint COS model")+ # ,subtitle = "integrated",
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank())+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

map_plot <- plot_grid(
  simu_plot,sci_plot,
  int_Yi_plot,int_Dj_plot,
  nrow = 2,
  ncol = 2,
  align="hv",
  labels = c("A","B","C","D"))

ggsave("../../paper_reallocation/images/Map_multi_square.png",width = 9, height = 9)


## N zone simulations
#--------------------
load(paste0(results_file,"results/n_zone/Results_n_zone.RData"))
Results_plot <- Results_2 %>%
  filter(converge == 0) %>%
  # filter(simu_type == "Unsampled Rectangles") %>%  # and the one that will be plotted
  filter(Model %in% c("Integrated model","Point-level model"))

# And those that will be plotted
Results_plot$obs <- NA
Results_plot$obs[which(Results_plot$aggreg_obs == T)] <- "Joint COS model"
Results_plot$obs[which(Results_plot$aggreg_obs == F)] <- "Two-step approach"
Results_plot$obs[which(Results_plot$Estimation_model == 2)] <- "Point-level model"
Results_plot$obs <- factor(Results_plot$obs,levels = c("Two-step approach","Joint COS model","Point-level model"))

Results_plot <- full_join(Results_plot,Results_sci_plot)
Results_plot$n_zone[which(is.na(Results_plot$n_zone))] <- " "

Results_plot$n_zone <- factor(Results_plot$n_zone,levels = c("1","3","5"," "))


# MSPE
mspe_nzone_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = as.factor(n_zone),
                   y = mspe,
                   fill = obs),outliers = F)+
  scale_fill_manual(breaks = c("Two-step approach","Joint COS model","Point-level model"), 
                    values=c("#F8766D", "#00BA38", "#619CFF"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        legend.position = "none",plot.title = element_text(hjust = 0.5,face = "bold"))+
  xlab("")+ylab("MSPE")+
  ylim(0,NA)+
  ggtitle("Several zone scenarios",subtitle = " ")


# Covariate effect
beta_nzone_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = as.factor(n_zone),
                   y = beta1_est,
                   fill = obs),outliers = F)+
  scale_fill_manual(breaks = c("Two-step approach","Joint COS model","Point-level model"), 
                    values=c("#F8766D", "#00BA38", "#619CFF"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        legend.position = "none")+
  xlab("")+ylab(TeX("\\beta"))+
  ylim(0,NA)+
  geom_hline(yintercept = 2,color="red",linetype="dashed")


## zero inflation simulations
#----------------------------
load(paste0(results_file,"results/q1/Results_q1.RData"))
Results_plot <- Results_2 %>%
  filter(converge == 0) %>%
  # filter(simu_type == "Unsampled Rectangles") %>%  # and the one that will be plotted
  filter(Model %in% c("Integrated model","Point-level model"))

# And those that will be plotted
Results_plot$obs <- NA
Results_plot$obs[which(Results_plot$aggreg_obs == T)] <- "Joint COS model"
Results_plot$obs[which(Results_plot$aggreg_obs == F)] <- "Two-step approach"
Results_plot$obs[which(Results_plot$Estimation_model == 2)] <- "Point-level model"
Results_plot$obs <- factor(Results_plot$obs,levels = c("Two-step approach","Joint COS model","Point-level model"))


# Add proportion of zero in data frame
# probability of getting a zero: exp(-exp(xi)*exp(intercept))
# xi=0 : nearly 0% of zero
# xi=-1 : 7% of zero
# xi=-2 : 37% of zero
# xi=-3 : 70% of zero


Results_plot$zero_prop <- NA
Results_plot$zero_prop[which(Results_plot$q1 == 0)] <- "0 %"
Results_plot$zero_prop[which(Results_plot$q1 == -1)] <- "7 %"
Results_plot$zero_prop[which(Results_plot$q1 == -2)] <- "37 %"
Results_plot$zero_prop[which(Results_plot$q1 == -3)] <- "70 %"

Results_plot <- full_join(Results_plot,Results_sci_plot)
Results_plot$zero_prop[which(is.na(Results_plot$zero_prop))] <- " "

Results_plot$zero_prop <- factor(Results_plot$zero_prop,levels = c("0 %","7 %","37 %","70 %"," "))

# MSPE
mspe_0s_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = zero_prop,
                   y = mspe,
                   fill = obs),outliers = F)+
  scale_fill_manual(breaks = c("Two-step approach","Joint COS model","Point-level model"), 
                    values=c("#F8766D", "#00BA38", "#619CFF"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,face = "bold"))+
  xlab("")+ylab("MSPE")+
  ylim(0,NA)+
  ggtitle("Zero-inflation scenarios",subtitle = " ")


# Covariate effect
beta_0s_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = zero_prop,
                   y = beta1_est,
                   fill = obs),outliers = F)+
  scale_fill_manual(breaks = c("Two-step approach","Joint COS model","Point-level model"), 
                    values=c("#F8766D", "#00BA38", "#619CFF"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        legend.position = "none")+
  xlab("")+ylab(TeX("\\beta"))+
  ylim(0,NA)+
  geom_hline(yintercept = 2,color="red",linetype="dashed")

# Proportion of zeros
prop0_plot <- ggplot()+
  geom_boxplot(data = Results_2,
               aes(x = as.factor(q1),
                   y = prop_zero,
                   fill = obs),outliers = F)+
  scale_fill_manual(breaks = c("Two-step approach","Joint COS model","Point-level model"), 
                    values=c("#F8766D", "#00BA38", "#619CFF"))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1)+
  xlab("")+ylab("MSPE")+
  ylim(0,NA)

add_scenratios_plot <- cowplot::plot_grid(mspe_nzone_plot,mspe_0s_plot,
                                          beta_nzone_plot,beta_0s_plot,ncol=2,align = "hv")

legend <- ggpubr::as_ggplot(ggpubr::get_legend(mspe_nzone_plot+
                                                 theme(legend.position = "right",
                                                       legend.text = element_text(size = 12)),position="right"))

add_scenratios_plot <- cowplot::plot_grid(add_scenratios_plot,legend,ncol=1,rel_heights = c(0.9,0.1))

ggsave("../../paper_reallocation/images/boxplot_2.png",height = 8,width = 8)

## Table of copnvergence
Results_conv <- Results_2
Results_conv$one <- 1

# add columns
Results_conv$obs <- NA
Results_conv$obs[which(Results_conv$aggreg_obs == T)] <- "Joint COS model"
Results_conv$obs[which(Results_conv$aggreg_obs == F)] <- "Two-step approach"
Results_conv$obs[which(Results_conv$Estimation_model == 2)] <- "Point-level model"
Results_conv$obs <- factor(Results_conv$obs,levels = c("Two-step approach","Joint COS model","Point-level model"))

Results_conv$zero_prop <- NA
Results_conv$zero_prop[which(Results_conv$q1 == 0)] <- "0 %"
Results_conv$zero_prop[which(Results_conv$q1 == -1)] <- "7 %"
Results_conv$zero_prop[which(Results_conv$q1 == -2)] <- "37 %"
Results_conv$zero_prop[which(Results_conv$q1 == -3)] <- "70 %"

Results_conv$zero_prop <- factor(Results_conv$zero_prop,levels = c("0 %","7 %","37 %","70 %"))

table_df <- doBy::summaryBy(converge+one~
                  obs+
                  zero_prop,
                data=Results_conv,
                FUN=sum) %>%
  mutate(perc_convergence = round((1 - converge.sum / one.sum)*100,digits = 3)) %>%
  dplyr::select(-converge.sum,-one.sum)

table_df %>%
  gt() %>%
  tab_options(
    summary_row.background.color = "#ACEACE80",
    grand_summary_row.background.color = "#990000",
    row_group.background.color = "#FFEFDB80",
    heading.background.color = "#EFFBFC",
    column_labels.background.color = "#EFFBFC",
    stub.background.color = "#EFFBFC",
    table.font.color = "#323232",
    table_body.hlines.color = "#989898",
    table_body.border.top.color = "#989898",
    heading.border.bottom.color = "#989898",
    row_group.border.top.color = "#989898",
    row_group.border.bottom.style = "none",
    stub.border.style = "dashed",
    stub.border.color = "#989898",
    stub.border.width = "1px",
    summary_row.border.color = "#989898",
    table.width = "80%",
    column_labels.font.weight = "bold"
  )

knitr::kable(
  table_df,"latex")

