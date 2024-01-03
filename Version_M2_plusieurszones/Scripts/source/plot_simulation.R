###################
## Simulation plots
###################

## N zone simulations
#--------------------
load(paste0(results_file,"results/n_zone/Results_n_zone.RData"))
Results_plot <- Results_2 %>%
  filter(converge == 0) %>%
  # filter(simu_type == "Unsampled Rectangles") %>%  # and the one that will be plotted
  filter(Model %in% c("Integrated model","Scientific model"))

# And those that will be plotted
Results_plot$obs <- NA
Results_plot$obs[which(Results_plot$aggreg_obs == T)] <- "COS model"
Results_plot$obs[which(Results_plot$aggreg_obs == F)] <- "GeoCatch model"
Results_plot$obs[which(Results_plot$Estimation_model == 2)] <- "Scientific model"
Results_plot$obs <- factor(Results_plot$obs,levels = c("GeoCatch model","COS model","Scientific model"))

# MSPE
mspe_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = mspe,
                   fill = as.factor(n_zone)))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  xlab("")+ylab("MSPE")+
  ylim(0,NA)


# Covariate effect
beta_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = beta1_est,
                   fill = as.factor(n_zone)))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1)+
  xlab("")+ylab("MSPE")+
  ylim(0,NA)


## zero inflation simulations
#----------------------------
load(paste0(results_file,"results/n_zone/Results_n_zone.RData"))
Results_plot <- Results_2 %>%
  # filter(converge == 0) %>%
  # filter(simu_type == "Unsampled Rectangles") %>%  # and the one that will be plotted
  filter(Model %in% c("Integrated model","Scientific model"))

# And those that will be plotted
Results_plot$obs <- NA
Results_plot$obs[which(Results_plot$aggreg_obs == T)] <- "COS model"
Results_plot$obs[which(Results_plot$aggreg_obs == F)] <- "GeoCatch model"
Results_plot$obs[which(Results_plot$Estimation_model == 2)] <- "Scientific model"
Results_plot$obs <- factor(Results_plot$obs,levels = c("GeoCatch model","COS model","Scientific model"))

# MSPE
mspe_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = mspe,
                   fill = as.factor(q1)))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1)+
  xlab("")+ylab("MSPE")+
  ylim(0,NA)


# Covariate effect
beta_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = beta1_est,
                   fill = as.factor(q1)))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1)+
  xlab("")+ylab("MSPE")+
  ylim(0,NA)



#-------------------------------------------------------------------------------------------------------------------------


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
Results_plot$obs[which(Results_plot$aggreg_obs == T)] <- "COS model"
Results_plot$obs[which(Results_plot$aggreg_obs == F)] <- "GeoCatch model"
Results_plot$obs[which(Results_plot$Estimation_model == 2)] <- "Scientific model"
Results_plot$obs <- factor(Results_plot$obs,levels = c("GeoCatch model","COS model","Scientific model"))

# MSPE
mspe_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = mspe,
                   fill = obs))+
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
                   fill = obs))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  geom_hline(yintercept = 2,color="red",linetype="dashed")+
  xlab("")+ylab(TeX("$\\beta \\, _{S}$"))

# Range values
range_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = Range,
                   fill = obs))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  geom_hline(yintercept = 0.6,color="red",linetype="dashed")+
  xlab("")+ylim(0,2.5)

# Variance values
variance_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = sigma_com_est,
                   fill = obs))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  geom_hline(yintercept = 1,color="red",linetype="dashed")+
  xlab("")+
  scale_fill_manual(breaks = c("Yri","Dj","Scientific model"), 
                    values=c("#F8766D", "#00BA38", "#619CFF"))+
  ylim(0,NA)

# zero-inflation parameter
q1_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = q1_com_est,
                   fill = obs))+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  geom_hline(yintercept = -1,color="red",linetype="dashed")+
  xlab("")+
  scale_fill_manual(breaks = c("Yri","Dj","Scientific model"), 
                    values=c("#F8766D", "#00BA38", "#619CFF"))

legend <- ggpubr::as_ggplot(ggpubr::get_legend(mspe_plot,position="right"))

boxplot_1 <- plot_grid(mspe_plot,beta_plot,legend,align = "hv",ncol = 3)

ggsave("../../paper_reallocation/images/boxplot_1.png",height = 4,width = 8)

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

## Model outputs
load("draft/paper_3/res/sci_df.RData")
load("draft/paper_3/res/int_Yi_df.RData")
load("draft/paper_3/res/int_Dj_df.RData")

max_val <- max(simu_df$S_x/sum(simu_df$S_x),
               sci_df$S_x/sum(sci_df$S_x),
               int_Yi_df$S_x/sum(int_Yi_df$S_x),
               com_Yi_df$S_x/sum(com_Yi_df$S_x),
               int_Dj_df$S_x/sum(int_Dj_df$S_x),
               com_Dj_df$S_x/sum(com_Dj_df$S_x))

simu_plot <- ggplot(simu_df)+
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
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
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Scientific model")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank())+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

int_Yi_plot <- ggplot(int_Yi_df)+
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("GeoCatch model")+ # integrated model
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank())+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

int_Dj_plot <- ggplot(int_Dj_df)+
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("COS model")+ # ,subtitle = "integrated",
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
