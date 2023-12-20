###################
## Simulation plots
###################

## Simulation boxplots
#---------------------

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


## Simulation maps
#-----------------

load(file=paste0(results_file,"results/real_simu/simu.data.RData"))
S_x <- simu.data$S_x
index_i <- simu.data$index_i
y_sci_i <- simu.data$y_sci_i
y_i <- simu.data$y_i
boats_i <- simu.data$boats_i
cov_x <- simu.data$cov_x
index_sci_i <- simu.data$index_sci_i

load(file=paste0(results_file,"results/real_simu/no.realloc_int_df.RData"))
res_no.realloc_int <- fit_IM_res
load(file=paste0(results_file,"results/real_simu/sci_df.RData"))
res_sci <- fit_IM_res
load(file=paste0(results_file,"results/real_simu/no.realloc_com_df.RData"))
res_no.realloc_com <- fit_IM_res
load(file=paste0(results_file,"results/real_simu/realloc_int_df.RData"))
res_realloc_int <- fit_IM_res
load(file=paste0(results_file,"results/real_simu/realloc_com_df.RData"))
res_realloc_com <- fit_IM_res

samp_sci <- data.frame(cell = index_sci_i) %>%
  dplyr::group_by(cell) %>%
  dplyr::summarise(n_sci=n()) %>%
  inner_join(loc_x[,c("cell","x","y")])

samp_com <- data.frame(cell = index_i) %>%
  dplyr::group_by(cell) %>%
  dplyr::summarise(n_com=n()) %>%
  inner_join(loc_x[,c("cell","x","y")])


#######
## Maps
#######

## Sampling locations/effort
library(viridis)
pal <- inferno(100)
sampling_plot <- ggplot()+
  geom_point(data=samp_com,aes(x=x,y=y,col=n_com),shape=15,size=2)+
  geom_point(data=samp_sci,aes(x=x,y=y),col="red",size=1)+
  # scale_color_gradientn(colours = pal) + 
  scale_color_distiller(palette = "RdBu",limits=c(0,NA))+
  ggtitle("Sample points",subtitle = "Red dots: scientific samples, scale color: commercial samples")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5,size=10))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")

# BrBG, PiYG, PRGn, PuOr, , RdGy, RdYlBu, RdYlGn, Spectral

## Simulation
simu_df <- data.frame(loc_x,S_x=S_x)
simu_plot <- ggplot(simu_df)+
  geom_point(aes(x=x,y=y,col=S_x),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
  ggtitle("Simulation")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")

## Scientific
sci_df <- data.frame(loc_x,S_x=res_sci$Report$S_x)
sci_plot <- ggplot(sci_df)+
  geom_point(aes(x=x,y=y,col=S_x),shape=15,size=1.5)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
  ggtitle("Scientific model")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")

## No reallocation - integrated
no.realloc_int_df <- data.frame(loc_x,S_x=res_no.realloc_int$Report$S_x)
no.realloc_int_plot <- ggplot(no.realloc_int_df)+
  geom_point(aes(x=x,y=y,col=S_x),shape=15,size=1.5)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
  ggtitle("Integrated model",
          subtitle = "No reallocation in estimation - Likelihood: Yi")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5,size=10))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")

## No reallocation - commercial
no.realloc_com_df <- data.frame(loc_x,S_x=res_no.realloc_com$Report$S_x)
no.realloc_com_plot <- ggplot(no.realloc_com_df)+
  geom_point(aes(x=x,y=y,col=S_x),shape=15,size=1.5)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
  ggtitle("Commercial model",
          subtitle = "No reallocation in estimation - Likelihood: Yi")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5,size=10))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")

## Reallocation - integrated
realloc_int_df <- data.frame(loc_x,S_x=res_realloc_int$Report$S_x)
realloc_int_plot <- ggplot(realloc_int_df)+
  geom_point(aes(x=x,y=y,col=S_x),shape=15,size=1.5)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
  ggtitle("Integrated model",
          subtitle = "Reallocation in estimation - Likelihood: Dj")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5,size=10))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")

## Reallocation - commercial
realloc_com_df <- data.frame(loc_x,S_x=res_realloc_com$Report$S_x)
realloc_com_plot <- ggplot(realloc_com_df)+
  geom_point(aes(x=x,y=y,col=S_x),shape=15,size=1.5)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
  ggtitle("Commercial model",
          subtitle = "Reallocation in estimation - Likelihood: Dj")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5,size=10))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")

full_plot <- plot_grid(sampling_plot,simu_plot,sci_plot,
                       NULL,no.realloc_int_plot,no.realloc_com_plot,
                       NULL,realloc_int_plot,realloc_com_plot,
                       ncol = 3,nrow = 3,align="hv")

ggsave(file="images/real_simu/simu_real_maps.png",width = 15, height = 15)
