###################
## Simulation plots
###################

## Simulation boxplots
#---------------------

load("draft/paper_3/res/Results_full_multi_square.RData")

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
