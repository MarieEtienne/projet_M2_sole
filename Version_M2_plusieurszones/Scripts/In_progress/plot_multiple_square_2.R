##############################################################
## Multiple square simulations - Plot single simulation ouputs
##############################################################

compute_profile <- F
save_data <- F

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

# ## To build list input
# list_input <- list(loc_x=loc_x,
#                    S_x=S_x,
#                    c_com_x=c_com_x,
#                    y_i=y_i,
#                    y_i2=y_i2,
#                    index_i=index_i,
#                    y_sci_i=y_sci_i,
#                    index_sci_i=index_sci_i,
#                    boats_i=boats_i,
#                    cov_x_2=cov_x_2,
#                    ref_level=ref_level,
#                    spde=spde,
#                    mesh=mesh,
#                    n_cells=n_cells,
#                    sampling=sampling)
# save(data=list_input,file="draft/paper_3/res/list_input.RData")


## Simulated latent field
simu_df <- data.frame(loc_x,S_x=S_x)

if(save_data == T) save(data=simu_df,file="draft/paper_3/res/simu_df.RData")

## Simulated data
sci.data_df <- data.frame(cell = index_sci_i,
                          y_sci = y_sci_i) %>%
  inner_join(loc_x[,c("x","y","cell","ICESNAME")])

com.data_df <- data.frame(cell = index_i,
                          y_com = y_i2) %>%
  inner_join(loc_x[,c("x","y","cell","ICESNAME")]) %>%
  group_by(x,y,cell,ICESNAME) %>%
  dplyr::summarise(n=n())

## Scientific data
fit_IM_res <- fit_IM(Estimation_model_i = 2,
                     Samp_process = 0,
                     EM = "fix_b",
                     TmbFile = "Scripts/",
                     ignore.uncertainty = F,
                     c_com_x = c_com_x,
                     y_com_i = y_i2,
                     index_com_i = index_i,
                     y_sci_i = y_sci_i,
                     index_sci_i = index_sci_i,
                     aggreg_obs=T,
                     boats_number = boats_i,
                     Cov_x = as.matrix(cov_x_2), # NULL, # 
                     ref_level = ref_level,
                     lf_param = "RE", # lf_param,
                     spde=spde,
                     mesh=mesh,
                     n_cells=n_cells,
                     cov_est = T,
                     Params_step.est=NULL,
                     Map_step.est=NULL,
                     ObsM=F,
                     y_ObsM_i=NULL,
                     index_ObsM_i=NULL,
                     sampling = sampling,
                     landings = T,
                     quadratic_cov = F)

if(save_data == T) save(data=fit_IM_res,file="draft/paper_3/res/fit_IM_res_sci.RData")

sci_df <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)

if(compute_profile){
  sci_prof <- TMB::tmbprofile(fit_IM_res$Obj,lincomb = c(0,1,0,0,0,0))
  sci_prof_df <- sci_prof %>%
    mutate(Estimation_model = 2) %>%
    mutate(value_2 = (value - mean(value)) / sd(value) ) %>%
    mutate(lkl = " ")
}

## Integrated model - Yi
fit_IM_res <- fit_IM(Estimation_model_i = 1,
                     Samp_process = 0,
                     EM = "fix_b",
                     TmbFile = "Scripts/",
                     ignore.uncertainty = F,
                     c_com_x = c_com_x,
                     y_com_i = y_i2,
                     index_com_i = index_i,
                     y_sci_i = y_sci_i,
                     index_sci_i = index_sci_i,
                     aggreg_obs=F,
                     boats_number = boats_i,
                     Cov_x = as.matrix(cov_x_2), # NULL, # 
                     ref_level = ref_level,
                     lf_param = "RE", # lf_param,
                     spde=spde,
                     mesh=mesh,
                     n_cells=n_cells,
                     cov_est = T,
                     Params_step.est=NULL,
                     Map_step.est=NULL,
                     ObsM=F,
                     y_ObsM_i=NULL,
                     index_ObsM_i=NULL,
                     sampling = sampling,
                     landings = T,
                     quadratic_cov = F)

if(save_data == T) save(data=fit_IM_res,file="draft/paper_3/res/fit_IM_res_int_Yi.RData")

int_Yi_df <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)

if(compute_profile){
  int_Yi_prof <- TMB::tmbprofile(fit_IM_res$Obj,lincomb = c(0,1,0,0,0,0,0,0,0))
  int_Yi_prof_df <- int_Yi_prof %>%
    mutate(Estimation_model = 1) %>%
    mutate(value_2 = (value - mean(value)) / sd(value) ) %>%
    mutate(lkl = "Yi")
}

## Commercial model - Yi
fit_IM_res <- fit_IM(Estimation_model_i = 3,
                     Samp_process = 0,
                     EM = "fix_b",
                     TmbFile = "Scripts/",
                     ignore.uncertainty = F,
                     c_com_x = c_com_x,
                     y_com_i = y_i2,
                     index_com_i = index_i,
                     y_sci_i = y_sci_i,
                     index_sci_i = index_sci_i,
                     aggreg_obs=F,
                     boats_number = boats_i,
                     Cov_x = as.matrix(cov_x_2), # NULL, # 
                     ref_level = ref_level,
                     lf_param = "RE", # lf_param,
                     spde=spde,
                     mesh=mesh,
                     n_cells=n_cells,
                     cov_est = T,
                     Params_step.est=NULL,
                     Map_step.est=NULL,
                     ObsM=F,
                     y_ObsM_i=NULL,
                     index_ObsM_i=NULL,
                     sampling = sampling,
                     landings = T,
                     quadratic_cov = F)

if(save_data == T) save(data=fit_IM_res,file="draft/paper_3/res/fit_IM_res_com_Yi.RData")

com_Yi_df <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)

if(compute_profile){
  com_Yi_prof <- TMB::tmbprofile(fit_IM_res$Obj,lincomb = c(0,1,0,0,0,0))
  com_Yi_prof_df <- com_Yi_prof %>%
    mutate(Estimation_model = 3) %>%
    mutate(value_2 = (value - mean(value)) / sd(value) ) %>%
    mutate(lkl = "Yi")
}

## Integrated model - Dj
fit_IM_res <- fit_IM(Estimation_model_i = 1,
                     Samp_process = 0,
                     EM = "fix_b",
                     TmbFile = "Scripts/",
                     ignore.uncertainty = F,
                     c_com_x = c_com_x,
                     y_com_i = y_i2,
                     index_com_i = index_i,
                     y_sci_i = y_sci_i,
                     index_sci_i = index_sci_i,
                     aggreg_obs=T,
                     boats_number = boats_i,
                     Cov_x = as.matrix(cov_x_2), # NULL, # 
                     ref_level = ref_level,
                     lf_param = "RE", # lf_param,
                     spde=spde,
                     mesh=mesh,
                     n_cells=n_cells,
                     cov_est = T,
                     Params_step.est=NULL,
                     Map_step.est=NULL,
                     ObsM=F,
                     y_ObsM_i=NULL,
                     index_ObsM_i=NULL,
                     sampling = sampling,
                     landings = T,
                     quadratic_cov = F)

if(save_data == T) save(data=fit_IM_res,file="draft/paper_3/res/fit_IM_res_int_Dj.RData")

int_Dj_df <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)

if(compute_profile){
  int_Dj_prof <- TMB::tmbprofile(fit_IM_res$Obj,lincomb = c(0,1,0,0,0,0,0,0,0))
  int_Dj_prof_df <- int_Dj_prof %>%
    mutate(Estimation_model = 1) %>%
    mutate(value_2 = (value - mean(value)) / sd(value) ) %>%
    mutate(lkl = "Dj")
}


## Commercial model - Dj
fit_IM_res <- fit_IM(Estimation_model_i = 3,
                     Samp_process = 0,
                     EM = "fix_b",
                     TmbFile = "Scripts/",
                     ignore.uncertainty = F,
                     c_com_x = c_com_x,
                     y_com_i = y_i2,
                     index_com_i = index_i,
                     y_sci_i = y_sci_i,
                     index_sci_i = index_sci_i,
                     aggreg_obs=T,
                     boats_number = boats_i,
                     Cov_x = as.matrix(cov_x_2), # NULL, # 
                     ref_level = ref_level,
                     lf_param = "RE", # lf_param,
                     spde=spde,
                     mesh=mesh,
                     n_cells=n_cells,
                     cov_est = T,
                     Params_step.est=NULL,
                     Map_step.est=NULL,
                     ObsM=F,
                     y_ObsM_i=NULL,
                     index_ObsM_i=NULL,
                     sampling = sampling,
                     landings = T,
                     quadratic_cov = F)

if(save_data == T) save(data=fit_IM_res,file="draft/paper_3/res/fit_IM_res_com_Dj.RData")

com_Dj_df <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)

if(compute_profile){
  com_Dj_prof <- TMB::tmbprofile(fit_IM_res$Obj,lincomb = c(0,1,0,0,0,0))
  com_Dj_prof_df <- com_Dj_prof %>%
    mutate(Estimation_model = 3) %>%
    mutate(value_2 = (value - mean(value)) / sd(value) ) %>%
    mutate(lkl = "Dj")
}



############
## Plot maps
############


max_val <- max(simu_df$S_x/sum(simu_df$S_x),
               sci_df$S_x/sum(sci_df$S_x),
               int_Yi_df$S_x/sum(int_Yi_df$S_x),
               com_Yi_df$S_x/sum(com_Yi_df$S_x),
               int_Dj_df$S_x/sum(int_Dj_df$S_x),
               com_Dj_df$S_x/sum(com_Dj_df$S_x))

simu_df_2 <- simu_df %>%
  mutate(Ab = sum(S_x))

data_plot <- ggplot()+
  geom_point(data=simu_df_2,aes(x=x,y=y,col=S_x/Ab),shape=15,size=2,alpha=0.125)+
  geom_point(data=simu_df_2[which(simu_df_2$ICESNAME %in% com.data_df$ICESNAME),],
             aes(x=x,y=y,col=S_x/Ab),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val))+
  # geom_point(data=com.data_df,aes(x=x,y=y),col="grey",shape=15,size=2,alpha=0.45)+
  geom_point(data=sci.data_df,aes(x=x,y=y),col="red",size=1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "none"
  )+
  geom_sf(data=ICES_rect_2,alpha=0)+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

simu_plot <- ggplot(simu_df)+
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Simulated latent field")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

sci_plot <- ggplot(sci_df)+
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Scientific model")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

int_Yi_plot <- ggplot(int_Yi_df)+
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Integrated model - Yr")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

com_Yi_plot <- ggplot(com_Yi_df)+
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Commercial model - Yr")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

int_Dj_plot <- ggplot(int_Dj_df)+
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Integrated model - Dk")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

com_Dj_plot <- ggplot(com_Dj_df)+
  geom_point(aes(x=x,y=y,col=S_x/sum(S_x)),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,max_val+0.01*max_val))+
  ggtitle("Commercial model - Dk")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")

legend <- as_ggplot(cowplot::get_legend(
  com_Dj_plot+
    theme(legend.position="right",
          legend.title = element_blank()
    )))

map_plot <- plot_grid(
  data_plot,simu_plot,sci_plot,
  legend,int_Yi_plot,com_Yi_plot,
  NULL,int_Dj_plot,com_Dj_plot,
  nrow = 3,
  ncol = 3,
  align="hv",
  labels = c("A","B","C","","D","E","","F","G")
)

ggsave("draft/paper_3/images/Map_multi_square.png",width = 12, height = 12)



## Likelihood profile
sci_prof_df_1 <- sci_prof_df
sci_prof_df_2 <- sci_prof_df
sci_prof_df_1$lkl <- "Dj"
sci_prof_df_2$lkl <- "Yi"

sci_prof_df_1$value <- sci_prof_df_1$value_2 - min(sci_prof_df_1$value_2)
sci_prof_df_2$value <- sci_prof_df_2$value_2 - min(sci_prof_df_2$value_2)
int_Yi_prof_df$value <- int_Yi_prof_df$value_2 - min(int_Yi_prof_df$value_2)
com_Yi_prof_df$value <- com_Yi_prof_df$value_2 - min(com_Yi_prof_df$value_2)
int_Dj_prof_df$value <- int_Dj_prof_df$value_2 - min(int_Dj_prof_df$value_2)
com_Dj_prof_df$value <- com_Dj_prof_df$value_2 - min(com_Dj_prof_df$value_2)

lkl_prof_df <- rbind(sci_prof_df_1,
                     sci_prof_df_2,
                     int_Yi_prof_df,
                     com_Yi_prof_df,
                     int_Dj_prof_df,
                     com_Dj_prof_df)

lkl_prof_df$Model <- NA
lkl_prof_df$Model[which(lkl_prof_df$Estimation_model == 1)] <- "Integrated"
lkl_prof_df$Model[which(lkl_prof_df$Estimation_model == 2)] <- "Scientific"
lkl_prof_df$Model[which(lkl_prof_df$Estimation_model == 3)] <- "Commercial"

lkl_prof_df$lkl <- factor(lkl_prof_df$lkl,levels = c("Yi","Dj"))

if(save_data == T) save(data=fit_IM_res,file="draft/paper_3/res/lkl_prof_beta.RData")

ggplot(lkl_prof_df)+
  geom_line(aes(x=parameter,y=-value,col=Model),size=1)+
  theme()+
  facet_wrap(.~lkl)+
  theme_bw()+
  ylab("Rescaled log-likelihood values")+
  xlab("Parameter value (beta)")+
  geom_vline(xintercept = 2,col="darkgrey",size=1,linetype="dashed")



