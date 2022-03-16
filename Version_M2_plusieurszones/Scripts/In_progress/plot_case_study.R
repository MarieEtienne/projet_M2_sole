
#######
## Maps
#######
## Fitted data
# ggplot()+
#   geom_sf(data=points_sf_2,alpha=0.5,size=0.5)+
#   geom_point(data=HH_obsmer_2,aes(x=x,y=y),col="orange")+
#   geom_point(data=ptsurvey_wgs84_2,aes(x=x,y=y),col="red")+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=st_geometry(st_as_sf(survey_layer_3)),
#           size=1.5,col="darkgrey",
#           alpha=0)+
#   geom_sf(data = mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48),
#            expand = FALSE)+
#   scale_color_distiller(palette="Spectral")+
#   theme(panel.background = element_rect(fill = "skyblue"))+
#   xlab("")+ylab("")

load("results/case_study/no.realloc_int_df.RData")
# save(data=fit_IM_res,file="results/case_study/no.realloc_int_df.RData")

no.realloc_df <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)
no.realloc_plot <- ggplot(no.realloc_df)+
  geom_point(aes(x=x,y=y,col=S_x),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
  ggtitle("Integrated model",subtitle = "Commercial likelihood on Yi")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")

# load("results/case_study/realloc_int_df.RData")
# 
# realloc_df <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)
# realloc_plot <- ggplot(realloc_df)+
#   geom_point(aes(x=x,y=y,col=S_x),shape=15,size=2)+
#   scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
#   ggtitle("Integrated model",subtitle = "Commercial likelihood on Dj")+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5))+
#   geom_sf(data = mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
#   xlab("")+ylab("")

load("results/case_study/realloc_int_df_rest.RData")
# save(data=fit_IM_res,file="results/case_study/realloc_int_df_rest.RData")

realloc_rest_df <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)
realloc_rest_plot <- ggplot(realloc_rest_df)+
  geom_point(aes(x=x,y=y,col=S_x),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
  ggtitle("Integrated model",subtitle = "Commercial likelihood on Dj - step estimation")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")


load("results/case_study/sci_df.RData")
# save(data=fit_IM_res,file="results/case_study/sci_df.RData")

sci_df <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)
sci_plot <- ggplot(sci_df)+
  geom_point(aes(x=x,y=y,col=S_x),shape=15,size=2)+
  scale_color_distiller(palette = "Spectral",limits=c(0,NA))+
  ggtitle("Scientific model")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_sf(data = mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)+
  xlab("")+ylab("")

case_study_plot <- plot_grid(sci_plot,no.realloc_plot,realloc_rest_plot,ncol=3,align="hv")

ggsave(file="images/case_study/case_study_maps.png",width = 9, height = 3)


############
## Estimates
############
file_vec <- c("no.realloc_int_df","sci_df","realloc_int_df_rest") # ,"realloc_int_df"
for(file_i in file_vec){
  
  print(file_i)
  load(paste0("results/case_study/",file_i,".RData"))
  
  if(fit_IM_res$Converge == 0){
    par_est <- fit_IM_res$SD$value[which(str_detect(names(fit_IM_res$SD$value),"beta_j|MargSD|Range|q1_com|Sigma_com|q1_sci|Sigma_sci|k_com"))]
    sd_est <- fit_IM_res$SD$sd[which(str_detect(names(fit_IM_res$SD$value),"beta_j|MargSD|Range|q1_com|Sigma_com|q1_sci|Sigma_sci|k_com"))]
    x_axis <- 1:length(par_est)
  }else{
    par_est <- rep(0,length(par_est))
    sd_est <- rep(0,length(par_est))
  }
  
  if(str_detect(file_i,"int")) Model <- "Integrated"
  if(str_detect(file_i,"com")) Model <- "Commercial"
  if(str_detect(file_i,"sci")) Model <- "Scientific"
  
  if(str_detect(file_i,"no.realloc")) lkl <- "Yi"
  if(!str_detect(file_i,"no.realloc")) lkl <- "Dj"
  if(str_detect(file_i,"_rest")) lkl <- "Dj - r.est"
  if(str_detect(file_i,"sci")) lkl <- "Scientific"
  
  
  names(par_est)[which(str_detect(names(par_est),"beta_j"))[1]] <- "intercept"
  names(par_est)[which(str_detect(names(par_est),"beta_j"))] <- colnames(fit_IM_res$Data$Cov_xj)[-1]
  
  est_par_df <- data.frame(par_names = names(par_est),
                           par_val = par_est,
                           sd = sd_est,
                           Model = Model,
                           lkl = lkl)
  
  if(file_i == file_vec[1]){
    est_par_df_full <- est_par_df
  }else{
    est_par_df_full <- rbind(est_par_df_full,est_par_df)
  }
  
}

est_par_df_full_2 <- est_par_df_full %>%
  mutate(CI.inf = par_val - 1.96 * sd,
         CI.sup = par_val + 1.96 * sd) %>%
  mutate(CI.inf = ifelse(CI.inf < -2.5,-Inf,CI.inf),
         CI.sup = ifelse(CI.sup > 5,Inf,CI.sup))

ref_capt <- "k_com"
est_par_df_full_2$Model <- factor(est_par_df_full_2$Model,levels = c("Scientific","Integrated","Integrated_rest"))
est_par_df_full_2$par_names <- factor(est_par_df_full_2$par_names,levels = rev(c("intercept",colnames(fit_IM_res$Data$Cov_xj),"MargSD","Range","q1_sci","Sigma_sci","q1_com","Sigma_com",ref_capt)))

est_par_df_full_2$lkl[which(est_par_df_full_2$lkl == "Yi")] <- "Yi"
est_par_df_full_2$lkl[which(est_par_df_full_2$lkl == "Dj - r.est")] <- "Dj"
est_par_df_full_2$lkl[which(est_par_df_full_2$lkl == "Scientific")] <- "Scientific model"

est_par_df_full_3 <- est_par_df_full_2 # %>%
  # filter(par_names %in% c("intercept",
  #                         "bathy2",
  #                         "substr_Sand_Coarse_substrate",
  #                         "substr_Rock",
  #                         "substr_Mud_sediment",
  #                         "Sigma_com",
  #                         "q1_com"))


est_par_df_full_3$lkl <- factor(est_par_df_full_3$lkl,levels = c("Dj","Yi","Scientific model"))
par_plot <- ggplot(est_par_df_full_3, aes(y=par_val, x=par_names))+
  geom_point(
    aes(color = lkl),
    position = position_dodge(0.5),
    size=2
  )+
  geom_errorbar(aes(ymin = CI.inf, ymax = CI.sup, color = lkl),
                position = position_dodge(0.5),width=0.4
  )+
  scale_color_manual(breaks = c("Scientific model",
                                "Yi",
                                "Dj"),
                     values=c("black","chartreuse3","red"))+
  geom_hline(yintercept=0, linetype="dashed", color = "red",alpha=0.4)+
  xlab("")+ylab("")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "right",
        aspect.ratio = 1)+
  coord_flip()# +
  # ylim(-2.5,5)+
  # facet_wrap(.~Model,ncol = 1)

ggsave(file="images/case_study/par_plot.png",width = 10, height = 5)

