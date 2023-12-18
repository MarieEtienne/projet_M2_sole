
library(rnaturalearth)
library(sf)
library(mapdata)
mapBase <- map("worldHires", fill = T, plot = F)
mapBase <- st_as_sf(mapBase) %>% filter(ID %in% c("France","Spain"))
ref_capt <- "k_com"

load(file="results/real_simu/simu.data.RData")
S_x <- simu.data$S_x
index_i <- simu.data$index_i
y_sci_i <- simu.data$y_sci_i
y_i <- simu.data$y_i
boats_i <- simu.data$boats_i
cov_x <- simu.data$cov_x
index_sci_i <- simu.data$index_sci_i

load(file="results/real_simu/no.realloc_int_df.RData")
res_no.realloc_int <- fit_IM_res
load(file="results/real_simu/sci_df.RData")
res_sci <- fit_IM_res
load(file="results/real_simu/no.realloc_com_df.RData")
res_no.realloc_com <- fit_IM_res
load(file="results/real_simu/realloc_int_df.RData")
res_realloc_int <- fit_IM_res
load(file="results/real_simu/realloc_com_df.RData")
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

############
## Estimates
############
file_vec <- c("no.realloc_int_df","sci_df","no.realloc_com_df","realloc_int_df","realloc_com_df")
for(file_i in file_vec){
        
        print(file_i)
        load(paste0("results/real_simu/",file_i,".RData"))
        
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
        if(str_detect(file_i,"sci")) lkl <- "Scientific"
        
        
        names(par_est)[which(str_detect(names(par_est),"beta_j"))[1]] <- "intercept"
        
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

est_par_df_full <- est_par_df_full %>%
        mutate(CI.inf = par_val - 1.96 * sd,
               CI.sup = par_val + 1.96 * sd) %>%
        mutate(CI.inf = ifelse(CI.inf < -2.5,-Inf,CI.inf),
               CI.sup = ifelse(CI.sup > 5,Inf,CI.sup))

true_par_df <- data.frame(par_names = c("intercept","beta_j","MargSD","Range","q1_sci","Sigma_sci","q1_com","Sigma_com",ref_capt),
                          par_val_true = c(intercept,beta1,SD_delta,range_delta,q1_sci,exp(logSigma_sci),q1,SD_obs,1))

est_par_df_full_2 <- est_par_df_full %>%
        filter(Model != "Commercial") %>%
        inner_join(true_par_df) %>%
        mutate(CI.inf.rel = CI.inf/par_val_true,
               CI.sup.rel = CI.sup/par_val_true,
               par_val.rel = par_val/par_val_true,
               par_val_true.rel = 1)

est_par_df_full_2$Model <- factor(est_par_df_full_2$Model,levels = c("Scientific","Integrated"))
est_par_df_full_2$par_names <- factor(est_par_df_full_2$par_names,levels = rev(c("intercept","beta_j","MargSD","Range","q1_sci","Sigma_sci","q1_com","Sigma_com",ref_capt)))

est_par_df_full_2$lkl[which(est_par_df_full_2$lkl == "Yi")] <- "Integrated model - Commercial likelihood on Yi"
est_par_df_full_2$lkl[which(est_par_df_full_2$lkl == "Dj")] <- "Integrated model - Commercial likelihood on Dj"
est_par_df_full_2$lkl[which(est_par_df_full_2$lkl == "Scientific")] <- "Scientific model"

par_plot <- ggplot(est_par_df_full_2, aes(y=par_val, x=par_names))+
        geom_point(
                aes(color = lkl),
                position = position_dodge(0.5),
                size=2
        )+
        geom_errorbar(aes(ymin = CI.inf, ymax = CI.sup, color = lkl),
                      position = position_dodge(0.5),width=0.4
                      )+
        scale_color_manual(breaks = c("Scientific model","Integrated model - Commercial likelihood on Yi", "Integrated model - Commercial likelihood on Dj"),
                           values=c("black","chartreuse3","cadetblue3"))+
        geom_hline(yintercept=0, linetype="dashed", color = "red",alpha=0.4)+
        xlab("")+ylab("")+
        theme_bw()+
        theme(legend.title = element_blank(),
              legend.position = "right",
              aspect.ratio = 1)+
        coord_flip()+
        # ylim(-2.5,5)+
        # facet_wrap(.~Model,ncol = 1)+
        geom_point(aes(y=par_val_true, x=par_names),col="red",size=2)

ggsave(file="images/real_simu/par_plot.png",width = 10, height = 5)

