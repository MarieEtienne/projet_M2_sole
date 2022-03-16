################################################
## Plot outputs from mutliple square simulations
################################################
## B. Alglave


#----------------------------------------
## Plot latent field and covariate effect
#----------------------------------------
if(plot_outputs){
  
  if(sampling == "simulate"){
    
    ## Plot simulations outputs
    plot(x=log(S_x),
         y=log(Report$S_x),
         main="Estimated S_x vs. simulated S_x",
         xlab="Simulated",
         ylab="Estimated")
    mtext(text = paste0("Est_mod: ",Data_source[Estimation_model_i],
                        " Realloc: ",aggreg_obs,
                        " Converge: ",ifelse(opt$convergence==0,T,F)),
          cex=0.75)
    
    for(i in 1:2){
      if(i==1){
        S_plot <- S_x
        main_title <- "Simulated S_x"
        breaks_S <- quantile(S_plot,probs = seq(0,1,length.out=101))
        breaks_S <- unique(breaks_S)
        S_plot.ref <- S_x
        S_interval <- cut(S_plot,breaks = 100,include.lowest = T)
        pal_col <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(100))
        # pal_col <- rev(heat.colors(100))
      }else if(i==2){
        S_plot <- Report$S_x
        main_title <- "Estimated S_x"
        breaks_S <- quantile(S_plot,probs = seq(0,1,length.out=101))
        breaks_S <- unique(breaks_S)
        S_plot.ref <- S_x
        S_interval <- cut(S_plot,breaks = 100,include.lowest = T)
        pal_col <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(100))
      }
      S_interval <- cut(S_plot,breaks = breaks_S, include.lowest = T)
      S_xint <- data.frame(S_interval = S_interval, cell=1:n_cells)
      col_S_xint <- data.frame(S_interval=levels(S_interval),col=pal_col)
      color_df <- merge(x = S_xint, y = col_S_xint, by = "S_interval")
      color_df <- color_df[order(color_df$cell),]
      plot(x=loc_x$x,y=loc_x$y,
           col=as.character(color_df$col),
           pch = 15,cex=0.5,
           asp = 1,
           main = main_title,xlab = "x",ylab="y")
      pal_col2 <- pal_col[c(1,seq(10,100,by=10))]
      lgd_ = rep(NA, length(pal_col2))
      lgd_[c(1,5,10)] = (c(0,round(median(S_plot.ref),digits = 3),round(max(S_plot.ref),digits = 3)))
      legend(x = 17.5, y = 25,
             legend = lgd_,
             fill = pal_col2,
             border = NA,
             y.intersp = 0.5,
             cex = 1, text.font = 2,
             bg="white")
    }
    
  }else if(sampling == "from_tacsateflalo" & fit_IM_res$Converge == 0){
    
    S_plot <- Report$S_x
    breaks_S <- quantile(S_plot,probs = seq(0,1,length.out=101))
    breaks_S <- unique(breaks_S)
    S_interval <- cut(S_plot,breaks = 100,include.lowest = T)
    pal_col <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(breaks_S)-1))
    main_title <- paste0("Est model: ",Estimation_model_i)
    
    S_interval <- cut(S_plot,breaks = breaks_S, include.lowest = T)
    S_xint <- data.frame(S_interval = S_interval, cell=1:n_cells)
    col_S_xint <- data.frame(S_interval=levels(S_interval),col=pal_col)
    color_df <- merge(x = S_xint, y = col_S_xint, by = "S_interval")
    color_df <- color_df[order(color_df$cell),]
    plot(x=loc_x$x,y=loc_x$y,
         col=as.character(color_df$col),
         pch = 15,cex=0.5,
         asp = 1,
         main = main_title,xlab = "x",ylab="y")
    pal_col2 <- pal_col[c(1,seq(10,100,by=10))]
    lgd_ = rep(NA, length(pal_col2))
    lgd_[c(1,5,10)] = (c(0,round(median(S_plot),digits = 3),round(max(S_plot),digits = 3)))
    legend(x = 17.5, y = 25,
           legend = lgd_,
           fill = pal_col2,
           border = NA,
           y.intersp = 0.5,
           cex = 1, text.font = 2,
           bg="white")
  }
  
  
  if(fit_IM_res$Converge == 0){
    par_est <- fit_IM_res$SD$value[which(str_detect(names(fit_IM_res$SD$value),"beta_j|MargSD|Range|q1_com|Sigma_com|q1_sci|Sigma_sci|k_com"))]
    sd_est <- fit_IM_res$SD$sd[which(str_detect(names(fit_IM_res$SD$value),"beta_j|MargSD|Range|q1_com|Sigma_com|q1_sci|Sigma_sci|k_com"))]
    x_axis <- 1:length(par_est)
  }else{
    par_est <- rep(0,length(par_est))
    sd_est <- rep(0,length(par_est))
  }
  
  plot(x=x_axis,
       y=par_est,
       ylim=range(c(par_est - 2, par_est + 2)),
       xaxt = "n",
       pch=19, xlab="", ylab="Parameters values")
  
  # hack: we draw arrows but with very special "arrowheads"
  arrows(x_axis, par_est - 1.96 * sd_est,
         x_axis, par_est + 1.96 * sd_est,
         length=0.05, angle=90, code=3)
  
  names(par_est)[which(str_detect(names(par_est),"beta_j"))[1]] <- "intercept"
  
  axis(1,
       at=1:(length(par_est)),
       labels=names(par_est),
       tick=T,las=2,cex.axis=0.75)
  
  abline(h=0, col="skyblue", lwd=2, lty=2)
  
  if(sampling == "simulate" & fit_IM_res$Converge == 0){
    
    ref_capt <- "k_com"
    
    if(Estimation_model_i == 1){
      par_name <- c("intercept","beta_j","MargSD","Range","q1_com","Sigma_com","q1_sci","Sigma_sci",ref_capt)
      par_value <- c(intercept,beta1,SD_delta,range_delta,q1,SD_obs,q1_sci,exp(logSigma_sci),1)
    }
    
    if(Estimation_model_i == 2){
      par_name <- c("intercept","beta_j","MargSD","Range","q1_sci","Sigma_sci")
      par_value <- c(intercept,beta1,SD_delta,range_delta,q1_sci,exp(logSigma_sci))
    }
    
    
    if(Estimation_model_i == 3){
      par_name <- c("intercept","beta_j","MargSD","Range","q1_com","Sigma_com")
      par_value <- c(intercept,beta1,SD_delta,range_delta,q1,SD_obs)
    }
    
    names(par_value) <- par_name
    par_value <- par_value[names(par_est)]    
    
    points(x=x_axis,y=par_value,col="red",pch=19)
    
  }
  
}

#-------------------
## Consistency check
#-------------------
if(Estimation_model_i == 1){
  
  obj_int <- Obj
  SD_int <- SD
  obj_int$fn()
  pl_int <- as.list(SD_int,"Est")
  
}

if(Estimation_model_i == 2){
  
  obj_sci <- Obj
  SD_sci <- SD
  obj_sci$fn()
  pl_sci <- as.list(SD_sci,"Est")
  par_sci <- SD_sci$par.fixed
  f_sci <- as.numeric(obj_sci$fn(par_sci))
  
  ## Fixed
  par_int <- SD_int$par.fixed # extractBaseParameters(obj_sci, pl_sci, pl_int)
  par_int <- par_int[which(names(par_int) %in% names(par_sci))]
  f_int <- as.numeric(obj_sci$fn(par_int))
  def_free <- sum(par_sci != 0)
  fixed <- 1 - pchisq( 2 * (f_int-f_sci), df=def_free )
  
  ## Random
  # par_int.all <- extractBaseParameters(obj_sci, pl_sci, pl_int, all=TRUE)
  par_sci.all <- obj_sci$env$last.par.best
  par_int.all <- obj_int$env$last.par.best
  par_int.all <- par_int.all[which(names(par_int.all) %in% names(par_sci.all))]
  f_sci.all <- obj_sci$env$f(par_sci.all) #Best evaluated parameters
  f_int.all <- obj_sci$env$f(par_int.all)
  def_free.all <- def_free + length(obj_sci$env$random)
  random <- 1 - pchisq( 2 * (f_int.all-f_sci.all), df=def_free.all )
  
}

if(plot_outputs & consistency_check == T & Estimation_model_i == 3){
  
  x11(width = 20,height = 10)
  par(mfrow = c(1,2))
  plot(x=log(SD_sci$value[which(names(SD_sci$value)=="S_x")]),
       y=log(SD_int$value[which(names(SD_int$value)=="S_x")]),
       xlab="log(scientific predictions)",ylab="log(integrated predictions)",
       xlim = log(range(c(SD_sci$value[which(names(SD_sci$value)=="S_x")]))),
       ylim = log(range(c(SD_int$value[which(names(SD_int$value)=="S_x")]))))
  abline(a = 0,b = 1)
  
  mtext(text = paste0("Consistency check p-value \n",
                      "fixed = ",format(signif(fixed,digits = 3),scientific = T),
                      " | random = ",format(signif(random,digits = 3),scientific = T)),
        side = 3,cex = 0.75)
  
  plot(x=SD_sci$par.fixed,
       y=SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))],
       xlab="Scientific",ylab="Integrated",
       xlim = c(min(c(SD_sci$par.fixed,SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))]))-1,
                max(c(SD_sci$par.fixed,SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))])))+1,
       ylim = c(min(c(SD_sci$par.fixed,SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))]))-1,
                max(c(SD_sci$par.fixed,SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))]))+1)
  )
  
  abline(a = 0,b = 1)
  text(x=SD_sci$par.fixed,
       y=SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))],
       labels=names(SD_int$par.fixed)[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))],
       pos=1)
  
}


# test <- cbind(loc_x[,c("x","y")],
#               # S_est = fit_IM_res$Report$S_p,
#               S_sim = S_x)
#
# simu_plot <- ggplot(test)+
#   geom_point(aes(x=x,y=y,col=S_sim))+
#   scale_color_distiller(palette = "Spectral")+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Simulation")
#
# est_plot <- ggplot(test)+
#   geom_point(aes(x=x,y=y,col=S_est))+
#   scale_color_distiller(palette = "Spectral")+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Estimation")
#
# x11()
# plot_grid(effort_plot,seq_plot,simu_plot,est_plot)



# data.points_df <- cbind(loc_x[,c("x","y")],c_com_x)
# effort_plot <- ggplot()+
#   geom_point(data = data.points_df,aes(x=x,y=y,col=c_com_x))+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Sampling effort")+
#   scale_color_distiller(palette = "Spectral")
#
# sample_df <- data.frame(cell=index_i,
#                         seq=boats_i,
#                         catch_no.realloc=y_i,
#                         catch_realloc=y_i2) %>%
#   inner_join(loc_x[,c("cell","x","y")]) %>%
#   filter(seq < 10)
# 
# seq_plot <- ggplot()+
#   geom_point(data = sample_df,aes(x=x,y=y,col=factor(seq)))+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Fishing sequence")

# extractBaseParameters <- function(obj1, pl1, pl2, all=FALSE) {
#   npar1 <- sapply(pl1,length) #Number of survey parameters by component
#   plnew <- Map(head, pl2, npar1) #Extract survey parameters from second fit
#   applyMap <- function(parameters, map) {
#     param.map <- lapply(names(map),
#                         function(nam)
#                         {
#                           TMB:::updateMap(parameters[[nam]], map[[nam]])
#                         })
#     parameters[names(map)] <- param.map
#     parameters
#   }
#   par2 <- unlist(applyMap(plnew, obj1$env$map))
#   if (!all) par2 <- par2[-obj1$env$random]
#   par2
# }



# ## Plot fishing sequence catch values
# seq_plot <- ggplot()+
#   geom_point(data = sample_df,aes(x=x,y=y,col=factor(seq)))+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   xlab("")+ylab("")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust=0.5),
#         legend.position = "none",
#         panel.background = element_rect(fill="skyblue"))
# 
# catch.no.realloc_plot <- ggplot()+
#   geom_point(data = test,aes(x=x,y=y,col=S_sim),alpha=0.2,shape=15)+
#   # geom_point(data = sample_df,aes(x=x,y=y,fill=catch_no.realloc),col="black",shape=21)+
#   scale_color_distiller(palette="Spectral")+
#   scale_fill_distiller(palette="Spectral")+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Exact catches")+xlab("")+ylab("")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust=0.5),legend.position = "none")
# 
# catch.realloc_plot <- ggplot()+
#   geom_point(data = test,aes(x=x,y=y,col=S_sim),alpha=0.2,shape=15)+
#   geom_point(data = sample_df,aes(x=x,y=y,fill=catch_realloc),col="black",shape=21)+
#   scale_color_distiller(palette="Spectral")+
#   scale_fill_distiller(palette="Spectral")+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Reallocated catches")+xlab("")+ylab("")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust=0.5),legend.position = "none")
# 
# plot_grid(NULL,seq_plot,catch.no.realloc_plot,catch.realloc_plot)
