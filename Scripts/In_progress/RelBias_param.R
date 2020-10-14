############################
## Load full simulation data
############################
#' @title RelBias_param_Results_loop : produce boxplots of relative bias of models parameters (or just bias when true value == 0 because in this case it is not possible to divide by the true value)
#' 
#' @param file : file of Results_loop data
#' 
#' file = "E:/outputs/Simulation/com_x_sci_data_10-2020-04-28_14_51_11_Alpha2_gamma-Results_loop.RData"
RelBias_param <- function(Results_loop){
  
  # Load results and separate simulation according to model that was used
  counter_table <- Results %>%
    filter(Convergence == 0) %>%
    dplyr::select(counter,Data_source,b_true)
  counter <- 0
  # parameter extraction loop each model
  for(lev in levels(factor(c("scientific_commercial","scientific_only","commercial_only")))){ # counter_table$Data_source
    count_plot <- 0
    x11()
    print(lev)
    par(mfrow = c(4,4) , mar = c(5.1, 4.1, 4.1, 2.1))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, lev, 
         cex = 1.1, col = "black")
    
    counter_table_f <- counter_table %>%
      filter(Data_source == lev)
    
    # if(lev == "scientific_only") counter_table_f <- counter_table_f[sample(x = c(1:nrow(counter_table_f)),size = 100,replace = F),] # sample size may have an effect on dispersion so I don't take all simu but only 100
    
    # dataframes
    RelBias_latent_linpred <- c()
    RelBias_fb_linpred <- c()
    RelBias_obs_param <- c()
    RelBias_spatial <- c()
    for(i in counter_table_f$counter){
      if(! is.null(List_param[[i]])){
        print(i)
        
        # load estimated parameter
        estim <- List_param[[i]]$SD$par.fixed # - 1 : because in TMB the first replicate is ine position before the first TRUE (which corresponds to the second replicate)
        if(! is.null(estim)){
          while(TRUE %in% duplicated(names(estim))){
            min_duplicat <- min(which(duplicated(names(estim)) == T)) - 1
            names(estim)[which(names(estim) == names(estim[min_duplicat]))] <- paste0(names(estim[min_duplicat]),which(names(estim) == names(estim[min_duplicat])) - min_duplicat)
          }
          
          # linear predictor of latent field
          latent_param_est <- estim[which(str_detect(names(estim),"beta_j") == T)]
          latent_param_true <- c(beta0 = List_param[[i]]$data.info$beta0,beta = List_param[[i]]$data.info$beta)
          names(latent_param_true) <- str_replace(names(latent_param_true),"beta","beta_j")
          RelBias_latent_linpred <- as.data.frame(rbind(RelBias_latent_linpred,c(counter = i,if_else(latent_param_true > 0,(latent_param_est - latent_param_true)/latent_param_true,(latent_param_est - latent_param_true)))))
          names(RelBias_latent_linpred) <- c("counter",names(latent_param_true))
          
          
          # linear predictor of fishing behavior
          if(TRUE %in% str_detect(names(estim),"beta_k")){
            fb_param_est <- estim[which(str_detect(names(estim),"beta_k") == T | names(estim) == "b" | names(estim) == "par_b")]
            # fb_param_est[which(names(fb_param_est) == "b" | names(fb_param_est) == "par_b")] <- exp(fb_param_est[which(names(fb_param_est) == "b" | names(fb_param_est) == "par_b")])
            if("b" %in% names(estim) | "par_b" %in% names(estim)){
              fb_param_true <- c(List_param[[i]]$data.info$beta0_fb,List_param[[i]]$data.info$beta_fb,b = counter_table_f$b_true[which(i == counter_table_f$counter)])
            } else{
              fb_param_true <- c(List_param[[i]]$data.info$beta0_fb,List_param[[i]]$data.info$beta_fb)
            }
            names(fb_param_true) <- names(fb_param_est)
            RelBias_fb_linpred <- as.data.frame(rbind(RelBias_fb_linpred,c(counter = i,if_else(fb_param_true != 0, (fb_param_est - fb_param_true)/fb_param_true, fb_param_est-fb_param_true))))
            names(RelBias_fb_linpred) <- c("counter",names(fb_param_true))
          }
          
          # observation equation parameter
          obs_param_est <- estim[which(str_detect(names(estim),"q1_") | str_detect(names(estim),"q2_") | str_detect(names(estim),"logSigma_") == T)]
          obs_param_true <- c(q1_com = List_param[[i]]$data.info$q1_com,q1_sci = List_param[[i]]$data.info$q1_sci,logSigma_com = List_param[[i]]$data.info$logSigma_com,logSigma_sci = List_param[[i]]$data.info$logSigma_sci)
          obs_param_true <- obs_param_true[which(names(obs_param_true) %in% names(obs_param_est) == T)]
          obs_param_est <- obs_param_est[which(names(obs_param_est) %in% names(obs_param_true))]
          obs_param_est <- obs_param_est[order(names(obs_param_est))] # order parameters so that obs_param_estand obs_param_true are in the same order 
          obs_param_true <- obs_param_true[order(names(obs_param_true))]
          RelBias_obs_param <- as.data.frame(rbind(RelBias_obs_param,c(counter = i,if_else(obs_param_true != 0, (obs_param_est - obs_param_true)/obs_param_true,obs_param_est - obs_param_true))))
          names(RelBias_obs_param) <- c("counter",names(obs_param_true))
        }
        
        # Spatial process parameter
        Spatial_param_est <- c(List_param[[i]]$Report$Range_S,List_param[[i]]$Report$MargSD_S)
        Spatial_param_est <- unique(Spatial_param_est)
        names(Spatial_param_est) <- c(rep("range",length(unique(List_param[[i]]$Report$Range_S))),
                                      rep("SD",length(List_param[[i]]$Report$MargSD_S)))
        Spatial_param_true <- c(range = List_param[[i]]$data.info$range,SD_delta = List_param[[i]]$data.info$SD_delta,SD_eta = List_param[[i]]$data.info$SD_eta)
        Spatial_param_true <- c(Spatial_param_true[which(names(Spatial_param_true) == "range")],
                                Spatial_param_true[which(names(Spatial_param_true) == "SD_delta" | names(Spatial_param_true) == "SD_eta")])
        
        Spatial_param_true <- c(rep(Spatial_param_true[which(str_detect("range",names(Spatial_param_true)))],length(str_which("range",names(Spatial_param_est)))),
                                  Spatial_param_true[which(names(Spatial_param_true) == "SD_delta")],
                                  Spatial_param_true[which(names(Spatial_param_true) == "SD_eta")])
        # 
        # if(lev == "scientific_only"){
        #   # Spatial_param_true <- c(Spatial_param_true[which(names(Spatial_param_true) == "range")],  # order parameters so that Spatial_param_est and Spatial_param_true are in the same order (in model "commercial_scientific" or "commercial only" the first SD is delta and the second is eta) 
        #   #                         Spatial_param_true[which(names(Spatial_param_true) == "SD_delta")])
        #   Spatial_param_est <- Spatial_param_est[1:2]
        # }
        
        RelBias_spatial <- as.data.frame(rbind(RelBias_spatial,c(counter = i,if_else(Spatial_param_true != 0,(Spatial_param_est - Spatial_param_true)/Spatial_param_true,(Spatial_param_est - Spatial_param_true)))))
        names(RelBias_spatial) <- c("counter",names(Spatial_param_true))
        
      }
    }
    
    # RelBias_latent_linpred <- RelBias_latent_linpred %>% dplyr::select(counter, beta_j0, beta_j1)
    RelBias_spatial <- RelBias_spatial %>% dplyr::select(-SD_eta)
    
    ## plots
    RelBias_latent_linpred <- inner_join(RelBias_latent_linpred,counter_table_f,by = "counter")
    for(name_param in names(RelBias_latent_linpred)[which(str_detect(names(RelBias_latent_linpred),"beta_j"))]){
      if(count_plot > 16) x11()
      count_plot <- count_plot+1
      boxplot(RelBias_latent_linpred[,name_param] ~ if_else(is.na(RelBias_latent_linpred$b_true) == T,0,RelBias_latent_linpred$b_true),
              main = name_param,
              xlab ="",ylab="Relative Bias")
      abline(h=0, lty = 2)
      
    }
    
    if(! is.null(RelBias_fb_linpred)){
      RelBias_fb_linpred  <- inner_join(RelBias_fb_linpred,counter_table_f,by = "counter")
      for(name_param in names(RelBias_fb_linpred)[which(str_detect(names(RelBias_fb_linpred),"beta_k") | names(RelBias_fb_linpred) == "b")]){
        if(count_plot > 16) x11()
        count_plot <- count_plot+1
        boxplot(RelBias_fb_linpred[,name_param] ~ if_else(is.na(RelBias_latent_linpred$b_true) == T,0,RelBias_latent_linpred$b_true),
                main = name_param,
                xlab ="",ylab="Relative Bias")
        abline(h=0, lty = 2)

      }
    }

    RelBias_obs_param  <- inner_join(RelBias_obs_param,counter_table_f,by = "counter")
    for(name_param in names(RelBias_obs_param)[which(str_detect(names(RelBias_obs_param),"q1_") == T | str_detect(names(RelBias_obs_param),"logSigma_") == T)]){
      if(count_plot > 16) x11()
      count_plot <- count_plot+1
      boxplot(RelBias_obs_param[,name_param] ~ if_else(is.na(RelBias_latent_linpred$b_true) == T,0,RelBias_latent_linpred$b_true),
              main = name_param,
              xlab ="",ylab="Relative Bias")
      abline(h=0, lty = 2)
    }
    
    # if(lev == "scientific_only") RelBias_spatial <- dplyr::select(RelBias_spatial,counter,range,SD_delta)
    
    if(length(which(names(RelBias_spatial) == "range")) == 2) names(RelBias_spatial)[which(names(RelBias_spatial) == "range")] <- paste0("range",c("_delta","_eta"))
    if(length(which(names(RelBias_spatial) == "range")) > 2) names(RelBias_spatial)[which(names(RelBias_spatial) == "range")] <- paste0("range.",seq(1:length(which(names(RelBias_spatial) == "range"))))
    
    RelBias_spatial  <- inner_join(RelBias_spatial,counter_table_f,by = "counter")
    for(name_param in names(RelBias_spatial)[which(str_detect(names(RelBias_spatial),"range") | str_detect(names(RelBias_spatial),"SD_") | str_detect(names(RelBias_spatial),"scale"))]){
      if(count_plot > 16) x11()
      count_plot <- count_plot+1
      boxplot(RelBias_spatial[,name_param] ~ if_else(is.na(RelBias_latent_linpred$b_true) == T,0,RelBias_latent_linpred$b_true),
      main = name_param,
      xlab ="",ylab="Relative Bias")
      abline(h=0, lty = 2)
    }
  }

}
