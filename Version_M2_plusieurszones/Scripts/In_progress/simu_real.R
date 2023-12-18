## specific loops
counter <- 1201
simu_file <- paste0("/media/balglave/Elements/results/q1/")
load(paste0(simu_file,"/Results.RData"))

i0 = 70

for(q1 in c(-1,-3,-5)){
# for(n_zone in c(1,3,5)){
  
  for(i in i0:100){
    
    set.seed( i ) # for figures : i = 2
    
    cov_x <- sim_GF_Matern(as.matrix(loc_x[,c("x","y")]), nu, range_cov, SD_cov^2)[[1]]$y
    delta_x <- sim_GF_Matern(as.matrix(loc_x[,c("x","y")]), nu, range_delta, SD_delta^2)[[1]]$y
    
    S_x <- exp(intercept + cov_x %*% beta1 + delta_x)
    vec <- cov_x %*% beta1
    loc_x_2 <- cbind(loc_x,cov=cov_x,delta=delta_x,S_x=S_x)
    
    # cov_plot <- ggplot(loc_x_2)+
    #   geom_point(aes(x=x,y=y,col=cov_x))+
    #   scale_color_distiller(palette = "Spectral")+
    #   ggtitle("Covariate")
    # 
    # delta_plot <- ggplot(loc_x_2)+
    #   geom_point(aes(x=x,y=y,col=delta_x))+
    #   scale_color_distiller(palette = "Spectral")+
    #   ggtitle("Random effect")
    # 
    # lf_plot <- ggplot(loc_x_2)+
    #   geom_point(aes(x=x,y=y,col=S_x))+
    #   scale_color_distiller(palette = "Spectral")+
    #   ggtitle("Latent field")
    # 
    # plot_grid(cov_plot,delta_plot,lf_plot,nrow = 1)
    
    
    #---------------------------------------------------------------------------------
    ##-------------------------------- Simulate data ---------------------------------
    #---------------------------------------------------------------------------------
    
    if(sampling == "simulate"){
      
      source("Scripts/source/simulate_real.R")
      ObsM <- F
      y_ObsM_i <- NULL
      index_ObsM_i <- NULL
      
    }
    
    #-----------------------------------------------------------------------------------
    ##----------------------------------- Fit Model ------------------------------------ 
    #-----------------------------------------------------------------------------------
    
    logSigma_com <- log(0)+0.1 # to initialize 
    b <- 0
    SP_est <- 0 # Do not account for sampling process
    b_est <- 0 # Do not account for preferential sampling
    eta_est <- 0 # Do not account for sampling processes other than preferential sampling
    lf_param <- "cov"
    lf_param_num <- ifelse(lf_param=="cov",0,1)
    obs_mod <- 2
    consistency_check <- T
    
    ## Mesh
    mesh <- inla.mesh.create( loc_x[,c('x','y')] )
    spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]  # define sparse matrices for parametrisation of precision matrix
    n_cells <- nrow(loc_x)
    Data_source <- c("Integrated","Scientific","Commercial")
    Obj_step.est <- NULL
    step.est <- F
    quadratic_cov <- F
    plot_outputs <- F
    
    if(sampling == "from_tacsateflalo"){
      if(quadratic_cov == T){
        cov_x_2 <- cov_x
        bathy_c <- cov_x[,"bathy"]
        bathy_c[which(bathy_c < -200)] <- -110
        bathy_c <- (bathy_c - mean(bathy_c))/sd(bathy_c)
        cov_x_2$bathy <- bathy_c
        cov_x_2$bathy2 <- (bathy_c)^2
        # cov_x_2 <- cov_x_2[,c("bathy","bathy2")] # "bathy"
        cov_x_2 <- cov_x_2[,c("substr_Sand_Coarse_substrate","substr_Rock","substr_Mud_sediment")]
        cov_x_2 <- cov_x_2 %>%
          mutate(substr_Sand_Coarse_substrate = ifelse(substr_Rock==1,1,substr_Sand_Coarse_substrate)) %>%
          dplyr::select(-substr_Rock)
      }else{
        cov_x_2 <- cov_x[,c("substr_Sand_Coarse_substrate","substr_Rock","substr_Mud_sediment")]
        cov_x_2 <- cov_x_2 %>%
          mutate(substr_Sand_Coarse_substrate = ifelse(substr_Rock==1,1,substr_Sand_Coarse_substrate)) %>%
          dplyr::select(-substr_Rock)
      }
    }else if(sampling == "simulate"){
      
      cov_x_2 <- cov_x
      
    }
    
    for(aggreg_obs in c(T)){
      
      # for(aggreg_obs in c(F,T)){
      
      if(plot_outputs){
        x11(width = 20,height = 15)
        par(mfrow = c(3,4))
      }
      
      for(Estimation_model_i in c(1,3)){
        
        print(paste0("sim: ",i," | counter: ",counter," | aggreg_obs :",aggreg_obs," | Estimation_model:", Estimation_model_i))
        
        # =T
        # Estimation_model_i=1
        do_step.est <- F
        Params_step.est <- NULL
        Map_step.est <- NULL
        step.est <- 0
        if(step.est==0 & do_step.est) aggreg_obs <- F
        
        ## Fit model
        source("Scripts/function/fit_IM.R")
        TmbFile = "Scripts/"
        TMB::compile(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.cpp"),"-O1 -g",DLLFLAGS="")
        dyn.load( dynlib(paste0("Scripts/inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )
        
        fit_IM_res <- fit_IM(Estimation_model_i = Estimation_model_i, # = 1
                             Samp_process = 0,
                             EM = "fix_b",
                             TmbFile = "Scripts/",
                             ignore.uncertainty = F,
                             c_com_x = c_com_x,
                             y_com_i = y_i2,
                             index_com_i = index_i,
                             y_sci_i = y_sci_i,
                             index_sci_i = index_sci_i,
                             aggreg_obs=aggreg_obs,
                             boats_number = boats_i,
                             Cov_x = as.matrix(cov_x_2), # NULL, # 
                             ref_level = ref_level,
                             lf_param = "RE", # lf_param,
                             spde=spde,
                             mesh=mesh,
                             n_cells=n_cells,
                             cov_est = T,
                             Params_step.est=Params_step.est,
                             Map_step.est=Map_step.est,
                             ObsM=F,
                             y_ObsM_i=y_ObsM_i,
                             index_ObsM_i=index_ObsM_i,
                             sampling = sampling,
                             landings = T,
                             quadratic_cov = quadratic_cov)
        
        Report <- fit_IM_res$Report
        Opt <- fit_IM_res$Opt
        Obj <- fit_IM_res$Obj
        SD <- fit_IM_res$SD
        
        source("Scripts/source/step_estimation.R")
        
        # if(Estimation_model_i == 1 & aggreg_obs == F & sampling == "simulate") save(file="results/real_simu/no.realloc_int_df.RData",data=fit_IM_res)
        # if(Estimation_model_i == 2 & sampling == "simulate") save(file="results/real_simu/sci_df.RData",data=fit_IM_res)
        # if(Estimation_model_i == 3 & aggreg_obs == F & sampling == "simulate") save(file="results/real_simu/no.realloc_com_df.RData",data=fit_IM_res)
        # if(Estimation_model_i == 1 & aggreg_obs == T & sampling == "simulate") save(file="results/real_simu/realloc_int_df.RData",data=fit_IM_res)
        # if(Estimation_model_i == 3 & aggreg_obs == T & sampling == "simulate") save(file="results/real_simu/realloc_com_df.RData",data=fit_IM_res)
        
        # if(Estimation_model_i == 1 & aggreg_obs == F & sampling == "from_tacsateflalo") save(file="results/case_study/no.realloc_int_df.RData",data=fit_IM_res)
        # if(Estimation_model_i == 1 & aggreg_obs == T & sampling == "from_tacsateflalo") save(file="results/case_study/realloc_int_df.RData",data=fit_IM_res)
        # if(Estimation_model_i == 2 & sampling == "from_tacsateflalo") save(file="results/case_study/sci_df.RData",data=fit_IM_res)
        
        source("Scripts/source/plot_outputs.R")
        
        if(sampling == "simulate"){
          
          source("Scripts/source/save_outputs.R")
          
        }
        
        counter <- counter + 1
        
      }
    }
  }
  
}
