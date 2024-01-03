i0 <- i + 1
for(q1 in c(-3)){
  
  for(i in i0:100){
    
    set.seed( i ) # for figures : i = 2
    
    cov_x <- sim_GF_Matern(as.matrix(loc_x[,c("x","y")]), nu, range_cov, SD_cov^2)[[1]]$y
    delta_x <- sim_GF_Matern(as.matrix(loc_x[,c("x","y")]), nu, range_delta, SD_delta^2)[[1]]$y
    
    S_x <- exp(intercept + cov_x %*% beta1 + delta_x)
    vec <- cov_x %*% beta1
    loc_x_2 <- cbind(loc_x,cov=cov_x,delta=delta_x,S_x=S_x)
    
    
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
    if(sampling == "simulate"){
      
      cov_x_2 <- cov_x
      
    }
    
    for(aggreg_obs in c(F,T)){
      
      if(plot_outputs){
        x11(width = 20,height = 15)
        par(mfrow = c(3,4))
      }
      
      
      for(Estimation_model_i in 1){
        
        skip_to_next <- FALSE
        
        
        print(paste0("sim: ",i," | counter: ",counter," | aggreg_obs :",aggreg_obs," | Estimation_model:", Estimation_model_i))
        
        # Step estimation actualization
        if(step.est==0 & do_step.est){
          
          Params_step.est <- NULL
          Map_step.est <- NULL
          aggreg_obs <- F
          
        }else if(step.est>0 & do_step.est){
          
          Params_step.est <- fit_IM_res$Params_step.est
          Map_step.est <- fit_IM_res$Map_step.est
          aggreg_obs <- F
          
        }
        
        if(initialize_params){
          
          Params_step.est = list("beta_j"=c(intercept,beta1)+rnorm(1,0,0.0001), 
                                 "beta_k"=0, # intercept of fishing intensity
                                 "par_b"=0, # link between abundance and sampling intensity
                                 "logSigma_com"=log(SD_obs)+rnorm(1,0,0.0001),
                                 "logSigma_sci"=logSigma_sci+rnorm(1,0,0.0001),
                                 "q1_sci"=0+rnorm(1,0,0.0001),
                                 "q1_com"=q1+rnorm(1,0,0.0001),
                                 "k_com" = 1+rnorm(1,0,0.0001),
                                 "k_sci" = 1,
                                 "deltainput_x"=rep(0,mesh$n),
                                 "logtau"=c(0)+rnorm(1,0,0.0001),
                                 "logkappa"=c(log(sqrt(8)/range_delta))+rnorm(1,0,0.0001))
          
          
        }
        
        
        # Fit model
        tryCatch(
          fit_IM_res <- fit_IM(Estimation_model_i = Estimation_model_i, # = 1
                               Samp_process = 0,
                               EM = "fix_b",
                               TmbFile = "Scripts/",
                               ignore.uncertainty = T,
                               c_com_x = c_com_x,
                               y_com_i = y_i2,
                               index_com_i = index_i,
                               y_sci_i = y_sci_i,
                               index_sci_i = index_sci_i,
                               aggreg_obs=aggreg_obs,
                               boats_number = boats_i,
                               Cov_x = as.matrix(cov_x_2),
                               ref_level = ref_level,
                               lf_param = "RE",
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
          , error = function(e) { skip_to_next <<- TRUE})
        
        if(skip_to_next) { next }
        
        
        Report <- fit_IM_res$Report
        Opt <- fit_IM_res$Opt
        Obj <- fit_IM_res$Obj
        SD <- fit_IM_res$SD
        
        source("Scripts/source/step_estimation.R")
        
        source("Scripts/source/plot_outputs.R")
        
        if(sampling == "simulate"){
          
          source("Scripts/source/save_outputs.R")
          
        }
        
        counter <- counter + 1
        
      }
      
    }
  }
}
   
