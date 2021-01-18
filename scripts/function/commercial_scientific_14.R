########################################################
## Simulation loop function (simu_commercial_scientific)
########################################################

#' @title simu_commercial_scientific()
#' 
#' @param Results dataframe summarising simulation/estimation
#' @param simu_file simulation file to save 'Results' and outputs
#' 
#' ## Simulations scenarios
#' # Latent field
#' @param grid_dim grid dimension
#' @param n_cells cell number
#' @param latent_fields_simu if == "simulate" : simulate data, if == "load_from_RData" : load already simulated data, to have the same data between different simulation loops (here `latent_field` and `scientific_data`)
#' @param latent_field already simulated data for the abundance distribution equation
#' @param scientific_data already simulated data for the scientific data
#' @param beta0 intercept of the abundance distribution equation
#' @param beta  fixed parameters of the abundance distribution equation
#' @param simu_tool simulation functions ("RandomFields" or "MaternCov")
#' @param SpatialScale spatial scale parameter. Must be provided for the `RandomFields` package
#' @param range range parameter. Must be provided for the "MaternCov" package.
#' @param nu nu parameter. Must be provided for the "MaternCov" package.
#' @param SD_x Standard deviation for the simulation of covariate related to the abundance distribution equation.
#' @param SD_delta Standard deviation for the random effect delta (abundance distribution).
#' @param SD_xfb Standard deviation for the simulation of covariate related to the sampling process.
#' @param SD_eta Standard deviation for the random effect eta (sampling process).
#' 
#' 
#' # Commercial sampling
#' @param beta0_fb intercept of the sampling process
#' @param beta_fb fixed parameters of the sampling process
#' 
#' # Observations
#' @param DataObs Observation model (1 : zero-inflated gamma, 2 : zero-inflated lognormal, 3 : lognormal)
#' @param zero.infl_model Zero-inflated model parametrisation
#' # Scientific
#' @param n_samp_sci number number of scientific samples
#' @param logSigma_sci observation error for the scientific data (log scale)
#' @param q1_sci parameter related to the zero-inflation for scientific data
#' @param q2_sci relative catchability of the scientific data
#' @param sci_sampling type of sampling plan for scientific data ("random_stratified";"fixed")
#' @param n_strate parameter for the size of the strata (n=9)
#' # Commercial
#' @param n_samp_com number of commercial samples
#' @param logSigma_com observation error for the scientific data (log scale)
#' @param q1_com parameter related to the zero-inflation for commercial data
#' @param q2_com relative catchability of the commercial data
#' @param b_set preferential sampling parameter levels
#' @param zone_without_fishing if TRUE, part of the fishing area is not sampled, if FALSE, part of the fishing area is sampled
#' 
## Model configuration
#' @param Data_source Alternative for the estimation model ("scientific_commercial" : both data are fitted to the model, "scientific_only" : only scientific data is fitted to the model, "commercial_only" : only commercial data is fitted to the model)
#' @param Samp_process If 1 : the sampling process contribute to the likelihood, else it doesn't
#' @param EM if "est_b" : b is estimated, if "fix_b" : b is fixed to 0
#' @param weights_com weight related to the commercial data in the likelihood
#' @param commercial_obs if 1 : commercial data contribute to likelihood
#' @param catchability_random if T, catchability are considered as random effects
#' @param b_constraint constraint on b (if b==1 : b is positive, b==2 : b is free)
#' 
#' # models and SPDE objects configuration
#' @param Spatial_model (DEPRECATED) Spatial model estimation could be "SPDE_GMRF" or "ICAR"
#' @param Use_REML (DEPRECATED) Use restricted (or residual) maximum likelihood method (TRUE/FALSE)
#' @param Alpha parameter related to the smoothness of the gaussian fields (default == 2)
#' @param RandomSeed 
#' 
#' # TMB model version = "com_data"
#' @param Version estimation model version
#' @param TmbFile file for the C++ template
#' @param ignore.uncertainty if TRUE, ignore uncertainty in estimation
#' 
#' ## Loop indices
#' @param counter index for the dataframe
#' @param i simulation index
#' @param n_sim ending iteration
#' 
#' @return Results : summary dataframe of simulation/estimation
#' @return List_param : List of model outputs
#' @return counter : index to browse Results and List_param
#' 
#' @author B. Alglave

# S -- true abundance
# Y -- sampled abundance
# lambda -- intensity of sampling process
# delta -- affects density
# eta -- affects sampling probability

# For comments on previous version see "Rmd/Model/notebook_model"
# For an example ==> run the main script

simu_commercial_scientific <- function(Results,
                                       simu_file,
                                       grid_dim,
                                       n_cells,
                                       latent_fields_simu,
                                       latent_field,
                                       scientific_data,
                                       beta0,
                                       beta,
                                       simu_tool,
                                       no_eta,
                                       SpatialScale,
                                       SpatialScale_eta,
                                       range,
                                       range_eta,
                                       nu,
                                       SD_x,
                                       SD_delta,
                                       SD_xfb,
                                       SD_eta,
                                       beta0_fb,
                                       beta_fb,
                                       DataObs,
                                       zero.infl_model,
                                       nb_par_zinfl,
                                       n_samp_sci,
                                       logSigma_sci,
                                       q1_sci,
                                       q2_sci,
                                       sci_sampling,
                                       n_strate,
                                       n_samp_com,
                                       logSigma_com,
                                       q1_com,
                                       q2_com,
                                       b_set,
                                       zone_without_fishing,
                                       Data_source,
                                       Samp_process,
                                       EM,
                                       weights_com,
                                       commercial_obs,
                                       catchability_random,
                                       b_constraint,
                                       Spatial_model,
                                       Use_REML,
                                       Alpha,
                                       RandomSeed,
                                       Version,
                                       TmbFile,
                                       ignore.uncertainty,
                                       counter,
                                       i,
                                       n_sim){
    
  ################
  ## Simulate data
  ################
  
  set.seed( RandomSeed + i ) # for figures : i = 2
  
  #---------------
  # Construct grid
  #---------------
  
  # create grid and cells + strata
  n_cells <- grid_dim['x'] * grid_dim['y']
  loc_x = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
  loc_x = cbind(loc_x,cell = 1:n_cells)
  diff.dim_y_x <- grid_dim['y'] - grid_dim['x']
  loc_x %>%
    mutate(strata_1 = ifelse(x <= n_strate & y > grid_dim['y'] - n_strate , 1, 0)) %>%
    mutate(strata_4 = ifelse(x >= grid_dim['x'] - n_strate & y <= n_strate + diff.dim_y_x, 1, 0), 1, 0) %>%
    mutate(strata_3 = ifelse((x < grid_dim['x'] - n_strate & y <= grid_dim['y'] - n_strate & x + y <= grid_dim['y'] + 1),1,0)) %>%
    mutate(strata_2 = ifelse(strata_1 == 0 & strata_4 == 0 & strata_3 == 0, 1, 0)) %>%
    dplyr::select(x,y,cell,strata_1,strata_2,strata_3,strata_4) -> loc_x
  
  # plot strata
  loc_x %>%
    tidyr::pivot_longer(cols = starts_with("strata"),names_to = "strata") %>%
    data.frame() %>% filter(value > 0) -> loc_x_plot

  ggplot(loc_x_plot)+
  geom_point(aes(x,y,col=strata)) + theme_bw()
  
  
  #---------------
  #  Latent field
  #---------------
  # simulate or load data for abundance distribution
  simu_latent_field_outputs <- simu_latent_field(loc_x,
                                                 latent_fields_simu,
                                                 latent_field,
                                                 scientific_data,
                                                 beta0,
                                                 beta,
                                                 simu_tool,
                                                 range,
                                                 nu,
                                                 SD_x,
                                                 SD_delta,
                                                 SD_eta)
  if(latent_fields_simu == "simulate"){
    
    Cov_x <- simu_latent_field_outputs$Cov_x
    Strue_x <- simu_latent_field_outputs$Strue_x
    beta <- simu_latent_field_outputs$beta
    delta_x <- simu_latent_field_outputs$delta_x
    xfb_x <- simu_latent_field_outputs$xfb_x
    eta_x <- simu_latent_field_outputs$eta_x
    
  }else if(latent_fields_simu == "load_from_RData"){
    
    Cov_x <- latent_field[[i]]$Cov_x
    Strue_x <- latent_field[[i]]$Strue_x
    beta <- latent_field[[i]]$beta
    delta_x <- latent_field[[i]]$delta_x
    xfb_x <- latent_field[[i]]$xfb_x
    eta_x <- latent_field[[i]]$eta_x
    
  }
  
  #-----------------
  #  Scientific data
  #-----------------
  # simulate or load data for scientific data
  
  if(latent_fields_simu == "simulate"){
    
    simu_scientific_data_outputs <- simu_scientific_data(loc_x,
                                                         grid_dim,
                                                         Strue_x,
                                                         zero.infl_model,
                                                         n_samp_sci,
                                                         logSigma_sci,
                                                         q1_sci,
                                                         q2_sci,
                                                         n_strate,
                                                         scientific_data)
    
    index_sci_i <- simu_scientific_data_outputs$index_sci_i
    c_sci_x <- simu_scientific_data_outputs$c_sci_x
    y_sci_i <- simu_scientific_data_outputs$y_sci_i
    
  }else if(latent_fields_simu == "load_from_RData"){
    
    index_sci_i <- scientific_data[[i]]$index_sci_i
    c_sci_x <- scientific_data[[i]]$c_sci_x
    y_sci_i <- scientific_data[[i]]$y_sci_i
    
  }
  
  # loop on preferential sampling levels
  for(b in b_set){
    
    cat(paste("counter ",counter,"Simulation ",i, "b ",b ,"\n"))
    
    #-----------------
    #  Commercial data
    #-----------------
    
    simu_commercial_data_outputs <- simu_commercial_data(loc_x,
                                                         grid_dim,
                                                         n_cells,
                                                         Strue_x,
                                                         beta0_fb,
                                                         beta_fb,
                                                         xfb_x,
                                                         eta_x,
                                                         zero.infl_model,
                                                         n_samp_com,
                                                         logSigma_com,
                                                         q1_com,
                                                         q2_com,
                                                         b)
    
    index_com_i <-  simu_commercial_data_outputs$index_com_i
    y_com_i <- simu_commercial_data_outputs$y_com_i
    b_com_i <- simu_commercial_data_outputs$b_com_i
    VE_i <- simu_commercial_data_outputs$b_com_i
    c_com_x <- simu_commercial_data_outputs$c_com_x
    
    print(paste0("% of pos. values : ",length((y_com_i[which(y_com_i > 0)]))/length(y_com_i)))
    
    ############
    ## Fit model
    ############
    # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
    mesh = inla.mesh.create( loc_x[,c('x','y')] )
    spde <- (inla.spde2.matern(mesh, alpha=Alpha)$param.inla)[c("M0","M1","M2")]  # define sparse matrices for parametrisation of precision matrix
    
    # Loop on alternative configuration models
    for(Estimation_model in Data_source){
      Estimation_model_i <- which(Data_source == Estimation_model)
      cat(paste("counter ",counter," | Simulation ",i, " | b ",b ," | zone_without_fishing ",zone_without_fishing," | EM ",EM," | Estimation_model ",Estimation_model," | n_samp_com ",n_samp_com," | b_constraint ",b_constraint,"\n"))
      
      ############
      ## Fit Model
      ############
      skip_to_next <- F
      tryCatch({fit_IM_res <- fit_IM(Estimation_model_i,
                                     Samp_process,
                                     EM,
                                     TmbFile,
                                     ignore.uncertainty,
                                     c_com_x,
                                     y_com_i,
                                     index_com_i,
                                     y_sci_i,
                                     index_sci_i,
                                     Cov_x)
      
      SD <- fit_IM_res$SD
      Report <- fit_IM_res$Report
      Opt <- fit_IM_res$Opt
      Converge <- fit_IM_res$Converge
      
      # source("Scripts/Simulation/format_outputs.R")
      # format_outputs_res <- format_outputs()
      # Results <- format_outputs_res$Results
      # List_param <- format_outputs_res$List_param
      # 
      #######################
      ## Fill outputs objects
      #######################
      
      if(Converge == 0){
        Results[counter,"N_est"]=SD$unbiased$value['total_abundance']
        Results[counter,"SD_N"]=SD$sd[which(names(SD$value)=="total_abundance")]
      }
      
      # All data and info
      list_simu.data.info <- list(grid_dim=grid_dim,
                                  latent_fields_simu=latent_fields_simu,
                                  beta0=beta0,
                                  beta=beta,
                                  simu_tool=simu_tool,
                                  SpatialScale=SpatialScale,
                                  range=range,nu=nu,
                                  SD_x=SD_x,
                                  SD_delta=SD_delta,
                                  SD_xfb=SD_xfb,
                                  SD_eta=SD_eta,
                                  beta0_fb=beta0_fb,
                                  beta_fb=beta_fb,
                                  DataObs=DataObs,
                                  zero.infl_model=zero.infl_model,
                                  n_samp_sci=n_samp_sci,
                                  logSigma_sci=logSigma_sci,
                                  q1_sci=q1_sci,
                                  q2_sci=q2_sci,
                                  sci_sampling=sci_sampling,
                                  n_strate=n_strate,
                                  n_samp_com=n_samp_com,
                                  logSigma_com=logSigma_com,
                                  q1_com=q1_com,
                                  q2_com=q2_com,
                                  b_set=b_set,
                                  zone_without_fishing=zone_without_fishing,
                                  Data_source=Data_source,
                                  Samp_process=Samp_process,
                                  EM=EM,
                                  weights_com=weights_com,
                                  commercial_obs=commercial_obs,
                                  b_constraint=b_constraint,
                                  Spatial_model=Spatial_model,
                                  Use_REML=Use_REML,
                                  Alpha=Alpha,
                                  RandomSeed=RandomSeed,
                                  Version=Version,
                                  counter=counter,
                                  i=i,
                                  n_sim=n_sim,
                                  Strue_x=Strue_x,
                                  delta_x=delta_x,
                                  xfb_x=xfb_x,
                                  eta_x=eta_x,
                                  index_sci_i=index_sci_i,
                                  c_sci_x=c_sci_x,
                                  y_sci_i=y_sci_i,
                                  b=b,
                                  index_com_i=index_com_i,
                                  y_com_i=y_com_i,
                                  b_com_i=b_com_i,
                                  c_com_x=c_com_x)
      
      
      # Full outputs list
      List_param <- list(data.info = list_simu.data.info,
                         Opt_par = Opt,
                         SD = SD,
                         Report = Report)


      
      Results[counter,"N_true"]=sum(Strue_x)
      Results[counter,"LogLik"]=-Opt$objective
      # Results[counter,"BiasBetaj1"]=Opt$par[1]-beta0
      # Results[counter,"BiasBetaj2"]=Opt$par[2]-beta
      Results[counter,"MSPE_S"] = sum((Strue_x - Report$S_x)^2)/n_cells
      Results[counter,"MSPE_S"] = sum((Strue_x - Report$S_x)^2)/n_cells
      
      MSPE_S_2_df <- cbind(loc_x,Strue_x,Report$S_x) %>%
        dplyr::mutate(S_x = Report$S_x) %>%
        filter(x < 9 & y < 9)

      Results[counter,"MSPE_S_2"] <- sum((MSPE_S_2_df$Strue_x - MSPE_S_2_df$S_x)^2)/(6*6)

      Results[counter,"Alpha"]=Alpha
      
      if(Samp_process == 1 & Estimation_model != "scientific_only"){
        Results[counter,"type_b"]=EM
        Results[counter,"b_est"] = SD$unbiased$value['b']
      }
      
      if(Estimation_model == "scientific_commercial" | Estimation_model == "scientific_commercial_q_est"){
        
        Results[counter,"Sigma_com"]=Report$Sigma_com
        Results[counter,"Sigma_sci"]=Report$Sigma_sci
        Results[counter,"Bias_Sigma_com"]=Report$Sigma_com - exp(logSigma_com)
        Results[counter,"Bias_Sigma_sci"]=Report$Sigma_sci - exp(logSigma_sci)
        Results[counter,"Sigma_com_true"]=exp(logSigma_com)
        Results[counter,"Sigma_sci_true"]=exp(logSigma_sci)
        
      }else if(Estimation_model == "scientific_only"){
        Results[counter,"Sigma_sci"]=Report$Sigma_sci
        Results[counter,"Bias_Sigma_sci"]=Report$Sigma_sci - exp(logSigma_sci)
        Results[counter,"Sigma_sci_true"]=exp(logSigma_sci)
        # Results[counter,"sci_sampling"]=sci_sampling
        Results[counter,"n_samp_sci"]=n_samp_sci
        
        
      }else if(Estimation_model == "commercial_only"){
        Results[counter,"b_true"]= b
        Results[counter,"Sigma_com"]=Report$Sigma_com
        Results[counter,"Bias_Sigma_com"]=Report$Sigma_com - exp(logSigma_com)
        Results[counter,"Sigma_com_true"]=exp(logSigma_com)
        Results[counter,"n_samp_com"]=n_samp_com
        Results[counter,"b_constraint"]=b_constraint
        Results[counter,"zone_without_fishing"]=zone_without_fishing
        
      }
      
      if (Estimation_model == "scientific_commercial_q_est"){
        Results[counter,"k"]=Report$k
        
      }
      
      Results[counter,"Convergence"]=Converge
      
      # save data
      if(!file.exists(simu_file)) dir.create(simu_file)
      save(List_param, file = paste0(simu_file,"/List_param_",counter,".RData"))
      save(Results, file = paste0(simu_file,"/Results.RData"))
      }, error = function(e) { skip_to_next <<- TRUE})
      # Fill Results --> summary of simulation loops
      Results[counter,"counter"]=counter
      Results[counter,"sim"]=i
      Results[counter,"Data_source"]=Estimation_model
      Results[counter,"ObsMod"]=DataObs
      Results[counter,"b_true"]= b
      Results[counter,"sci_sampling"]=sci_sampling
      Results[counter,"n_samp_com"]=n_samp_com
      Results[counter,"n_samp_sci"]=n_samp_sci
      Results[counter,"b_constraint"]=b_constraint
      Results[counter,"zone_without_fishing"]=zone_without_fishing
      
      if(skip_to_next) { next }     
      
      counter <- counter + 1
    }
  }
  res <- list(Results,List_param,counter)
  return(res)
}
