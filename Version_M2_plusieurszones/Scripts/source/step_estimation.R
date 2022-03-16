#################################
## Step estimation initialization
#################################
## B. Alglave

if(Estimation_model_i == 1){
  
  Obj_step.est <- Obj
  init_val <- Obj_step.est$env$last.par.best
  Params <- fit_IM_res$Params
  Map <- fit_IM_res$Map
  
  #-------------------
  ## PARAMETERS UPDATE
  #-------------------
  ## Observation parameters
  Params$logSigma_com <- init_val[which(names(init_val) == "logSigma_com")]+0.00001
  Params$logSigma_sci <- init_val[which(names(init_val) == "logSigma_sci")]
  Params$q1_sci <- init_val[which(names(init_val) == "q1_sci")]
  Params$q1_com <- init_val[which(names(init_val) == "q1_com")]
  
  ## Covariates
  Params$beta_j[which(!is.na(Map$beta_j))] <- init_val[which(names(init_val) == "beta_j")]
  
  ## Spatial random effect structure
  if("logtau" %in% names(init_val)) Params$logtau <- init_val[which(names(init_val) == "logtau")]
  if("logkappa" %in% names(init_val)) Params$logkappa <- init_val[which(names(init_val) == "logkappa")]
  Params$deltainput_x <- init_val[which(names(init_val) == "deltainput_x")]
  
  if(step.est>0) aggreg_obs <- T
  
  #---------------------------------------
  # COVARIATES - 1st PARAMETER TO LET FREE
  #---------------------------------------
  Map[["beta_j"]] <- as.character(Map[["beta_j"]])
  if(step.est == 1) Map[["beta_j"]][1] <- factor(NA)
  if(step.est > 1) Map$beta_j[1] <- 1
  Map[["beta_j"]] <- factor(Map[["beta_j"]])
  
  #-----------------------------------------
  # CATCHABILITY - 2nd PARAMETER TO LET FREE
  #-----------------------------------------
  if("k_com" %in% names(init_val)) Params$k_com <- init_val[which(names(init_val) == "k_com")]
  if("k_sci" %in% names(init_val)) Params$k_sci <- init_val[which(names(init_val) == "k_sci")]
  
  if(step.est %in% c(1,2)){
    Map[["k_com"]] <- factor(NA)
    Map[["k_sci"]] <- factor(NA)
  }
  
  if(step.est > 2){
    # Map[["k_com"]] <- NULL
    # Map[["k_sci"]] <- NULL
  }
  
  #-------------------------------------------------
  # CATCHABILITY - 3rd AND 4th PARAMETER TO LET FREE
  #-------------------------------------------------
  if(step.est %in% c(1,2,3)) Map[["logkappa"]] <- factor(NA)
  if(step.est > 3) Map[["logkappa"]] <- NULL
  if(step.est %in% c(1,2,3,4)) Map[["logtau"]] <- factor(NA)
  if(step.est > 4) Map[["logtau"]] <- NULL
  
  if("par_b" %in% names(init_val)) Params$par_b <- init_val[which(names(init_val) == "par_b")]
  
  Params_step.est <- Params
  Map_step.est <- Map
  step.est <- step.est + 1
  
}


# ##-------------------------------------------------------------------------
# model_input <- list(Estimation_model_i=1, # = 1
#                     Samp_process = 0,
#                     EM = "fix_b",
#                     TmbFile = "Scripts/",
#                     ignore.uncertainty = F,
#                     c_com_x = c_com_x,
#                     y_i2 = y_i2,
#                     index_i = index_i,
#                     y_sci_i = y_sci_i,
#                     index_sci_i = index_sci_i,
#                     aggreg_obs=F,
#                     boats_i = boats_i,
#                     cov_x_2 = as.matrix(cov_x_2), # NULL, #
#                     ref_level = ref_level,
#                     lf_param = "RE", # lf_param,
#                     spde=spde,
#                     mesh=mesh,
#                     n_cells=n_cells,
#                     cov_est = T,
#                     Params_step.est=Params_step.est,
#                     Map_step.est=Map_step.est,
#                     ObsM=T,
#                     y_ObsM_i=y_ObsM_i,
#                     index_ObsM_i=index_ObsM_i,
#                     sampling = sampling,
#                     landings = T,
#                     quadratic_cov = quadratic_cov)
# save(file="results/model_input.RData",data=model_input)
# ##-------------------------------------------------------------------------
# 
# load("results/model_input.RData")
# 
# c_com_x <- model_input$c_com_x
# y_i2 <- model_input$y_i2
# index_i <- model_input$index_i
# y_sci_i <- model_input$y_sci_i
# index_sci_i <- model_input$index_sci_i
# boats_i <- model_input$boats_i
# cov_x_2 <- model_input$cov_x_2
# ref_level <- model_input$ref_level
# spde <- model_input$spde
# mesh <- model_input$mesh
# n_cells <- model_input$n_cells
# spde <- model_input$spde
# y_ObsM_i <- model_input$y_ObsM_i
# index_ObsM_i <- model_input$index_ObsM_i
# sampling <- model_input$sampling
# quadratic_cov <- model_input$quadratic_cov

