############
## Fit model
############

#' @title fit_IM()
#' 
## Model configuration
#' @param Estimation_model_i index of the estimation model (scientific only, commercial only, integrated)
#' @param Samp_process If 1 : the sampling process contribute to the likelihood, else it doesn't
#' @param EM if "est_b" : b is estimated, if "fix_b" : b is fixed to 0
#' 
#' # TMB model version = "com_data"
#' @param Version estimation model version
#' @param TmbFile file for the C++ template
#' @param ignore.uncertainty if TRUE, ignore uncertainty in estimation
#' 
#' # Data/model Inputs
#' 
#' #ce sont les outputs de de commercial data
#' @param c_com_x number of sample in each cell (line) for each fleet (column)
#' @param y_com_i catch data
#' @param index_com_i sampled cell for commercial observation
#' 
#' #ce sont les outputs de scientific data
#' @param y_sci_i scientific observation
#' @param index_sci_i sampled cell for scientific observation
#' 
#' #output de latent field
#' @param Cov_x covariate for species distribution
#' 
#' @return SD 
#' @return Opt 
#' @return Report
#' @return Converge
#' @return Data
#' @return Params
#' @return Map


fit_IM <- function(Estimation_model_i = 1,
                   Samp_process = 1,
                   EM = "est_b",
                   TmbFile,
                   ignore.uncertainty,
                   c_com_x,
                   y_com_i,
                   index_com_i,
                   y_sci_i,
                   index_sci_i,
                   Cov_x){

  #----------------------
  ## Shape models outputs
  #----------------------
  # Option_vec : model configuration
  # Data : list of data (catch data, covariates, Option_vec, spde objects, etc.)
  # Params : list of parameters
  # Map : parameters that are fixed to a specific value
  
  #configuration du modele considéré pour l'estimation 
  Options_vec = c( 'Prior'=0, # (DEPRECATED)
                   'Alpha'=2, # (DEPRECATED)
                   'IncludeDelta'=1, # (DEPRECATED)
                   'IncludeEta'=1, # (DEPRECATED)
                   'SE'=1, # bias correction for S
                   'DataSource' = Estimation_model_i,
                   'DataObs' = 2,  # (DEPRECATED)
                   'SamplingProcess' = Samp_process , # 1 : sampling process is activated, else : it is ignored
                   'zero.infl_model' = 2,  # (DEPRECATED)
                   'commercial_obs' = 1, # (DEPRECATED)
                   'b_constraint' = 2, # (DEPRECATED)
                   'catchability_random' = 0)  # (DEPRECATED)
  
  ## Data & Params
  Map = list() #liste vide
  
  #Pour chacune des trois configurations de modele (1,2 ou 3 ie les deux, scienti
  #only ou commer only) on définit les données d'entrée, les paramètres à estimer 
  if(Estimation_model_i == 1){   # Integrated model (scientific_commercial)
    
    #on prend bien en argument les sorties des fonctions commerciales et scientifiques
    Data = list( "Options_vec"=Options_vec,
                 "c_com_x"=c_com_x,
                 "y_com_i"=y_com_i,
                 "y_sci_i"=y_sci_i,
                 "index_com_i"=index_com_i-1,
                 "index_sci_i"=index_sci_i-1,
                 "Cov_xj"=cbind(1,Cov_x),
                 "Cov_xk"=array(1,c(list(nrow(Cov_x),1))),
                 "q2_sci" =  1,
                 "q2_com" = 1 )
    
    #on initialise tous les parametres a estimer
    Params = list("beta_j"=rep(0,ncol(Data$Cov_xj)), # linear predictor for abundance 
                  "beta_k"=0, # intercept of fishin intensity
                  "par_b"=0, # link between abundance and sampling intensity
                  "logSigma_com"=log(1),
                  "logSigma_sci"=log(1),
                  "q1_sci"=0,
                  "q1_com"=0,
                  "k_com" = 1,
                  "k_sci" = 1)
    
    Map[["k_com"]] <- seq(1:(length(Params$k_com))) # first k is for scientific data
    Map[["k_com"]][1] <- NA # reference level is the first fleet
    Map[["k_com"]] <- factor(Map[["k_com"]])
    
  }else if(Estimation_model_i == 2){ # scientific model (scientific_only)
    
    Data = list( "Options_vec"=Options_vec,
                 "y_sci_i"=y_sci_i,
                 "index_sci_i"=index_sci_i-1,
                 "Cov_xj"=cbind(1,Cov_x),
                 "q2_sci" = 1)
    
    Params = list("beta_j"=rep(0,ncol(Data$Cov_xj)), # linear predictor for abundance 
                  "logSigma_sci"=log(1),
                  "q1_sci"=0,
                  "k_sci" = 1)
    
    Map[["k_sci"]] <- factor(NA)
    
    
  }else if(Estimation_model_i == 3){ # commercial model (commercial_only)
    
    Data = list( "Options_vec"=Options_vec,
                 "c_com_x"=c_com_x,
                 "y_com_i"=y_com_i,
                 "index_com_i"=index_com_i-1,
                 "Cov_xj"=cbind(1,Cov_x),
                 "Cov_xk"=array(1,c(list(nrow(Cov_x),1))),
                 "q2_sci" =  1,
                 "q2_com" = 1)
    
    Params = list("beta_j"=rep(0,ncol(Data$Cov_xj)), # linear predictor for abundance 
                  "beta_k"=0, # intercept of fishin intensity
                  "par_b"=0, # link between abundance and sampling intensity
                  "logSigma_com"=log(1),
                  "q1_com"=0,
                  "k_com" = 1)
    
    Map[["k_com"]] <- seq(1:(length(Params$k_com))) # first k is for scientific data
    Map[["k_com"]][1] <- NA
    Map[["k_com"]] <- factor(Map[["k_com"]])
    
  }
  
  # fix reference level for latent field covariate
  map_beta_j <- seq(1:ncol(Data$Cov_xj))
  if(exists("ref_level")) map_beta_j[which(colnames(Data$Cov_xj) %in% ref_level)] <- NA
  Map[["beta_j"]] = factor(map_beta_j)
  
  # no linkeage between sampling process and biomass field
  if( EM=="fix_b" ) Map[["par_b"]] <- factor(rep(NA,length(Params$par_b)))
  
  #-----------
  ## Run model
  #-----------
  
  Start_time = Sys.time()
  # library(TMB)
  # TMB::compile(paste0(TmbFile,"inst/executables/",Version,"_scientific_commercial.cpp"),"-O1 -g",DLLFLAGS="")
  Obj = MakeADFun( data=Data, parameters=Params, map = Map, silent = TRUE,hessian = T)
  Obj$fn( Obj$par )
  
  # ## Likelihood profile
  # prof <- tmbprofile(Obj,"logSigma_sci")
  # plot(prof)
  # confint(prof)
  
  # Sys.setenv(PATH="%PATH%;C:/Rtools/mingw64/bin;c:/Rtools/bin")
  

  # # For debugging
  # fixwinpath <- function() {
  #   PATH <- Sys.getenv("PATH")
  #   PATH <- paste0(R.home(), "/bin/x64;", PATH)
  #   PATH <- paste0("c:/rtools40/mingw64/bin;", PATH)
  #   Sys.setenv(PATH=PATH)
  # }
  # fixwinpath()
  # shell("where g++")
  # shell("where gdb")
  # shell("where Rterm")
  # source("Scripts/function/MakeADFun_windows_debug.R")
  # library(TMB)
  # TMB::compile(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.cpp"),"-O1 -g",DLLFLAGS="")
  # dyn.load( dynlib(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )
  # MakeADFun_windows_debug(cpp_name = paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple"),  data=Data, parameters=Params)
  # TMB::gdbsource(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.R"),interactive = T) ## Non-interactive
  # dyn.unload( dynlib(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )

  # # Hessian
  #   hessian.fixed <- optimHess(Obj$par,Obj$fn,Obj$gr)
  # chol(hessian.fixed)
  # which(diag(hessian.fixed)[which(names(diag(hessian.fixed)) == "k_com")]<0)
  
  # Run
  #Lower = -Inf
  #Upper = Inf
  Lower = -50  #getting a few cases where -Inf,Inf bounds result in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  Report = Obj$report()
  
  # https://www.stats.bris.ac.uk/R/web/packages/glmmTMB/vignettes/troubleshooting.html
  
  # # If one random effect can be neglected --> rerun the model without this random effect
  # if( any(Report$MargSD_S<0.001) & Opt$convergence == 0){
  #   Which = which(Report$MargSD_S<0.001)
  #   
  #   # fix negligible random effect
  #   Map[["logtau_S"]] <- c(1:length(Report$MargSD_S))
  #   Map[["logtau_S"]][Which] <- NA
  #   Map[["logtau_S"]] <- factor(Map[["logtau_S"]])
  #   
  #   # Map[["logkappa_S"]] <- c(1:length(Report$MargSD_S))
  #   # Map[["logkappa_S"]][Which] <- NA
  #   # Map[["logkappa_S"]] <- factor(Map[["logkappa_S"]])
  #   if(length(Which)==2){
  #     Map[["logkappa_S"]] = factor( c(NA,NA) )
  #   }
  #   
  #   if( any(Which==1) ){
  #     Map[["deltainput_x"]] = factor( rep(NA,length(Params[["deltainput_x"]])) )
  #     Params[["deltainput_x"]][] = 0
  #   }
  #   
  #   if( any(Which==2) & Estimation_model_i != "scientific_only" ){
  #     Map[["etainput_x"]] = factor( rep(NA,length(Params[["etainput_x"]])) )
  #     Params[["etainput_x"]][] = 0
  #   }
  #   
  #   Data$Options_vec[Which+2] = 0
  #   
  #   print(c("Map ",names(Map)))
  #   print(c("Params ",names(Params)))
  #   
  #   # Re-run
  #   if( length(Which)!=2 ) Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE)
  #   if( length(Which)==2 ) Obj = MakeADFun( data=Data, parameters=Params, random=NULL, map=Map, silent=TRUE)
  #   Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
  #   Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  #   
  # }
  
  Converge=Opt$convergence
  
  
  #-------------------------
  ## Compute standard errors
  #-------------------------
  
  # SD  --> very long with catchability and I got NANs
  if(Converge==0){
    Report = Obj$report()
    SD = sdreport( Obj,ignore.parm.uncertainty = F)
    SD$unbiased$value = c("total_abundance"=Report$total_abundance)

  }else{SD = NULL}
  Opt[["run_time"]] = Sys.time()-Start_time
  
  res <- list(SD = SD, Opt = Opt, Report = Report,  Converge = Converge,Data = Data, Params = Params, Map = Map)
  return(res)
}
