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

#C’est bien d’aller regarder le modèle codé en C++

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
  
  #Le Option_vec sert à la configuration du modele considéré pour l'estimation 
  #il y a pleins de choses qui sont "deprecated" : ce sont des éléments qui ne sont
  #pas pris en compte (car le modèle est simplifié ici)
  
  #ATTENTION : mieux vaut ne pas toucher au Options_vec (car relié au code C++)
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
  Map = list() #liste vide : correspond aux paramètres fixés dans l'estimation
  
  #Pour chacune des trois configurations de modele (1,2 ou 3 ie les deux, scienti
  #only ou commer only) on définit les données d'entrée, les paramètres à estimer 
  if(Estimation_model_i == 1){   # Integrated model (scientific_commercial)
    
    #on prend bien en argument les sorties des fonctions commerciales et scientifiques
    Data = list( "Options_vec"=Options_vec,
                 "c_com_x"=c_com_x,
                 "y_com_i"=y_com_i,
                 "y_sci_i"=y_sci_i,
                 "index_com_i"=index_com_i-1,
                 "index_sci_i"=index_sci_i-1, #les -1 permettent de décaler les indices de 1 car en C++ c’est comme en Python, les objets sont numérotés à partir de 0 (et pas à partir de 1 comme dans R)
                 "Cov_xj"=cbind(1,Cov_x),
                 "Cov_xk"=array(1,c(list(nrow(Cov_x),1))),#Cov_xk : pour faire varier b dans l'espace : nous on s'en fiche !
                 "q2_sci" =  1,
                 "q2_com" = 1 )
    
    #les des parametres du modèle a estimer (on les initialise ici)
    Params = list("beta_j"=rep(0,ncol(Data$Cov_xj)), # linear predictor for abundance 
                  "beta_k"=0, # intercept of fishin intensity
                  "par_b"=0, # link between abundance and sampling intensity
                  "logSigma_com"=log(1),
                  "logSigma_sci"=log(1),
                  "q1_sci"=0,
                  "q1_com"=0,
                  "k_com" = 1,
                  "k_sci" = 1)
    #beta_k : effets des covariables qui jouent sur l’échantillonnage, mais ici on n’en prend pas en compte (cf formule 4 du papier de Baptiste)
    #k_com = 1 #capturabilite relative
    
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
  
  #on cree cet objet pour le TMB qui va contenir tout ce dont on a besoin pour 
  #l'estimation
  #il n'est pas necessaire que l'on comprenne cette fonction
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
  #permet de regarder les paramètres problématiques quand l’algorithme n’arrive pas
  #à converger, quand il y a des pbs pour estimer certains paramètres
  
  #   hessian.fixed <- optimHess(Obj$par,Obj$fn,Obj$gr)
  # chol(hessian.fixed)
  # which(diag(hessian.fixed)[which(names(diag(hessian.fixed)) == "k_com")]<0)
  
  
  # Run
  #Lower = -Inf
  #Upper = Inf
  #Lower et Upper permettent de fixer les bornes inférieures et supérieures
  #pour l'estimation des paramètres
  #gamme de valeurs où notre algorithme va chercher à estimer nos paramètres, 
  #c’est la même pour tous les paramètres
  #on peut mettre -inf / inf, c’est pareil, juste là on rajoute une contrainte
  
  Lower = -50  #getting a few cases where -Inf,Inf bounds result in nlminb failure (NaN gradient)
  Upper = 50
  
  
  #La fonction nlminb est l’algorithme qui permet d’optimiser la vraisemblance, 
  #cet algorithme est fait en TMB, avec des méthodes de différentiation automatique
  #pour la descente de gradient
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, iter.max=100000000))
  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  
  #Avec cet objet Report, on va chercher dans les sorties de l'algo qui a convergé 
  #ce qui nous intéresse : valeur des paramètres, distribution du champ latent, 
  #lambda etc
  #Dans le fichier C++ on spécifie les paramètres que l’on veut sortir, 
  #ça correspond à ce $report()
  #Donc ça vaut le coup d’aller voir le C++ pour savoir ce qu’on ressort, 
  #et potentiellement sortir d’autres choses si ça nous intéresse
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
  
  
  #est ce que l'algo a convergé : 0 c'est bon, sinon non
  Converge=Opt$convergence
  
  
  #-------------------------
  ## Compute standard errors
  #-------------------------
  
  # SD  --> very long with catchability and I got NANs
  
  #dans ce if (on rentre dans la boucle si l'algorithme a convergé),
  #on calcule la matrice de variance covariance de notre modèle, on récupère les 
  #écarts types des paramètres/grandeurs du modèle
  if(Converge==0){
    Report = Obj$report()
    SD = sdreport( Obj,ignore.parm.uncertainty = F)
    SD$unbiased$value = c("total_abundance"=Report$total_abundance)
    
  }else{SD = NULL}
  Opt[["run_time"]] = Sys.time()-Start_time
  
  res <- list(SD = SD, Opt = Opt, Report = Report,  Converge = Converge,Data = Data, Params = Params, Map = Map)
  return(res)
}

#map: parametres fixés dans l'estimation