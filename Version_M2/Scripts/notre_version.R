##### 1. Initialisation des paramètres

source("Scripts/function/load_packages.R")
simu_name = "test"


## Latent field

# Grid dimension
grid_side_x <- 25
grid_side_y <- 25
grid_dim = c("x"=grid_side_x, "y"=grid_side_y)
n_cells = grid_dim[1]*grid_dim[2]

# Intercepts and covariates of the latent field
beta0 = 2 # intercept
beta = c(0,0,0,0) # covariates
# parameters covariates vector : 4 modalités de la covariable catégorielle

# Correlation structure of the covariate
SD_x = 0.5 # a combien on le met ? 1 ou 0.1 ? cf notre rmd simulation


## Commercial sampling process

# intercept
beta0_fb <- 2 #intercept dans la formule du lambda


## Commercial data

# Number of commercial samples
n_samp_com = 3000
# observation error 
logSigma_com = log(1) # A quoi ca correspond? 
# zero-inflation parameter
q1_com <- 1
# Relative catchability
q2_com <- 1

# Levels of preferential sampling
b <- 1 


## Model configuration

# Data sources feeding the model
Data_source = c("scientific_commercial","scientific_only","commercial_only")

# b estimated or b fixed to 0 in estimation
EM_set = c("fix_b","est_b")
EM = EM_set[2]

# if 1, sampling process is accounted for in estimation
Samp_process = 1

# TMB model version
TmbFile = "Scripts/"

# if F, compute uncertainty for parameters estimates
ignore.uncertainty = T

## Loop indices

# Index for Results
counter <- 1 
# index of the first iteration
i0 <- 1
# Number of simulation
n_sim = 1

RandomSeed = 123456

## Results dataframe

n_cov = 5
colnames_Results <- c("counter","sim","b_true","Data_source","type_b","Alpha","b_est","ObsMod","Sigma_com_true","Sigma_sci_true","Sigma_com","Sigma_sci","n_samp_com","n_samp_sci","N_true","N_est","SD_N","Convergence","LogLik","MSPE_S","k")
Results = data.frame(matrix(NA,1,length(colnames_Results)))
colnames(Results)=colnames_Results

## List for simulated parameters and parameters estimates
List_param <- list()



##### 2. On définit les fonctions qui nous interessent

# When simu crashes --> param to re-run the loop.
# Before re-runing the loop load the last Results_loop list (Results/Simu/)
restart_after_crash = F
if(restart_after_crash == T){
  counter <- max(Results$counter,na.rm=T)
  i0 <- max(Results$sim,na.rm=T) + 1
  n_sim <- n_sim
  Results <- Results
  List_param <- Results_loop$List_param
  simu_file <- "results/com_x_sci_data_14/com_x_sci_data_14-2020-08-29_11_59_26_Lognormal_Nsamp"
}

#ou est ce que latent_fields_simu, latent_field, simu_tool ont été fixés ??
simu_latent_field <- function(loc_x,
                              latent_fields_simu,
                              latent_field,
                              beta0,
                              beta,
                              simu_tool,
                              SD_x){
  
  # Create design matrix for covariates
  Cov_x <- as.matrix(loc_x[,which(str_detect(colnames(loc_x),"strata"))])
  
  # Total abundance
  Strue_x = exp(beta0 + as.matrix(Cov_x) %*% beta)
  res <- list(Cov_x = Cov_x,
              Strue_x = Strue_x,
              beta = beta)
  return(res)
}


# ou est ce que beta_fb, xfb_x, zero.infl_model? 
simu_commercial_data <- function(loc_x,
                                 grid_dim,
                                 n_cells,
                                 Strue_x,
                                 beta0_fb,
                                 beta_fb,
                                 xfb_x,
                                 zero.infl_model,
                                 n_samp_com,
                                 logSigma_com,
                                 q1_com,
                                 q2_com,
                                 b){
  
  #initialisation
  index_com_i_2 = c() #index_com_i : sampled cell for each observation
  y_com_i_2 = c() # y_com_i : catch data
  b_com_i_2 = c() #b_com_i : fleet related to each observation
  c_com_x_2 = c() #c_com_x : number of sample in each cell (line) for each fleet (column)
  
  Strue_fb <- log(Strue_x) # covariate related to abundance in the sampling intensity process
  
  #on rajoute une nouvelle colonne (Strue_fb) à loc_x : il contient désormais les
  #colonnes x,y,cell, 4 col de strat binaires, et la col Strue_fb
  loc_x <- cbind(loc_x,Strue_fb) # add log-scaled abundance to design matrix
  
  # grid to simulate point process : grille remplie de NA de 25 lignes et 25 col
  #qui va contenir les points de peche commerciaux
  Int_ps <- matrix(ncol = grid_dim['x'],nrow = grid_dim['y'])
  
  # sampling intensity in each cell : Int_ps va contenir les intensites lambda
  #pour chaque point
  for(k in 1:n_cells){
    Int_ps[loc_x$y[k],loc_x$x[k]] <- exp(beta0_fb + b*Strue_fb[k])
  }
  
  #on simule le processus de poisson non homogene
  win_df <- owin(xrange=c(0.5,grid_dim['x'] + 0.5), yrange=c(0.5,grid_dim['y'] + 0.5),ystep  = 1,xstep=1) # parameter of the window to simulate point process (fishing shoots)
  X <- rpoint(n_samp_com,as.im(Int_ps,win_df)) # simulate point process
  
  #plot(X) #representation des points de peche commerciaux simules dans la grille
  
  ## Aggregate data / discretise (use Raster, same as VMS data)
  #------------------------------------------------------------
  
  #definition de la grille raster
  grid <- raster(extent(c(0.5,grid_dim['x']+0.5,0.5,grid_dim['y']+0.5))) # create grid
  res(grid) <- 1 # resolution of the grid
  # Transform this raster into a polygon 
  gridpolygon <- rasterToPolygons(grid)
  coord_grid <- as.data.frame(coordinates(gridpolygon)) #combinaisons de v1/v2 : 
  #tous les couples de coordonnées de la grilles
  colnames(coord_grid)[1] <- "x"
  colnames(coord_grid)[2] <- "y"
  coord_grid #les colonnes de coord_grid sont desormais appelées x et y 
  
  
  # # check that loc_x and gridpolygon have same coordinate system
  # library(ggrepel)
  # ggplot(loc_x)+geom_point(aes(x,y))+geom_text(aes(x,y,label = cell), size = 3.5)
  
  #on rajoute la col cell a coord_grid qui contient les numéros de cellules 
  #correspondant a chaque couple (x,y)
  coord_grid <- inner_join(coord_grid,loc_x[,c("x","y","cell")],by = c("x","y"))
  gridpolygon$layer <- coord_grid$cell #on rentre la col des cellules dans le raster
  
  X_1 <- SpatialPoints(coords=cbind(X$x,X$y))
  
  # intersect of grid and VMS data
  X_2 <- over(X_1, gridpolygon) #3000 lignes pour les 3000 points de peche commerciales
  #et une colonne layer (numero de la cellule correspondant a chacun des 3000 points de peche ? )
  colnames(X_2)[1] <- "layer"
  
  #creation de X3 : nombre de peches commerciales pour chacune des 625 cellules 
  X_2 %>%
    dplyr::count(layer) %>%
    arrange(layer) %>%
    data.frame()  -> X_3
  
  
  gridpolygon@data <- full_join(gridpolygon@data,X_3,by=c("layer"))
  
  #creation de c_com_x : nombre de peches commerciales pour chacune des 625 cellules
  #rangé dans l'ordre croissant
  gridpolygon@data %>%
    arrange(layer) %>%
    mutate(n = ifelse(is.na(n), 0,n)) %>%
    dplyr::select(n) %>%
    c() %>% unlist() %>% as.matrix() -> c_com_x
  
  
  ## For each shot simulate data (presence/absence and positive values)
  #--------------------------------------------------------------------
  
  # Generate a vector where number of cell are repeated as many times as the cell was sampled by fishers
  # vecteur des cellules : chaque cellule apparait autant de fois que de tentatives 
  #de peches qu'elle a recu, les cellules sont rangées dans l'ordre croissant
  #par exemple il y a eu 5 tentatives de peche dans la cellule 1 : le vecteur
  #commence par 4 uns puis 2 deux par ex 
  index_com_i <- do.call(c,lapply(1:length(c_com_x),function(k){
    n<-c_com_x[k]   # not clean --> call fishing_shots which is outside the function
    titi <- rep(k,n)
  }))
  
  
  #calcul de la valeur de capture des 3000 tentatives de peche commerciales 
  y_com_i <- do.call(c,lapply(1:length(index_com_i), function(j){
    exp_catch <- q2_com * Strue_x[index_com_i[j]] # expected catch #mu cf formule papier baptiste
    
    # proba of encounter : probabilité d'une peche fructueuse (cf papier baptiste)
    prob_encount <- 1-exp(- exp(q1_com[1]) * exp_catch) # c'est pour ca qu'il y a aucune peche ratée
    abs_com_i <- rbinom(1,1,prob_encount) #vaut 1 : que du succès pour chacun des 3000 points de peche
    
    #c'est l'équivalent de notre log normale pour calculer la quantité péchée
    #expliquer car cela ne correspond pas à la formule du papier 
    if(abs_com_i>0){
      # y_com_i <- rlnorm(1,meanlog = log(exp_catch), sdlog = exp(logSigma_com))
      # shape <- exp(logSigma_com)^-2
      # y_com_i <-  rgamma(1,shape=shape,scale= exp_catch * exp(logSigma_com)^2)
      y_com_i <- exp(rnorm(1,mean = log(exp_catch) - exp(logSigma_com)^2/2, sd = exp(logSigma_com)))
    }else{
      y_com_i <- 0
    }
  })) # q2_com * Strue_x : expected catch     |      q1_com : parameter linking number of zero and density for scientific data
  
  
  # fleet target factor
  #i est le numéro de de la boucle, le num de la simulation, varie entre 1 et 100
  #vecteur qui contient 3000 fois i
  b_com_i <- rep(i,length(y_com_i))
  
  #index_com_i_2 contient les index_com_i (vecteur des cellules : chaque cellule 
  #apparait autant de fois que de tentatives 
  #de peches qu'elle a recu, les cellules sont rangées dans l'ordre croissant
  #par exemple il y a eu 4 tentatives de peche dans la cellule 1 : le vecteur
  #commence par 4 uns) de l'ITERATION i
  index_com_i_2 = c(index_com_i_2,index_com_i) 
  
  #contient les valeurs de peche des 3000 tentatives de peche 
  #de l'iteration i
  y_com_i_2 = c(y_com_i_2,y_com_i)
  capturetotale <- sum(y_com_i_2)
  y_com_i_2_new <- rep(capturetotale/3000, 3000)
  y_com_i_2_new
  
  #vecteur de 3000 * i
  b_com_i_2 = c(b_com_i_2,b_com_i)
  
  #nombre de peches commerciales pour chacune des 625 cellules
  #rangé dans l'ordre croissant pour l'itération i
  c_com_x_2 = cbind(c_com_x_2,c_com_x)
  
  b_com_i_2 <- factor(b_com_i_2)
  
  res <- list(index_com_i = index_com_i_2, y_com_i = y_com_i_2_new, b_com_i = b_com_i_2, c_com_x = c_com_x_2)
  return(res)
}


simu_commercial_scientific <- function(Results,
                                       simu_file,
                                       grid_dim,
                                       n_cells,
                                       beta0,
                                       beta,
                                       range,
                                       nu,
                                       SD_x,
                                       SD_delta,
                                       SD_eta,
                                       n_samp_sci,
                                       logSigma_sci,
                                       q1_sci,
                                       q2_sci,
                                       n_strate,
                                       n_samp_com,
                                       logSigma_com,
                                       q1_com,
                                       q2_com,
                                       b_set,
                                       Data_source,
                                       Samp_process,
                                       EM,
                                       RandomSeed,
                                       Version,
                                       TmbFile,
                                       ignore.uncertainty,
                                       counter,
                                       i,
                                       n_sim)
{
  ################
  ## Simulate data
  ################
  
  set.seed( RandomSeed + i ) # for figures : i = 2
  
  #---------------
  # Construct grid
  #---------------
  
  # create grid and cells + strata
  n_cells <- grid_dim['x'] * grid_dim['y'] # 10*10 = 100
  loc_x = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y']) #toutes les 
  #combinaisons possibles de (x,y)
  loc_x = cbind(loc_x,cell = 1:n_cells) #chaque combinaison possible est une cellule
  #il les numerote de 1 à 625
  diff.dim_y_x <- grid_dim['y'] - grid_dim['x'] #10-10 =0
  loc_x %>%
    mutate(strata_1 = ifelse(x <= n_strate & y > grid_dim['y'] - n_strate , 1, 0)) %>%
    mutate(strata_4 = ifelse(x >= grid_dim['x'] - n_strate & y <= n_strate + diff.dim_y_x, 1, 0), 1, 0) %>%
    mutate(strata_3 = ifelse((x < grid_dim['x'] - n_strate & y <= grid_dim['y'] - n_strate & x + y <= grid_dim['y'] + 1),1,0)) %>%
    mutate(strata_2 = ifelse(strata_1 == 0 & strata_4 == 0 & strata_3 == 0, 1, 0)) %>%
    dplyr::select(x,y,cell,strata_1,strata_2,strata_3,strata_4) -> loc_x
  #a ce stae, loc_x contient une colonne x, une colonne y, une colonne cellule (jusqu'a 100)
  #et 4 colonnes de strates codées de facon binaires : si le point de coordonnées
  #(x,y) est dans la strate valeur 1, s'il est pas dans la strate valeur 0
  
  # # plot strata : representation graphique des 4 strates dans la grille
  #  loc_x %>%
  #    tidyr::pivot_longer(cols = starts_with("strata"),names_to = "strata") %>%
  #    data.frame() %>% filter(value > 0) -> loc_x_plot
  # 
  #  ggplot(loc_x_plot)+
  #  geom_point(aes(x,y,col=strata)) + theme_bw()
  
  
  #---------------
  #  Latent field
  #---------------
  # simulate or load data for abundance distribution
  #l'out put de cette fonction est une liste qui contient Cov_x (matrix for 
  #covariates levels / values), Strue_x (Strue_x Latent field values), 
  #beta (fixed parameters for the abundance distribution equation), 
  #delta_x (random effect for the abundance distribution equation), 
  #eta_x (random effect for the sampling process equation)
  simu_latent_field_outputs <- simu_latent_field(loc_x,
                                                 latent_fields_simu,
                                                 latent_field,
                                                 scientific_data,
                                                 beta0,
                                                 beta,
                                                 simu_tool,
                                                 SD_x)
  
  Cov_x <- simu_latent_field_outputs$Cov_x #Cov_x contient une colonne cov_x
  #correspondant à la covariable continue et les 4 colonnes de strates codées
  #binaires pour l'appartenance du point a la strate ou non
  Strue_x <- simu_latent_field_outputs$Strue_x #contient 625 lignes et une colonne
  #correspondant aux valeurs du champ latent en chacun des points
  beta <- simu_latent_field_outputs$beta # contient le coefficient beta de la 
  #var continue et les coeff des 4 modalités/strates de la covariable catégorielle
  # delta_x <- simu_latent_field_outputs$delta_x
  # eta_x <- simu_latent_field_outputs$eta_x
  
  ##debut de l'ancienne boucle sur b_set
    
  #-----------------
  #  Commercial data
  #-----------------
  ##l'out put de cette fonction est une liste qui contient :
  # index_com_i (cells sampled by the scientific survey),
  # y_com_i (vector filled with 0/1.
  # b_com_i(scientific observations)
  # c_com_x
  simu_commercial_data_outputs <- simu_commercial_data(
    loc_x,
    grid_dim,
    n_cells,
    Strue_x,
    beta0_fb,
    beta_fb,
    xfb_x,
    zero.infl_model,
    n_samp_com,
    logSigma_com,
    q1_com,
    q2_com,
    b)
  
  index_com_i <-  simu_commercial_data_outputs$index_com_i
  y_com_i <- simu_commercial_data_outputs$y_com_i
  c_com_x <- simu_commercial_data_outputs$c_com_x
  
  # print(paste0("% of pos. values : ",length((y_com_i[which(y_com_i > 0)]))/length(y_com_i)))
  
  ############
  ## Fit model
  ############
  
  ##CETTE BOUCLE SERAIT A ENLEVER vu qu'on a que commercial : a voir comment
  #apres avoir compris fitIM
  # Loop on alternative configuration models
  for (Estimation_model in Data_source) {
    #Estimation_model_i prend a chaque simulation (de 1 à 100) alternativement
    #les trois configurations de data source (que scientifique, que commercial
    #ou les deux)
    
    # Estimation_model_i vaut 1 (commercial scientifique),2 (scientifique only)
    #3 (commercial only)
    Estimation_model_i <- which(Data_source == Estimation_model)
    cat(
      paste(
        "counter ",
        counter,
        " | Simulation ",
        i,
        " | b ",
        b ,
        " | EM ",
        EM,
        " | Estimation_model ",
        Estimation_model,
        " | n_samp_com ",
        n_samp_com,
        "\n"
      )
    )
    
    ############
    ## Fit Model
    ############
    fit_IM_res <- fit_IM(
      Estimation_model_i,
      Samp_process,
      EM,
      TmbFile,
      ignore.uncertainty,
      c_com_x,
      y_com_i,
      index_com_i,
      y_sci_i,
      index_sci_i,
      as.matrix(Cov_x[, 1])
    )
    
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
    
    if (Converge == 0) {
      Results[counter, "N_est"] = SD$unbiased$value['total_abundance']
      Results[counter, "SD_N"] = SD$sd[which(names(SD$value) == "total_abundance")]
    }
    
    # All data and info
    list_simu.data.info <- list(
      grid_dim = grid_dim,
      # latent_fields_simu=latent_fields_simu,
      beta0 = beta0,
      beta = beta,
      range = range,
      nu = nu,
      SD_x = SD_x,
      SD_delta = SD_delta,
      n_samp_sci = n_samp_sci,
      logSigma_sci = logSigma_sci,
      q1_sci = q1_sci,
      q2_sci = q2_sci,
      n_strate = n_strate,
      n_samp_com = n_samp_com,
      logSigma_com = logSigma_com,
      q1_com = q1_com,
      q2_com = q2_com,
      b_set = b_set,
      Data_source = Data_source,
      Samp_process = Samp_process,
      EM = EM,
      RandomSeed = RandomSeed,
      counter = counter,
      i = i,
      n_sim = n_sim,
      # delta_x=delta_x,
      # eta_x=eta_x,
      Strue_x = Strue_x,
      index_sci_i = index_sci_i,
      c_sci_x = c_sci_x,
      y_sci_i = y_sci_i,
      b = b,
      index_com_i = index_com_i,
      y_com_i = y_com_i,
      # b_com_i=b_com_i,
      c_com_x = c_com_x)
    
    
    # Full outputs list
    List_param <- list(
      data.info = list_simu.data.info,
      Opt_par = Opt,
      SD = SD,
      Report = Report)
    
    
    # Fill Results --> summary of simulation loops
    Results[counter, "counter"] = counter
    Results[counter, "sim"] = i
    Results[counter, "Data_source"] = Estimation_model
    # Results[counter,"ObsMod"]=DataObs
    
    Results[counter, "N_true"] = sum(Strue_x)
    Results[counter, "Convergence"] = Converge
    Results[counter, "LogLik"] = -Opt$objective
    Results[counter, "MSPE_S"] = sum((Strue_x - Report$S_x) ^ 2) / n_cells
    Results[counter, "MSPE_S"] = sum((Strue_x - Report$S_x) ^ 2) / n_cells
    
    MSPE_S_2_df <- cbind(loc_x, Strue_x, Report$S_x) %>%
      dplyr::mutate(S_x = Report$S_x) %>%
      filter(x < 9 & y < 9)
    
    Results[counter, "MSPE_S_2"] <-
      sum((MSPE_S_2_df$Strue_x - MSPE_S_2_df$S_x) ^ 2) / (6 * 6)
    
    # Results[counter,"Alpha"]=Alpha
    
    if (Samp_process == 1 &
        Estimation_model != "scientific_only") {
      Results[counter, "type_b"] = EM
      Results[counter, "b_est"] = Report$par_b
    }
    
    if (Estimation_model == "scientific_commercial" | Estimation_model == "scientific_commercial_q_est") {
      Results[counter, "b_true"] = b
      Results[counter, "Sigma_com"] = Report$Sigma_com
      Results[counter, "Sigma_sci"] = Report$Sigma_sci
      Results[counter, "Bias_Sigma_com"] = Report$Sigma_com - exp(logSigma_com)
      Results[counter, "Bias_Sigma_sci"] = Report$Sigma_sci - exp(logSigma_sci)
      Results[counter, "Sigma_com_true"] = exp(logSigma_com)
      Results[counter, "Sigma_sci_true"] = exp(logSigma_sci)
      # Results[counter,"sci_sampling"]=sci_sampling
      Results[counter, "n_samp_com"] = n_samp_com
      Results[counter, "n_samp_sci"] = n_samp_sci
    } else if (Estimation_model == "scientific_only") {
      Results[counter, "Sigma_sci"] = Report$Sigma_sci
      Results[counter, "Bias_Sigma_sci"] = Report$Sigma_sci - exp(logSigma_sci)
      Results[counter, "Sigma_sci_true"] = exp(logSigma_sci)
      # Results[counter,"sci_sampling"]=sci_sampling
      Results[counter, "n_samp_sci"] = n_samp_sci
    } else if (Estimation_model == "commercial_only") {
      Results[counter, "b_true"] = b
      Results[counter, "Sigma_com"] = Report$Sigma_com
      Results[counter, "Bias_Sigma_com"] = Report$Sigma_com - exp(logSigma_com)
      Results[counter, "Sigma_com_true"] = exp(logSigma_com)
      Results[counter, "n_samp_com"] = n_samp_com
    }
    
    if (Estimation_model == "scientific_commercial_q_est") {
      Results[counter, "k"] = Report$k
    }
    
    # save data
    if (!dir.exists(simu_file)){
      dir.create(simu_file)
      save(List_param,
         file = paste0(simu_file, "/List_param_", counter, ".RData"))
      save(Results, file = paste0(simu_file, "/Results.RData"))
      counter <- counter + 1
    }
  }
  res <- list(Results,List_param,counter)
  return(res)
}


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
  
  # Run
  #Lower = -Inf
  #Upper = Inf
  Lower = -50  #getting a few cases where -Inf,Inf bounds result in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  Report = Obj$report()
  
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




#### Maintenant on fait tourner


####### Compile TMB model ########

TMB::compile(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.cpp"),"-O1 -g",DLLFLAGS="")


# file name for savinf outputs
Start_time.tot = Sys.time()
Start_time.tot_2 <- str_replace_all(Start_time.tot, " ", "_")
Start_time.tot_2 <- str_replace_all(Start_time.tot_2, ":", "_")
simu_file <- paste0("results/com_x_sci_data_14_scientific_commercial_simple-",Start_time.tot_2,"_",simu_name,"/")

# load TMB model
dyn.load( dynlib(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )

for(i in i0:n_sim){
  
  res <- simu_commercial_scientific(Results,
                                    simu_file,
                                    grid_dim,
                                    n_cells,
                                    beta0,
                                    beta,
                                    range,
                                    nu,
                                    SD_x,
                                    SD_delta,
                                    SD_eta,
                                    n_samp_sci,
                                    logSigma_sci,
                                    q1_sci,
                                    q2_sci,
                                    n_strate,
                                    n_samp_com,
                                    logSigma_com,
                                    q1_com,
                                    q2_com,
                                    b_set,
                                    Data_source,
                                    Samp_process,
                                    EM,
                                    RandomSeed,
                                    Version,
                                    TmbFile,
                                    ignore.uncertainty,
                                    counter,
                                    i,
                                    n_sim)
  
  Results <- res[[1]]
  List_param <- res[[2]]
  counter <- res[[3]]
}

dyn.unload( dynlib(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )
