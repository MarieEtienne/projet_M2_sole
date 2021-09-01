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
#' @param beta0 intercept of the abundance distribution equation
#' @param beta  fixed parameters of the abundance distribution equation
#' @param range range parameter. Must be provided for the "MaternCov" package.
#' @param nu nu parameter. Must be provided for the "MaternCov" package.
#' @param SD_x Standard deviation for the simulation of covariate related to the abundance distribution equation.
#' 
#' 
#' # Observations
#' # Scientific
#' @param n_samp_sci number number of scientific samples
#' @param logSigma_sci observation error for the scientific data (log scale)
#' @param q1_sci parameter related to the zero-inflation for scientific data
#' @param q2_sci relative catchability of the scientific data
#' # Commercial
#' @param n_samp_com number of commercial samples
#' @param logSigma_com observation error for the scientific data (log scale)
#' @param q1_com parameter related to the zero-inflation for commercial data
#' @param q2_com relative catchability of the commercial data
#' @param b_set preferential sampling parameter levels
#' 
## Model configuration
#' @param Data_source Alternative for the estimation model ("scientific_commercial" : both data are fitted to the model, "scientific_only" : only scientific data is fitted to the model, "commercial_only" : only commercial data is fitted to the model)
#' @param Samp_process If 1 : the sampling process contribute to the likelihood, else it doesn't
#' @param EM if "est_b" : b is estimated, if "fix_b" : b is fixed to 0
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
#' @param RandomSeed 
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
                                       run_datarmor,
                                       data.res_folder,
                                       r_folder,
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
                                       b,
                                       Data_source,
                                       Samp_process,
                                       EM,
                                       aggreg_obs,
                                       RandomSeed,
                                       Version,
                                       TmbFile,
                                       ignore.uncertainty,
                                       counter,
                                       i,
                                       n_sim,
                                       k,
                                       xlim,
                                       ylim,
                                       sequencesdepeche,
                                       zonespersequence,
                                       taillezone,
                                       i0,
                                       i1,
                                       n_fact,
                                       cluster_nb,
                                       n_nodes){

  ################
  ## Simulate data
  ################

  set.seed( RandomSeed + i )
  
  #---------------
  # Construct grid
  #---------------
  
  # create grid and cells + strata
  n_cells <- grid_dim['x'] * grid_dim['y'] # 25*25 = 625
  loc_x = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y']) #toutes les 
  # Combinaisons possibles de (x,y)
  loc_x = cbind(loc_x,cell = 1:n_cells) # Chaque combinaison possible est une cellule, numerotée de 1 à 625
  diff.dim_y_x <- grid_dim['y'] - grid_dim['x'] #25-25 =0
  loc_x %>%
    mutate(strata_1 = ifelse(x <= n_strate & y > grid_dim['y'] - n_strate , 1, 0)) %>%
    mutate(strata_4 = ifelse(x >= grid_dim['x'] - n_strate & y <= n_strate + diff.dim_y_x, 1, 0), 1, 0) %>%
    mutate(strata_3 = ifelse((x < grid_dim['x'] - n_strate & y <= grid_dim['y'] - n_strate & x + y <= grid_dim['y'] + 1),1,0)) %>%
    mutate(strata_2 = ifelse(strata_1 == 0 & strata_4 == 0 & strata_3 == 0, 1, 0)) %>%
    dplyr::select(x,y,cell,strata_1,strata_2,strata_3,strata_4) -> loc_x
  # A ce stade, loc_x contient une colonne x, une colonne y, une colonne cellule (jusqu'a 625)
  # et 4 colonnes de strates codées de facon binaires : si le point de coordonnées
  # (x,y) est dans la strate valeur 1, s'il est pas dans la strate valeur 0
  
  # On a donc dès lors un plan d'échantillonage stratifié, le nombre de points d'échantillonage
  # dans chaque strate est proportionnel à la taille de la strate 
  # Ce plan d'échantillonage stratifié présente une structure en escalier
  # (c'est le choix de baptiste, représentation assez réaliste)
  # On peut visualiser ce plan d'échantillonage avec le plot ci dessous
  
  # # plot_strates : representation graphique des 4 strates dans la grille 
  # loc_x %>%
  #   tidyr::pivot_longer(cols = starts_with("strata"),names_to = "strata") %>%
  #   data.frame() %>% filter(value > 0) -> loc_x_plot
  # 
  # plot_strates = ggplot(loc_x_plot) + geom_point(aes(x,y,col=strata), size=4) + theme_bw() +
  #   labs(title = paste0("Strates de la zone etudiee, pour i = ", i, " et k = ", k))

  
  #---------------
  #  Latent field
  #---------------
  # simulate or load data for abundance distribution
  # L'out put de cette fonction est une liste qui contient :
  # Cov_x (matrix for  covariates levels / values),
  # Strue_x (Strue_x Latent field values), 
  # beta (fixed parameters for the abundance distribution equation), 
  # delta_x (random effect for the abundance distribution equation), 
  # eta_x (random effect for the sampling process equation)
  validation = 0
  wronglatentfield = 1
  
  while (validation == 0)
  {
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
    Strue_x <- simu_latent_field_outputs$Strue_x
    # Contient 625 lignes et une colonne correspondant aux valeurs du champ latent en chacun des points
    if (max(Strue_x) < 400 & sd(Strue_x) > 10){
      validation = 1
    }else{
      wronglatentfield = wronglatentfield + 1
      set.seed( RandomSeed - wronglatentfield )
    }
  }
  

  
  Cov_x <- simu_latent_field_outputs$Cov_x
  # Cov_x contient une colonne cov_x
  # correspondant à la covariable continue et les 4 colonnes de strates codées
  # binaires pour l'appartenance du point a la strate ou non
  beta <- simu_latent_field_outputs$beta
  # Contient le coefficient beta de lavar continue et les coeff des 4 modalités/strates de la
  # covariable catégorielle
  # delta_x <- simu_latent_field_outputs$delta_x
  # eta_x <- simu_latent_field_outputs$eta_x
  
  # Représentation graphique du champ latent
  Strue_x_2 = as.data.frame(cbind(x=loc_x$x, y=loc_x$y, Champ_latent_reel=as.vector(Strue_x)))
  # plot_latentfield = ggplot(Strue_x_2) + geom_point(aes(x,y,col=Champ_latent), size=4) + theme_bw() +
  #   scale_color_gradient2(midpoint = mean(Strue_x_2$Champ_latent), low = "#E6F2FC", mid = "#62B4FC",
  #                         high = "#02182C", space = "Lab" ) +
  #   labs(title = paste0("Representation du champ latent, pour i = ", i, " et k = ", k))
  
  
  #-----------------
  #  Scientific data
  #-----------------
  # On n'a plus besoin des données scientifiques pour estimer le modèle
  # Cela ne nous empeche pas de les garder, ca ne prend pas de temps de calcul important
  # simulate or load data for scientific data
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
  
  # # Représentation graphique des points de prelevement scientifique
  # pointsdepeche_sci = as.data.frame(cbind(index_sci_i, y_sci_i))
  # pointsdepeche_sci$x = rep(0, length(pointsdepeche_sci$index_sci_i))
  # pointsdepeche_sci$y = rep(0, length(pointsdepeche_sci$index_sci_i))
  # for (j in 1:length(pointsdepeche_sci$index_sci_i))
  # {
  #   pointsdepeche_sci[j, "x"] = loc_x$x[which(loc_x[, "cell"] == pointsdepeche_sci[j, "index_sci_i"])]
  #   pointsdepeche_sci[j, "y"] = loc_x$y[which(loc_x[, "cell"] == pointsdepeche_sci[j, "index_sci_i"])]
  # }
  # plot_scientificdata = ggplot(pointsdepeche_sci) + geom_point(aes(x=x, y=y, col=y_sci_i), size=4) +
  #   labs(color= "Quantite pechee") + theme_bw() +
  #   scale_color_gradient2(midpoint = mean(pointsdepeche_sci$y_sci_i), low = "#E6F2FC", mid = "#62B4FC",
  #                         high = "#02182C", space = "Lab" ) +
  #   labs(title = paste0("Quantites pechees pour les peches scientifiques, pour i = ", i, " et k = ", k))
  
  
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
                                                       zero.infl_model,
                                                       n_samp_com,
                                                       logSigma_com,
                                                       q1_com,
                                                       q2_com,
                                                       b,
                                                       sequencesdepeche,
                                                       zonespersequence,
                                                       taillezone)
  
  index_com_i <-  simu_commercial_data_outputs$index_com_i
  y_com_i <- simu_commercial_data_outputs$y_com_i
  c_com_x <- simu_commercial_data_outputs$c_com_x
  boats_number <- simu_commercial_data_outputs$boats_number
  x_com <- simu_commercial_data_outputs$x_com
  y_com <- simu_commercial_data_outputs$y_com
  centres <- simu_commercial_data_outputs$centres
  peche_com_old <- simu_commercial_data_outputs$peche_com_old
  y_com_i_ref <- y_com_i

  # par(mfrow = c(3, 1))
  # plot(Strue_x[index_com_i],(y_com_i_ref),xlab="S_x",ylab="commercial obs. (true)")
  if (k>0){ # Si on fait de la reallocation uniforme (cad k = 1)
    y_com_i <- commercial_reallocation_uniforme(k, xlim, ylim, y_com_i, n_samp_com,
                                                index_com_i, loc_x, sequencesdepeche, boats_number)
  }
  plot(Strue_x[index_com_i],(y_com_i), xlab="S_x",ylab="commercial obs. (reallocated)")
  plot(x = y_com_i_ref,y= y_com_i,col = boats_number,xlab="true commercial obs.",ylab="reallocated commercial obs.")
  
  
  # Représentation graphique des donnees commerciales (4 graphes)
  
  # 1. Centres des zones de peche
  # plot_centres = ggplot(centres) + geom_point(aes(x=x, y=y, col=as.factor(boats)), size=4) +
  #   labs(color= "Bateau") +
  #   theme_bw() +
  #   theme(axis.title = element_blank(),legend.position='none')
  
  # 2. Points de peche
  # plot_pointsdepechecomperboat = ggplot(peche_com_old) + geom_point(aes(x=x, y=y, col=as.factor(boats)), size=4) +
  #   labs(color= "Bateau") +
  #   theme_bw() +
  #   theme(axis.title = element_blank(),legend.position='none')
  
  # 3. Points de peche dans les cellules
  pointsdepeche_com = as.data.frame(cbind(x = x_com, y = y_com, boats = boats_number, ncell = index_com_i, y_com_i = y_com_i))
  # plot_pointsdepechecell = ggplot(pointsdepeche_com) +
  #   geom_point(aes(x=x, y=y, col=as.factor(boats)), size=4) +
  #   labs(color= "Bateau") +
  #   theme_bw() +
  #   theme(axis.title = element_blank(),legend.position='none')
  
  # 4. Quantite pechee
  # if (k==0){
  #   plot_pointsdepecheqte = ggplot(pointsdepeche_com) + geom_point(aes(x=x, y=y, col=y_com_i), size=4) +
  #     labs(color= "Quantite pechee") +
  #     theme_bw() +
  #     scale_color_gradient2(midpoint = mean(pointsdepeche_com$y_com_i), low = "#E6F2FC", mid = "#62B4FC",
  #                           high = "#02182C", space = "Lab" ) + 
  #     theme(axis.title = element_blank(),legend.position='none')
  # } else{
  #   plot_pointsdepecheqte = ggplot(pointsdepeche_com) + geom_point(aes(x=x, y=y, col=as.factor(round(y_com_i, 1))), size=4) +
  #     labs(color= "Quantite pechee") + theme_bw() +
  #     theme(axis.title = element_blank(),legend.position='none')
  # }
  
  # print(paste0("% of pos. values : ",length((y_com_i[which(y_com_i > 0)]))/length(y_com_i)))
  
  
  ############
  ## Fit model
  ############
  
  # Loop on alternative configuration models
  for(Estimation_model in Data_source[3]){ # On ne regarde que le modèle commercial
    
    Estimation_model_i <- which(Data_source == Estimation_model)
    cat(paste("counter ",counter," | Simulation ",i, " | b ",b ," | EM ",EM," | Estimation_model ",Estimation_model," | n_samp_com ",n_samp_com,"\n"))
    
    ############
    ## Fit Model
    ############
    fit_IM_res <- fit_IM(Estimation_model_i,
                         Samp_process,
                         EM,
                         TmbFile,
                         ignore.uncertainty,
                         c_com_x,
                         y_com_i,
                         index_com_i,
                         y_sci_i,
                         index_sci_i,
                         aggreg_obs,
                         boats_number,
                         as.matrix(Cov_x[,1]))
    
    SD <- fit_IM_res$SD
    Report <- fit_IM_res$Report
    Opt <- fit_IM_res$Opt
    Converge <- fit_IM_res$Converge
    champ_latent_estime <- Report$S_x
    
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
    # list_simu.data.info contient toutes les informations sur la simulation
    list_simu.data.info <- list(grid_dim=grid_dim,
                                # latent_fields_simu=latent_fields_simu,
                                beta0=beta0,
                                beta=beta,
                                range=range,nu=nu,
                                SD_x=SD_x,
                                SD_delta=SD_delta,
                                n_samp_sci=n_samp_sci,
                                logSigma_sci=logSigma_sci,
                                q1_sci=q1_sci,
                                q2_sci=q2_sci,
                                n_strate=n_strate,
                                n_samp_com=n_samp_com,
                                logSigma_com=logSigma_com,
                                q1_com=q1_com,
                                q2_com=q2_com,
                                b_set=b_set,
                                Data_source=Data_source,
                                Samp_process=Samp_process,
                                EM=EM,
                                RandomSeed=RandomSeed,
                                counter=counter,
                                i=i,
                                n_sim=n_sim,
                                # delta_x=delta_x,
                                # eta_x=eta_x,
                                Strue_x=Strue_x,
                                index_sci_i=index_sci_i,
                                c_sci_x=c_sci_x,
                                y_sci_i=y_sci_i,
                                b=b,
                                index_com_i=index_com_i,
                                y_com_i=y_com_i,
                                # b_com_i=b_com_i,
                                c_com_x=c_com_x,
                                i0=i0,
                                i1=i1,
                                n_fact=n_fact,
                                cluster_nb=cluster_nb,
                                n_nodes=n_nodes)
    
    
    # Full outputs list
    # List_param : données d’entrée + données de sortie
    # Les valeurs du champ latent pour chaque simulation sont disponibles dans 
    # List_param (et non dans Results)
    List_param <- list(data.info = list_simu.data.info,
                       Opt_par = Opt,
                       SD = SD,
                       Report = Report)
    
    
    # Fill Results --> summary of simulation loops
    # Results : seulement ce qui nous intéresse
    # Ce Results permet de créer rapidement les plots des métriques de performance
    # MSPE  biais sur l’abondance totale, biais sur b
    
    Results[counter,"counter"]=(i0-1)*n_fact+counter
    Results[counter,"sim"]=i
    Results[counter,"Data_source"]=Estimation_model
    # Results[counter,"ObsMod"]=DataObs
    
    Results[counter,"N_true"]=sum(Strue_x)
    Results[counter,"Convergence"]=Converge
    Results[counter,"LogLik"]=-Opt$objective
    Results[counter,"MSPE_S"] = sum((Strue_x - Report$S_x)^2)/n_cells
    
    # Results[counter,"Alpha"]=Alpha
    
    if(Samp_process == 1 & Estimation_model != "scientific_only"){
      Results[counter,"type_b"]=EM
      Results[counter,"b_est"]=Report$par_b
    }
    
    if(Estimation_model == "scientific_commercial" | Estimation_model == "scientific_commercial_q_est"){
      Results[counter,"b_true"]= b
      Results[counter,"Sigma_com"]=Report$Sigma_com
      Results[counter,"Sigma_sci"]=Report$Sigma_sci
      Results[counter,"Bias_Sigma_com"]=Report$Sigma_com - exp(logSigma_com)
      Results[counter,"Bias_Sigma_sci"]=Report$Sigma_sci - exp(logSigma_sci)
      Results[counter,"Sigma_com_true"]=exp(logSigma_com)
      Results[counter,"Sigma_sci_true"]=exp(logSigma_sci)
      # Results[counter,"sci_sampling"]=sci_sampling
      Results[counter,"n_samp_com"]=n_samp_com
      Results[counter,"n_samp_sci"]=n_samp_sci
      
      
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
      
    }
    
    if (Estimation_model == "scientific_commercial_q_est"){
      Results[counter,"k"]=Report$k
      
    }
    
    Results[counter, "x"] = x
    Results[counter, "reallocation"] = k
    
    Results[counter,"core"]=cluster_nb
    
    # save data
    if(run_datarmor) setwd(data.res_folder)
    
    if(! dir.exists(simu_file)) dir.create(simu_file,recursive = T)
    save(List_param, file = paste0(simu_file,"/List_param_",(i0-1)*n_fact+counter,".RData"))
    if(n_nodes == 1) save(Results, file = paste0(simu_file,"/Results.RData"))
    if(n_nodes > 1) save(Results, file = paste0(simu_file,"/Results",cluster_nb,".RData"))
    
    if(run_datarmor) setwd(r_folder)
    
    counter <- counter + 1
  }
  
  # map_results :
  # Liste de 4 dataframes : Strue_x_2, centres, peche_com_old et pointsdepeche_com
  # Ces dataframes serviront à construire les 5 maps associées à chaque simulation :
  # champ latent, centres de peche, points de peche, points de peche par cellule, quantite pechee par point de peche
  Strue_x_2 = cbind(Strue_x_2, "Champ_latent_estime" = champ_latent_estime)
  map_results = list(Strue_x_2, centres, peche_com_old, pointsdepeche_com)
  # On sauvegarde cette liste dans un workspace
  if(! dir.exists(simu_file)) dir.create(simu_file,recursive = T)
  path = paste0(simu_file,"/map_results_i",i,"_b",b,"_k",k,"_x",x,".RData")
  save(map_results, file = path)
  
  res <- list(Results, List_param, counter)
  return(res)
}
