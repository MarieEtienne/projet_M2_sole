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
                                       n_sim,
                                       k,
                                       xlim,
                                       ylim,
                                       nboat){
  #on se fiche de l'argument version ici, on pourrait le supprimer
    
  ################
  ## Simulate data
  ################
  
  set.seed( RandomSeed + i ) # for figures : i = 2
  
  #---------------
  # Construct grid
  #---------------
  
  # create grid and cells + strata
  n_cells <- grid_dim['x'] * grid_dim['y'] # 25*25 = 625
  loc_x = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y']) #toutes les 
  #combinaisons possibles de (x,y)
  loc_x = cbind(loc_x,cell = 1:n_cells) #chaque combinaison possible est une cellule
  #il les numerote de 1 à 625
  diff.dim_y_x <- grid_dim['y'] - grid_dim['x'] #25-25 =0
  loc_x %>%
    mutate(strata_1 = ifelse(x <= n_strate & y > grid_dim['y'] - n_strate , 1, 0)) %>%
    mutate(strata_4 = ifelse(x >= grid_dim['x'] - n_strate & y <= n_strate + diff.dim_y_x, 1, 0), 1, 0) %>%
    mutate(strata_3 = ifelse((x < grid_dim['x'] - n_strate & y <= grid_dim['y'] - n_strate & x + y <= grid_dim['y'] + 1),1,0)) %>%
    mutate(strata_2 = ifelse(strata_1 == 0 & strata_4 == 0 & strata_3 == 0, 1, 0)) %>%
    dplyr::select(x,y,cell,strata_1,strata_2,strata_3,strata_4) -> loc_x
  #a ce stade, loc_x contient une colonne x, une colonne y, une colonne cellule (jusqu'a 625)
  #et 4 colonnes de strates codées de facon binaires : si le point de coordonnées
  #(x,y) est dans la strate valeur 1, s'il est pas dans la strate valeur 0
  
  #on a donc dès lors un plan d'échantillonage stratifié, le nombre de points d'échantillonage
  #dans chaque strate est proportionnel à la taille de la strate 
  #ce plan d'échantillonage stratifié présente une structure en esclier (c'est le choix de baptiste, représentation assez réaliste)
  #on peut visualiser ce plan d'échantillonage avec le plot ci dessous
  
 # # plot strata : representation graphique des 4 strates dans la grille 
  loc_x %>%
    tidyr::pivot_longer(cols = starts_with("strata"),names_to = "strata") %>%
    data.frame() %>% filter(value > 0) -> loc_x_plot

  plot_strates = ggplot(loc_x_plot) + geom_point(aes(x,y,col=strata), size=4) + theme_bw() +
    labs(title = paste0("Strates de la zone etudiee, pour i = ", i, " et k = ", k))

  
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
                                                 range,
                                                 nu,
                                                 SD_x,
                                                 SD_delta,
                                                 SD_eta)
  
  Cov_x <- simu_latent_field_outputs$Cov_x #Cov_x contient une colonne cov_x
  #correspondant à la covariable continue et les 4 colonnes de strates codées
  #binaires pour l'appartenance du point a la strate ou non
  Strue_x <- simu_latent_field_outputs$Strue_x #contient 625 lignes et une colonne
  #correspondant aux valeurs du champ latent en chacun des points
  beta <- simu_latent_field_outputs$beta # contient le coefficient beta de la 
  #var continue et les coeff des 4 modalités/strates de la covariable catégorielle
  # delta_x <- simu_latent_field_outputs$delta_x
  # eta_x <- simu_latent_field_outputs$eta_x
  
  # Représentation graphique du champ latent
  Strue_x_2 = as.data.frame(cbind(x=loc_x$x, y=loc_x$y, Champ_latent=as.vector(Strue_x)))
  plot_latentfield = ggplot(Strue_x_2) + geom_point(aes(x,y,col=Champ_latent), size=4) + theme_bw() +
    scale_color_gradient2(midpoint = mean(Strue_x_2$Champ_latent), low = "#E6F2FC", mid = "#62B4FC",
                          high = "#02182C", space = "Lab" ) +
    labs(title = paste0("Representation du champ latent, pour i = ", i, " et k = ", k))
  
    
  #-----------------
  #  Scientific data
  #-----------------
  # simulate or load data for scientific data
  #l'out put de cette fonction est une liste qui contient :
  #index_sci_i (cells sampled by the scientific survey), 
  #c_sci_x (vector filled with 0/1. 
    #If 1 : the cell has been sampled. If 0 :  the cell has not been sampled.), 
  #y_sci_i (scientific observations)
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
  
  # Représentation graphique des points de prelevement scientifique
  pointsdepeche_sci = as.data.frame(cbind(index_sci_i, y_sci_i))
  pointsdepeche_sci$x = rep(0, length(pointsdepeche_sci$index_sci_i))
  pointsdepeche_sci$y = rep(0, length(pointsdepeche_sci$index_sci_i))
  for (j in 1:length(pointsdepeche_sci$index_sci_i))
  {
    pointsdepeche_sci[j, "x"] = loc_x$x[which(loc_x[, "cell"] == pointsdepeche_sci[j, "index_sci_i"])]
    pointsdepeche_sci[j, "y"] = loc_x$y[which(loc_x[, "cell"] == pointsdepeche_sci[j, "index_sci_i"])]
  }
  plot_scientificdata = ggplot(pointsdepeche_sci) + geom_point(aes(x=x, y=y, col=y_sci_i), size=4) +
    labs(color= "Quantite pechee") + theme_bw() +
    scale_color_gradient2(midpoint = mean(pointsdepeche_sci$y_sci_i), low = "#E6F2FC", mid = "#62B4FC",
                          high = "#02182C", space = "Lab" ) +
    labs(title = paste0("Quantites pechees pour les peches scientifiques, pour i = ", i, " et k = ", k))

  # loop on preferential sampling levels
  #3 valeurs de b (0,1 et 3) ont été initialisées dans simu_main_script
  for(b in b_set){
    
    #-----------------
    #  Commercial data
    #-----------------
    ##l'out put de cette fonction est une liste qui contient :
    # index_com_i (cells sampled by the scientific survey), 
    # y_com_i (vector filled with 0/1. 
    # b_com_i(scientific observations)
    # c_com_x
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
                                                         nboat)
    
    index_com_i <-  simu_commercial_data_outputs$index_com_i
    y_com_i <- simu_commercial_data_outputs$y_com_i
    c_com_x <- simu_commercial_data_outputs$c_com_x
    boats_number <- simu_commercial_data_outputs$boats_number
    x_com <- simu_commercial_data_outputs$x_com
    y_com <- simu_commercial_data_outputs$y_com
    plot_centres <- simu_commercial_data_outputs$plot_centres
    plot_pointsdepechecomperboat <- simu_commercial_data_outputs$plot_pointsdepechecomperboat
    
    if (k>0){ # Si on fait de la reallocation uniforme (cad k = 1)
      y_com_i <- commercial_reallocation_uniforme(k, xlim, ylim, y_com_i, n_samp_com,
                                                  index_com_i, loc_x, nboat, boats_number)
    }
    
    
    # Représentation graphique des points de peche, donnees commerciales
    pointsdepeche_com = as.data.frame(cbind(x = x_com, y = y_com, boats = boats_number, ncell = index_com_i, y_com_i = y_com_i))
    plot_pointsdepechecell = ggplot(pointsdepeche_com) +
      geom_point(aes(x=x, y=y, col=as.factor(boats)), size=4) +
      labs(color= "Bateau") +
      theme_bw() +
      theme(axis.title = element_blank(),legend.position='none')
    if (b==b_set[1]){
      plot_pointsdepechecell <- plot_pointsdepechecell +
        facet_wrap(~paste0("b = ", b_set[1]))
    }else if (b==b_set[2]){
      plot_pointsdepechecell <- plot_pointsdepechecell +
        facet_wrap(~paste0("b = ", b_set[2]))
    } else {
      plot_pointsdepechecell <- plot_pointsdepechecell +
        facet_wrap(~paste0("b = ", b_set[3]))
    }
    
    if (k==0){
      plot_pointsdepecheqte = ggplot(pointsdepeche_com) + geom_point(aes(x=x, y=y, col=y_com_i), size=4) +
        labs(color= "Quantite pechee") +
        theme_bw() +
        scale_color_gradient2(midpoint = mean(pointsdepeche_com$y_com_i), low = "#E6F2FC", mid = "#62B4FC",
                              high = "#02182C", space = "Lab" ) + 
        theme(axis.title = element_blank(),legend.position='none')
      if (b==b_set[1]){
        plot_pointsdepecheqte <- plot_pointsdepecheqte +
          facet_wrap(~paste0("b = ", b_set[1]))
      }else if (b==b_set[2]){
        plot_pointsdepecheqte <- plot_pointsdepecheqte +
          facet_wrap(~paste0("b = ", b_set[1]))
      } else {
        plot_pointsdepecheqte <- plot_pointsdepecheqte +
          facet_wrap(~paste0("b = ", b_set[1]))
      }
        
    } else{
      plot_pointsdepecheqte = ggplot(pointsdepeche_com) + geom_point(aes(x=x, y=y, col=as.factor(round(y_com_i, 1))), size=4) +
        labs(color= "Quantite pechee") + theme_bw() +
        theme(axis.title = element_blank(),legend.position='none')
      if (b==b_set[1]){
        plot_pointsdepecheqte <- plot_pointsdepecheqte +
          facet_wrap(~c("b = ", b_set[1]))
      }else if (b==b_set[2]){
        plot_pointsdepecheqte <- plot_pointsdepecheqte +
          facet_wrap(~c("b = ", b_set[2]))
      } else {
        plot_pointsdepecheqte <- plot_pointsdepecheqte +
          facet_wrap(~c("b = ", b_set[3]))
      }
    }
    
    if (b==b_set[1]) {
      plot_centres_b0 = plot_centres
      plot_pointsdepechecomperboat_b0 = plot_pointsdepechecomperboat
      plot_pointsdepechecell_b0 = plot_pointsdepechecell
      plot_pointsdepecheqte_b0 = plot_pointsdepecheqte
    } else if (b==b_set[2]) {
      plot_centres_b1 = plot_centres
      plot_pointsdepechecomperboat_b1 = plot_pointsdepechecomperboat
      plot_pointsdepechecell_b1 = plot_pointsdepechecell
      plot_pointsdepecheqte_b1 = plot_pointsdepecheqte
    } else {
      plot_centres_b3 = plot_centres
      plot_pointsdepechecomperboat_b3 = plot_pointsdepechecomperboat
      plot_pointsdepechecell_b3 = plot_pointsdepechecell
      plot_pointsdepecheqte_b3 = plot_pointsdepecheqte
    }
    
    # print(paste0("% of pos. values : ",length((y_com_i[which(y_com_i > 0)]))/length(y_com_i)))
    
    ############
    ## Fit model
    ############
    
    # Loop on alternative configuration models
    for(Estimation_model in Data_source){
      
      #Estimation_model_i prend a chaque simulation (de 1 à 100) alternativement
      #les trois configurations de data source (que scientifique, que commercial
      #ou les deux)
      
      # Estimation_model_i vaut 1 (commercial scientifique),2 (scientifique only)
      #3 (commercial only)
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
                           as.matrix(Cov_x[,1]))
      
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
                                  c_com_x=c_com_x)
      
      
      # Full outputs list
      # List_param : données d’entrée + données de sortie
      #les valeurs du champ latent pour chaque simulation sont disponnibles dans 
      #List_param (et non dans Results)
      List_param <- list(data.info = list_simu.data.info,
                         Opt_par = Opt,
                         SD = SD,
                         Report = Report)


      # Fill Results --> summary of simulation loops
      #Results : seulement ce qui nous intéresse
      #Ce Results permet de créer rapidement les plots des métriques de performance
      #(→ MSPE, biais sur l’abondance totale, biais sur b)
      Results[counter,"counter"]=counter
      Results[counter,"sim"]=i
      Results[counter,"Data_source"]=Estimation_model
      # Results[counter,"ObsMod"]=DataObs
      
      Results[counter,"N_true"]=sum(Strue_x)
      Results[counter,"Convergence"]=Converge
      Results[counter,"LogLik"]=-Opt$objective
      Results[counter,"MSPE_S"] = sum((Strue_x - Report$S_x)^2)/n_cells
      Results[counter,"MSPE_S"] = sum((Strue_x - Report$S_x)^2)/n_cells
      
      MSPE_S_2_df <- cbind(loc_x,Strue_x,Report$S_x) %>%
        dplyr::mutate(S_x = Report$S_x) %>%
        filter(x < 9 & y < 9)

      Results[counter,"MSPE_S_2"] <- sum((MSPE_S_2_df$Strue_x - MSPE_S_2_df$S_x)^2)/(6*6)

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
      
      # save data
      if(! dir.exists(simu_file)) dir.create(simu_file)
      save(List_param, file = paste0(simu_file,"/List_param_",counter,".RData"))
      save(Results, file = paste0(simu_file,"/Results.RData"))

      counter <- counter + 1
    }
  }
  
  plot_centres = ggarrange(plot_centres_b0, plot_centres_b1, plot_centres_b3, ncol=3,
                           common.legend=TRUE, legend="right")
  plot_centres <- annotate_figure(plot_centres,
                  top = text_grob(paste0("Centres des zones de peche pour les differentes valeurs de b et pour i = ", i, " et k = ", k, ", pour ", nboat, " bateaux"), face = "bold", size = 20),
                  left = text_grob("y",
                                   rot = 90,
                                   size = 11),
                  bottom = text_grob("x", size = 11))
  
  plot_pointsdepechecomperboat = ggarrange(plot_pointsdepechecomperboat_b0,
                                           plot_pointsdepechecomperboat_b1,
                                           plot_pointsdepechecomperboat_b3,
                                           ncol=3,
                                           common.legend=TRUE, legend="right")
  plot_pointsdepechecomperboat <- annotate_figure(plot_pointsdepechecomperboat,
                                  top = text_grob(paste0("Points de peches commerciaux par bateaux pour les differentes valeurs de b et pour i = ", i, " et k = ", k, ", pour ", nboat, " bateaux"), face = "bold", size = 20),
                                  left = text_grob("y",
                                                   rot = 90,
                                                   size = 11),
                                  bottom = text_grob("x", size = 11))
  
  plot_pointsdepechecell = ggarrange(plot_pointsdepechecell_b0, plot_pointsdepechecell_b1,
                                     plot_pointsdepechecell_b3,
                                     ncol=3,
                                     common.legend=TRUE, legend="right")
  plot_pointsdepechecell <- annotate_figure(plot_pointsdepechecell,
                                                  top = text_grob(paste0("Cellules ou il y a de la peche par bateaux pour les differentes valeurs de b et pour i = ", i, " et k = ", k, ", pour ", nboat, " bateaux"), face = "bold", size = 20),
                                                  left = text_grob("y",
                                                                   rot = 90,
                                                                   size = 11),
                                                  bottom = text_grob("x", size = 11))
  
  plot_pointsdepecheqte = ggarrange(plot_pointsdepecheqte_b0, plot_pointsdepecheqte_b1,
                                    plot_pointsdepecheqte_b3, ncol=3,
                                    common.legend=TRUE, legend="right")
  if (k==0){
    plot_pointsdepecheqte <- annotate_figure(plot_pointsdepecheqte,
                                            top = text_grob(paste0("Quantite pechee par cellule pour les differentes valeurs de b et pour i = ", i, " et k = ", k, ", pour ", nboat, " bateaux"), face = "bold", size = 20),
                                            left = text_grob("y",
                                                             rot = 90,
                                                             size = 11),
                                            bottom = text_grob("x", size = 11))
  }else{
    plot_pointsdepecheqte <- annotate_figure(plot_pointsdepecheqte,
                                             top = text_grob(paste0("Quantite pechee par cellule pour les differentes valeurs de b et pour i = ", i, " et k = ", k, ", pour ", nboat, " bateaux"), face = "bold", size = 20),
                                             left = text_grob("y",
                                                              rot = 90,
                                                              size = 11),
                                             bottom = text_grob("x", size = 11))
  }

  # On sauve les 7 graphes produits pour ce couple de valeurs (k, i)
  if(! dir.exists(simu_file_plots)) dir.create(simu_file_plots)
  
  path = paste0(simu_file_plots, "strates_i",i,"_k",k,"_boat",nboat,".png")
  png(file = path, width = 600, height = 500)
  plot(plot_strates)
  dev.off()
  
  path = paste0(simu_file_plots, "latentfield_i",i,"_k",k,"_boat",nboat,".png")
  png(file = path, width = 600, height = 500)
  plot(plot_latentfield)
  dev.off()
  
  path = paste0(simu_file_plots, "scientificdata_i",i,"_k",k,"_boat",nboat,".png")
  png(file = path, width = 600, height = 500)
  plot(plot_scientificdata)
  dev.off()
  
  path = paste0(simu_file_plots, "centres_i",i,"_k",k,"_boat",nboat,".png")
  png(file = path, width = 1500, height = 500)
  plot(plot_centres)
  dev.off()
  
  path = paste0(simu_file_plots, "pecheperboat_i",i,"_k",k,"_boat",nboat,".png")
  png(file = path, width = 1500, height = 500)
  plot(plot_pointsdepechecomperboat)
  dev.off()
  
  path = paste0(simu_file_plots, "pechecell_i",i,"_k",k,"_boat",nboat,".png")
  png(file = path, width = 1500, height = 500)
  plot(plot_pointsdepechecell)
  dev.off()
  
  path = paste0(simu_file_plots, "pecheqte_i",i,"_k",k,"_boat",nboat,".png")
  png(file = path, width = 1500, height = 500)
  plot(plot_pointsdepecheqte)
  dev.off()
  
  res <- list(Results,List_param,counter)
  return(res)
}
