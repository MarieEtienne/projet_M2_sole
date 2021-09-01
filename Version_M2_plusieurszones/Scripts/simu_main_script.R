####################
## Simulation script
####################
# B. Alglave based on Conn et al. (2017)

# source("Scripts/function/load_packages.R")
library(parallel)

# Simulation name
simu_name = "SimuTest"

# Parallel configuration
n_nodes <- 8
n_sim <- 8
run_datarmor <- F
cl <- makeCluster(n_nodes, outfile = "")
for(cluster_nb in 1:n_nodes){
  clusterExport(cl[cluster_nb], c('cluster_nb','n_sim','n_nodes','simu_name','run_datarmor'))
}


x_parallel1 <- clusterEvalQ(cl,{
  
  library(dplyr)
  library(TMB)
  library(stringr)
  library(spatstat)
  library(raster)

  # folder name for codes and data
  r_folder <- "/home1/datahome/balglave/projet_M2_sole/Version_M2_plusieurszones"
  data.res_folder <- "/home1/datawork/balglave/reallocation"
  
  if(run_datarmor) setwd(data.res_folder)
  
  load("results/exploratory_analysis/reallocation/nb_sequence_vs_nb_pings.RData")
  load("results/exploratory_analysis/reallocation/simulation_scenarios.RData")
  
  # file name for savinf outputs
  Start_time.tot = round(Sys.time(), units = "mins")
  Start_time.tot_2 <- str_replace_all(Start_time.tot, " ", "_")
  Start_time.tot_2 <- str_replace_all(Start_time.tot_2, ":", "_")
  if(n_nodes == 1) simu_file <- paste0("results/com_x_sci_data_14_scientific_commercial_simple-",Start_time.tot_2,"_",simu_name,"/")
  if(n_nodes > 1) simu_file <- paste0("results/com_x_sci_data_14_scientific_commercial_simple-",Start_time.tot_2,"_",simu_name,"/",cluster_nb)
  
  if(run_datarmor) setwd(r_folder)
  
  #-----------------------------------------------------------
  #-------------- Simulation/estimation settings -------------
  #-----------------------------------------------------------
  
  #---------------------
  # Simulation scenarios
  #---------------------
  
  ###
  ## Latent field
  ###
  
  # Grid dimension
  grid_side_x <- 25
  grid_side_y <- 25
  grid_dim = c("x"=grid_side_x, "y"=grid_side_y)
  n_cells = grid_dim[1]*grid_dim[2]
  
  # Intercepts and covariates of the latent field
  beta0 = 3 # intercept
  beta = c(0,0,0,0,0) # covariates
  # parameters covariates vector : 1st param : continuous variable, 4 nexts : strata_1, strata_2, strata_3, strata_4
  # La premiere covariable est continue et a une structure spatiale
  # Afin de simuler cette corrélation spatiale on a besoin de deux paramètres :
  # - la portée de cette corrélation (nu) (+ la portée est grande et + la distance sur 
  #   laquelle les points sont corrélées entre eux est importante, et donc plus la 
  #   distance est grande, plus il faut "s'éloigner" pour avoir des points indépendants)
  # - la variance marginale de ce processus (écart type SD_x) (cad à quel point ce signal varie)
  
  
  # Correlation structure of the continuous covariate 
  simu_tool <- "MaternCov" # 2 way to simulate latent fields : "RandomFields", "MaternCov"
  range = sqrt(prod(grid_dim))/(5)*2 # --> MaternCov
  nu = 1  #portée de la covariable continue
  SD_x = 2 # ecart type/ variance marginale du processus de la var continue
  # Initialement on avait SD_x = 0.5, la on passe a SD_x = 2 pour avoir plus de variabilité dans le champ latent
  
  #effets spatiaux aléatoires
  SD_delta = SD_eta= 0.5
  
  
  ###
  ## Commercial sampling process
  ###
  
  # intercept pour l'intensite lambda
  beta0_fb = 2 
  
  ###
  ## Observation process
  ###
  
  ## Scientific data
  # Number of scientific data
  n_samp_sci = 50
  # observation error 
  logSigma_sci = log(1)
  # zero-inflation parameter
  q1_sci <- 1
  # Relative catchability
  q2_sci <- 1
  # parameter controling size of the strata
  n_strate <- 9
  
  ## Commercial data
  # Number of commercial samples
  # n_samp_com = 3000
  n_samp_com = 3000
  # observation error 
  logSigma_com = log(0.1)
  # zero-inflation parameter
  q1_com <- -3
  # Relative catchability
  q2_com <- 1
  # Levels of preferential sampling
  b_set = c(0,1,3)
  
  
  
  #--------------------
  # Model configuration
  #--------------------
  
  # Data sources feeding the model
  Data_source = c("scientific_commercial","scientific_only","commercial_only")
  # Data_source = c("commercial_only") # On ne fait plus que le modele commercial
  
  # b estimated or b fixed to 0 in estimation
  EM_set = c("fix_b","est_b")
  EM = EM_set[2]
  
  aggreg_obs = F
  
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
  # i0 <- 1
  # Number of simulation
  # n_sim = 100
  
  RandomSeed = 123456

  ##------------------------------------------------------------
  ##---------------------- Scenarios ---------------------------
  ##------------------------------------------------------------
  # Reallocation uniforme ?
  reallocation = c(0,1)
  # Reallocation = 0 si pas de reallocation et structuration, = 1 si reallocation et structuration
  
  # Nombre de bateaux ?
  n_samp_com_vect = c(50, 100, 1000)
  sequencesdepeche_vect = c(1, 10, 100)
  zonespersequence_vect = c(1, 3, 5)
  NPZ <- expand.grid(sequencesdepeche = sequencesdepeche_vect,
                     zonespersequence = zonespersequence_vect)
  NPZ <- cbind(n_samp_com=n_samp_com_vect,
               NPZ)
  
  
  # PZ = list(c(sequencesdepeche_vect[1], zonespersequence_vect[1]),
  #           c(sequencesdepeche_vect[1], zonespersequence_vect[2]),
  #           c(sequencesdepeche_vect[1], zonespersequence_vect[3]),
  #           c(sequencesdepeche_vect[2], zonespersequence_vect[1]),
  #           c(sequencesdepeche_vect[2], zonespersequence_vect[2]),
  #           c(sequencesdepeche_vect[2], zonespersequence_vect[3]),
  #           c(sequencesdepeche_vect[3], zonespersequence_vect[1]),
  #           c(sequencesdepeche_vect[3], zonespersequence_vect[2]),
  #           c(sequencesdepeche_vect[3], zonespersequence_vect[3]),
  #           c(sequencesdepeche_vect[4], zonespersequence_vect[4]))
  # x=1 : P=2 et Z=1
  # x=2 : P=2 et Z=3
  # x=3 : P=2 et Z=5
  # x=4 : P=5 et Z=1
  # x=5 : P=5 et Z=3
  # x=6 : P=5 et Z=5
  # x=7 : P=10 et Z=1
  # x=8 : P=10 et Z=3
  # x=9 : P=10 et Z=5
  # x=10 : P=1 et Z=1
  
  # ggplot(data=n.seq.rect_n.pings.rect,aes(x=n_seq.rect,y=n_pings.rect))+
  #   # geom_density_2d(data=n.seq.rect_n.pings.rect,aes(x=log10(n_seq.rect),y=log10(n_pings.rect)))
  #   stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  #   scale_fill_continuous(low="lavenderblush", high="red")+
  #   # geom_jitter(alpha=0.5, size = 1.1)+
  #   geom_point(alpha = 0.25,shape=16,size=0.25)+
  #   theme_bw()+theme(legend.position = "none")+
  #   geom_point(data = scenarios_df_3,aes(x=round(n_seq.rect_discr),
  #                                        y=round(n_pings.rect_discr),
  #                                        col=factor(tau)),size=3)+
  #   scale_y_continuous(trans = "log10")+
  #   scale_x_continuous(trans = "log10")+
  #   ylab("Number of pings per ICES rectangle")+
  #   xlab("Number of sequence per rectangle")
  
  
  
  ################
  # source Scripts
  ################
  ## Simulation loop function : simu_commercial_scientific()
  source("Scripts/function/commercial_scientific_14.R")
  ## Simulate Matérn field : sim_GF_Matern()
  source("Scripts/function/sim_GF_Matern.R")
  ## Fit model : fit_IM()
  source("Scripts/function/fit_IM.R")
  ## Simulate latent field : simu_latent_field()
  source("Scripts/function/simu_latent_field.R")
  ## Simulate scientific data : simu_scientific_data()
  source("Scripts/function/simu_scientific_data.R")
  ## Simulate commercial data : simu_commercial_data()
  source("Scripts/function/simu_commercial_data.R")
  ## Reallocation uniforme des peches commerciales
  source("Scripts/function/commercial_reallocation_uniforme.R")
  ## Function for plotting elative bias of abundance, of b and MSPE : Plot_Results()
  source("Scripts/function/plot_simu.R")
  ## Les plots des resultats
  source("Scripts/function/plot_results.R")
  
  
  #-------------------------------------------------------------
  #-------------- Dataframes and model compilation -------------
  #-------------------------------------------------------------
  
  ## Results dataframe
  n_cov = 5
  colnames_Results <- c("counter","sim","b_true","Data_source","type_b","Alpha","b_est","ObsMod","Sigma_com_true","Sigma_sci_true","Sigma_com","Sigma_sci","n_samp_com","n_samp_sci","N_true","N_est","SD_N","Convergence","LogLik","MSPE_S","k", "x","reallocation")
  
  Results = data.frame(matrix(NA,1,length(colnames_Results)))
  colnames(Results)=colnames_Results
  
  ## list for simulated parameters and parameters estimates
  List_param <- list()
  
  
  ####### Compile TMB model ########
  
  # https://kaskr.github.io/adcomp/Introduction.html : for TMB details
  # https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html 
  
  TMB::compile(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.cpp"),"-O1 -g",DLLFLAGS="")
  
  # If problems with the PATH
  # Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # if "C:/Rtools/bin" is not in the PATH
  
  
  #-------------------------------------------------------------
  #------------------- Simulation/Estimation loop --------------
  #-------------------------------------------------------------
  
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

  n_fact <- length(b_set) * nrow(NPZ) * length(reallocation)  # 3 level of b, 10 combination of p and z
  i0 <- (cluster_nb - 1) * n_sim / n_nodes + 1
  i1 <- (cluster_nb) * n_sim / n_nodes
  
  
  # load TMB model
  dyn.load( dynlib(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )
  # dyn.unload( dynlib(paste0(TmbFile,"inst/executables/",Version,"_scientific_commercial") ) )
  
  ## loop
  # 1er niveau de boucle : on repete chaque simulation n_sim fois pour avoir de la variabilité
  for(i in i0:n_sim){
    # 2ème niveau de boucle : on teste chaque valeur de b
    for (b in b_set){
      # 3ème niveau de boucle : on teste sans structuration et sans reallocation, avec structuration et sans
      # reallocation, avec structuration et avec reallocation
      for (k in reallocation){
        # 4ème niveau de boucle : on teste tous les couples P*Z
        for (x in 1:nrow(NPZ)){
          
          n_samp_com = NPZ$n_samp_com[x]
          sequencesdepeche = NPZ$sequencesdepeche[x]
          zonespersequence = NPZ$zonespersequence[x]
          
          if (x!=10){
            taillezone = 3
          }else{
            taillezone = 25
          }
          
          
          res <- simu_commercial_scientific(Results,
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
                                            n_nodes)
          
          Results <- res[[1]]
          List_param <- res[[2]]
          counter <- res[[3]]
          
          
        }
      }
    }
  }

dyn.unload( dynlib(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )

})


# Close cluster
stopCluster(cl)


# #-------------------------------------------------------------
# #------------------------ Plot results -----------------------
# #-------------------------------------------------------------
# 
# # D'abord on définit comment accéder aux résultats de la simulation qui nous intéresse
# path_base = "~/ACO/3A/Projet_ingenieur/Projet_soles/GitHub_Baptiste/projet_M2_sole/Version_M2_plusieurszones/results/"
# simu_file = "com_x_sci_data_14_scientific_commercial_simple-2021-01-26_12_22_12_SimuComplete/"
# 
# 
# # On va commencer par générer les différentes map
# # Ces maps seront stockées dans le dossier results_map
# 
# # 1er niveau de boucle : on repete chaque simulation n_sim fois pour avoir de la variabilité
# for(i in i0:10){
#   # 2ème niveau de boucle : on teste chaque valeur de b
#   # On mute ce niveau de boucle car finalement on va varier les b à l'intérieur de la fonction generation_map
#   # for (b in b_set){
#     # 3ème niveau de boucle : on teste avec réallocation et sans réallocation
#     # On mute ce niveau de boucle car finalement on va varier les k à l'intérieur de la fonction generation_map
#     # for (k in reallocation){
#       # 4ème niveau de boucle : on teste tous les couples P*Z
#       for (x in 1:length(PZ)){
#         generation_map(path_base, simu_file, i, b_set, reallocation, x)
#       }
#     # }
#   # }
# }
# 
# 
# # Maintenant, on va créer les différents graphes représentant les métriques de performance
# # 1er niveau de boucle : on va générer 3 graphes par couple PZ, donc par x
# # 1er graphe : biais de l'abondance
# # 2ème graphe : biais de b
# # 3ème graphe : MSPE
# for (couplePZ in 1:(length(PZ)-1)){ # On fait tous les cas sauf le dernier, qui est la situation de reference
#   generation_metriques(path_base, simu_file, couplePZ, b_set)
# }
# 
# 
# 
# # On genere 3 listes :
# # 1. Liste des biais de b :
# # Elle contient 6 tableaux : b = 0, 1, 3 et k = 0, 1
# # Chaque tableau contient le biais median de b pour les 9 situations ZP
# # 2. Liste des biais de N :
# # Idem
# # 3. Liste des MSPE :
# # Idem
# 
# 
# ## On crée les listes
# 
# # On crée la liste des biais de b
# 
# biaisb_b0_k0 = matrix(ncol = 3, nrow = 3)
# biaisb_b0_k0 = as.data.frame(biaisb_b0_k0)
# colnames(biaisb_b0_k0) = sequencesdepeche_vect[-4] # Les colonnes correspondent aux differentes valeurs de P
# rownames(biaisb_b0_k0) = zonespersequence_vect[-4] # Les colonnes correspondent aux differentes valeurs de Z
# 
# biaisb_b1_k0 = matrix(ncol = 3, nrow = 3)
# biaisb_b1_k0 = as.data.frame(biaisb_b1_k0)
# colnames(biaisb_b1_k0) = sequencesdepeche_vect[-4]
# rownames(biaisb_b1_k0) = zonespersequence_vect[-4]
# 
# biaisb_b3_k0 = matrix(ncol = 3, nrow = 3)
# biaisb_b3_k0 = as.data.frame(biaisb_b3_k0)
# colnames(biaisb_b3_k0) = sequencesdepeche_vect[-4]
# rownames(biaisb_b3_k0) = zonespersequence_vect[-4]
# 
# biaisb_b0_k1 = matrix(ncol = 3, nrow = 3)
# biaisb_b0_k1 = as.data.frame(biaisb_b0_k1)
# colnames(biaisb_b0_k1) = sequencesdepeche_vect[-4]
# rownames(biaisb_b0_k1) = zonespersequence_vect[-4]
# 
# biaisb_b1_k1 = matrix(ncol = 3, nrow = 3)
# biaisb_b1_k1 = as.data.frame(biaisb_b1_k1)
# colnames(biaisb_b1_k1) = sequencesdepeche_vect[-4]
# rownames(biaisb_b1_k1) = zonespersequence_vect[-4]
# 
# biaisb_b3_k1 = matrix(ncol = 3, nrow = 3)
# biaisb_b3_k1 = as.data.frame(biaisb_b3_k1)
# colnames(biaisb_b3_k1) = sequencesdepeche_vect[-4]
# rownames(biaisb_b3_k1) = zonespersequence_vect[-4]
# 
# biaisb = list(biaisb_b0_k0, biaisb_b1_k0, biaisb_b3_k0,
#               biaisb_b0_k1, biaisb_b1_k1, biaisb_b3_k1)
# 
# # On crée la liste des biais de N
# 
# biaisN_b0_k0 = matrix(ncol = 3, nrow = 3)
# biaisN_b0_k0 = as.data.frame(biaisN_b0_k0)
# colnames(biaisN_b0_k0) = sequencesdepeche_vect[-4] # Les colonnes correspondent aux differentes valeurs de P
# rownames(biaisN_b0_k0) = zonespersequence_vect[-4] # Les colonnes correspondent aux differentes valeurs de Z
# 
# biaisN_b1_k0 = matrix(ncol = 3, nrow = 3)
# biaisN_b1_k0 = as.data.frame(biaisN_b1_k0)
# colnames(biaisN_b1_k0) = sequencesdepeche_vect[-4]
# rownames(biaisN_b1_k0) = zonespersequence_vect[-4]
# 
# biaisN_b3_k0 = matrix(ncol = 3, nrow = 3)
# biaisN_b3_k0 = as.data.frame(biaisN_b3_k0)
# colnames(biaisN_b3_k0) = sequencesdepeche_vect[-4]
# rownames(biaisN_b3_k0) = zonespersequence_vect[-4]
# 
# biaisN_b0_k1 = matrix(ncol = 3, nrow = 3)
# biaisN_b0_k1 = as.data.frame(biaisN_b0_k1)
# colnames(biaisN_b0_k1) = sequencesdepeche_vect[-4]
# rownames(biaisN_b0_k1) = zonespersequence_vect[-4]
# 
# biaisN_b1_k1 = matrix(ncol = 3, nrow = 3)
# biaisN_b1_k1 = as.data.frame(biaisN_b1_k1)
# colnames(biaisN_b1_k1) = sequencesdepeche_vect[-4]
# rownames(biaisN_b1_k1) = zonespersequence_vect[-4]
# 
# biaisN_b3_k1 = matrix(ncol = 3, nrow = 3)
# biaisN_b3_k1 = as.data.frame(biaisN_b3_k1)
# colnames(biaisN_b3_k1) = sequencesdepeche_vect[-4]
# rownames(biaisN_b3_k1) = zonespersequence_vect[-4]
# 
# biaisN = list(biaisN_b0_k0, biaisN_b1_k0, biaisN_b3_k0,
#               biaisN_b0_k1, biaisN_b1_k1, biaisN_b3_k1)
# 
# # On crée la liste des MSPE
# 
# MSPE_b0_k0 = matrix(ncol = 3, nrow = 3)
# MSPE_b0_k0 = as.data.frame(MSPE_b0_k0)
# colnames(MSPE_b0_k0) = sequencesdepeche_vect[-4] # Les colonnes correspondent aux differentes valeurs de P
# rownames(MSPE_b0_k0) = zonespersequence_vect[-4] # Les colonnes correspondent aux differentes valeurs de Z
# 
# MSPE_b1_k0 = matrix(ncol = 3, nrow = 3)
# MSPE_b1_k0 = as.data.frame(MSPE_b1_k0)
# colnames(MSPE_b1_k0) = sequencesdepeche_vect[-4]
# rownames(MSPE_b1_k0) = zonespersequence_vect[-4]
# 
# MSPE_b3_k0 = matrix(ncol = 3, nrow = 3)
# MSPE_b3_k0 = as.data.frame(MSPE_b3_k0)
# colnames(MSPE_b3_k0) = sequencesdepeche_vect[-4]
# rownames(MSPE_b3_k0) = zonespersequence_vect[-4]
# 
# MSPE_b0_k1 = matrix(ncol = 3, nrow = 3)
# MSPE_b0_k1 = as.data.frame(MSPE_b0_k1)
# colnames(MSPE_b0_k1) = sequencesdepeche_vect[-4]
# rownames(MSPE_b0_k1) = zonespersequence_vect[-4]
# 
# MSPE_b1_k1 = matrix(ncol = 3, nrow = 3)
# MSPE_b1_k1 = as.data.frame(MSPE_b1_k1)
# colnames(MSPE_b1_k1) = sequencesdepeche_vect[-4]
# rownames(MSPE_b1_k1) = zonespersequence_vect[-4]
# 
# MSPE_b3_k1 = matrix(ncol = 3, nrow = 3)
# MSPE_b3_k1 = as.data.frame(MSPE_b3_k1)
# colnames(MSPE_b3_k1) = sequencesdepeche_vect[-4]
# rownames(MSPE_b3_k1) = zonespersequence_vect[-4]
# 
# MSPE = list(MSPE_b0_k0, MSPE_b1_k0, MSPE_b3_k0,
#             MSPE_b0_k1, MSPE_b1_k1, MSPE_b3_k1)
# 
# 
# ## On remplit les différents tableaux
# 
# for (table in (1:6)){
#   if (table == 1){
#     b = b_set[1]
#     k = reallocation[1]
#   } else if (table == 2){
#     b = b_set[2]
#     k = reallocation[1]
#   } else if (table == 3) {
#     b = b_set[3]
#     k = reallocation[1]
#   } else if (table == 4){
#     b = b_set[1]
#     k = reallocation[2]
#   } else if (table == 5){
#     b = b_set[2]
#     k = reallocation[2]
#   } else {
#     b = b_set[3]
#     k = reallocation[2]
#   }
#   
#   res = indicateurs(path_base, simu_file, table, b, k, biaisb, biaisN, MSPE)
#   biaisb = res[[1]]
#   biaisN = res[[2]]
#   MSPE = res[[3]]
# }
