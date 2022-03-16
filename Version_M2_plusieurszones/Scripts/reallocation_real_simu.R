##################################################
## Reallociation - Build simulations on real data
##################################################
## B. Alglave

memory.size(1e5)

library(cowplot)
library(dplyr)
library(ggplot2)
library(raster)
library(rnaturalearth)
library(RNetCDF)
library(sf)
library(spatstat)
library(stringr)
library(TMB)
library(tidyr)

set.seed(2) # 3

folder_phd_codes <- "C:/R_projects/phd_zfh_baptiste_alglave"
folder_project <- "C:/R_projects/projet_M2_sole/Version_M2_plusieurszones"

setwd(folder_phd_codes)

##--------------
## Load raw data
##--------------
## Map info
domain_shapefile <- "ORHAGO" # "ORHAGO", "EVHOE"
grid_projection <- "+proj=longlat +datum=WGS84"
mapBase <- ne_countries(scale = "medium", returnclass = "sf")
grid_xmin <- -6 # Coordinates of study area
grid_xmax <- 0
grid_ymin <- 42
grid_ymax <- 48
grid_limit <- extent(c(grid_xmin,grid_xmax,grid_ymin,grid_ymax))
Resol_map <- 0.05

# Species name
species_to_plot <- "Solea_solea"
# species_to_plot <- "Merluccius_merluccius"

# Time step
if(species_to_plot == "Solea_solea"){
  year_y <- 2017
  month_m <- 11
  quarter_q <- 4
}else if(species_to_plot == "Merluccius_merluccius"){
  year_y <- 2018
  month_m <- 11
  quarter_q <- 4
}

## data file name
sampling <- "simulate" # "simulate": simulate sampling, "from_tacsateflalo": take data points from tacsateflalo
simu_name <- "SimuFullArea"

# file name for saving outputs
Start_time.tot = round(Sys.time(), units = "mins")
Start_time.tot_2 <- str_replace_all(Start_time.tot, " ", "_")
Start_time.tot_2 <- str_replace_all(Start_time.tot_2, ":", "_")
simu_file <- paste0("results/severa_rect-",Start_time.tot_2,"_",simu_name,"/")

commercial_data_folder <- 'data/raw/vms_logbooks/bob_cs/MergeTACSAT_EFLALO/'
ObsMer_file <- 'data/raw/SIH/obsmer'
ObsM <- T


## Load ICES rectangle shapefile
ICES_rect <- st_read("data/raw/shapefile/ices_rect_stats/ICES_Statistical_Rectangles_Eco.shp")
st_crs(ICES_rect) <- grid_projection
ICES_rect_2 <- ICES_rect %>%
  filter(Ecoregion == "Bay of Biscay and the Iberian Coast")

setwd(folder_project)

## Load domain
source("Scripts/source/load_domain.R")

## Add other covariates to datapoint
source("Scripts/function/extract_covariates.R")
cov_file <- "data/raw/Envir.Data/" # file of covariates within 'data' folder
aggreg_substrate = T

setwd(folder_phd_codes)

## Substrate
datapoint_cov <- extract_substrate_rasterCelticBoB(cov_file,datapoint,datapoint_5,aggreg_substrate)

## Bathymetry
datapoint_cov <- extract_bathy_rasterCelticBoB(cov_file,datapoint,datapoint_cov)

## Copernicus database (SST)
year_month <- paste0(year_y,"-",month_m)
copernicus.phys.var <- "thetao"  # "bottomT", "so"
name_var <- "SST"
datapoint_cov <- extract_copernicus.phys_rasterCelticBoB(cov_file,datapoint,datapoint_cov,copernicus.phys.var,name_var,year_month,grid_projection)

## Copernicus database (bottom temperature)
copernicus.phys.var <- "bottomT"  # "bottomT", "so"
name_var <- "bottomT"
datapoint_cov <- extract_copernicus.phys_rasterCelticBoB(cov_file,datapoint,datapoint_cov,copernicus.phys.var,name_var,year_month,grid_projection)

## Copernicus database (salinity)
copernicus.phys.var <- "so"  # "bottomT", "so"
name_var <- "salinity"
datapoint_cov <- extract_copernicus.phys_rasterCelticBoB(cov_file,datapoint,datapoint_cov,copernicus.phys.var,name_var,year_month,grid_projection)

## Copernicus database (chlorophylle a)
copernicus.bio.var <- "chl"
name_var <- "chl.a"
datapoint_cov <- extract_copernicus.bio_rasterCelticBoB(cov_file,datapoint,datapoint_cov,copernicus.bio.var,name_var,year_month,grid_projection)
useless_levels <- c("NA","Seabed")
normalize_covariates <- T

## Reshape datapoint_cov
latent_field_cov <- c("bathy") # ,"substr","SST","salinity", 
discret_cov <- c("aggreg_substrate","strata") # "bathy"
if(species_to_plot == "Solea_solea") bathy_breaks <- c(-Inf,-50,Inf)
if(species_to_plot == "Merluccius_merluccius") bathy_breaks <- c(-Inf,-100,Inf)
if(species_to_plot == "Solea_solea") ref_level <- c("strata_zs","bathy_(-Inf,-50]","substr_Sand_Coarse_substrate")
if(species_to_plot == "Merluccius_merluccius") ref_level <- c("strata_zs","bathy_(-Inf,-100]","substr_Sand_Coarse_substrate")
setwd(folder_project)
reshape_cov_lf_results <- reshape_cov_lf(datapoint_cov,bathy_breaks,latent_field_cov,useless_levels,discret_cov)
datapoint_cov <- reshape_cov_lf_results$loc_x %>%
  as.data.frame

loc_x <- datapoint_cov %>% mutate(cell = 1:nrow(datapoint_5))
loc_x_sf <- st_as_sf(loc_x, coords = c("x", "y"))
st_crs(loc_x_sf) <- grid_projection

## Cross with statistical rectangle
loc_x_sf <- loc_x_sf[st_intersects(loc_x_sf,ICES_rect_2) %>% lengths > 0,]
loc_x_sf <- st_join(loc_x_sf,ICES_rect_2[,"ICESNAME"])
loc_x_sf$strata <- datapoint_5$strata
loc_x <- loc_x_sf %>%
  as.data.frame %>%
  mutate(x=loc_x$x,
         y=loc_x$y)

## Load tacsateflalo and obsmer
source("Scripts/source/load_commercial_data.R")

## Initialize loop
#-----------------
counter <- 1
i0 <- 1

## Restart after crash
restart_after_crash = F
if(restart_after_crash == T){
  
  simu_file <- "results/severa_rect-2022-01-12_13_03_00_SimuUnsampledRect/"
  load(file=paste0(simu_file,"Results.RData"))
  counter <- max(Results$counter,na.rm=T) + 1
  i0 <- max(Results$sim,na.rm=T) + 1
  
}

##------------------------------------------------------------------------------------
##-------------------- Simulate latent field on the domain ---------------------------
##------------------------------------------------------------------------------------
source("Scripts/function/sim_GF_Matern.R")

## Simulation configuration
#--------------------------
## Latent field
intercept <- 2 # intercept of the latent field
beta1 <- 2 # effect of the covariate
nu <- 1
range_cov <- 1.5 # range of the covariate
SD_cov <- 0.5 # marginal variance of the covariate
range_delta <- 0.6 # range of random effect delta
SD_delta <- 1 # marginal variance of random effect delta

## Commercial data
q1 <- -1 # zero-inflation parameter
SD_obs <- 1 # observation error

n_seq_com <- 300 # number of declarations
n_ping_per_seq <- 10 # number of fishing pings per fishing sequence
# n_samp_com = n_seq_com * n_ping_per_seq # number of fishing points
zonesize <- 0.2 # size of the fishin zones
n_zone <- 1 # c(1,3,5) # number of zone visited during a fishing sequence
sample_full_rect <- F # fishing sequence is limited to a smaller zone than a statistical rectangle
k_sim <- 1 # 0: do not reallocate catch, 1: reallocate catch
unsampled_rect <- T # if T, some statistical rectangles are not sampled

## Scientific data
q1_sci <- 0 # zero-inflation parameter
Sigma_sci <- exp(-0.2) # observation error
logSigma_sci <- log(Sigma_sci)
n_samp_sci <- 100 # number of scientific points

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
  library(INLA)
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
  }
  
  
  for(aggreg_obs in c(F,T)){
    
    if(plot_outputs){
      x11(width = 20,height = 15)
      par(mfrow = c(3,4))
    }
    
    for(Estimation_model_i in c(1:3)){
      
      print(paste0("sim: ",i," | counter: ",counter," | aggreg_obs :",aggreg_obs," | Estimation_model:", Estimation_model_i))
      
      # step.est=T
      # Estimation_model_i=1
      
      Params_step.est <- NULL
      Map_step.est <- NULL
      step.est <- 0
      if(step.est==0) aggreg_obs <- F
      
      ## Fit model
      source("Scripts/function/fit_IM.R")
      TmbFile = "Scripts/"
      TMB::compile(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.cpp"),"-O1 -g",DLLFLAGS="")
      dyn.load( dynlib(paste0("Scripts/inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )
      
      fit_IM_res <- fit_IM(Estimation_model_i = 1, # = 1
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

