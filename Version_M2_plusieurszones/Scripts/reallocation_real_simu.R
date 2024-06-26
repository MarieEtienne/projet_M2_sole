##################################################
## Reallociation - Build simulations on real data
##################################################
## B. Alglave

## Load packages
library(cowplot)
library(dplyr)
library(ggplot2)
library(gt)
library(HistogramTools)
library(INLA)
library(latex2exp)
library(nngeo)
library(raster)
library(rnaturalearth)
library(RNetCDF)
library(sf)
library(spatstat)
library(stringr)
library(TMB)
library(tidyr)

set.seed(2)

## Folder to switch
folder_phd_codes <- "/media/balglave/Elements/backup_phd/phd_zfh_baptiste_alglave"
folder_project <- "/home/balglave/Desktop/Research/projet_M2_sole/Version_M2_plusieurszones"


##-----------------------------------------------------------------------------------
##------------------------ Load raw data and spatial domain -------------------------
##-----------------------------------------------------------------------------------

setwd(folder_phd_codes)

## Map info
domain_shapefile <- "ORHAGO" # "ORHAGO", "EVHOE"
grid_projection <- "+proj=longlat +datum=WGS84"
mapBase <- ne_countries(scale = "large", returnclass = "sf") %>% 
  filter(admin %in% c("France","Spain"))
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
simu_file <- paste0("/media/balglave/Elements/results/severa_rect-",Start_time.tot_2,"_",simu_name,"/")

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
foCatEu6_filter <- "OTT_DEF_>=70_0|OTB_CEP_>=70_0|OTB_DEF_>=70_0"
source("Scripts/source/load_commercial_data.R")


# ## Make figure manuscript
# source("Scripts/In_progress/figure_intro_manuscript.R")


##----------------------------------------------------------------------------------------------
##-------------------- Simulation parameterization and model fitting ---------------------------
##----------------------------------------------------------------------------------------------

## Gaussian field function
source("Scripts/function/sim_GF_Matern.R")

## Fitting function
source("Scripts/function/fit_IM.R")
TmbFile = "Scripts/"
TMB::compile(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.cpp"),"-O1 -g",DLLFLAGS="")
dyn.load( dynlib(paste0("Scripts/inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )
consistency_check <- F

## Mesh
mesh <- inla.mesh.create(loc = as.matrix( loc_x[,c('x','y')] ))
spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]  # define sparse matrices for parametrisation of precision matrix
n_cells <- nrow(loc_x)

## Latent field
intercept <- 2 # intercept of the latent field
beta1 <- 2 # effect of the covariate
nu <- 1
range_cov <- 1.5 # range of the covariate
SD_cov <- 0.5 # marginal variance of the covariate
range_delta <- 0.6 # range of random effect delta
SD_delta <- 1 # marginal variance of random effect delta

## Commercial data
q1 <- 0 # zero-inflation parameter
SD_obs <- 1 # observation error

n_seq_com <- 300 # number of declarations
n_ping_per_seq <- 10 # number of fishing pings per fishing sequence
# n_samp_com = n_seq_com * n_ping_per_seq # number of fishing points
zonesize <- 0.1 # size of the fishin zones
n_zone <- 1 # c(1,3,5) # number of zone visited during a fishing sequence
sample_full_rect <- F # fishing sequence is limited to a smaller zone than a statistical rectangle
k_sim <- 1 # 0: do not reallocate catch, 1: reallocate catch
unsampled_rect <- T # if T, some statistical rectangles are not sampled

## Scientific data
q1_sci <- 0 # zero-inflation parameter
Sigma_sci <- exp(-0.2) # observation error
logSigma_sci <- log(Sigma_sci)
n_samp_sci <- 100 # number of scientific points

## Data sources that feed the model
Data_source <- c("Integrated","Scientific","Commercial")



## If fit real data
##-----------------

quadratic_cov <- F # allow for quadratic form in the effect of the covariate

# Shape covariate data frame
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

##-------------------------------------------------------------------
##-------------------- Loop and file settings -----------------------
##-------------------------------------------------------------------
i0 <- 1 # simulation starting value
counter <- 1 # counter of simulation
simu_file <- paste0("/media/balglave/Elements/results/q1/") # results file
plot_outputs <- F # plot simulation outputs

## Step estimation
do_step.est <- F
Obj_step.est <- NULL
Params_step.est <- NULL
Map_step.est <- NULL
step.est <- F
initialize_params <- T

## Restart after crash
restart_after_crash = F
if(restart_after_crash == T){
  
  simu_file <- "results/severa_rect-2022-01-12_13_03_00_SimuUnsampledRect/"
  load(file=paste0(simu_file,"Results.RData"))
  counter <- max(Results$counter,na.rm=T) + 1
  i0 <- max(Results$sim,na.rm=T) + 1
  
}

##------------------------------------------------------------------------
##-------------------------- Simulation loop -----------------------------
##------------------------------------------------------------------------
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
  
  
  for(Estimation_model_i in c(1:3)){
    
    for(aggreg_obs in c(F,T)){
    
    if(plot_outputs){
      x11(width = 20,height = 15)
      par(mfrow = c(3,4))
    }
    
    
      
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
        
        Params_step.est = list("beta_j"=c(intercept,beta1), 
                               "beta_k"=0, # intercept of fishing intensity
                               "par_b"=0, # link between abundance and sampling intensity
                               "logSigma_com"=log(SD_obs),
                               "logSigma_sci"=logSigma_sci,
                               "q1_sci"=0,
                               "q1_com"=q1,
                               "k_com" = 1,
                               "k_sci" = 1,
                               "deltainput_x"=rep(0,mesh$n),
                               "logtau"=c(0),
                               "logkappa"=c(log(sqrt(8)/range_delta)))
        
        
      }
      
      # Fit model
      fit_IM_res <- fit_IM(Estimation_model_i = Estimation_model_i, # = 1
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
                           landings = T,
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
                           quadratic_cov = quadratic_cov)
      
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


##------------------------------------------------------------------------------------
##----------------------------------- Make figures -----------------------------------
##------------------------------------------------------------------------------------

results_file <- "/media/balglave/Elements/backup_phd/projet_M2_sole/Version_M2_plusieurszones/"

## Simulation plots
source("Scripts/source/plot_simulation.R")


## Case study plots
source("Scripts/source/plot_case_study.R")

if(Estimation_model_i == 2){S_x_sci <- fit_IM_res$Report$S_x;fit_IM_res_sci <- fit_IM_res}
if(Estimation_model_i == 1 & aggreg_obs == F){S_x_pred <- fit_IM_res$Report$S_x;fit_IM_res_Two_step <- fit_IM_res}
if(Estimation_model_i == 1 & aggreg_obs == T){S_x_pred_2 <- fit_IM_res$Report$S_x;fit_IM_res_COS <- fit_IM_res}
test1 <- data.frame(loc_x[,c("layer","x","y")],
                    S_x = S_x,
                    S_x_sci = S_x_sci,
                    S_x_TwoStep = S_x_pred,
                    S_x_COS = S_x_pred_2)  %>% 
   pivot_longer(S_x:S_x_COS) %>% 
  rename(model = name) %>% 
  mutate(model = ifelse(model == "S_x","Simulation",model)) %>% 
  mutate(model = ifelse(model == "S_x_sci","Scientific model",model)) %>% 
  mutate(model = ifelse(model == "S_x_TwoStep","Two-step approach",model)) %>% 
  mutate(model = ifelse(model == "S_x_COS","Joint COS model",model))
test1$model <- factor(test1$model,c("Simulation","Scientific model","Two-step approach","Joint COS model"))


ggplot(test1)+
  geom_point(aes(x = x, y = y, col = value),size=0.5)+
  scale_color_distiller(palette = "Spectral")+
  geom_sf(data=mapBase)+
  facet_wrap(.~ model)+
  theme_classic()+
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5, linetype = "solid"))+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)
  
test2 <- data.frame(loc_x[,c("layer","x","y")],
                   rel_bias_sci = (S_x_sci - S_x)/S_x,
                   rel_bias_COS = (S_x_pred_2 - S_x)/S_x,
                   rel_bias_Twostep = (S_x_pred - S_x)/S_x) %>% 
  pivot_longer(rel_bias_sci:rel_bias_Twostep) %>% 
  rename(model = name) %>% 
  mutate(model = ifelse(model == "rel_bias_sci","Scientific model",model)) %>% 
  mutate(model = ifelse(model == "rel_bias_Twostep","Two-step approach",model)) %>% 
  mutate(model = ifelse(model == "rel_bias_COS","Joint COS model",model))
test2$model <- factor(test2$model,c("Scientific model","Two-step approach","Joint COS model"))

ggplot(test2)+
  geom_point(aes(x = x, y = y, col =  value))+
  scale_color_gradient2(low = "blue",
                        midpoint = 0,
                        mid = "white",
                        high = "red",
                        space="Lab")+
  facet_wrap(.~ model)+
  geom_sf(data=mapBase)+
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5, linetype = "solid"))+
  coord_sf(xlim = c(-6,0), ylim = c(43,48+0.25), expand = FALSE)

