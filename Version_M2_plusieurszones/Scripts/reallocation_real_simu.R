##################################################
## Reallociation - Build simulations on real data
##################################################
## B. Alglave


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
sampling <- "from_tacsateflalo" # "simulate": simulate sampling, "from_tacsateflalo": take data points from tacsateflalo
commercial_data_folder <- 'data/raw/vms_logbooks/bob_cs/MergeTACSAT_EFLALO/'
ObsMer_file <- 'data/raw/SIH/obsmer'
ObsM <- T

## If simulate data points
# n_samp_com <- 5000
n_seq_com <- 300
n_ping_per_seq <- 10

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
latent_field_cov <- c("substr","bathy") # ,"SST","salinity", 
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

##------------------------------------------------------------------------------------
##-------------------- Simulate latent field on the domain ---------------------------
##------------------------------------------------------------------------------------
source("Scripts/function/sim_GF_Matern.R")
intercept <- 2
beta1 <- 2
nu <- 1
range_cov <- 1.5
SD_cov <- 0.5
range_delta <- 0.6
SD_delta <- 1

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


#-----------------------------------------------------------------------------------
##-------------------------------- Commercial data ---------------------------------
#-----------------------------------------------------------------------------------

## From tacsatEflalo and ObsMer
#------------------------------
if(sampling == "from_tacsateflalo"){
  
  ## Load tacsatEflalo
  setwd(folder_phd_codes)
  filter_with_Eflalo <- F
  source("r/real_data/load_data/extract_commercial_data.R")
  
  ## Load ObsMer
  HH_obsmer <- read.csv(paste0(ObsMer_file,"/HH_OBSMER_VIII_VII_2005_2019_230320.csv"),sep = ";",header = T)
  HL_obsmer <- read.csv(file = paste0(ObsMer_file,"/HL_OBSMER_VIII_VII_2005_2019_230320.csv"),sep = ";",header = T)
  SL_obsmer <- read.csv(file = paste0(ObsMer_file,"/SL_OBSMER_VIII_VII_2005_2019_230320.csv"),sep = ";",header = T)
  TR_obsmer <- read.csv(file = paste0(ObsMer_file,"/TR_OBSMER_VIII_VII_2005_2019_230320.csv"),sep = ";",header = T)
  
  ## Load scientific data
  if(species_to_plot == "Solea_solea") load("data/tidy/fitting_intermediate/bob_cs/Solea_solea/OTB/with_size/mature_survey_data.RData")
  if(species_to_plot == "Merluccius_merluccius") load("data/tidy/fitting_intermediate/bob_cs/Merluccius_merluccius/OTB/with_size/mature_survey_data.RData")
  
  
  
  setwd(folder_project)
  
  ## 'VMS x logbooks'
  #------------------
  source("Scripts/source/shape_VMS_logbooks.R")
  
  ## ObsMer data
  #-------------
  source("Scripts/source/shape_ObsMer.R")
  
  ## Scientific data
  #-----------------
  source("Scripts/source/shape_scientific.R")
  
}

## Simulate data points
#----------------------
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
quadratic_cov = T

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
  cov_x_2 <- cov_x
}


# if(sampling == "simulate"){
#   simu.data <- list(S_x=S_x,
#                     index_i=index_i,
#                     y_sci_i=y_sci_i,
#                     y_i=y_i,
#                     boats_i=boats_i,
#                     cov_x=cov_x,
#                     index_sci_i=index_sci_i,
#                     intercept=intercept,
#                     beta1=beta1,
#                     SD_delta=SD_delta,
#                     range_delta=range_delta,
#                     q1_com=q1,
#                     Sigma_com=SD_obs,
#                     q1_sci=q1_sci,
#                     Sigma_sci=exp(logSigma_sci),
#                     k_com=1)
# 
#   save(data=simu.data,file="results/real_simu/simu.data.RData")
# }

for(aggreg_obs in c(F,T)){
  
  x11(width = 20,height = 15)
  par(mfrow = c(3,4))
  
  for(Estimation_model_i in c(1:3)){
    
    Params_step.est <- NULL
    Map_step.est <- NULL
    
    ## Fit model
    source("Scripts/function/fit_IM.R")
    TmbFile = "Scripts/"
    TMB::compile(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.cpp"),"-O1 -g",DLLFLAGS="")
    dyn.load( dynlib(paste0("Scripts/inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )
    fit_IM_res <- fit_IM(Estimation_model_i = 1, #  = 1
                         Samp_process = 0,
                         EM = "fix_b",
                         TmbFile = "Scripts/",
                         ignore.uncertainty = F,
                         c_com_x = c_com_x,
                         y_com_i = y_i2,
                         index_com_i = index_i,
                         y_sci_i = y_sci_i,
                         index_sci_i = index_sci_i,
                         aggreg_obs=T,
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
                         ObsM=T,
                         y_ObsM_i=y_ObsM_i,
                         index_ObsM_i=index_ObsM_i,
                         sampling = sampling,
                         landings = T,
                         quadratic_cov = quadratic_cov)
    
    Report <- fit_IM_res$Report
    opt <- fit_IM_res$Opt
    Obj <- fit_IM_res$Obj
    SD <- fit_IM_res$SD
    
    if(Estimation_model_i == 1 & step.est == T){
      
      Obj_step.est <- Obj
      init_val <- Obj_step.est$env$last.par.best
      Params <- fit_IM_res$Params
      Map <- fit_IM_res$Map
      
      ## Covariates
      Params$beta_j[which(!is.na(Map$beta_j))] <- init_val[which(names(init_val) == "beta_j")]
      Map$beta_j <- as.character(Map$beta_j)
      
      Map[["beta_j"]] <- as.character(Map[["beta_j"]])
      # Map$beta_j[1] <- 1
      Map[["beta_j"]][1] <- factor(NA)
      Map[["beta_j"]] <- factor(Map[["beta_j"]])
      
      
      ## Observation parameters
      Params$logSigma_com <- init_val[which(names(init_val) == "logSigma_com")]
      Params$logSigma_sci <- init_val[which(names(init_val) == "logSigma_sci")]
      Params$q1_sci <- init_val[which(names(init_val) == "q1_sci")]
      Params$q1_com <- init_val[which(names(init_val) == "q1_com")]
      
      if("k_com" %in% names(init_val)) Params$k_com <- init_val[which(names(init_val) == "k_com")]
      if("k_sci" %in% names(init_val)) Params$k_sci <- init_val[which(names(init_val) == "k_sci")]
      Map[["k_com"]] <- factor(NA)
      # Map[["k_com"]] <- NULL
      
      
      ## Spatial random effect structure
      Params$logtau <- init_val[which(names(init_val) == "logtau")]
      Params$logkappa <- init_val[which(names(init_val) == "logkappa")]
      Params$deltainput_x <- init_val[which(names(init_val) == "deltainput_x")]
      
      Map[["logkappa"]] <- factor(NA)
      # Map[["logkappa"]] <- NULL
      Map[["logtau"]] <- factor(NA)
      # Map[["logtau"]] <- NULL
      
      if("par_b" %in% names(init_val)) Params$par_b <- init_val[which(names(init_val) == "par_b")]
      
      Params_step.est <- Params
      Map_step.est <- Map
      
    }
    
    
    # if(Estimation_model_i == 1 & aggreg_obs == F & sampling == "simulate") save(file="results/real_simu/no.realloc_int_df.RData",data=fit_IM_res)
    # if(Estimation_model_i == 2 & sampling == "simulate") save(file="results/real_simu/sci_df.RData",data=fit_IM_res)
    # if(Estimation_model_i == 3 & aggreg_obs == F & sampling == "simulate") save(file="results/real_simu/no.realloc_com_df.RData",data=fit_IM_res)
    # if(Estimation_model_i == 1 & aggreg_obs == T & sampling == "simulate") save(file="results/real_simu/realloc_int_df.RData",data=fit_IM_res)
    # if(Estimation_model_i == 3 & aggreg_obs == T & sampling == "simulate") save(file="results/real_simu/realloc_com_df.RData",data=fit_IM_res)

    if(Estimation_model_i == 1 & aggreg_obs == F & sampling == "from_tacsateflalo") save(file="results/case_study/no.realloc_int_df.RData",data=fit_IM_res)
    if(Estimation_model_i == 1 & aggreg_obs == T & sampling == "from_tacsateflalo") save(file="results/case_study/realloc_int_df.RData",data=fit_IM_res)
    if(Estimation_model_i == 2 & sampling == "from_tacsateflalo") save(file="results/case_study/sci_df.RData",data=fit_IM_res)
    
    if(sampling == "simulate"){
      
      ## Plot simulations outputs
      plot(x=log(S_x),
           y=log(Report$S_x),
           main="Estimated S_x vs. simulated S_x",
           xlab="Simulated",
           ylab="Estimated")
      mtext(text = paste0("Est_mod: ",Data_source[Estimation_model_i],
                          " Realloc: ",aggreg_obs,
                          " Converge: ",ifelse(opt$convergence==0,T,F)),
            cex=0.75)
      
      for(i in 1:2){
        if(i==1){
          S_plot <- S_x
          main_title <- "Simulated S_x"
          breaks_S <- quantile(S_plot,probs = seq(0,1,length.out=101))
          breaks_S <- unique(breaks_S)
          S_plot.ref <- S_x
          S_interval <- cut(S_plot,breaks = 100,include.lowest = T)
          pal_col <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(100))
          # pal_col <- rev(heat.colors(100))
        }else if(i==2){
          S_plot <- Report$S_x
          main_title <- "Estimated S_x"
          breaks_S <- quantile(S_plot,probs = seq(0,1,length.out=101))
          breaks_S <- unique(breaks_S)
          S_plot.ref <- S_x
          S_interval <- cut(S_plot,breaks = 100,include.lowest = T)
          pal_col <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(100))
        }
        S_interval <- cut(S_plot,breaks = breaks_S, include.lowest = T)
        S_xint <- data.frame(S_interval = S_interval, cell=1:n_cells)
        col_S_xint <- data.frame(S_interval=levels(S_interval),col=pal_col)
        color_df <- merge(x = S_xint, y = col_S_xint, by = "S_interval")
        color_df <- color_df[order(color_df$cell),]
        plot(x=loc_x$x,y=loc_x$y,
             col=as.character(color_df$col),
             pch = 15,cex=0.5,
             asp = 1,
             main = main_title,xlab = "x",ylab="y")
        pal_col2 <- pal_col[c(1,seq(10,100,by=10))]
        lgd_ = rep(NA, length(pal_col2))
        lgd_[c(1,5,10)] = (c(0,round(median(S_plot.ref),digits = 3),round(max(S_plot.ref),digits = 3)))
        legend(x = 17.5, y = 25,
               legend = lgd_,
               fill = pal_col2,
               border = NA,
               y.intersp = 0.5,
               cex = 1, text.font = 2,
               bg="white")
      }
      
    }else if(sampling == "from_tacsateflalo" & fit_IM_res$Converge == 0){
      
      S_plot <- Report$S_x
      breaks_S <- quantile(S_plot,probs = seq(0,1,length.out=101))
      breaks_S <- unique(breaks_S)
      S_interval <- cut(S_plot,breaks = 100,include.lowest = T)
      pal_col <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(breaks_S)-1))
      main_title <- paste0("Est model: ",Estimation_model_i)
      
      S_interval <- cut(S_plot,breaks = breaks_S, include.lowest = T)
      S_xint <- data.frame(S_interval = S_interval, cell=1:n_cells)
      col_S_xint <- data.frame(S_interval=levels(S_interval),col=pal_col)
      color_df <- merge(x = S_xint, y = col_S_xint, by = "S_interval")
      color_df <- color_df[order(color_df$cell),]
      plot(x=loc_x$x,y=loc_x$y,
           col=as.character(color_df$col),
           pch = 15,cex=0.5,
           asp = 1,
           main = main_title,xlab = "x",ylab="y")
      pal_col2 <- pal_col[c(1,seq(10,100,by=10))]
      lgd_ = rep(NA, length(pal_col2))
      lgd_[c(1,5,10)] = (c(0,round(median(S_plot),digits = 3),round(max(S_plot),digits = 3)))
      legend(x = 17.5, y = 25,
             legend = lgd_,
             fill = pal_col2,
             border = NA,
             y.intersp = 0.5,
             cex = 1, text.font = 2,
             bg="white")
    }
    
    
    if(fit_IM_res$Converge == 0){
      par_est <- fit_IM_res$SD$value[which(str_detect(names(fit_IM_res$SD$value),"beta_j|MargSD|Range|q1_com|Sigma_com|q1_sci|Sigma_sci|k_com"))]
      sd_est <- fit_IM_res$SD$sd[which(str_detect(names(fit_IM_res$SD$value),"beta_j|MargSD|Range|q1_com|Sigma_com|q1_sci|Sigma_sci|k_com"))]
      x_axis <- 1:length(par_est)
    }else{
      par_est <- rep(0,length(par_est))
      sd_est <- rep(0,length(par_est))
    }
    
    plot(x=x_axis,
         y=par_est,
         ylim=range(c(par_est - 2, par_est + 2)),
         xaxt = "n",
         pch=19, xlab="", ylab="Parameters values")
    
    # hack: we draw arrows but with very special "arrowheads"
    arrows(x_axis, par_est - 1.96 * sd_est,
           x_axis, par_est + 1.96 * sd_est,
           length=0.05, angle=90, code=3)
    
    names(par_est)[which(str_detect(names(par_est),"beta_j"))[1]] <- "intercept"
    
    axis(1,
         at=1:(length(par_est)),
         labels=names(par_est),
         tick=T,las=2,cex.axis=0.75)
    
    abline(h=0, col="skyblue", lwd=2, lty=2)
    
    if(sampling == "simulate" & fit_IM_res$Converge == 0){
      
      ref_capt <- "k_com"
      
      if(Estimation_model_i == 1){
        par_name <- c("intercept","beta_j","MargSD","Range","q1_com","Sigma_com","q1_sci","Sigma_sci",ref_capt)
        par_value <- c(intercept,beta1,SD_delta,range_delta,q1,SD_obs,q1_sci,exp(logSigma_sci),1)
      }
      
      if(Estimation_model_i == 2){
        par_name <- c("intercept","beta_j","MargSD","Range","q1_sci","Sigma_sci")
        par_value <- c(intercept,beta1,SD_delta,range_delta,q1_sci,exp(logSigma_sci))
      }
      
      
      if(Estimation_model_i == 3){
        par_name <- c("intercept","beta_j","MargSD","Range","q1_com","Sigma_com")
        par_value <- c(intercept,beta1,SD_delta,range_delta,q1,SD_obs)
      }
      
      names(par_value) <- par_name
      par_value <- par_value[names(par_est)]    
      
      points(x=x_axis,y=par_value,col="red",pch=19)
      
    }
    
    ## Concistency check
    if(Estimation_model_i == 1){
      
      obj_int <- Obj
      SD_int <- SD
      obj_int$fn()
      pl_int <- as.list(SD_int,"Est")
      
    }
    
    if(Estimation_model_i == 2){
      
      obj_sci <- Obj
      SD_sci <- SD
      obj_sci$fn()
      pl_sci <- as.list(SD_sci,"Est")
      par_sci <- SD_sci$par.fixed
      f_sci <- as.numeric(obj_sci$fn(par_sci))
      
      ## Fixed
      par_int <- SD_int$par.fixed # extractBaseParameters(obj_sci, pl_sci, pl_int)
      par_int <- par_int[which(names(par_int) %in% names(par_sci))]
      f_int <- as.numeric(obj_sci$fn(par_int))
      def_free <- sum(par_sci != 0)
      fixed <- 1 - pchisq( 2 * (f_int-f_sci), df=def_free )
      
      ## Random
      # par_int.all <- extractBaseParameters(obj_sci, pl_sci, pl_int, all=TRUE)
      par_sci.all <- obj_sci$env$last.par.best
      par_int.all <- obj_int$env$last.par.best
      par_int.all <- par_int.all[which(names(par_int.all) %in% names(par_sci.all))]
      f_sci.all <- obj_sci$env$f(par_sci.all) #Best evaluated parameters
      f_int.all <- obj_sci$env$f(par_int.all)
      def_free.all <- def_free + length(obj_sci$env$random)
      random <- 1 - pchisq( 2 * (f_int.all-f_sci.all), df=def_free.all )
      
    }
    
    if(consistency_check == T & Estimation_model_i == 3){
      
      x11(width = 20,height = 10)
      par(mfrow = c(1,2))
      plot(x=log(SD_sci$value[which(names(SD_sci$value)=="S_x")]),
           y=log(SD_int$value[which(names(SD_int$value)=="S_x")]),
           xlab="log(scientific predictions)",ylab="log(integrated predictions)",
           xlim = log(range(c(SD_sci$value[which(names(SD_sci$value)=="S_x")]))),
           ylim = log(range(c(SD_int$value[which(names(SD_int$value)=="S_x")]))))
      abline(a = 0,b = 1)
      
      mtext(text = paste0("Consistency check p-value \n",
                          "fixed = ",format(signif(fixed,digits = 3),scientific = T),
                          " | random = ",format(signif(random,digits = 3),scientific = T)),
            side = 3,cex = 0.75)
      
      plot(x=SD_sci$par.fixed,
           y=SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))],
           xlab="Scientific",ylab="Integrated",
           xlim = c(min(c(SD_sci$par.fixed,SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))]))-1,
                    max(c(SD_sci$par.fixed,SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))])))+1,
           ylim = c(min(c(SD_sci$par.fixed,SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))]))-1,
                    max(c(SD_sci$par.fixed,SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))]))+1)
      )
      
      abline(a = 0,b = 1)
      text(x=SD_sci$par.fixed,
           y=SD_int$par.fixed[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))],
           labels=names(SD_int$par.fixed)[which(names(SD_int$par.fixed) %in% names(SD_sci$par.fixed))],
           pos=1)
      
    }
  }
}

# test <- cbind(loc_x[,c("x","y")],
#               # S_est = fit_IM_res$Report$S_p,
#               S_sim = S_x)
#
# simu_plot <- ggplot(test)+
#   geom_point(aes(x=x,y=y,col=S_sim))+
#   scale_color_distiller(palette = "Spectral")+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Simulation")
#
# est_plot <- ggplot(test)+
#   geom_point(aes(x=x,y=y,col=S_est))+
#   scale_color_distiller(palette = "Spectral")+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Estimation")
#
# x11()
# plot_grid(effort_plot,seq_plot,simu_plot,est_plot)



# data.points_df <- cbind(loc_x[,c("x","y")],c_com_x)
# effort_plot <- ggplot()+
#   geom_point(data = data.points_df,aes(x=x,y=y,col=c_com_x))+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Sampling effort")+
#   scale_color_distiller(palette = "Spectral")
#
# sample_df <- data.frame(cell=index_i,
#                         seq=boats_i,
#                         catch_no.realloc=y_i,
#                         catch_realloc=y_i2) %>%
#   inner_join(loc_x[,c("cell","x","y")]) %>%
#   filter(seq < 10)
# 
# seq_plot <- ggplot()+
#   geom_point(data = sample_df,aes(x=x,y=y,col=factor(seq)))+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Fishing sequence")

# extractBaseParameters <- function(obj1, pl1, pl2, all=FALSE) {
#   npar1 <- sapply(pl1,length) #Number of survey parameters by component
#   plnew <- Map(head, pl2, npar1) #Extract survey parameters from second fit
#   applyMap <- function(parameters, map) {
#     param.map <- lapply(names(map),
#                         function(nam)
#                         {
#                           TMB:::updateMap(parameters[[nam]], map[[nam]])
#                         })
#     parameters[names(map)] <- param.map
#     parameters
#   }
#   par2 <- unlist(applyMap(plnew, obj1$env$map))
#   if (!all) par2 <- par2[-obj1$env$random]
#   par2
# }



# ## Plot fishing sequence catch values
# seq_plot <- ggplot()+
#   geom_point(data = sample_df,aes(x=x,y=y,col=factor(seq)))+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   xlab("")+ylab("")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust=0.5),
#         legend.position = "none",
#         panel.background = element_rect(fill="skyblue"))
# 
# catch.no.realloc_plot <- ggplot()+
#   geom_point(data = test,aes(x=x,y=y,col=S_sim),alpha=0.2,shape=15)+
#   # geom_point(data = sample_df,aes(x=x,y=y,fill=catch_no.realloc),col="black",shape=21)+
#   scale_color_distiller(palette="Spectral")+
#   scale_fill_distiller(palette="Spectral")+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Exact catches")+xlab("")+ylab("")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust=0.5),legend.position = "none")
# 
# catch.realloc_plot <- ggplot()+
#   geom_point(data = test,aes(x=x,y=y,col=S_sim),alpha=0.2,shape=15)+
#   geom_point(data = sample_df,aes(x=x,y=y,fill=catch_realloc),col="black",shape=21)+
#   scale_color_distiller(palette="Spectral")+
#   scale_fill_distiller(palette="Spectral")+
#   geom_sf(data=st_geometry(ICES_rect_2),alpha=0)+
#   geom_sf(data=mapBase)+
#   coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
#   ggtitle("Reallocated catches")+xlab("")+ylab("")+
#   theme_classic()+
#   theme(plot.title = element_text(hjust=0.5),legend.position = "none")
# 
# plot_grid(NULL,seq_plot,catch.no.realloc_plot,catch.realloc_plot)
