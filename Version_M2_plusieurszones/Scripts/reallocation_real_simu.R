##################################################
## Reallociation - Build simulations on real data
##################################################
## B. Alglave

memory.size(1e+9)

library(cowplot)
library(dplyr)
library(ggplot2)
library(raster)
library(rnaturalearth)
library(sf)
library(spatstat)
library(stringr)
library(TMB)

folder_data <- "C:/R_projects/phd_zfh_baptiste_alglave"
folder_project <- "C:/R_projects/projet_M2_sole/Version_M2_plusieurszones"

setwd(folder_data)

##--------------
## Load raw data
##--------------
## Map info
grid_projection <- "+proj=longlat +datum=WGS84"
mapBase <- ne_countries(scale = "medium", returnclass = "sf")
grid_xmin <- -6 # Coordinates of study area 
grid_xmax <- 0
grid_ymin <- 42
grid_ymax <- 48
grid_limit <- extent(c(grid_xmin,grid_xmax,grid_ymin,grid_ymax))
Resol_map <- 0.05

## Time interval
year_y <- 2018
month_m <- 11

## data file name
sampling <- "simulate" # "simulate": simulate sampling, "from_tacsateflalo": take data points from tacsateflalo
commercial_data_folder <- 'data/raw/vms_logbooks/bob_cs/MergeTACSAT_EFLALO/'

## If simulate data points
# n_samp_com <- 5000
n_seq_com <- 500
n_ping_per_seq <- 10

## Load ICES rectangle shapefile
ICES_rect <- st_read("data/raw/shapefile/ices_rect_stats/ICES_Statistical_Rectangles_Eco.shp")
st_crs(ICES_rect) <- grid_projection
ICES_rect_2 <- ICES_rect %>%
  filter(Ecoregion == "Bay of Biscay and the Iberian Coast")


## Load domain
strata_file <- "data/raw/shapefile/QGIS_ORHAGO/shape/" 

zone_centre_large <- sf::st_read(paste0(strata_file,"zone centre au large.shp"))
zone_centre_cote <-sf::st_read(paste0(strata_file,"zone centre cote.shp"))
zone_nord <- sf::st_read(paste0(strata_file,"zone nord.shp"))
zone_sud <- sf::st_read(paste0(strata_file,"zone sud.shp"))

setwd(folder_project)

## Orhago polygon
survey_layer <- st_union(zone_centre_large,zone_centre_cote)
survey_layer <- st_union(survey_layer,zone_nord)
survey_layer <- st_union(survey_layer,zone_sud)
survey_layer_2 <- nngeo::st_remove_holes(survey_layer)
survey_layer_3 <- sf::as_Spatial(survey_layer_2)
crs(survey_layer_3) <- grid_projection
study_domain_2 <- survey_layer_3

## Matrices to plot maps
grid <- raster(grid_limit)
res(grid) <- Resol_map
crs(grid) <- grid_projection
gridpolygon <- rasterToPolygons(grid)
gridpolygon$layer <- c(1:length(gridpolygon$layer))

raster_to_point <- rasterToPoints(grid)
datapoint <- SpatialPointsDataFrame(coords=raster_to_point,
                                    data=data.frame(layer = 1:nrow(raster_to_point)),
                                    proj4string=crs(grid_projection))

crs(datapoint) <- crs(study_domain_2) <- grid_projection
datapoint_2 <- over(datapoint,study_domain_2)
datapoint_2 <- cbind(datapoint_2,layer = datapoint@data$layer)
datapoint_2 <- datapoint_2[which(!is.na(datapoint_2[,1])),]
## Add coordinates to gridpoint
tmp <- cbind(datapoint@data, coordinates(datapoint))
colnames(tmp) <- c(names(datapoint@data), "x", "y")
datapoint_3 <- inner_join(datapoint_2,tmp) %>%
  dplyr::select(layer,x,y)

## Add strata to datapoint
survey_layer_zcl <- as_Spatial(zone_centre_large)
survey_layer_zcc <- as_Spatial(zone_centre_cote)
survey_layer_zn <- as_Spatial(zone_nord)
survey_layer_zs <- as_Spatial(zone_sud)

survey_layer_zcl@data <- cbind(survey_layer_zcl@data,strata = "zcl")
survey_layer_zcc@data <- cbind(survey_layer_zcc@data,strata = "zcc")
survey_layer_zn@data <- cbind(survey_layer_zn@data,strata = "zn")
survey_layer_zs@data <- cbind(survey_layer_zs@data,strata = "zs")

# cross sp datapoint and strata
datapoint_zcl <- over(datapoint,survey_layer_zcl)
datapoint_zcl_2 <- cbind(datapoint_zcl,layer = datapoint@data$layer)
datapoint_zcl_3 <- filter(datapoint_zcl_2,!is.na(id))

datapoint_zcc <- over(datapoint,survey_layer_zcc)
datapoint_zcc_2 <- cbind(datapoint_zcc,layer = datapoint@data$layer)
datapoint_zcc_3 <- filter(datapoint_zcc_2,!is.na(id))

datapoint_zn <- over(datapoint,survey_layer_zn)
datapoint_zn_2 <- cbind(datapoint_zn,layer = datapoint@data$layer)
datapoint_zn_3 <- filter(datapoint_zn_2,!is.na(id))

datapoint_zs <- over(datapoint,survey_layer_zs)
datapoint_zs_2 <- cbind(datapoint_zs,layer = datapoint@data$layer)
datapoint_zs_3 <- filter(datapoint_zs_2,!is.na(id))

colnames(datapoint_zcl_3) <- c("id","surface","strata","layer")
colnames(datapoint_zcc_3) <- c("id","surface","strata","layer")
colnames(datapoint_zn_3) <- c("id","surface","strata","layer")
colnames(datapoint_zs_3) <- c("id","surface","strata","layer")

layer_strata <- rbind(datapoint_zcl_3,datapoint_zcc_3)
layer_strata <- rbind(layer_strata,datapoint_zn_3)
layer_strata <- rbind(layer_strata,datapoint_zs_3)

# Add coordinates to datapoint
tmp <- cbind(datapoint@data, coordinates(datapoint))
colnames(tmp) <- c(names(datapoint@data), "long", "lati")
datapoint_4 <- inner_join(datapoint_3 , tmp)

# Add strata
datapoint_5 <- inner_join(datapoint_4 , layer_strata[,c("strata","layer")],by="layer")

loc_x <- datapoint_5 %>% mutate(cell = 1:nrow(datapoint_5))
loc_x_sf <- st_as_sf(loc_x, coords = c("x", "y"))
st_crs(loc_x_sf) <- grid_projection

## Cross with statistical rectangle
loc_x_sf <- loc_x_sf[st_intersects(loc_x_sf,ICES_rect_2) %>% lengths > 0,]
loc_x_sf <- st_join(loc_x_sf,ICES_rect_2)
loc_x <- loc_x_sf %>%
  data.frame %>%
  dplyr::select(layer,cell,ICESNAME,strata,-geometry) %>%
  mutate(x=loc_x$x,
         y=loc_x$y)

##------------------------------------------------------------------------------------
##-------------------- Simulate latent field on the domain ---------------------------
##------------------------------------------------------------------------------------
source("Scripts/function/sim_GF_Matern.R")
intercept <- 2
beta1 <- runif(1,1,2)
nu <- 1
range_cov <- 1
SD_cov <- 0.5
range_delta <- 0.5
SD_delta <- 1

cov_x <- sim_GF_Matern(as.matrix(loc_x[,c("x","y")]), nu, range_cov, SD_cov^2)[[1]]$y
delta_x <- sim_GF_Matern(as.matrix(loc_x[,c("x","y")]), nu, range_delta, SD_delta^2)[[1]]$y

S_x <- exp(intercept + cov_x %*% beta1 + delta_x)

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

## From tacsatEflalo
#-------------------
if(sampling == "from_tacsateflalo"){
  
  ## Load tacsatEflalo
  setwd(folder_data)
  filter_with_Eflalo <- F
  source("r/real_data/load_data/extract_commercial_data.R")
  setwd(folder_project)
  
  ## Real dataset
  #--------------
  tacsatEflalo_2 <- tacsatEflalo %>%
    filter(SI_LATI < 48 &  SI_LONG > -6 & SI_STATE == 1) %>%
    filter(str_detect(LE_MET_level6,"OTB|OTT|PTM|GTR|OTM|PS"))
  
  # convert the points to an sf object
  points_sf <- st_as_sf(tacsatEflalo_2, coords = c("SI_LONG", "SI_LATI"))
  
  # set CRS for the points to be the same as shapefile
  st_crs(points_sf) <- st_crs(grid_projection)
  
  # cross with statistical rectangles
  points_sf <- points_sf[st_intersects(points_sf,ICES_rect_2) %>% lengths > 0,]
  points_sf <- st_join(points_sf,ICES_rect_2)
  
  points_sf <- points_sf %>%
    mutate(seq_id = paste0(VE_REF,"_",SI_DATE,"_",LE_GEAR,"_",ICESNAME))
  
  # ## Number of pings per fishing sequence
  # points_sf  <- points_sf %>% filter(SI_STATE == 1)
  # points_df <- points_df %>% filter(SI_STATE == 1)
  # 
  # nb_ping_seq <- as.data.frame(table(points_sf$seq_id)) %>%
  #   group_by(Freq) %>%
  #   dplyr::summarise(n = n()) %>%
  #   dplyr::rename(ping_per_seq = Freq) %>%
  #   filter(ping_per_seq <= 24) %>%
  #   mutate(freq = n / sum(n))
  # 
  # ggplot(nb_ping_seq,aes(x=ping_per_seq,y=freq))+
  #   geom_point()+
  #   geom_line()+
  #   ylim(0,NA)
  # save(data=nb_ping_seq,file="results/exploratory_analysis/reallocation/nb_ping_seq.RData")
  
  nb_ping_gear <- as.data.frame(table(points_sf$VE_REF,points_sf$LE_GEAR)) %>%
    filter(Freq != 0)
  
  vessel_vec <- c("103","157914","14868","151271") # OTB_DEF, GTR, OTT, PTM
  
  points_sf_2 <- points_sf %>%
    # filter(VE_REF %in% vessel_vec) %>%
    mutate(year = as.numeric(str_sub(SI_DATE,7,10)),
           month = as.numeric(str_sub(SI_DATE,4,5))) %>%
    filter(str_detect(LE_MET_level6,"OTB_DEF_>=70_0")) %>%
    filter(month == 11 & year == 2018)
  
  ## cross datapoints with grid
  #----------------------------
  gridpolygon_sf <- st_as_sf(gridpolygon)
  points_sf_3 <- points_sf_2[st_intersects(points_sf_2,gridpolygon_sf) %>% lengths > 0,]
  points_sf_3 <- st_join(points_sf_3,gridpolygon_sf)
  points_sf_4 <- points_sf_3 %>%
    dplyr::select(SI_DATE,LE_MET_level6,LE_KG_Solea_solea,ICESNAME,seq_id,layer)
  points_df <- points_sf_4 %>% data.frame
  points_df <- inner_join(points_df,loc_x)
  
  ## Create vector of the number of fishing locations per cell
  #-----------------------------------------------------------
  c_com_x <- points_df %>%
    dplyr::count(cell) %>%
    full_join(loc_x) %>%
    arrange(cell) %>%
    data.frame() %>% 
    mutate(n = ifelse(is.na(n), 0,n)) %>%
    dplyr::select(n) %>%
    c() %>% unlist() %>%
    as.matrix()

  index_i <- points_df$cell
  boats_i <- as.numeric(factor(points_df$seq_id))
  
}

## Simulate data points
#----------------------
if(sampling == "simulate"){

  area_df <- loc_x %>%
    group_by(ICESNAME) %>%
    dplyr::summarise(area=n()) %>%
    mutate(perc_area = area/(sum(area))) %>%
    mutate(ICESNAME = factor(ICESNAME)) %>%
    arrange(ICESNAME)
  
  ## Sequential sampling (not that there is no preferential sampling)
  # First select a statistical rectangle
  sample_full_rect <- F
  zonesize <- 0.2
  n_zone <- 1 # c(1,3,5)
  samp_df <- do.call(rbind,lapply(1:n_seq_com, function(j){
    
    ## Select a statistical rectangle
    rect_stat <- sample(x = area_df$ICESNAME,
                        size = 1,
                        replace=T,
                        prob =  area_df$perc_area)
    
    loc_x_rect <- loc_x %>%
      filter(ICESNAME == rect_stat)
    
    if(sample_full_rect == T){
      
      ## Sample within the statistical rectangle
      samp_cell <- sample(x = loc_x_rect$cell,
                          size = n_ping_per_seq,
                          replace=T)
      
    }else{
      
      samp_cell <- do.call(c,lapply(1:n_zone, function(j){
        ## Sample within the statistical rectangle
        center <- sample(x = loc_x_rect$cell,
                         size = 1,
                         replace=T)
        
        loc_x_zone <- loc_x_rect[which(loc_x_rect$x < loc_x_rect$x[which(loc_x_rect$cell == center)] + zonesize &
                                         loc_x_rect$x > loc_x_rect$x[which(loc_x_rect$cell == center)] - zonesize &
                                         loc_x_rect$y < loc_x_rect$y[which(loc_x_rect$cell == center)] + zonesize &
                                         loc_x_rect$y > loc_x_rect$y[which(loc_x_rect$cell == center)] - zonesize),]
        
        ## Sample within the statistical rectangle
        samp_cell <- sample(x = loc_x_zone$cell,
                            size = n_ping_per_seq/n_zone,
                            replace=T)
      }))
      
      
    }
    
    seq_id <- rep(j,n_ping_per_seq)
    df <- data.frame(cell=samp_cell,seq_id=seq_id,rect_stat)
    return(df)
    
  }))
  
  index_i <- samp_df$cell
  boats_i <- samp_df$seq_id
  
  c_com_x <- samp_df %>% 
    dplyr::count(cell) %>%
    full_join(loc_x) %>%
    arrange(cell) %>%
    data.frame() %>% 
    mutate(n = ifelse(is.na(n), 0,n)) %>%
    dplyr::select(n) %>%
    c() %>% unlist() %>%
    as.matrix()
  
}


## Simulate catch observations (zero-inflated lognormal distribution)
q1 <- -1
SD_obs <- 0.5
y_i <- do.call(c,lapply(1:length(index_i), function(j){
  exp_catch <- S_x[index_i[j]]
  prob_encount <- 1-exp(- exp(q1) * exp_catch)
  abs_i <- rbinom(1,1,prob_encount)
  if(abs_i>0){
    y_i <- exp(rnorm(1,mean = log(exp_catch/prob_encount) - SD_obs^2/2, sd = SD_obs))
  }else{
    y_i <- 0
  }
}))

## Reallocation
k_sim <- 1
if(k_sim==1){ ## Reallocate data
  mean_y_i <- aggregate(x = y_i, 
                        by = list(unique.values = boats_i),
                        FUN = mean)[,2]
  
  sum_y_i <- aggregate(x = y_i, 
                       by = list(unique.values = boats_i),
                       FUN = sum)[,2]
  
  y_i2 <- rep(NA,length(y_i))
  y_i2 <- do.call(c,lapply(1:length(y_i), function(k){
    y_i2[k] <- mean_y_i[boats_i[k]]
  }))
}else{ ## Do not reallocate data
  y_i2 <- y_i
}

# # check:
# x11()
# par(mfrow = c(3,1))
# plot((S_x[index_i]),log(y_i))
# plot((y_i2),(y_i))
# plot((S_x[index_i]),log(y_i2))

#-----------------------------------------------------------------------------------
##-------------------------------- Scientific data ---------------------------------
#-----------------------------------------------------------------------------------

index_sci_i <- c()
n_cells <- nrow(loc_x)
n_samp_sci <- 100
nb_hauls_strata <- loc_x %>% 
  mutate(value=1) %>%
  group_by(strata) %>%
  dplyr::summarise(value = sum(value)) %>%
    mutate(hauls = round(value/n_cells*n_samp_sci))

# positions of hauls (straitified random sampling)
index_sci_i <- do.call(c,lapply(1:nrow(nb_hauls_strata),function(j){
  index_sci_i <- c(index_sci_i,sample(loc_x$cell[which(loc_x$strata == nb_hauls_strata$strata[j])], size=nb_hauls_strata$hauls[which(nb_hauls_strata$strata == nb_hauls_strata$strata[j])],replace=FALSE))
}))

# c_sci_x = ifelse(1:n_cells %in% index_sci_i, 1, 0) # shots for scientific data
# # check sampling of sci data
# verif_samp_sci <- cbind(loc_x,c_sci_x)
# simu_plot_1 <- ggplot()+
#   geom_point(data = verif_samp_sci,aes(x=x,y=y,col=strata)) +
#   geom_point(data = verif_samp_sci[which(verif_samp_sci$c_sci_x > 0),],aes(x=x,y=y)) +
#   theme_bw()

# scientific catch data (zero-inflated lognormal)
q1_sci <- 0
logSigma_sci <- log(0.5)
y_sci_i <- do.call(c,lapply(1:length(index_sci_i), function(j){
  exp_catch <- S_x[index_sci_i[j]]
  prob_encount <- 1-exp(- exp(q1_sci) * exp_catch)
  
  abs_sci_i <- rbinom(1,1,prob_encount)
  if(abs_sci_i>0){
    y_sci_i <-  exp(rnorm(1,mean = log(exp_catch)-exp(logSigma_sci)^2/2, sd = exp(logSigma_sci)))
  }else{
    y_sci_i <- 0 
  }
  return(y_sci_i)
  
}))



#-----------------------------------------------------------------------------------
##----------------------------------- Fit Model ------------------------------------ 
#-----------------------------------------------------------------------------------

logSigma_com <- log(SD_obs)+0.1 # to initialize 
b <- 0
SP_est <- 0 # Do not account for sampling process
b_est <- 0 # Do not account for preferential sampling
eta_est <- 0 # Do not account for sampling processes other than preferential sampling
lf_param <- "cov"
lf_param_num <- ifelse(lf_param=="cov",0,1)
obs_mod <- 2

## Mesh
library(INLA)
mesh <- inla.mesh.create( loc_x[,c('x','y')] )
spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]  # define sparse matrices for parametrisation of precision matrix
n_cells <- nrow(loc_x)



for(aggreg_obs in c(F,T)){
  
  x11(width = 20,height = 10)
  par(mfrow = c(2,4))
  
  for(Estimation_model_i in c(1,3)){
  
    
    ## Fit model
    source("Scripts/function/fit_IM.R")
    TmbFile = "Scripts/"
    TMB::compile(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.cpp"),"-O1 -g",DLLFLAGS="")
    dyn.load( dynlib(paste0("Scripts/inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )
    fit_IM_res <- fit_IM(Estimation_model_i = Estimation_model_i,
                         Samp_process = 0,
                         EM = "fix_b",
                         TmbFile = "Scripts/",
                         ignore.uncertainty = F,
                         c_com_x = c_com_x,
                         y_com_i = y_i2,
                         index_com_i = index_i,
                         y_sci_i = y_sci_i,
                         index_sci_i = index_sci_i,
                         aggreg_obs = aggreg_obs,
                         boats_number = boats_i,
                         Cov_x = as.matrix(cov_x),
                         lf_param = "RE",# lf_param,
                         spde=spde,
                         mesh=mesh)
    
    Report <- fit_IM_res$Report
    opt <- fit_IM_res$Opt
    
    ## Plot simulations outputs
    plot(x=log(S_x),
         y=log(Report$S_x),
         main="Estimated S_x vs. simulated S_x",
         xlab="Simulated",
         ylab="Estimated")
    mtext(text = paste0("Est_mod: ",ifelse(Estimation_model_i==1,"Integrated","Commercial"),
                        " Realloc: ",aggreg_obs,
                        " Converge: ",ifelse(opt$convergence==0,T,F)),
          cex=0.75)
    
    for(i in 1:2){
      if(i==1){
        S_plot <- S_x
        main_title <- "Simulated S_x"
        breaks_S <- quantile(S_plot,probs = seq(0,1,length.out=101))
        S_plot.ref <- S_x
        S_interval <- cut(S_plot,breaks = 100,include.lowest = T)
        pal_col <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(100))
        # pal_col <- rev(heat.colors(100))
      }else if(i==2){
        S_plot <- Report$S_x
        main_title <- "Estimated S_x"
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
    
    par_est <- fit_IM_res$SD$value[which(str_detect(names(fit_IM_res$SD$value),"beta_j|MargSD|Range|q1_com|Sigma_com|q1_sci|Sigma_sci"))]
    sd_est <- fit_IM_res$SD$sd[which(str_detect(names(fit_IM_res$SD$value),"beta_j|MargSD|Range|q1_com|Sigma_com|q1_sci|Sigma_sci"))]
    x_axis <- 1:length(par_est)
    
    plot(x=x_axis,
         y=par_est,
         ylim=range(c(par_est - 2, par_est + 2)),
         xaxt = "n",
         pch=19, xlab="", ylab="Parameters values",
    )
    
    # hack: we draw arrows but with very special "arrowheads"
    arrows(x_axis, par_est - 1.96 * sd_est,
           x_axis, par_est + 1.96 * sd_est,
           length=0.05, angle=90, code=3)
    axis(1,
         at=1:(length(par_est)),
         labels=names(par_est),
         tick=T,las=2,cex.axis=0.75)
    
    if(Estimation_model_i == 1){
      par_name <- c("beta_j","beta_j","MargSD","Range","q1_com","Sigma_com","q1_sci","Sigma_sci")
      par_value <- c(intercept,beta1,SD_delta,range_delta,q1,SD_obs,q1_sci,exp(logSigma_sci))
    }     

    if(Estimation_model_i == 3){
      par_name <- c("beta_j","beta_j","MargSD","Range","q1_com","Sigma_com")
      par_value <- c(intercept,beta1,SD_delta,range_delta,q1,SD_obs)
    }
    
    names(par_value) <- par_name
    par_value <- par_value[names(par_est)]    
    
    points(x=x_axis,y=par_value,col="red",pch=19)
    
  }
  
}


test <- cbind(loc_x[,c("x","y")],
              S_est = fit_IM_res$Report$S_x,
              S_sim = S_x)

simu_plot <- ggplot(test)+
  geom_point(aes(x=x,y=y,col=S_sim))+
  scale_color_distiller(palette = "Spectral")+
  geom_sf(data=mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  ggtitle("Simulation")

est_plot <- ggplot(test)+
  geom_point(aes(x=x,y=y,col=S_est))+
  scale_color_distiller(palette = "Spectral")+
  geom_sf(data=mapBase)+
  coord_sf(xlim = c(-6,0), ylim = c(43,48), expand = FALSE)+
  ggtitle("Estimation")

x11()
plot_grid(effort_plot,seq_plot,simu_plot,est_plot)



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
