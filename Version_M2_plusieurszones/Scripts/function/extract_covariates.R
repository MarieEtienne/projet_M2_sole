############################################
## Extract covariates from shapefiles/raster
############################################

#---------
# Sediment
#---------
#' @title extract_substrate_rasterCelticBoB()
#' @param cov_file file where are located covariates
#' @param gridpoint grid points
#' @param gridpoint_cov grid points with covariates
#' @param aggreg_substrate if T, substrate are aggregated
#' 
#' @return gridpoint_cov gridpoint with substrate levels

# cov_file <- "E:/data2/Envir.Data/"
# aggreg_substrate = T
extract_substrate_rasterCelticBoB <- function(cov_file,gridpoint,gridpoint_cov,aggreg_substrate){
  # Extract from raster and project in the grid
  str_name<-paste0(cov_file,'Substrate.Habitat.Atl.Bob.tif')
  imported_raster=raster(str_name)
  substrate <- raster::extract(stack(imported_raster), coordinates(gridpoint))
  substratepoint.sp <- gridpoint
  substratepoint.sp@data <- cbind(substratepoint.sp@data,substrate)
  colnames(substratepoint.sp@data)[which(colnames(substratepoint.sp@data)=="Substrate.Habitat.Atl.Bob")] <- "Substrate.Code"
  # spplot(substratepoint.sp,"Substrate.Code",cex = 0.01)
  
  # numbers --> sediment name
  Substrate.names <- c("Fine mud",
                       "Rock or other hard substrata",
                       "Sandy mud or Muddy sand",
                       "Sand",
                       "Muddy sand",
                       "Seabed",
                       "Mixed sediment",
                       "Sandy mud",
                       "Coarse substrate",
                       "Sediment",
                       "Fine mud or Sandy mud or Muddy sand")
  
  Substrate.Code <- c(1:length(Substrate.names))
  Substrate.Code_df <-  data.frame(Substrate.names = Substrate.names,
                                   Substrate.Code = Substrate.Code)
  
  substratepoint.sp@data <- left_join(substratepoint.sp@data,Substrate.Code_df,by = "Substrate.Code")
  substratepoint.sp@data$Substrate.names <- as.character(substratepoint.sp@data$Substrate.names)
  
  # Aggregate type of sediment
  if(aggreg_substrate == T){
    substratepoint.sp@data <- substratepoint.sp@data %>%
      mutate(aggreg_substrate = if_else(str_detect(Substrate.names,"Sediment") | str_detect(Substrate.names,"sediment") | str_detect(Substrate.names,"Mud") | str_detect(Substrate.names,"mud"),"Mud_sediment",Substrate.names)) %>%
      mutate(aggreg_substrate = if_else(str_detect(Substrate.names,"Rock") ,"Rock",aggreg_substrate)) %>%
      mutate(aggreg_substrate = if_else(str_detect(Substrate.names,"Coarse") | str_detect(Substrate.names,"Sand"),"Sand_Coarse_substrate",aggreg_substrate))
  }else{
    substratepoint.sp@data <- substratepoint.sp@data %>%
      mutate(aggreg_substrate = Substrate.names)
  }
  substratepoint.sp@data$Substrate.names <- as.factor(substratepoint.sp@data$Substrate.names)
  substratepoint.sp@data$aggreg_substrate <- as.factor(substratepoint.sp@data$aggreg_substrate)
  
  
  crs(substratepoint.sp) <- grid_projection
  
  gridpoint_cov <- merge(gridpoint_cov, substratepoint.sp@data[,c("layer","aggreg_substrate")], by = c("layer"),all.x=T)
  
  test <- merge(gridpoint_cov, substratepoint.sp@data[,c("layer","aggreg_substrate")], by = c("layer"),all.x=T)
  
  return(gridpoint_cov)
}

# # Other version for sediment (faster to plot full data) ----------------------------------
# str_name<-'E:/data2/Envir.Data/Substrate.Habitat.Atl.Bob.tif'
# imported_raster=raster(str_name)
# Substrate.sp <- as(imported_raster, 'SpatialPointsDataFrame')
# 
# gridpoint_cov <- merge(gridpoint_cov , substratepoint.sp@data[,c("layer","aggreg_substrate")], by = "layer",all.x=T)
# colnames(Substrate.sp@data)[1] <- "Substrate.Code"
# Substrate.sp@data <- inner_join(Substrate.sp@data,Substrate.Code_df,by = "Substrate.Code")
# Substrate.sp@data$Substrate.names <- as.character(Substrate.sp@data$Substrate.names)
# Substrate.sp@data <- Substrate.sp@data %>%
#   mutate(aggreg_substrate = if_else(str_detect(Substrate.names,"Sediment") | str_detect(Substrate.names,"sediment") | str_detect(Substrate.names,"Mud") | str_detect(Substrate.names,"mud"),"Mud_sediment",Substrate.names)) %>%
#   mutate(aggreg_substrate = if_else(str_detect(Substrate.names,"Coarse") | str_detect(Substrate.names,"Rock") ,"Rock_Coarse_substrate",aggreg_substrate)) %>%
#   mutate(aggreg_substrate = if_else(str_detect(Substrate.names,"Sand"),"Sand",aggreg_substrate))
# 
# Substrate.sp@data$Substrate.names <- as.factor(Substrate.sp@data$Substrate.names)
# Substrate.sp@data$aggreg_substrate <- as.factor(Substrate.sp@data$aggreg_substrate)
# crs(Substrate.sp) <- grid_projection
# 
# substrpoint_1 <- over(gridpoint,Substrate.sp)
# substrpoint_2 <- cbind(gridpoint@data, substrpoint_1)
# substrpoint_3 <- filter(substrpoint_2,!is.na(Substrate.Code))
# substrpoint_4 <- merge(substrpoint_3 , tmp, all.x=T)
# substrpoint_5 <- dplyr::select(substrpoint_4,layer,long,lati,aggreg_substrate)
# 
# gridpoint_cov <- merge(gridpoint_cov, substrpoint_5[,c("layer","aggreg_substrate")], by = "layer",all.x=T)
# spplot(Substrate.sp,"Substrate.names",cex = 0.1)
# 
# # ----------------------------------

#--------
# Bathy
#--------

#' @title extract_bathy_rasterCelticBoB()
#' @param cov_file file where are located covariates
#' @param gridpoint grid points
#' @param gridpoint_cov grid points with covariates
#' 
#' @return gridpoint_cov gridpoint with bathymetry values

extract_bathy_rasterCelticBoB <- function(cov_file,gridpoint,gridpoint_cov){
  str_name<-paste0(cov_file,'Bathy.BayBiscay.grd')
  imported_raster=raster(str_name)
  bathy <- raster::extract(stack(imported_raster), coordinates(gridpoint))
  colnames(bathy) <- "bathy"
  bathypoint.sp <- gridpoint
  bathypoint.sp@data <- cbind(bathypoint.sp@data,bathy)
  # spplot(bathypoint.sp,"bathy",cex = 0.01)
  gridpoint_cov <- merge(gridpoint_cov , bathypoint.sp@data, by = "layer",all.x=T)
  return(gridpoint_cov)
}

#-----------------
## Copernicus data
#-----------------

## Physical data

# year_month <- "2018-11"
# copernicus.phys.var <- "thetao"  # "bottomT", "so"
# name_var <- "SST"

#' @title extract_copernicus.phys_rasterCelticBoB()
#' @param cov_file file where are located covariates
#' @param gridpoint grid points
#' @param gridpoint_cov grid points with covariates
#' @param copernicus.phys.var covariate to select in copernicus database
#' @param name_var name to give to this covariate
#' @param year_month month to select
#' @param grid_projection grid projection
#' 
#' @return gridpoint_cov gridpoint with physical covariates values

extract_copernicus.phys_rasterCelticBoB <- function(cov_file,gridpoint,gridpoint_cov,copernicus.phys.var,name_var,year_month,grid_projection){
  # read variable in nc table
  nc_in <- open.nc(paste0(cov_file,"CopernicusData/dataset-ibi-reanalysis-phys-005-002-monthly_1576082594192_0m.nc"))
  # print.nc(nc_in) # variable
  nc_lon <- var.get.nc(nc_in, "longitude")
  nc_lat <- var.get.nc(nc_in, "latitude")
  nc_var <- var.get.nc(nc_in, copernicus.phys.var,unpack=TRUE)
  nc_t <- var.get.nc(nc_in, "time")
  nc_t <- as.POSIXct(nc_t*60*60, origin="1950-01-01",tz="GMT", format="%Y-%m-%d")
  close.nc(nc_in)
  
  # assuming a regular grid (which by the way appears not entirely true for this data set):
  gt <- GridTopology(cellcentre.offset = c(nc_lon[1], nc_lat[1]),
                     cellsize = c((tail(nc_lon, 1) - nc_lon[1])/(length(nc_lon) - 1), (tail(nc_lat, 1) - nc_lat[1])/(length(nc_lat) - 1)),
                     cells.dim = c(length(nc_lon), length(nc_lat)))
  sg <- SpatialGrid(gt, CRS(grid_projection))
  sgdf <- SpatialGridDataFrame(sg, data.frame(var = as.vector(nc_var[,,which(str_detect(nc_t,year_month))])))
  raster_var <- flip(raster(sgdf),direction = 2)
  # plot(flip(raster(sgdf),direction = 2))
  
  extract_raster_var <- raster::extract(stack(raster_var), coordinates(gridpoint))
  copernicuspoint.sp <- gridpoint
  copernicuspoint.sp@data <- as.data.frame(cbind(copernicuspoint.sp@data[,"layer"],extract_raster_var))
  colnames(copernicuspoint.sp@data)[1] <- "layer"
  colnames(copernicuspoint.sp@data)[2] <- name_var
  # spplot(copernicuspoint.sp,name_var,cex = 0.01)
  
  gridpoint_cov <- merge(gridpoint_cov, copernicuspoint.sp@data, by = "layer",all.x=T)
  return(gridpoint_cov)
}


## Biological data

#' @title extract_copernicus.bio_rasterCelticBoB()
#' @param cov_file file where are located covariates
#' @param gridpoint grid points
#' @param gridpoint_cov grid points with covariates
#' @param copernicus.phys.var covariate to select in copernicus database
#' @param name_var name to give to this covariate
#' @param year_month month to select
#' @param grid_projection grid projection
#' 
#' @return gridpoint_cov gridpoint with physical covariates values


# year_month <- "2018-11"
# copernicus.bio.var <- "chl"
# name_var <- "chl.a"
extract_copernicus.bio_rasterCelticBoB <- function(cov_file,gridpoint,gridpoint_cov,copernicus.phys.var,name_var,year_month,grid_projection){
  # read variable in nc table
  nc_in <- open.nc(paste0(cov_file,"CopernicusData/dataset-ibi-reanalysis-bio-005-003-monthly_1576085356161_0m.nc"))
  # print.nc(nc_in) # variable
  nc_lon <- var.get.nc(nc_in, "longitude")
  nc_lat <- var.get.nc(nc_in, "latitude")
  nc_var <- var.get.nc(nc_in, copernicus.bio.var,unpack=TRUE)
  nc_t <- var.get.nc(nc_in, "time")
  nc_t <- as.POSIXct(nc_t*60*60, origin="1950-01-01",tz="GMT", format="%Y-%m-%d")
  close.nc(nc_in)
  
  # assuming a regular grid (which by the way appears not entirely true for this data set):
  gt <- GridTopology(cellcentre.offset = c(nc_lon[1], nc_lat[1]),
                     cellsize = c((tail(nc_lon, 1) - nc_lon[1])/(length(nc_lon) - 1), (tail(nc_lat, 1) - nc_lat[1])/(length(nc_lat) - 1)),
                     cells.dim = c(length(nc_lon), length(nc_lat)))
  sg <- SpatialGrid(gt, CRS(grid_projection))
  sgdf <- SpatialGridDataFrame(sg, data.frame(var = as.vector(nc_var[,,which(str_detect(nc_t,year_month))])))
  raster_var <- flip(raster(sgdf),direction = 2)
  # plot(flip(raster(sgdf),direction = 2))
  
  extract_raster_var <- raster::extract(stack(raster_var), coordinates(gridpoint))
  copernicuspoint.sp <- gridpoint
  copernicuspoint.sp@data <- as.data.frame(cbind(copernicuspoint.sp@data[,"layer"],extract_raster_var))
  colnames(copernicuspoint.sp@data)[1] <- "layer"
  colnames(copernicuspoint.sp@data)[2] <- name_var
  # spplot(copernicuspoint.sp,name_var,cex = 0.01)
  
  gridpoint_cov <- merge(gridpoint_cov, copernicuspoint.sp@data, by = "layer",all.x=T)
  return(gridpoint_cov)
}


##############################
# Reshape covariates dataframe 
##############################
# reshape dataframe so that each line is a cell and column is a (discrete) factor level filled 
# with 0 (level covariate is not in this cell) and 1 (level covariate is in this cell)

#------------------------
# Latent field covariates
#------------------------
#' @title reshape_cov_lf()
#' 
#' @param gridpoint_cov grid points with covariates
#' @param bathy_breaks bathymetry interval limits
#' @param latent_field_cov covariates affecting abundance disribution
#' @param useless_levels factor levels that are not used in the data
#' 
#' @return loc_x : dataframe in correct format to fit the model

# latent_field_cov <- c("bathy", "substr", "strata")
# useless_levels <- c("NA","Seabed")
reshape_cov_lf <- function(gridpoint_cov,bathy_breaks,latent_field_cov,useless_levels,discret_cov=NULL){
  loc_x <- gridpoint_cov %>%
    mutate(cell = 1:nrow(gridpoint_cov))
  
  ## strata
  if("strata" %in% discret_cov){
    loc_x$aggreg_substrate <- as.character(loc_x$aggreg_substrate)
    loc_x <- mutate(loc_x,strata = paste0("strata_",strata))
    loc_x$aggreg_substrate <- factor(loc_x$aggreg_substrate)
    
    loc_x <- loc_x %>%
      mutate(value = 1) %>%
      pivot_wider(names_from = strata,
                  values_from = value,
                  values_fill = list(value = 0))
  }
  
  # bathy
  if("bathy" %in% discret_cov){
    loc_x <- loc_x %>%
      mutate(value = 1) %>%
      mutate(bathy = cut(bathy,breaks = bathy_breaks)) %>%
      mutate(bathy = paste0("bathy_",bathy)) %>%
      pivot_wider(names_from = bathy,
                  values_from = value,
                  values_fill = list(value = 0))
  }
  
  # sediment
  if("aggreg_substrate" %in% discret_cov){
    loc_x <- loc_x %>%
      mutate(aggreg_substrate = paste0("substr_",aggreg_substrate)) %>%
      mutate(value = 1) %>%
      pivot_wider(names_from = aggreg_substrate,
                  values_from = value,
                  values_fill = list(value = 0))
  }
  
  # filter useless covariates levels
  for(lev in useless_levels){
    loc_x <- loc_x[,which(str_detect(colnames(loc_x),lev,negate = T))]
  }
  
  list_res <- list(loc_x = loc_x)
  return(list_res)
}
