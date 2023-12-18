##############
## Load domain
##############

setwd(folder_phd_codes)

if(domain_shapefile == "ORHAGO"){
  
  strata_file <- "data/raw/shapefile/QGIS_ORHAGO/shape/" 
  
  zone_centre_large <- sf::st_read(paste0(strata_file,"zone centre au large.shp"))
  zone_centre_cote <-sf::st_read(paste0(strata_file,"zone centre cote.shp"))
  zone_nord <- sf::st_read(paste0(strata_file,"zone nord.shp"))
  zone_sud <- sf::st_read(paste0(strata_file,"zone sud.shp"))
  
  ## Orhago polygon
  colnames(zone_centre_large) <- c("id","Surface","geometry")
  colnames(zone_centre_cote) <- c("id","Surface","geometry")
  colnames(zone_nord) <- c("id","Surface","geometry")
  colnames(zone_sud) <- c("id","Surface","geometry")
  
  survey_layer <- st_union(rbind(zone_centre_large,zone_centre_cote,zone_nord,zone_sud))
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
                                      proj4string=CRS(grid_projection))
  
  crs(datapoint) <- grid_projection
  crs(study_domain_2) <- grid_projection
  datapoint_2 <- over(datapoint,study_domain_2)
  datapoint_2 <- cbind(datapoint_2,layer = datapoint@data$layer)
  datapoint_2 <- datapoint_2[which(!is.na(datapoint_2[,1])),]
  ## Add coordinates to gridpoint
  tmp <- cbind(datapoint@data, coordinates(datapoint))
  colnames(tmp) <- c(names(datapoint@data), "x", "y")
  datapoint_3 <- inner_join(as.data.frame(datapoint_2),tmp) %>%
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
  
}else if(domain_shapefile == "EVHOE"){
  
  strata_file <- "data/raw/shapefile/EVHOE/STRATES/Agreed_Strata_EVHOE_Polyg_WGS84.shp"
  limites <- sf::st_read(strata_file,crs=4326)
  survey_layer_3 <- sf::as_Spatial(limites)
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
  
  # Add coordinates to datapoint
  tmp <- cbind(datapoint@data, coordinates(datapoint))
  colnames(tmp) <- c(names(datapoint@data), "long", "lati")
  datapoint_4 <- inner_join(datapoint_3 , tmp)
  
  # Add strata
  datapoint_5 <- datapoint_4
  datapoint_5$strata <- datapoint_2$STRATE
  
}

setwd(folder_project)

