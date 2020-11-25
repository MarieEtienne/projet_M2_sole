###########################
## Simulate commercial data
###########################

#' @title simu_commercial_data()
#' 
#' @param loc_x dataframe with grid coordinate and strata
#' @param grid_dim grid dimension
#' @param n_cells number of cells in the grid
#' @param Strue_x spatial abundance data
#' 
#' # Commercial sampling
#' @param beta0_fb intercept of the sampling process
#' @param beta_fb fixed parameters of the sampling process
#' @param xfb_x covariate of the dampling process
#' 
#' # Observations
#' @param DataObs Observation model (1 : zero-inflated gamma, 2 : zero-inflated lognormal,3 : lognormal)
#' @param zero.infl_model Zero-inflated model parametrisation
#'
#' # Commercial
#' @param n_samp_com number of commercial samples
#' @param logSigma_com observation error for the scientific data (log scale)
#' @param q1_com parameter related to the zero-inflation for commercial data
#' @param q2_com relative catchability of the commercial data
#' @param b_fl dataframe or vector for b parameter values. There could be several fleets. In this case, fleet would be in columns and b values in lines.
#' @param zone_without_fishing if TRUE, part of the fishing area is not sampled, if FALSE, part of the fishing area is sampled
#' 
#' @return index_com_i : sampled cell for each observation
#' @return y_com_i : catch data
#' @return b_com_i : fleet related to each observation
#' @return c_com_x : number of sample in each cell (line) for each fleet (column)
#' 
#' @author B. Alglave


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
  
  index_com_i_2 = c()
  y_com_i_2 = c()
  b_com_i_2 = c()
  c_com_x_2 = c()
  
  Strue_fb <- log(Strue_x) # covariate related to abundance in the sampling intensity process
  
  loc_x <- cbind(loc_x,Strue_fb) # add log-scaled abundance to design matrix
  
  # grid to simulate point process  
  Int_ps <- matrix(ncol = grid_dim['x'],nrow = grid_dim['y'])
  
  # sampling intensity in each cell
  for(k in 1:n_cells){
    Int_ps[loc_x$y[k],loc_x$x[k]] <- exp(beta0_fb + b*Strue_fb[k])
  }
  
  win_df <- owin(xrange=c(0.5,grid_dim['x'] + 0.5), yrange=c(0.5,grid_dim['y'] + 0.5),ystep  = 1,xstep=1) # parameter of the window to simulate point process (fishing shoots)
  X <- rpoint(n_samp_com,as.im(Int_ps,win_df)) # simulate point process
  
  # plot(X)
  
  ## Aggregate data / discretise (use Raster, same as VMS data)
  #------------------------------------------------------------
  
  grid <- raster(extent(c(0.5,grid_dim['x']+0.5,0.5,grid_dim['y']+0.5))) # create grid
  res(grid) <- 1 # resolution of the grid
  # Transform this raster into a polygon 
  gridpolygon <- rasterToPolygons(grid)
  coord_grid <- as.data.frame(coordinates(gridpolygon))
  colnames(coord_grid)[1] <- "x"
  colnames(coord_grid)[2] <- "y"
  
  # # check that loc_x and gridpolygon have same coordinate system
  # library(ggrepel)
  # ggplot(loc_x)+geom_point(aes(x,y))+geom_text(aes(x,y,label = cell), size = 3.5)
  
  coord_grid <- inner_join(coord_grid,loc_x[,c("x","y","cell")],by = c("x","y"))
  gridpolygon$layer <- coord_grid$cell
  X_1 <- SpatialPoints(coords=cbind(X$x,X$y))
  # intersect of grid and VMS data
  X_2 <- over(X_1, gridpolygon)
  colnames(X_2)[1] <- "layer"
  X_2 %>%
    dplyr::count(layer) %>%
    arrange(layer) %>%
    data.frame()  -> X_3
  gridpolygon@data <- full_join(gridpolygon@data,X_3,by=c("layer"))
  gridpolygon@data %>%
    arrange(layer) %>%
    mutate(n = ifelse(is.na(n), 0,n)) %>%
    dplyr::select(n) %>%
    c() %>% unlist() %>% as.matrix() -> c_com_x
  
  # # Check point
  # library(prevR)
  # gridpolygon@data <- full_join(gridpolygon@data,loc_x,by=c("layer" = "cell"))
  # gridpolygon@data$c_Strue_x <- exp(gridpolygon@data$Strue_fb)
  # my.palette <- rev(prevR.colors.blue.inverse(50))[1:25] # rev(grep("^grey", colours(), value = TRUE))
  # simu_plot_2 <- spplot(gridpolygon, "c_Strue_x",col = "transparent", col.regions = my.palette) +
  #   layer(lpoints(X$x, X$y,col="black",cex = 0.1,pch=20))
  # plot(simu_plot_2)
  # trellis.device(device="png", filename=paste0("images/Simu/preferential_sampling_map_b_",b,".png"),height=300, width=300)
  # print(simu_plot_2)
  # dev.off()

  # grid.text("S(x)", x=unit(0.85, "npc"), y=unit(0.50, "npc"),
  #           gp = gpar(fontsize = 10, fontface = "bold"),
  #           rot=-90)
  
  ## For each shot simulate data (presence/absence and positive values)
  #--------------------------------------------------------------------
  
  # Generate a vector where number of cell are repeated as many times as the cell was sampled by fishers
  index_com_i <- do.call(c,lapply(1:length(c_com_x),function(k){
    n<-c_com_x[k]   # not clean --> call fishing_shots which is outside the function
    titi <- rep(k,n)
  }))
  
  
  # presence/absence data for commercial data
  q2_com <- 1 * q2_sci   ## commercial catchability
  y_com_i <- do.call(c,lapply(1:length(index_com_i), function(j){
    exp_catch <- q2_com * Strue_x[index_com_i[j]] # expected catch
    # proba of encounter
    prob_encount <- 1-exp(- exp(q1_com[1]) * exp_catch)
    abs_com_i <- rbinom(1,1,prob_encount)
    
    if(abs_com_i>0){
      y_com_i <- exp(rnorm(1,mean = log(exp_catch) - exp(logSigma_com)^2/2, sd = exp(logSigma_com)))
    }else{
      y_com_i <- 0
    }
  })) # q2_com * Strue_x : expected catch     |      q1_com : parameter linking number of zero and density for scientific data
  
  
  # fleet target factor
  b_com_i <- rep(i,length(y_com_i))
  
  # length((y_com_i[which(y_com_i > 0)]))/length(y_com_i)
  
  
  # # Plot stuff
  # verif_df_obs <- data.frame(cell = index_com_i,observation = y_com_i)
  # verif_df_S <- data.frame(cell = 1:n_cells,biom = Strue_x, counts = c_com_x)
  # verif_df <- inner_join(verif_df_obs,verif_df_S,by = "cell")
  # ggplot(data.frame(Abondance = Strue_x[index_sci_i],Captures_sci = y_sci_i),aes(x = Abondance,y = y_sci_i))+geom_point()
  # verif_1 <- ggplot(data.frame(Abondance = Strue_x,Nb_coups_peche = c_com_x),aes(x = log(Abondance),y = Nb_coups_peche))+geom_point()+theme_bw()
  # verif_2 <- ggplot(data.frame(Abondance = verif_df$biom,Captures_com = verif_df$observation),aes(x = Abondance,y = Captures_com))+geom_point()+theme_bw()
  # grid.arrange(verif_1,verif_2,ncol=2)
  # ggplot(verif_df,aes(x=observation)) +
  #   geom_histogram(aes(y=..density..),color="black", fill="white")
  
  index_com_i_2 = c(index_com_i_2,index_com_i) 
  y_com_i_2 = c(y_com_i_2,y_com_i)
  b_com_i_2 = c(b_com_i_2,b_com_i)
  c_com_x_2 = cbind(c_com_x_2,c_com_x)
  
  b_com_i_2 <- factor(b_com_i_2)
  
  res <- list(index_com_i = index_com_i_2, y_com_i = y_com_i_2, b_com_i = b_com_i_2, c_com_x = c_com_x_2)
  return(res)
}

# xfb_x <- latent_field[[i]]$xfb_x
# simu_commercial_data(loc_x,grid_dim,Strue_x,n_cells,beta0_fb,b,beta_fb,xfb_x,zone_without_fishing,n_samp_com,q1_com,q2_com,logSigma_com)
# 
# commercial_data <- list(list(),list())
# 
# for(i in 2:100){
#   print(i)
#   commercial_data <- list(commercial_data,list())
# }
# 
# # run the loop
# for(i in 1:100){
#   for(b in b_set){
#     
#     print(i)
#     set.seed( RandomSeed + i ) # i = 10
#     Strue_x <- latent_field[[i]]$Strue_x
#     xfb_x <- latent_field[[i]]$xfb_x
#     commercial_data[[i]][[which(b == b_set)]] <- simu_commercial_data(loc_x,grid_dim,Strue_x,n_cells,beta0_fb,b,beta_fb,xfb_x,zone_without_fishing,n_samp_com,q1_com,q2_com,logSigma_com)
#   }
# }
# save(data = latent_field,file = "C:/Users/test/Documents/GitHub/phd-baptiste_alglave/Results/Simu/simulated_data/scientific_data.RData")

