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
  
  #initialisation
  index_com_i_2 = c() #index_com_i : sampled cell for each observation
  y_com_i_2 = c() # y_com_i : catch data
  b_com_i_2 = c() #b_com_i : fleet related to each observation
  c_com_x_2 = c() #c_com_x : number of sample in each cell (line) for each fleet (column)
  
  Strue_fb <- log(Strue_x) # covariate related to abundance in the sampling intensity process
  
  #on rajoute une nouvelle colonne (Strue_fb) à loc_x : il contient désormais les
  #colonnes x,y,cell, 4 col de strat binaires, et la col Strue_fb
  loc_x <- cbind(loc_x,Strue_fb) # add log-scaled abundance to design matrix
  
  # grid to simulate point process : grille remplie de NA de 25 lignes et 25 col
  #qui va contenir les points de peche commerciaux
  Int_ps <- matrix(ncol = grid_dim['x'],nrow = grid_dim['y'])
  
  # sampling intensity in each cell : Int_ps va contenir les intensites lambda
  #pour chaque point
  for(k in 1:n_cells){
    Int_ps[loc_x$y[k],loc_x$x[k]] <- exp(beta0_fb + b*Strue_fb[k])
  }
  
  #on simule le processus de poisson non homogene
  win_df <- owin(xrange=c(0.5,grid_dim['x'] + 0.5), yrange=c(0.5,grid_dim['y'] + 0.5),ystep  = 1,xstep=1) # parameter of the window to simulate point process (fishing shoots)
  X <- rpoint(n_samp_com,as.im(Int_ps,win_df)) # simulate point process
  
  #ici, pour le processus de poisson homogène, baptiste n'utilise pas rpoispp mais
  #rpoint
  #rpoispp on ne fixe pas le nombre d'échantillon que l'on veut mais les paramètres
  #avec rpoint, on fixe en entrée le nb d'échantillons que l'on veut (3000 ici)
  #rpoint est donc plus adapté à notre situation
  
  #rpoint prend en argument une image (et non une fonction lambda comme le fait rpoispp)
  #On simule le processus ponctuel de poisson à partir de cette image sur une 
  #structure continue
  #mais nous on veut que ce soit défini sur la grille : 
  #pour cela, on discrétise tout, en passant par des raster
  
  
  #plot(X) #representation des points de peche commerciaux simules dans la grille
  
  ## Aggregate data / discretise (use Raster, same as VMS data)
  #------------------------------------------------------------
  
 
  #, on croise le processus spatial avec la grille
  
  
  #definition de la grille raster : On commence par créer un raster de même 
  #dimension que la grille
  grid <- raster(extent(c(0.5,grid_dim['x']+0.5,0.5,grid_dim['y']+0.5))) # create grid
  res(grid) <- 1 # resolution of the grid
  # Transform this raster into a polygon 
  gridpolygon <- rasterToPolygons(grid) #on transforme le raster en polygone
  coord_grid <- as.data.frame(coordinates(gridpolygon)) #combinaisons de v1/v2 : 
  #tous les couples de coordonnées de la grilles
  colnames(coord_grid)[1] <- "x"
  colnames(coord_grid)[2] <- "y"
  coord_grid #les colonnes de coord_grid sont desormais appelées x et y 
  
  
  # # check that loc_x and gridpolygon have same coordinate system
  # library(ggrepel)
  # ggplot(loc_x)+geom_point(aes(x,y))+geom_text(aes(x,y,label = cell), size = 3.5)

  #on rajoute la col cell a coord_grid qui contient les numéros de cellules 
  #correspondant a chaque couple (x,y) grace a une jointure
  coord_grid <- inner_join(coord_grid,loc_x[,c("x","y","cell")],by = c("x","y"))
  gridpolygon$layer <- coord_grid$cell #on rentre la col des cellules dans le raster
  
  X_1 <- SpatialPoints(coords=cbind(X$x,X$y)) #on transforme X (points du processus 
  #de poisson) en spatial point ie on recupere les coord
  
  # intersect of grid and VMS data
  X_2 <- over(X_1, gridpolygon) #3000 lignes pour les 3000 points de peche commerciales
  #et une colonne layer (numero de la cellule correspondant a chacun des 3000 points de peche ? )
  colnames(X_2)[1] <- "layer"
  
  #creation de X3 : nombre de peches commerciales pour chacune des 625 cellules 
  X_2 %>%
    dplyr::count(layer) %>%
    arrange(layer) %>%
    data.frame()  -> X_3

  #le @ permet d’accéder à un champ de l’objet, c’est un peu l’équivalent du $
  
  gridpolygon@data <- full_join(gridpolygon@data,X_3,by=c("layer"))
  
  #creation de c_com_x : nombre de peches commerciales pour chacune des 625 cellules
  #rangé dans l'ordre croissant
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
  # vecteur des cellules : chaque cellule apparait autant de fois que de tentatives 
  #de peches qu'elle a recu, les cellules sont rangées dans l'ordre croissant
  #par exemple il y a eu 5 tentatives de peche dans la cellule 1 : le vecteur
  #commence par 4 uns puis 2 deux par ex 
  index_com_i <- do.call(c,lapply(1:length(c_com_x),function(k){
    n<-c_com_x[k]   # not clean --> call fishing_shots which is outside the function
    titi <- rep(k,n)
  }))
  
  
  # presence/absence data for commercial data
  q2_com <- 1 * q2_sci   ## commercial catchability : la capturabilité commerciale
  #dépend de la capturabilité scientifique
  
  #calcul de la valeur de capture des 3000 tentatives de peche commerciales 
  y_com_i <- do.call(c,lapply(1:length(index_com_i), function(j){
    exp_catch <- q2_com * Strue_x[index_com_i[j]] # expected catch #mu cf formule papier baptiste
    
    # proba of encounter : probabilité d'une peche fructueuse (cf papier baptiste)
    prob_encount <- 1-exp(- exp(q1_com[1]) * exp_catch) # c'est pour ca qu'il y a aucune peche ratée
    abs_com_i <- rbinom(1,1,prob_encount) #vaut 1 : que du succès pour chacun des 3000 points de peche
    
    #c'est l'équivalent de notre log normale pour calculer la quantité péchée
    #expliquer car cela ne correspond pas à la formule du papier 
    if(abs_com_i>0){
      # y_com_i <- rlnorm(1,meanlog = log(exp_catch), sdlog = exp(logSigma_com))
      # shape <- exp(logSigma_com)^-2
      # y_com_i <-  rgamma(1,shape=shape,scale= exp_catch * exp(logSigma_com)^2)
      y_com_i <- exp(rnorm(1,mean = log(exp_catch) - exp(logSigma_com)^2/2, sd = exp(logSigma_com)))
    }else{
      y_com_i <- 0
    }
  })) 
  # q2_com * Strue_x : exp_catch (expected catch)
  #q1_com : parameter linking number of zero and density for scientific data (c'est la ksi)
  #ici c'est bien une loi lognormale, ce n'est juste pas modélisé exactement de la
  #meme facon que dans le papier de baptiste (d’où le fait qu’on n’ait pas exactement la même moyenne)
  
  # fleet target factor
  #i est le numéro de de la boucle, le num de la simulation, varie entre 1 et 100
  #vecteur qui contient 3000 fois i
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
  
  #index_com_i_2 contient les index_com_i (vecteur des cellules : chaque cellule 
  #apparait autant de fois que de tentatives 
  #de peches qu'elle a recu, les cellules sont rangées dans l'ordre croissant
  #par exemple il y a eu 5 tentatives de peche dans la cellule 1 : le vecteur
  #commence par 4 uns) de l'ITERATION i
  index_com_i_2 = c(index_com_i_2,index_com_i) 
  
  #contient les valeurs de peche des 3000 tentatives de peche 
  #de l'iteration i
  #ici, y n'a pas d'unité mais en réalité l'unité est des CPUE donc catch per unit 
  #effort donc kg pêchés / heure de pêche et on estime que cette unité est représentative de la biomasse

  y_com_i_2 = c(y_com_i_2,y_com_i)
  
  #vecteur de 3000 * i
  b_com_i_2 = c(b_com_i_2,b_com_i)
  
  #nombre de peches commerciales pour chacune des 625 cellules
  #rangé dans l'ordre croissant pour l'itération i
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

