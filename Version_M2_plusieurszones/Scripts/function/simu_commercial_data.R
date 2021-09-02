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
                                 b,
                                 sequencesdepeche,
                                 zonespersequence,
                                 taillezone){

  # Initialisation
  index_com_i_2 = c() # index_com_i : sampled cell for each observation
  y_com_i_2 = c() # y_com_i : catch data
  b_com_i_2 = c() # b_com_i : fleet related to each observation
  c_com_x_2 = c() # c_com_x : number of sample in each cell (line) for each fleet (column)
  
  Strue_fb <- log(Strue_x) # covariate related to abundance in the sampling intensity process
  
  # On rajoute une nouvelle colonne (Strue_fb) à loc_x : il contient désormais les
  # colonnes x,y,cell, 4 col de strat binaires, et la col Strue_fb
  loc_x <- cbind(loc_x,Strue_fb) # add log-scaled abundance to design matrix
  
  # grid to simulate point process
  # Grille remplie de NA de 25 lignes et 25 col, qui va contenir les points de peche commerciaux
  Int_ps <- matrix(ncol = grid_dim['x'],nrow = grid_dim['y'])
  
  # sampling intensity in each cell
  # Int_ps va contenir les intensites lambda pour chaque point
  for(q in 1:n_cells){
    Int_ps[loc_x$y[q],loc_x$x[q]] <- exp(beta0_fb + b*Strue_fb[q])
  }
  
  
  # On commence par définir les centres de peche
  # Ces points de centres de zone sont définis selon un processus de poisson inhomogène qui dépend du champ latent
  win_df <- owin(xrange=c(0.5,grid_dim['x'] + 0.5), yrange=c(0.5,grid_dim['y'] + 0.5),ystep  = 1,xstep=1)
  X <- rpoint(sequencesdepeche*zonespersequence,as.im(Int_ps,win_df)) # Simulation des points centraux
  # X correspond donc aux coordonnées des sequencesdepeche*zonespersequence "points centraux" des zones de peche de chaque bateau
  
  # Visualisation des sequencesdepeche*zonespersequence centres des zones de peche de chaque bateau
  boats = numeric(length=(sequencesdepeche*zonespersequence))
  boats[1:sequencesdepeche] = 1:sequencesdepeche
  if ((sequencesdepeche*zonespersequence) > sequencesdepeche){
    for (q in (sequencesdepeche+1):(sequencesdepeche*zonespersequence))
    {
      boats[q] = sample(x=(1:sequencesdepeche), size=1)
    }
  }
  
  centres = as.data.frame(cbind(x = X$x, y = X$y, boats = boats))
  # plot_centres = ggplot(centres) + geom_point(aes(x=x, y=y, col=as.factor(boats)), size=4) +
  #   labs(color= "Bateau") +
  #   theme_bw() +
  #   theme(axis.title = element_blank(),legend.position='none')

  
  # La zone de peche de chaque bateau = point central +- 2 sur y, point central +- 2 sur x
  
  # Pour chaque bateau on simule 150 / ( zonespersequence * sequencesdepeche) points de peche dans sa zone de peche
  boats_x = numeric()
  boats_y = numeric()
  # On va stocker les coordonnées de ces 150 points de peche dans ces deux vecteurs
  
  # On définit un vecteur qui correspond au numéro du bateau associé à chaque point de peche
  boats = numeric(0)
  
  # On considère successivement chacune des sequences de peche * zonespersequence, donc chacune des zones de peche
  for (centre in 1:(sequencesdepeche*zonespersequence))
  {
    # On simule le processus de Poisson inhomogene, avec a chaque fois les coordonnées de centre de peche
    # pour l'objet owin, n_samp_com/(sequencesdepeche*zonespersequence) pour le rpoint,
    # et on passe dans le rpoint Int_ps qui correspond
    # aux lambdas qui dépendent du champ latent
    win_df_boat <- owin(xrange=c(max(0.5, X$x[centre]-taillezone), min(25.5, X$x[centre]+taillezone)), yrange=c(max(0.5, X$y[centre]-taillezone), min(25.5, X$y[centre]+taillezone)), ystep=1, xstep=1)
    n_samp <- round(n_samp_com/(sequencesdepeche*zonespersequence)) # + 1 --> sinon il manque des échantillons
    X_boat <- rpoint(n_samp, as.im(Int_ps, win_df_boat))
    
    # if()
    
    # On ajoute les coordonnées de peche obtenues pour le bateau
    boats_x = c(boats_x, X_boat$x)
    boats_y = c(boats_y, X_boat$y)
    boats = c(boats, rep(centres$boats[centre], n_samp))
    
  }
  
  ## Si les chiffres tombent pas ronds, il manquera qql points de données
  ## si pas assez de données simulées
  ## --> on complète et les pings restant sont réalloués de facon aléatoire au bateaux
  ## Si trop de données
  ## --> on rééchantillone aléatoirement pour avoir le bon nombre de données
  if(length(boats) < n_samp_com){
    
    n_samp_comp <-  n_samp_com - length(boats)
    X_boat_comp <- rpoint(n_samp_comp, as.im(Int_ps, win_df_boat))
    boats_x = c(boats_x, X_boat_comp$x)
    boats_y = c(boats_y, X_boat_comp$y)
    boats = c(boats, sample(boats,n_samp_comp,replace = T))

  }else if(n_samp_com < length(boats)){
    
    n_samp_comp <-  length(boats) - n_samp_com
    sample_i <- sample(1:n_samp_com,n_samp_com)
    boats_x = boats_x[sample_i]
    boats_y = boats_y[sample_i]
    boats = boats[sample_i]
    
    
  }
  
  # # sample n_samp_com in boats_y, boats_x, boats
  # index_com <- sample(1:(length(boats_y)),n_samp_com,replace = FALSE)
  # boats_x <- boats_x[index_com]
  # boats_y <- boats_y[index_com]
  # boats <- boats[index_com]

  # On a donc :
  # boats_x les coordonnées sur l'axe des x des 150 points de peche
  # boats_y les coordonnées sur l'axe des y des 150 points de peche
  # boats les numéros des sequences de peche associés à chaque point de peche
  
  # Visualisation des points de peche de chaque bateau dans sa zone
  peche_com_old = as.data.frame(cbind(x=boats_x, y=boats_y, boats=boats))
  # plot_pointsdepechecomperboat = ggplot(peche_com_old) + geom_point(aes(x=x, y=y, col=as.factor(boats)), size=4) +
  #   labs(color= "Bateau") +
  #   theme_bw() + 
  #   theme(axis.title = element_blank(),legend.position='none')


  ## Aggregate data / discretise (use Raster, same as VMS data)
  #------------------------------------------------------------
  
  # On croise le processus spatial avec la grille
  
  # Definition de la grille raster : On commence par créer un raster de même dimension que la grille
  grid <- raster(extent(c(0.5,grid_dim['x']+0.5,0.5,grid_dim['y']+0.5))) # create grid
  res(grid) <- 1 # resolution of the grid
  # Transform this raster into a polygon 
  gridpolygon <- rasterToPolygons(grid) # On transforme le raster en polygone
  coord_grid <- as.data.frame(coordinates(gridpolygon))
  # Combinaisons de v1/v2 : tous les couples de coordonnées de la grilles
  colnames(coord_grid)[1] <- "x"
  colnames(coord_grid)[2] <- "y"
  # coord_grid # Les colonnes de coord_grid sont desormais appelées x et y 
  
  # On rajoute la col cell a coord_grid qui contient les numéros de cellules 
  # correspondant a chaque couple (x,y) grace a une jointure
  coord_grid <- inner_join(coord_grid,loc_x[,c("x","y","cell")],by = c("x","y"))
  gridpolygon$layer <- coord_grid$cell # On rentre la col des cellules dans le raster
  
  X_1 <- SpatialPoints(coords=cbind(boats_x,boats_y))
  # On transforme X (points du processus de Poisson) en spatial point ie on recupere les coord
  
  # intersect of grid and VMS data
  X_2 <- over(X_1, gridpolygon)
  # 3000 lignes pour les 3000 points de peche commerciales,
  # et une colonne layer (numero de la cellule correspondant a chacun des 3000 points de peche ?)
  colnames(X_2)[1] <- "layer"
  
  # Creation de X3 : nombre de peches commerciales pour chacune des 625 cellules 
  X_2 %>%
    dplyr::count(layer) %>%
    arrange(layer) %>%
    data.frame()  -> X_3
  
  # Le @ permet d’accéder à un champ de l’objet, c’est un peu l’équivalent du $
  
  gridpolygon@data <- full_join(gridpolygon@data,X_3,by=c("layer"))
  
  # Creation de c_com_x : nombre de peches commerciales pour chacune des 625 cellules, rangées dans l'ordre croissant
  gridpolygon@data %>%
    arrange(layer) %>%
    mutate(n = ifelse(is.na(n), 0,n)) %>%
    dplyr::select(n) %>%
    c() %>% unlist() %>% as.matrix() -> c_com_x
  

  ## For each shot simulate data (presence/absence and positive values)
  #--------------------------------------------------------------------
  
  # Generate a vector where number of cell are repeated as many times as the cell was sampled by fishers
  # Vecteur des cellules : chaque cellule apparait autant de fois que de tentatives 
  # de peches qu'elle a recu, les cellules sont rangées dans l'ordre croissant
  index_com_i <- do.call(c,lapply(1:length(c_com_x),function(q){
    n<-c_com_x[q]   # not clean --> call fishing_shots which is outside the function
    titi <- rep(q,n)
  }))
  
  
  # presence/absence data for commercial data
  q2_com <- 1 # commercial catchability
  
  # Calcul de la valeur de capture des tentatives de peche commerciales 
  y_com_i <- do.call(c,lapply(1:length(index_com_i), function(j){
    exp_catch <- q2_com * Strue_x[index_com_i[j]] # expected catch #mu --> cf formule papier baptiste
    # proba of encounter # probabilité d'une peche fructueuse (cf papier baptiste)
    prob_encount <- 1-exp(- exp(q1_com[1]) * exp_catch) # Aucune peche ratée
    abs_com_i <- rbinom(1,1,prob_encount) # Vaut 1 : que du succès pour chacun des points de peche
    
    # C'est l'équivalent de notre log normale pour calculer la quantité péchée
    # Cela ne correspond pas à la formule du papier 
    if(abs_com_i>0){
      # y_com_i <- rlnorm(1,meanlog = log(exp_catch), sdlog = exp(logSigma_com))
      # shape <- exp(logSigma_com)^-2
      # y_com_i <-  rgamma(1,shape=shape,scale= exp_catch * exp(logSigma_com)^2)
      y_com_i <- exp(rnorm(1,mean = log(exp_catch/prob_encount) - exp(logSigma_com)^2/2, sd = exp(logSigma_com)))
    }else{
      y_com_i <- 0
    }
  })) 
  # q2_com * Strue_x : exp_catch (expected catch)
  # q1_com : parameter linking number of zero and density for scientific data (c'est la ksi)
  # Ici c'est bien une loi lognormale, ce n'est juste pas modélisé exactement de la
  # meme facon que dans le papier de Baptiste (d’où le fait qu’on n’ait pas exactement la même moyenne)
  
  # fleet target factor
  # i est le numéro de de la boucle, le num de la simulation, varie entre 1 et 100
  # Vecteur qui contient 150 fois i
  b_com_i <- rep(i,length(y_com_i))
  
  # index_com_i_2 contient les index_com_i (vecteur des cellules : chaque cellule 
  # apparait autant de fois que de tentatives 
  # de peches qu'elle a recu, les cellules sont rangées dans l'ordre croissant
  # Par exemple il y a eu 5 tentatives de peche dans la cellule 1 : le vecteur
  # commence par 5 uns --> pour l'iteration i
  index_com_i_2 = c(index_com_i_2,index_com_i) 
  
  # Contient les valeurs de peche des 3000 tentatives de peche  de l'iteration i
  # Ici, y n'a pas d'unité mais en réalité l'unité est des CPUE donc catch per unit 
  # effort donc kg pêchés / heure de pêche et on estime que cette unité est représentative de la biomasse
  y_com_i_2 = c(y_com_i_2,y_com_i)
  
  # Vecteur de 3000 * i
  b_com_i_2 = c(b_com_i_2,b_com_i)
  
  # Nombre de peches commerciales pour chacune des 625 cellules, rangées dans l'ordre croissant pour l'itération i
  c_com_x_2 = cbind(c_com_x_2,c_com_x)
  
  b_com_i_2 <- factor(b_com_i_2)
  
  # # Visualisation des points de peche de chaque bateau dans sa zone et de la quantité pechee
  # peche_com_boat = peche_com_old
  # peche_com_boat$x = round(peche_com_boat$x)
  # peche_com_boat$y = round(peche_com_boat$y)
  # peche_com_boat$ncell = rep(0, length(peche_com_boat$boats))
  # for (j in 1:length(peche_com_boat$boats))
  # {
  #   peche_com_boat[j, "ncell"] = loc_x$cell[which(loc_x[, "x"] == peche_com_boat[i, "x"] & loc_x[, "y"] == peche_com_boat[j, "y"])]
  # }
  # qte_pechee = as.data.frame(cbind(ncell = index_com_i, y_com_i = y_com_i))
  # peche_com_boat = peche_com_boat[order(peche_com_boat$ncell),]
  # qte_pechee = qte_pechee[order(qte_pechee$ncell),]
  # peche_com_boat = cbind(peche_com_boat, y_com_i = qte_pechee$y_com_i)

  
  res <- list(index_com_i = index_com_i, y_com_i = y_com_i, b_com_i = b_com_i_2, c_com_x = c_com_x_2, boats_number = boats, x_com = boats_x, y_com = boats_y, centres = centres, peche_com_old = peche_com_old)
  return(res)
}


# test <- data.frame(y_com = y_com_i, cell = index_com_i, boats = boats, coord_x_com = boats_x, coord_y_com = boats_y)
# test_2 <- full_join(loc_x,test) # il faut qu'il y ait les valeurs du champ latent dans S_x
# ggplot(test_2)+
#   # geom_point(aes(x=x,y=y)) +
#   geom_point(aes(x=x,y=y,col=S_x),size=5,shape=15,alpha=0.2)+
#   geom_point(aes(x=coord_x_com,y=coord_y_com,size=y_com),col="red")

  
