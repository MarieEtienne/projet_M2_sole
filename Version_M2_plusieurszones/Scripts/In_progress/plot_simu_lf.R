#-------------------
# Build grid polygon
#-------------------
grid_xmin <- 0
grid_xmax <- grid_dim['x']
grid_ymin <- 0
grid_ymax <- grid_dim['y']
grid_limit <- extent(c(grid_xmin,grid_xmax,grid_ymin,grid_ymax))

grid <- raster(grid_limit)
res(grid) <- 1
gridpolygon <- rasterToPolygons(grid)
gridpolygon$layer <- c(1:length(gridpolygon$layer))
gridpolygon_sf <- st_as_sf(gridpolygon)

#------------
# Build loc_x
#------------
# create grid and cells + strata
n_cells <- grid_dim['x'] * grid_dim['y'] # 25*25 = 625
loc_x = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y']) #toutes les 
# Combinaisons possibles de (x,y)
loc_x = cbind(loc_x,cell = 1:n_cells) # Chaque combinaison possible est une cellule, numerotée de 1 à 625
diff.dim_y_x <- grid_dim['y'] - grid_dim['x'] #25-25 =0
loc_x %>%
  mutate(strata_1 = ifelse(x <= n_strate & y > grid_dim['y'] - n_strate , 1, 0)) %>%
  mutate(strata_4 = ifelse(x >= grid_dim['x'] - n_strate & y <= n_strate + diff.dim_y_x, 1, 0), 1, 0) %>%
  mutate(strata_3 = ifelse((x < grid_dim['x'] - n_strate & y <= grid_dim['y'] - n_strate & x + y <= grid_dim['y'] + 1),1,0)) %>%
  mutate(strata_2 = ifelse(strata_1 == 0 & strata_4 == 0 & strata_3 == 0, 1, 0)) %>%
  dplyr::select(x,y,cell,strata_1,strata_2,strata_3,strata_4) -> loc_x

#---------------
#  Latent field
#---------------
beta[1] = 2

# Simulate using Matérn covariance
delta_x = sim_GF_Matern(loc_x, nu, range, SD_delta^2)[[1]]$y

# Total abundance
Strue_x = (beta0 + (loc_x$x) / 10  * 1 + delta_x ) # + delta_x

gridpolygon_sf_2 <- cbind(gridpolygon_sf,Strue_x)

plot_1 <- ggplot(gridpolygon_sf_2)+
  geom_sf(aes(fill=Strue_x),alpha=0.9)+
  scale_fill_distiller(palette = "Spectral")+
  theme_void()+theme(legend.position = "none")





