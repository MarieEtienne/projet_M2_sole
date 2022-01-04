#########################
## Shape 'VMS x logbooks'
#########################

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

points_sf <- points_sf %>%
  # filter(VE_REF %in% vessel_vec) %>%
  mutate(year = as.numeric(str_sub(SI_DATE,7,10)),
         month = as.numeric(str_sub(SI_DATE,4,5)))

points_sf$quarter <- NA
points_sf$quarter[which(points_sf$month %in% c(1,2,3))] <- 1
points_sf$quarter[which(points_sf$month %in% c(4,5,6))] <- 2
points_sf$quarter[which(points_sf$month %in% c(7,8,9))] <- 3
points_sf$quarter[which(points_sf$month %in% c(10,11,12))] <- 4

points_sf_2 <- points_sf %>%
  filter(str_detect(LE_MET_level6,"OTB_CEP_>=70_0|OTB_DEF_>=70_0")) %>%
  filter(month == month_m & year == year_y)

## cross datapoints with grid
#----------------------------
gridpolygon_sf <- st_as_sf(gridpolygon)
points_sf_3 <- points_sf_2[st_intersects(points_sf_2,gridpolygon_sf) %>% lengths > 0,]
points_sf_3 <- st_join(points_sf_3,gridpolygon_sf)
colnames(points_sf_3)[which(colnames(points_sf_3) == paste0("LE_KG_",species_to_plot))] <- "LE_KG_spp"
points_sf_4 <- points_sf_3 %>%
  dplyr::select(SI_DATE,LE_MET_level6,LE_KG_spp,ICESNAME,seq_id,layer,INTV)
points_df <- points_sf_4 %>% data.frame
points_df <- inner_join(points_df,loc_x[,c("layer","cell","x","y","ICESNAME")])
points_df <- points_df %>%
  rename_at(vars(starts_with("LE_KG")),function(x){x <- "LE_KG"}) %>%
  group_by(layer,cell,x,y,ICESNAME,seq_id,LE_MET_level6) %>%
  dplyr::summarise(HF = sum(INTV)/60,
                   LE_KG = sum(LE_KG)) %>%
  mutate(CPUE = LE_KG/HF)

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

boats_i <- as.numeric(factor(points_df$seq_id))
boats_i <- boats_i
index_i <- points_df$cell
y_i2 <- points_df$CPUE

## Reallocation
mean_y_i <- aggregate(x = y_i2,
                      by = list(unique.values = boats_i),
                      FUN = mean)[,2]

sum_y_i <- aggregate(x = y_i2,
                     by = list(unique.values = boats_i),
                     FUN = sum)[,2]
