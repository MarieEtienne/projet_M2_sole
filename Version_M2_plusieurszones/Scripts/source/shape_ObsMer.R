###############
## Shape ObsMer
###############

HH_obsmer_2 <- HH_obsmer %>%
  mutate(y=(latIni+latFin)/2,x=(lonIni+lonFin)/2) %>%
  dplyr::select(-latIni,-latFin,-lonIni,-lonFin) %>%
  mutate(month=as.numeric(str_sub(date,6,7)))

HH_obsmer_2$quarter <- NA
HH_obsmer_2$quarter[which(HH_obsmer_2$month %in% c(1,2,3))] <- 1
HH_obsmer_2$quarter[which(HH_obsmer_2$month %in% c(4,5,6))] <- 2
HH_obsmer_2$quarter[which(HH_obsmer_2$month %in% c(7,8,9))] <- 3
HH_obsmer_2$quarter[which(HH_obsmer_2$month %in% c(10,11,12))] <- 4

HH_obsmer_2 <- HH_obsmer_2 %>%
  filter(str_detect(foCatEu6,"OTB_CEP_>=70_0|OTB_DEF_>=70_0")) %>%
  filter(year == year_y) %>%
  filter(y<48 & quarter == quarter_q) # month == month_m)

SL_obsmer_2 <- SL_obsmer %>%
  mutate(Species = str_replace_all(spp," ","_")) %>%
  filter(Species == species_to_plot)

ObsMer_df <- left_join(HH_obsmer_2,SL_obsmer_2) %>%
  dplyr::group_by_all() %>%
  dplyr::summarise(wt=sum(wt)) %>%
  mutate(Catch = ifelse(is.na(wt),0,wt)) %>%
  mutate(Catch = Catch/1000,
         foDur = foDur / 60) %>%
  mutate(CPUE = Catch/foDur)

# convert the points to an sf object
ObsMer_sf <- st_as_sf(ObsMer_df, coords = c("x", "y"))

# set CRS for the points to be the same as shapefile
st_crs(ObsMer_sf) <- st_crs(grid_projection)

ObsMer_sf_2 <- ObsMer_sf[st_intersects(ObsMer_sf,gridpolygon_sf) %>% lengths > 0,]
ObsMer_sf_2 <- st_join(ObsMer_sf_2,gridpolygon_sf)
colnames(ObsMer_sf_2)[which(colnames(ObsMer_sf_2) == paste0("LE_KG_",species_to_plot))] <- "LE_KG_spp"
ObsMer_sf_3 <- ObsMer_sf_2 %>%
  dplyr::select(date,year,time,foDur,foCatEu6,month,spp,quarter,wt,CPUE,layer)
ObsMer_df <- ObsMer_sf_3 %>% data.frame
ObsMer_df <- inner_join(ObsMer_df,loc_x[,c("layer","cell","x","y","ICESNAME")])

index_ObsM_i <- ObsMer_df$cell
y_ObsM_i <- ObsMer_df$CPUE
