########################
## Shape scientific data
########################
ptsurvey_wgs84_2 <- ptsurvey_wgs84 %>%
  group_by_at(setdiff(names(ptsurvey_wgs84),c("Sex","CatchWgt_mature","CatCatchWgt","scientificname","CatchWgt_spp"))) %>%
  dplyr::summarise(CatchWgt_spp = sum(CatchWgt_spp)) %>%
  dplyr::rename(x=long,y=lati) %>%
  ungroup %>%
  dplyr::select(-layer) %>%
  filter(Year == year_y)

ptsurvey_wgs84_sf <- st_as_sf(ptsurvey_wgs84_2,coords = c("x","y"),crs=grid_projection)
ptsurvey_wgs84_sf_2 <- ptsurvey_wgs84_sf[st_intersects(ptsurvey_wgs84_sf,gridpolygon_sf) %>% lengths > 0,]
ptsurvey_wgs84_sf_2 <- st_join(ptsurvey_wgs84_sf_2,gridpolygon_sf)
ptsurvey_wgs84_sf_3 <- ptsurvey_wgs84_sf_2 %>%
  data.frame
ptsurvey_wgs84_sf_4 <- inner_join(ptsurvey_wgs84_sf_3,loc_x[,c("layer","cell","x","y")])

y_sci_i <- ptsurvey_wgs84_sf_4$CatchWgt_spp
index_sci_i <- ptsurvey_wgs84_sf_4$cell

## Covariate matrix
cov_x <- loc_x %>%
  dplyr::select_at(vars(starts_with("bathy") | starts_with("substr")))
