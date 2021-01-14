commercial_reallocation_uniforme = function(k,
                                            xlim,
                                            ylim,
                                            y_com_i,
                                            n_samp_com,
                                            index_com_i,
                                            loc_x,
                                            sequencesdepeche,
                                            boats_number){
  
  y_com_i_new = cbind(y_com_i, boats_number)
  y_com_i_new = as.data.frame(y_com_i_new)
  y_com_i_new[, "pechemoyenne"] = rep(0, n_samp_com)
  y_com_i_new_2 = as.data.frame(y_com_i_new)
  y_com_i_new_2 = y_com_i_new_2 %>% dplyr::group_by(boats_number) %>% dplyr::summarise(y = mean(y_com_i))
  
  for (k in 1:n_samp_com)
  {
    value = which(y_com_i_new_2[, "boats_number"] == y_com_i_new[k, "boats_number"])
    y_com_i_new[k, "pechemoyenne"] = y_com_i_new_2[value, "y"]
  }
  
  y_com_i_new = as.vector(y_com_i_new[,"pechemoyenne"])
  
  return(y_com_i_new)
}
