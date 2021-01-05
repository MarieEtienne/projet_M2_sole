commercial_reallocation_uniforme = function(k,
                                            xlim,
                                            ylim,
                                            y_com_i,
                                            n_samp_com,
                                            index_com_i,
                                            loc_x,
                                            nboat){
  capturetotale <- sum(y_com_i)
  
  if (k==1){
    if (nboat==1){
      y_com_i_new <- rep(capturetotale/n_samp_com, n_samp_com)
    }
    
    else if (nboat > 1){
      # On transforme notre vecteur (3000 valeurs qui correspondent aux 3000 points de peche) en data.frame
      y_com_i_new = as.data.frame(y_com_i)
      # On ajoute une colonne qui correspond au bateau qui a réalisé la peche
      # On alloue les peches aux differents bateaux de façon aléatoire
      y_com_i_new[,"boat"] = sample(1:nboat, n_samp_com, replace=TRUE)
      # Pour chaque bateau, on calcule sa peche moyenne
      y_com_i_new_2 = y_com_i_new %>% dplyr::group_by(boat) %>% dplyr::summarise(y = mean(y_com_i))
      # On revient sur notre data.frame initial, et on remplace la peche de chaque bateau en chaque point par
      # la peche moyenne du bateau
      y_com_i_new[, "pechemoyenne"] = rep(0, n_samp_com)
      for (k in 1:n_samp_com){
        value = which(y_com_i_new_2[, "boat"] == y_com_i_new[k, "boat"])
        y_com_i_new[k, "pechemoyenne"] = y_com_i_new_2[value, "y"]
      }
      
      # On stocke le résultat final dans le vecteur que l'on va retourner
      y_com_i_new = as.vector(y_com_i_new[,"pechemoyenne"])
    }
  }
  
  else if (k==2){
    # D'abord on crée une table, y_com_i_new, qui contient :
    # - 3000 lignes correspondant aux 3000 points de pêche, aux 3000 échantillons commerciaux
    # - Pour chaque point de pêche, la quantité pêchée (issue de y_com_i)
    # - Pour chaque point de pêche, sa cellule et les coordonnées de celle-ci, x et y (issues de loc_x)
    index_com_i_new = as.data.frame(index_com_i)
    colnames(index_com_i_new) = "cell"
    y_com_i_new = as.data.frame(y_com_i)
    x = rep(0, n_samp_com)
    y = rep(0, n_samp_com)
    y_com_i_new = cbind(y_com_i_new, x, y, index_com_i_new)
    colnames(y_com_i_new) = c("y_com_i", "x", "y", "cell")
    for (i in 1:n_samp_com)
    {
      value = which(loc_x[, "cell"] == y_com_i_new[i, "cell"])
      y_com_i_new[i, "x"] = loc_x[value,"x"]
      y_com_i_new[i, "y"] = loc_x[value,"y"]
    }
    
    # On délimite les différentes zones
    y_com_i_new[, "zone"] = ifelse(y_com_i_new$x < xlim, 1, 2)
    
    if (nboat==1){
      # La moyenne des peches dans chaque zone
      mean_premierezone = mean(y_com_i_new[which(y_com_i_new[,"zone"] == 1), "y_com_i"])
      mean_deuxiemezone = mean(y_com_i_new[which(y_com_i_new[,"zone"] == 2), "y_com_i"])      
      
      # On re-attribue
      y_com_i_new$L_deuxzones = rep(0, nrow(y_com_i_new))
      y_com_i_new[which(y_com_i_new$zone == 1), "L_deuxzones"] = mean_premierezone
      y_com_i_new[which(y_com_i_new$zone == 2), "L_deuxzones"] = mean_deuxiemezone
      
      # On stocke le résultat final dans le vecteur que l'on va retourner
      y_com_i_new = as.vector(y_com_i_new$L_deuxzones)      
    }
    
    else if (nboat > 1){
      # On ajoute une colonne qui correspond au bateau qui a réalisé la peche
      # On alloue les peches aux differents bateaux de façon aléatoire
      y_com_i_new[,"boat"] = sample(1:nboat, n_samp_com, replace=TRUE)
      
      # Pour chaque bateau, on calcule sa peche moyenne par zone
      y_com_i_new_2 = y_com_i_new %>% dplyr::group_by(boat, zone) %>% dplyr::summarise(y = mean(y_com_i))
      
      # On revient sur notre data.frame initial, et on remplace la peche de chaque bateau en chaque point par
      # la peche moyenne du bateau dans la zone
      y_com_i_new[, "pechemoyenne"] = rep(0, n_samp_com)
      for (k in 1:n_samp_com){
        value = which(y_com_i_new_2[, "boat"] == y_com_i_new[k, "boat"] & y_com_i_new_2[, "zone"] == y_com_i_new[k, "zone"])
        y_com_i_new[k, "pechemoyenne"] = y_com_i_new_2[value, "y"]
      }      
      
      # On stocke le résultat final dans le vecteur que l'on va retourner
      y_com_i_new = as.vector(y_com_i_new[,"pechemoyenne"])
    }
  }
  
  else if (k==4){
    # D'abord on crée une table, y_com_i_new, qui contient :
    # - 3000 lignes correspondant aux 3000 points de pêche, aux 3000 échantillons commerciaux
    # - Pour chaque point de pêche, la quantité pêchée (issue de y_com_i)
    # - Pour chaque point de pêche, sa cellule et les coordonnées de celle-ci, x et y (issues de loc_x)
    
    # Code MPE a revoir
    #
    # dta_prov <-   data.frame(y = loc_x[index_com_i,"y"], x = loc_x[index_com_i,"x"] ) %>%
    #   mutate(n_cell = ifelse(x < xlim, 0, 1) + ifelse(y < ylim, 0, 1) * 2 + 1) 
    # 
    # dta_prov %>%
    #   group_by(n_cell) %>%
    #   dplyr::summarize(n = n())
    # 
    # data.frame(y = loc_x[index_com_i,"y"], x = loc_x[index_com_i,"x"], capt = y_com_i ) %>%
    #   mutate(n_cell = ifelse(x < xlim, 0, 1) + ifelse(y < ylim, 0, 1) * 2 + 1) %>%
    #   group_by(n_cell) %>%
    #   dplyr::summarize(capt_mean = mean(capt)) %>% inner_join(dta_prov) 
    #
    # Pour division en plusieurs zones voir switch plutot que ifelse
    # Sinon faire une fonction à qui on passe xlim et ylim et x et y et qui nous renvoie la zone
    
    index_com_i_new = as.data.frame(index_com_i)
    colnames(index_com_i_new) = "cell"
    y_com_i_new = as.data.frame(y_com_i)
    x = rep(0, n_samp_com)
  
    y = rep(0, n_samp_com)
    y_com_i_new = cbind(y_com_i_new, x, y, index_com_i_new)
    colnames(y_com_i_new) = c("y_com_i", "x", "y", "cell")
    for (i in 1:n_samp_com){
      value = which(loc_x[, "cell"] == y_com_i_new[i, "cell"])
      y_com_i_new[i, "x"] = loc_x[value,"x"]
      y_com_i_new[i, "y"] = loc_x[value,"y"]
    }
    
    # Ensuite on délimite les différentes zones
    y_com_i_new[, "zone"] = rep(0, n_samp_com)
    for (k in 1:n_samp_com){
      if (y_com_i_new[k, "x"] < xlim & y_com_i_new[k, "y"] < ylim){
        y_com_i_new[k, "zone"] = 1
      }
      else if (y_com_i_new[k, "x"] < xlim & y_com_i_new[k, "y"] >= ylim){
        y_com_i_new[k, "zone"] = 2
      }
      else if (y_com_i_new[k, "x"] >= xlim & y_com_i_new[k, "y"] < ylim){
        y_com_i_new[k, "zone"] = 3
      }
      else {
        y_com_i_new[k, "zone"] = 4
      }
    }
    
    if (nboat == 1){
      # La moyenne des peches dans chaque zone
      mean_premierezone = mean(y_com_i_new[which(y_com_i_new[,"zone"] == 1), "y_com_i"])
      mean_deuxiemezone = mean(y_com_i_new[which(y_com_i_new[,"zone"] == 2), "y_com_i"])
      mean_troisiemezone = mean(y_com_i_new[which(y_com_i_new[,"zone"] == 3), "y_com_i"])
      mean_quatriemezone = mean(y_com_i_new[which(y_com_i_new[,"zone"] == 4), "y_com_i"])
      
      # On re-attribue
      y_com_i_new$L_quatrezones = rep(0, nrow(y_com_i_new))
      y_com_i_new[which(y_com_i_new$zone == 1), "L_quatrezones"] = mean_premierezone
      y_com_i_new[which(y_com_i_new$zone == 2), "L_quatrezones"] = mean_deuxiemezone
      y_com_i_new[which(y_com_i_new$zone == 3), "L_quatrezones"] = mean_troisiemezone
      y_com_i_new[which(y_com_i_new$zone == 4), "L_quatrezones"] = mean_quatriemezone
      
      # On stocke le résultat final dans le vecteur que l'on va retourner
      y_com_i_new = as.vector(y_com_i_new$L_quatrezones)
    }
    
    else if (nboat > 1){
      # On ajoute une colonne qui correspond au bateau qui a réalisé la peche
      # On alloue les peches aux differents bateaux de façon aléatoire
      y_com_i_new[,"boat"] = sample(1:nboat, n_samp_com, replace=TRUE)
      
      # Pour chaque bateau, on calcule sa peche moyenne par zone
      y_com_i_new_2 = y_com_i_new %>% dplyr::group_by(boat, zone) %>% dplyr::summarise(y = mean(y_com_i))
      
      # On revient sur notre data.frame initial, et on remplace la peche de chaque bateau en chaque point par
      # la peche moyenne du bateau dans la zone
      y_com_i_new[, "pechemoyenne"] = rep(0, n_samp_com)
      for (k in 1:n_samp_com){
        value = which(y_com_i_new_2[, "boat"] == y_com_i_new[k, "boat"] & y_com_i_new_2[, "zone"] == y_com_i_new[k, "zone"])
        y_com_i_new[k, "pechemoyenne"] = y_com_i_new_2[value, "y"]
      }      
      
      # On stocke le résultat final dans le vecteur que l'on va retourner
      y_com_i_new = as.vector(y_com_i_new[,"pechemoyenne"])
    }
  }

  return(y_com_i_new)
}
