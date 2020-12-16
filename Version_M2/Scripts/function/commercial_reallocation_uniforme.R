commercial_reallocation_uniforme = function(reallocation,
                                            xlim,
                                            ylim,
                                            y_com_i,
                                            n_samp_com,
                                            index_com_i,
                                            loc_x){
  capturetotale <- sum(y_com_i)
  
  if (reallocation==1){
    y_com_i_new <- rep(capturetotale/n_samp_com, n_samp_com)
  }
  
  if (reallocation==2){
    # D'abord on crée une table, y_com_i_2_new, qui contient :
    # - 3000 lignes correspondant aux 3000 points de pêche, aux 3000 échantillons commerciaux
    # - Pour chaque point de pêche, la quantité pêchée (issue de y_com_i_2)
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
    
    # Ensuite on délimite les différentes zones
    premierezone = y_com_i_new[which(y_com_i_new$x < xlim),]
    deuxiemezone = y_com_i_new[which(y_com_i_new$x >= xlim),]
    
    # La moyenne des peches dans chaque zone
    mean_premierezone = mean(premierezone$y_com_i)
    mean_deuxiemezone = mean(deuxiemezone$y_com_i)
    
    # On re-attribue
    y_com_i_new$L_deuxzones = rep(0, nrow(y_com_i_new))
    y_com_i_new[which(y_com_i_new$x < xlim), "L_deuxzones"] = mean_premierezone
    y_com_i_new[which(y_com_i_new$x >= xlim), "L_deuxzones"] = mean_deuxiemezone

    # On stocke le résultat final dans le vecteur que l'on va retourner
    y_com_i_new = as.vector(y_com_i_new$L_deuxzones)
  }
  
  if (reallocation==4){
    # D'abord on crée une table, y_com_i_2_new, qui contient :
    # - 3000 lignes correspondant aux 3000 points de pêche, aux 3000 échantillons commerciaux
    # - Pour chaque point de pêche, la quantité pêchée (issue de y_com_i_2)
    # - Pour chaque point de pêche, sa cellule et les coordonnées de celle-ci, x et y (issues de loc_x)
    index_com_i_new = as.data.frame(index_com_i)
    colnames(index_com_i_new) = "cell"
    y_com_i_new = as.data.frame(y_com_i)
    x = rep(0, n_samp_com)
    y = rep(0, n_samp_com)
    y_com_i_new = cbind(y_com_i_new, x, y, index_com_i_new)
    colnames(y_com_i_new) = c("y_com_i_2", "x", "y", "cell")
    for (i in 1:n_samp_com)
    {
      value = which(loc_x[, "cell"] == y_com_i_new[i, "cell"])
      y_com_i_new[i, "x"] = loc_x[value,"x"]
      y_com_i_new[i, "y"] = loc_x[value,"y"]
    }
    
    # Ensuite on délimite les différentes zones
    premierezone = y_com_i_new[which(y_com_i_new$x < xlim & y_com_i_new$y < ylim),]
    deuxiemezone = y_com_i_new[which(y_com_i_new$x < xlim & y_com_i_new$y >= ylim),]
    troisiemezone = y_com_i_new[which(y_com_i_new$x >= xlim & y_com_i_new$y < ylim),]
    quatriemezone = y_com_i_new[which(y_com_i_new$x >= xlim & y_com_i_new$y >= ylim),]
    
    # La moyenne des peches dans chaque zone
    mean_premierezone = mean(premierezone$y_com_i)
    mean_deuxiemezone = mean(deuxiemezone$y_com_i)
    mean_troisiemezone = mean(troisiemezone$y_com_i)
    mean_quatriemezone = mean(quatriemezone$y_com_i)

    # On re-attribue
    y_com_i_new$L_quatrezones = rep(0, nrow(y_com_i_new))
    y_com_i_new[which(y_com_i_new$x < xlim & y_com_i_new$y < ylim), "L_quatrezones"] = mean_premierezone
    y_com_i_new[which(y_com_i_new$x < xlim & y_com_i_new$y >= ylim), "L_quatrezones"] = mean_deuxiemezone
    y_com_i_new[which(y_com_i_new$x >= xlim & y_com_i_new$y < ylim), "L_quatrezones"] = mean_troisiemezone
    y_com_i_new[which(y_com_i_new$x >= xlim & y_com_i_new$y >= ylim), "L_quatrezones"] = mean_quatriemezone

    # On stocke le résultat final dans le vecteur que l'on va retourner
    y_com_i_new = as.vector(y_com_i_new$L_quatrezones)
  }

  return(y_com_i_new)
}