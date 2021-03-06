generation_map = function(path_base,
                          simu_file,
                          i,
                          b_set,
                          reallocation,
                          x){
  ### On convertit les x en couple PZ
  if (x==1){
    sequencesdepeche=2
    zonespersequence=1
  } else if (x==2){
    sequencesdepeche=2
    zonespersequence=3
  } else if (x==3){
    sequencesdepeche=2
    zonespersequence=5
  } else if (x==4){
    sequencesdepeche=5
    zonespersequence=1
  } else if (x==5){
    sequencesdepeche=5
    zonespersequence=3
  } else if (x==6){
    sequencesdepeche=5
    zonespersequence=5
  } else if (x==7){
    sequencesdepeche=10
    zonespersequence=1
  } else if (x==8){
    sequencesdepeche=10
    zonespersequence=3
  } else if (x==9){
    sequencesdepeche=10
    zonespersequence=5
  } else{ # Cas ou x = 10 et donc on est dans le cas de reference
    sequencesdepeche = 1
    zonespersequence = 1
  }
  
  
  for (b in b_set){ # On va travailler successivement sur les differentes valeurs de b
    for (k in reallocation){ # On va travailler successivement sur le cas sans reallocation et avec reallocation
      ### On récupère les données
      
      path = paste0(path_base, simu_file, "map_results_i",i,"_b",b,"_k",k,"_x",x,".RData")
      load(path)
      Strue_x_2 = map_results[[1]]
      centres = map_results[[2]]
      peche_com_old = map_results[[3]]
      pointsdepeche_com = map_results[[4]]
      
      
      ### On crée les maps
      
      # 1ère map : champ latent reel
      plot_truelatentfield = ggplot(Strue_x_2) + geom_point(aes(x,y,col=Champ_latent_reel), size=4) + theme_bw() +
        scale_color_gradient2(midpoint = mean(Strue_x_2$Champ_latent_reel), low = "#E6F2FC", mid = "#62B4FC",
                              high = "#02182C", space = "Lab" ) +
        labs(color= "Champ latent") +
        xlim(0.5, 25.5) +
        ylim(0.5, 25.5) +
        labs(title = paste0("Representation du champ latent reel pour i = ", i))

      # 2ème map : champ latent estime
      
      # D'abord on refait le champ latent reel mais sans titre et avec un facet_wrap
      plot_latentfield = ggplot(Strue_x_2) + geom_point(aes(x,y,col=Champ_latent_reel), size=4) + theme_bw() +
        scale_color_gradient2(midpoint = mean(Strue_x_2$Champ_latent_reel), low = "#E6F2FC", mid = "#62B4FC",
                              high = "#02182C", space = "Lab" ) +
        labs(color= "Champ latent") +
        xlim(0.5, 25.5) +
        ylim(0.5, 25.5) +
        facet_wrap(~"Reel") +
        theme(axis.title = element_blank(),legend.position='none')
      
      # Puis le champ latent estime pour cette combinaison i*ZP*b*k
      plot_estlatentfield = ggplot(Strue_x_2) + geom_point(aes(x,y,col=Champ_latent_estime), size=4) + theme_bw() +
        scale_color_gradient2(midpoint = mean(Strue_x_2$Champ_latent_estime), low = "#E6F2FC", mid = "#62B4FC",
                              high = "#02182C", space = "Lab" ) +
        labs(color= "Champ latent") +
        xlim(0.5, 25.5) +
        ylim(0.5, 25.5) +
        theme(axis.title = element_blank(),legend.position='none')
      if (k==reallocation[1]){
        plot_estlatentfield <- plot_estlatentfield + facet_wrap(~"Estime sans reallocation")
      } else {
        plot_estlatentfield <- plot_estlatentfield + facet_wrap(~"Estime avec reallocation")
      }
      
      # 3ème map : centres des zones de peche
      plot_centres = ggplot(centres) + geom_point(aes(x=x, y=y, col=as.factor(boats)), size=4) +
        labs(color= "Sequence de peche") +
        theme_bw() +
        xlim(0.5, 25.5) +
        ylim(0.5, 25.5) +
        theme(axis.title = element_blank(),legend.position='none')
      if (b==b_set[1]){
        plot_centres <- plot_centres + facet_wrap(~paste0("b = ", b_set[1]))
      }else if (b==b_set[2]){
        plot_centres <- plot_centres + facet_wrap(~paste0("b = ", b_set[2]))
      } else {
        plot_centres <- plot_centres + facet_wrap(~paste0("b = ", b_set[3]))
      }
      
      # 4ème map : points de peche
      plot_pointsdepechecomperboat = ggplot(peche_com_old) +
        geom_point(aes(x=x, y=y, col=as.factor(boats)), size=4) +
        labs(color= "Sequence de peche") +
        theme_bw() +
        xlim(0.5, 25.5) +
        ylim(0.5, 25.5) +
        theme(axis.title = element_blank(),legend.position='none')
      if (b==b_set[1]){
        plot_pointsdepechecomperboat <- plot_pointsdepechecomperboat + facet_wrap(~paste0("b = ", b_set[1]))
      }else if (b==b_set[2]){
        plot_pointsdepechecomperboat <- plot_pointsdepechecomperboat + facet_wrap(~paste0("b = ", b_set[2]))
      } else {
        plot_pointsdepechecomperboat <- plot_pointsdepechecomperboat + facet_wrap(~paste0("b = ", b_set[3]))
      }
      
      # 5ème map : points de peche dans les cellules
      plot_pointsdepechecell = ggplot(pointsdepeche_com) +
        geom_point(aes(x=x, y=y, col=as.factor(boats)), size=4) +
        labs(color= "Sequence de peche") +
        theme_bw() +
        xlim(0.5, 25.5) +
        ylim(0.5, 25.5) +
        theme(axis.title = element_blank(),legend.position='none')
      if (b==b_set[1]){
        plot_pointsdepechecell <- plot_pointsdepechecell + facet_wrap(~paste0("b = ", b_set[1]))
      }else if (b==b_set[2]){
        plot_pointsdepechecell <- plot_pointsdepechecell + facet_wrap(~paste0("b = ", b_set[2]))
      } else {
        plot_pointsdepechecell <- plot_pointsdepechecell + facet_wrap(~paste0("b = ", b_set[3]))
      }
      
      # 6ème map : quantité péchée
      if (k==reallocation[1]){
        plot_pointsdepecheqte = ggplot(pointsdepeche_com) +
          geom_point(aes(x=x, y=y, col=y_com_i), size=4) +
          labs(color= "Quantite pechee") +
          theme_bw() +
          xlim(0.5, 25.5) +
          ylim(0.5, 25.5) +
          scale_color_gradient2(midpoint = mean(pointsdepeche_com$y_com_i),
                                low = "#E6F2FC", mid = "#62B4FC",
                                high = "#02182C", space = "Lab" ) +
          facet_wrap(~"Sans reallocation")
      } else {
        plot_pointsdepecheqte = ggplot(pointsdepeche_com) +
          geom_point(aes(x=x, y=y, col=as.factor(round(y_com_i, 1))), size=4) +
          labs(color= "Quantite pechee") +
          theme_bw() +
          xlim(0.5, 25.5) +
          ylim(0.5, 25.5) +
          facet_wrap(~"Avec reallocation")
      }
      
      
      ### Maintenant qu'on a crée les map, on les sauve dans des objets dépendant de b
      
      if (b==b_set[1]) {
        plot_centres_b0 = plot_centres
        plot_pointsdepechecomperboat_b0 = plot_pointsdepechecomperboat
        plot_pointsdepechecell_b0 = plot_pointsdepechecell
      } else if (b==b_set[2]) {
        plot_centres_b1 = plot_centres
        plot_pointsdepechecomperboat_b1 = plot_pointsdepechecomperboat
        plot_pointsdepechecell_b1 = plot_pointsdepechecell
      } else {
        plot_centres_b3 = plot_centres
        plot_pointsdepechecomperboat_b3 = plot_pointsdepechecomperboat
        plot_pointsdepechecell_b3 = plot_pointsdepechecell
      }
      
      ## Et dans des objets dépendant de k
      if (k==reallocation[1]){
        plot_estlatentfield_k0 = plot_estlatentfield
        plot_pointsdepecheqte_k0 = plot_pointsdepecheqte
      } else {
        plot_estlatentfield_k1 = plot_estlatentfield
        plot_pointsdepecheqte_k1 = plot_pointsdepecheqte
      }
    }
    
    ### On a fini de boucler sur les k, on ggarrange les 3 graphes necessaires pour le champ latent estime
    plot_estlatentfield = ggarrange(plot_latentfield, plot_estlatentfield_k0, plot_estlatentfield_k1,
                                    ncol = 3, common.legend=TRUE, legend="right")
    plot_estlatentfield <- annotate_figure(plot_estlatentfield,
                                    top = text_grob(paste0("Representation du champ latent estime pour i = ", i, ", b = ", b, " ; ", sequencesdepeche, " sequence(s) de peche et ", zonespersequence, " centre(s) par sequence en moyenne"), size = 15),
                                    left = text_grob("y",
                                                     rot = 90,
                                                     size = 11),
                                    bottom = text_grob("x", size = 11))
    
    # Et les 2 graphes necessaires pour les quantites pechees
    plot_pointsdepecheqte = ggarrange(plot_pointsdepecheqte_k0,
                                      plot_pointsdepecheqte_k1,
                                      ncol = 2)
    plot_pointsdepecheqte <- annotate_figure(plot_pointsdepecheqte,
                                           top = text_grob(paste0("Quantite pechee par cellule pour i = ", i, ", b = ", b, " ; ", sequencesdepeche, " sequence(s) de peche et ", zonespersequence, " centre(s) par sequence en moyenne"), size = 15),
                                           left = text_grob("y",
                                                            rot = 90,
                                                            size = 11),
                                           bottom = text_grob("x", size = 11))
    
    
    ## On sauvegarde les maps de champ latent estime
    # On le fait plus tot que les autres maps car celles ci dependent de b
    
    simu_file_plots <- paste0("results_map/", simu_file)
    if(! dir.exists(simu_file_plots)) dir.create(simu_file_plots)
    
    # On a un graphe par i*b*x
    path = paste0(simu_file_plots, "latentfieldest_i",i,"_b",b,"_x",x,".png")
    png(file = path, width = 1600, height = 500)
    plot(plot_estlatentfield)
    dev.off()

    path = paste0(simu_file_plots, "pecheqte_i",i,"_b",b,"_x",x,".png")
    png(file = path, width = 1200, height = 500)
    plot(plot_pointsdepecheqte)
    dev.off()
  }
  
  
  ### On ggarrange les differents objets pour diminuer le nombre de graphes sortis :
  # On crée un graphe pour toutes les valeurs de b
  
  # 3eme map : les centres
  plot_centres = ggarrange(plot_centres_b0, plot_centres_b1, plot_centres_b3, ncol=3,
                           common.legend=TRUE, legend="right")
  plot_centres <- annotate_figure(plot_centres,
                                  top = text_grob(paste0("Centres des zones de peche pour i = ", i, ", ", sequencesdepeche, " sequence(s) de peche et ", zonespersequence, " centre(s) par sequence en moyenne"), size = 15),
                                  left = text_grob("y",
                                                   rot = 90,
                                                   size = 11),
                                  bottom = text_grob("x", size = 11))
  
  # 4eme map : les points de peche  
  plot_pointsdepechecomperboat = ggarrange(plot_pointsdepechecomperboat_b0,
                                           plot_pointsdepechecomperboat_b1,
                                           plot_pointsdepechecomperboat_b3,
                                           ncol=3,
                                           common.legend=TRUE, legend="right")
  plot_pointsdepechecomperboat <- annotate_figure(plot_pointsdepechecomperboat,
                                                  top = text_grob(paste0("Points de peches pour i = ", i, ", ", sequencesdepeche, " sequence(s) de peche et ", zonespersequence, " centre(s) par sequence en moyenne"), size = 15),
                                                  left = text_grob("y",
                                                                   rot = 90,
                                                                   size = 11),
                                                  bottom = text_grob("x", size = 11))
  
  # 5eme map : les points de peche par cellule
  plot_pointsdepechecell = ggarrange(plot_pointsdepechecell_b0, plot_pointsdepechecell_b1,
                                     plot_pointsdepechecell_b3,
                                     ncol=3,
                                     common.legend=TRUE, legend="right")
  plot_pointsdepechecell <- annotate_figure(plot_pointsdepechecell,
                                            top = text_grob(paste0("Points de peche dans les cellules pour i = ", i, ", ", sequencesdepeche, " sequence(s) de peche et ", zonespersequence, " centre(s) par sequence en moyenne"), size = 15),
                                            left = text_grob("y",
                                                             rot = 90,
                                                             size = 11),
                                            bottom = text_grob("x", size = 11))
  
  
  ### On sauvegarde les maps crées dans results_map
  
  simu_file_plots <- paste0("results_map/", simu_file)
  if(! dir.exists(simu_file_plots)) dir.create(simu_file_plots)
  
  # Pour le champ latent (reel ou estime), on a un graphe / i
  if (x==1){
    path = paste0(simu_file_plots, "latentfieldtrue_i",i,".png")
    png(file = path, width = 600, height = 500)
    plot(plot_truelatentfield)
    dev.off()
  }
  
  # Pour les centres, les points de peche et les peches par cellule, on a un graphe / i*x
  # Pour ces 3 maps, toutes les valeurs de b sont representees sur le meme graphe
  path = paste0(simu_file_plots, "centres_i",i,"_x",x,".png")
  png(file = path, width = 1600, height = 500)
  plot(plot_centres)
  dev.off()
  
  path = paste0(simu_file_plots, "pecheperboat_i",i,"_x",x,".png")
  png(file = path, width = 1600, height = 500)
  plot(plot_pointsdepechecomperboat)
  dev.off()
  
  path = paste0(simu_file_plots, "pechecell_i",i,"_x",x,".png")
  png(file = path, width = 1600, height = 500)
  plot(plot_pointsdepechecell)
  dev.off()
}






generation_metriques = function(path_base,
                                simu_file,
                                couplePZ,
                                b_set){
  ### On convertit le couplePZ en vrai couple PZ
  
  if (couplePZ==1){
    sequencesdepeche=2
    zonespersequence=1
  } else if (couplePZ==2){
    sequencesdepeche=2
    zonespersequence=3
  } else if (couplePZ==3){
    sequencesdepeche=2
    zonespersequence=5
  } else if (couplePZ==4){
    sequencesdepeche=5
    zonespersequence=1
  } else if (couplePZ==5){
    sequencesdepeche=5
    zonespersequence=3
  } else if (couplePZ==6){
    sequencesdepeche=5
    zonespersequence=5
  } else if (couplePZ==7){
    sequencesdepeche=10
    zonespersequence=1
  } else if (couplePZ==8){
    sequencesdepeche=10
    zonespersequence=3
  } else if (couplePZ==9){
    sequencesdepeche=10
    zonespersequence=5
  }
  
  
  
  ### On recupère les données
  
  path = paste0(path_base, simu_file, "Results.RData")
  load(path)
  
  
  
  ### Tout d'abord, on construit le tableau de donnees de reference, celui pour x=10
  
  Results %>%
    filter(x == 10) -> Results_reference
  
  # Travail sur le tableau Results_reference (issu de plot_simu.R)
  Results_reference[,"OnBoundary"]=rep(0,nrow(Results_reference))
  Results_reference[which(abs(Results_reference[,"b_est"])>10),"OnBoundary"]=1
  Results_reference[,"Converge"]=Results_reference[,"OnBoundary"]+Results_reference[,"Convergence"]
  Results_reference[which(Results_reference[,"Converge"]==2),"Converge"]=1
  Converge_table = summaryBy(Converge~b_true+Data_source,data=Results_reference,FUN=sum)
  WhichFailed <- which(Results_reference[,"Convergence"]!=0)
  Results_cvg.failed <- Results_reference[WhichFailed,]
  WhichDone = which(Results_reference[,"Convergence"]==0)
  Results_reference=Results_reference[WhichDone,]
  Results_reference[,"RelBias_N"]=(Results_reference[,"N_est"]-Results_reference[,"N_true"])/Results_reference[,"N_true"]
  Results_reference[,"Bias_b"]=(Results_reference[,"b_est"]-Results_reference[,"b_true"]) / ifelse(Results_reference[,"b_true"] != 0,Results_reference[,"b_true"],1)
  Which_plot = which(Results_reference[,"type_b"] == "est_b")
  Results_reference_2 <- Results_reference
  Results_reference_2$b_est[which(is.na(Results_reference_2[,"b_est"]))] <- 0
  Results_reference_plot <- Results_reference_2[which(Results_reference_2[,"b_est"]<10),]
  Results_reference_plot$Data_source <- factor(Results_reference_plot$Data_source, levels = c("commercial_only"))
  # On ne s'intéresse qu'au cas sans réallocation pour la référence
  Results_reference_plot = Results_reference_plot %>% filter(reallocation == 0)
  
  
  
  ### Retour a notre x actuel
  
  # On va s'occuper des données qui concernent le couple PZ sur lequel on veut sortir les graphes
  Results %>%
    filter(x == couplePZ) -> Results
  

  ### Travail sur le tableau Results (issu de plot_simu.R)
  
  Results[,"OnBoundary"]=rep(0,nrow(Results))
  Results[which(abs(Results[,"b_est"])>10),"OnBoundary"]=1
  Results[,"Converge"]=Results[,"OnBoundary"]+Results[,"Convergence"]
  Results[which(Results[,"Converge"]==2),"Converge"]=1
  Converge_table = summaryBy(Converge~b_true+Data_source,data=Results,FUN=sum)
  WhichFailed <- which(Results[,"Convergence"]!=0)
  Results_cvg.failed <- Results[WhichFailed,]
  WhichDone = which(Results[,"Convergence"]==0)
  Results=Results[WhichDone,]
  Results[,"RelBias_N"]=(Results[,"N_est"]-Results[,"N_true"])/Results[,"N_true"]
  Results[,"Bias_b"]=(Results[,"b_est"]-Results[,"b_true"]) / ifelse(Results[,"b_true"] != 0,Results[,"b_true"],1)
  Which_plot = which(Results[,"type_b"] == "est_b")
  Results_2 <- Results
  Results_2$b_est[which(is.na(Results_2[,"b_est"]))] <- 0
  Results_plot <- Results_2[which(Results_2[,"b_est"]<10),]
  
  
  ### Plots (partiellement issus de plot_simu.R)

  Results_plot$Data_source <- factor(Results_plot$Data_source, levels = c("commercial_only"))
  
  Results_plot_k0 = Results_plot %>% filter(reallocation == 0)
  Results_plot_k1 = Results_plot %>% filter(reallocation == 1)
  
  
  ## Biais de l'abondance
  
  ymin = min(min(Results_plot_k0$RelBias_N), min(Results_plot_k1$RelBias_N), min(Results_reference_plot$RelBias_N))
  ymax = max(max(Results_plot_k0$RelBias_N), max(Results_plot_k1$RelBias_N), max(Results_reference_plot$RelBias_N))
  
  # Sans structuration
  plot_biaisabondance_ref <- ggplot() +
    geom_boxplot(data = Results_reference_plot,
                 aes(x = as.factor(b_true),
                     y = RelBias_N,
                     fill = Data_source)) +
    geom_hline(yintercept = 0,linetype="dashed") +
    theme_bw() +
    theme(axis.title = element_blank(), legend.position='none') +
    scale_fill_manual(breaks = c("commercial_only"),values=c("#E31D1D")) +
    facet_wrap(~"Sans structuration, sans reallocation") +
    ylim(ymin, ymax)
  
  # Sans reallocation
  plot_biaisabondance_k0 <- ggplot() +
    geom_boxplot(data = Results_plot_k0,
                 aes(x = as.factor(b_true),
                     y = RelBias_N,
                     fill = Data_source)) +
    geom_hline(yintercept = 0,linetype="dashed") +
    theme_bw() +
    theme(axis.title = element_blank(), legend.position='none') +
    scale_fill_manual(breaks = c("commercial_only"),values=c("#E31D1D")) +
    facet_wrap(~"Avec structuration, sans reallocation") +
    ylim(ymin, ymax)
  
  # Avec reallocation
  plot_biaisabondance_k1 <- ggplot() +
    geom_boxplot(data = Results_plot_k1,
                 aes(x = as.factor(b_true),
                     y = RelBias_N,
                     fill = Data_source)) +
    geom_hline(yintercept = 0,linetype="dashed") +
    theme_bw() +
    theme(axis.title = element_blank(), legend.position='none') +
    scale_fill_manual(breaks = c("commercial_only"),values=c("#E31D1D")) +
    facet_wrap(~"Avec structuration, avec reallocation") +
    ylim(ymin, ymax)
  
  # Graphe total
  plot_biaisabondance = ggarrange(plot_biaisabondance_ref,
                                  plot_biaisabondance_k0,
                                  plot_biaisabondance_k1,
                                  ncol=3)
  plot_biaisabondance = annotate_figure(plot_biaisabondance,
                                        top = text_grob(paste0("Biais de l'abondance ; structuration pour ", sequencesdepeche, " sequence(s) de peche et ", zonespersequence, " centre(s) par sequence en moyenne"), size = 12),
                                        left = text_grob("Biais de l'abondance",
                                                         rot = 90,
                                                         size = 11),
                                        bottom = text_grob("b", size = 11))
  
  
  ## Biais de b
  
  ymin = min(min(Results_plot_k0$Bias_b), min(Results_plot_k1$Bias_b), min(Results_reference_plot$Bias_b))
  ymax = max(max(Results_plot_k0$Bias_b), max(Results_plot_k1$Bias_b), max(Results_reference_plot$Bias_b))
  
  # Sans structuration
  plot_baisb_ref <- ggplot() +
    geom_boxplot(data = Results_reference_plot,
                 aes(x = as.factor(b_true),
                     y = Bias_b,
                     fill = Data_source)) +
    geom_hline(yintercept = 0,linetype="dashed") +
    theme_bw() +
    theme(axis.title = element_blank(),legend.position='none')+
    scale_fill_manual(breaks = c("commercial_only"), values=c("#F2690E")) +
    facet_wrap(~"Sans structuration, sans reallocation") +
    ylim(ymin, ymax)
  
  # Sans reallocation
  plot_baisb_k0 <- ggplot() +
    geom_boxplot(data = Results_plot_k0,
                 aes(x = as.factor(b_true),
                     y = Bias_b,
                     fill = Data_source)) +
    geom_hline(yintercept = 0,linetype="dashed") +
    theme_bw() +
    theme(axis.title = element_blank(),legend.position='none')+
    scale_fill_manual(breaks = c("commercial_only"), values=c("#F2690E")) +
    facet_wrap(~"Avec structuration, sans reallocation") +
    ylim(ymin, ymax)

  # Avec reallocation
  plot_baisb_k1 <- ggplot() +
    geom_boxplot(data = Results_plot_k1,
                 aes(x = as.factor(b_true),
                     y = Bias_b,
                     fill = Data_source)) +
    geom_hline(yintercept = 0,linetype="dashed") +
    theme_bw() +
    theme(axis.title = element_blank(),legend.position='none')+
    scale_fill_manual(breaks = c("commercial_only"), values=c("#F2690E")) +
    facet_wrap(~"Avec structuration, avec reallocation") +
    ylim(ymin, ymax)
  
  # Graphe total
  plot_baisb = ggarrange(plot_baisb_ref,
                         plot_baisb_k0,
                         plot_baisb_k1,
                         ncol=3)
  plot_baisb = annotate_figure(plot_baisb,
                               top = text_grob(paste0("Biais de b ; structuration pour ", sequencesdepeche, " sequence(s) de peche et ", zonespersequence, " centre(s) par sequence en moyenne"), size = 12),
                               left = text_grob("Biais de b",
                                                rot = 90,
                                                size = 11),
                               bottom = text_grob("b", size = 11))
  
  
  ## MSPE
  
  ymin = min(min(Results_plot_k0$MSPE_S), min(Results_plot_k1$MSPE_S), min(Results_reference_plot$MSPE_S))
  ymax = max(max(Results_plot_k0$MSPE_S), max(Results_plot_k1$MSPE_S), max(Results_reference_plot$MSPE_S))
  
  # Sans structuration
  plot_MSPE_ref <- ggplot() +
    geom_boxplot(data = Results_reference_plot,
                 aes(x = as.factor(b_true),
                     y = MSPE_S,
                     fill = Data_source)) +
    geom_hline(yintercept = 0,linetype="dashed") +
    theme_bw() +
    theme(axis.title = element_blank(),legend.position='none')+
    scale_fill_manual(breaks = c("commercial_only"), values=c("#FFA888")) +
    facet_wrap(~"Sans structuration, sans reallocation") +
    ylim(ymin, ymax)
  
  # Sans reallocation
  plot_MSPE_k0 <- ggplot() +
    geom_boxplot(data = Results_plot_k0,
                 aes(x = as.factor(b_true),
                     y = MSPE_S,
                     fill = Data_source)) +
    geom_hline(yintercept = 0,linetype="dashed") +
    theme_bw() +
    theme(axis.title = element_blank(),legend.position='none')+
    scale_fill_manual(breaks = c("commercial_only"), values=c("#FFA888")) +
    facet_wrap(~"Avec structuration, sans reallocation") +
    ylim(ymin, ymax)

  # Avec reallocation
  plot_MSPE_k1 <- ggplot() +
    geom_boxplot(data = Results_plot_k1,
                 aes(x = as.factor(b_true),
                     y = MSPE_S,
                     fill = Data_source)) +
    geom_hline(yintercept = 0,linetype="dashed") +
    theme_bw() +
    theme(axis.title = element_blank(),legend.position='none')+
    scale_fill_manual(breaks = c("commercial_only"), values=c("#FFA888")) +
    facet_wrap(~"Avec structuration, avec reallocation") +
    ylim(ymin, ymax)
  
  # Graphe total
  plot_MSPE = ggarrange(plot_MSPE_ref,
                        plot_MSPE_k0,
                        plot_MSPE_k1,
                        ncol=3)
  plot_MSPE = annotate_figure(plot_MSPE,
                              top = text_grob(paste0("MSPE (Mean Square Prediction Error) ; structuration pour ", sequencesdepeche, " sequence(s) de peche et ", zonespersequence, " centre(s) par sequence en moyenne"), size = 12),
                              left = text_grob("MSPE",
                                               rot = 90,
                                               size = 11),
                              bottom = text_grob("b", size = 11))
  
  
  ### Sauvegarde des plots dans results_metriques
  
  simu_file_metriques <- paste0("results_metriques/", simu_file)
  if(! dir.exists(simu_file_metriques)) dir.create(simu_file_metriques)
  
  # Biais de l'abondance
  path = paste0(simu_file_metriques, "biaisabondance_x",couplePZ,".png")
  png(file = path, width = 1000, height = 500)
  plot(plot_biaisabondance)
  dev.off()
  
  # Biais de b
  path = paste0(simu_file_metriques, "biaisb_x",couplePZ,".png")
  png(file = path, width = 1000, height = 500)
  plot(plot_baisb)
  dev.off()
  
  # MSPE
  path = paste0(simu_file_metriques, "MSPE_x",couplePZ,".png")
  png(file = path, width = 1000, height = 500)
  plot(plot_MSPE)
  dev.off()
}






indicateurs = function(path_base,
                       simu_file,
                       table,
                       b,
                       k,
                       biaisb,
                       biaisN,
                       MSPE){
  
  ### On recupère les données
  
  path = paste0(path_base, simu_file, "Results.RData")
  load(path)
  
  # On filtre pour avoir seulement les données du b et du k voulus
  
  Results %>%
    filter(b_true == b) -> Results
  if (k == reallocation[1]){
    Results = Results %>% filter(reallocation == 0)
  } else {
    Results = Results %>% filter(reallocation == 1)
  }
  
  ### Travail sur le tableau Results (issu de plot_simu.R)
  
  Results[,"OnBoundary"]=rep(0,nrow(Results))
  Results[which(abs(Results[,"b_est"])>10),"OnBoundary"]=1
  Results[,"Converge"]=Results[,"OnBoundary"]+Results[,"Convergence"]
  Results[which(Results[,"Converge"]==2),"Converge"]=1
  Converge_table = summaryBy(Converge~b_true+Data_source,data=Results,FUN=sum)
  WhichFailed <- which(Results[,"Convergence"]!=0)
  Results_cvg.failed <- Results[WhichFailed,]
  WhichDone = which(Results[,"Convergence"]==0)
  Results=Results[WhichDone,]
  Results[,"RelBias_N"]=(Results[,"N_est"]-Results[,"N_true"])/Results[,"N_true"]
  Results[,"Bias_b"]=(Results[,"b_est"]-Results[,"b_true"]) / ifelse(Results[,"b_true"] != 0,Results[,"b_true"],1)
  Which_plot = which(Results[,"type_b"] == "est_b")
  Results$b_est[which(is.na(Results[,"b_est"]))] <- 0
  Results <- Results[which(Results[,"b_est"]<10),]
  Results$Data_source <- factor(Results$Data_source, levels = c("commercial_only"))
  
  
  ### Calcul des indicateurs, remplissage des tableaux
  
  Results_x1 = Results %>% filter(x == 1)
  biaisb[[table]][1, 1] = median(Results_x1$Bias_b)
  biaisN[[table]][1, 1] = median(Results_x1$RelBias_N)
  MSPE[[table]][1, 1] = median(Results_x1$MSPE_S)
  
  Results_x2 = Results %>% filter(x == 2)
  biaisb[[table]][2, 1] = median(Results_x2$Bias_b)
  biaisN[[table]][2, 1] = median(Results_x2$RelBias_N)
  MSPE[[table]][2, 1] = median(Results_x2$MSPE_S)
  
  Results_x3 = Results %>% filter(x == 3)
  biaisb[[table]][3, 1] = median(Results_x3$Bias_b)
  biaisN[[table]][3, 1] = median(Results_x3$RelBias_N)
  MSPE[[table]][3, 1] = median(Results_x3$MSPE_S)
  
  Results_x4 = Results %>% filter(x == 4)
  biaisb[[table]][1, 2] = median(Results_x4$Bias_b)
  biaisN[[table]][1, 2] = median(Results_x4$RelBias_N)
  MSPE[[table]][1, 2] = median(Results_x4$MSPE_S)
  
  Results_x5 = Results %>% filter(x == 5)
  biaisb[[table]][2, 2] = median(Results_x5$Bias_b)
  biaisN[[table]][2, 2] = median(Results_x5$RelBias_N)
  MSPE[[table]][2, 2] = median(Results_x5$MSPE_S)
  
  Results_x6 = Results %>% filter(x == 6)
  biaisb[[table]][3, 2] = median(Results_x6$Bias_b)
  biaisN[[table]][3, 2] = median(Results_x6$RelBias_N)
  MSPE[[table]][3, 2] = median(Results_x6$MSPE_S)
  
  Results_x7 = Results %>% filter(x == 7)
  biaisb[[table]][1, 3] = median(Results_x7$Bias_b)
  biaisN[[table]][1, 3] = median(Results_x7$RelBias_N)
  MSPE[[table]][1, 3] = median(Results_x7$MSPE_S)
  
  Results_x8 = Results %>% filter(x == 8)
  biaisb[[table]][2, 3] = median(Results_x8$Bias_b)
  biaisN[[table]][2, 3] = median(Results_x8$RelBias_N)
  MSPE[[table]][2, 3] = median(Results_x8$MSPE_S)
  
  Results_x9 = Results %>% filter(x == 9)
  biaisb[[table]][3, 3] = median(Results_x9$Bias_b)
  biaisN[[table]][3, 3] = median(Results_x9$RelBias_N)
  MSPE[[table]][3, 3] = median(Results_x9$MSPE_S)
  
  res = list(biaisb, biaisN, MSPE)
  return(res)
}
