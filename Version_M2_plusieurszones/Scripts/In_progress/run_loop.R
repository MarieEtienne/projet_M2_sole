i0=19

for(i in i0:i1){
  
  for(aggreg_obs in aggreg_obs_moda){
    
    
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch({
      
      taillezone = 3
      
      # if (x!=10){
      #   taillezone = 3
      # }else{
      #   taillezone = 25
      # }
      
      
      res <- simu_commercial_scientific(Results,
                                        simu_file,
                                        run_datarmor,
                                        data.res_folder,
                                        r_folder,
                                        grid_dim,
                                        n_cells,
                                        beta0,
                                        beta,
                                        range,
                                        nu,
                                        SD_x,
                                        SD_delta,
                                        SD_eta,
                                        n_samp_sci,
                                        logSigma_sci,
                                        q1_sci,
                                        q2_sci,
                                        n_strate,
                                        n_samp_com,
                                        logSigma_com,
                                        q1_com,
                                        q2_com,
                                        b,
                                        Data_source,
                                        Samp_process,
                                        EM,
                                        aggreg_obs,
                                        RandomSeed,
                                        Version,
                                        TmbFile,
                                        ignore.uncertainty,
                                        counter,
                                        i,
                                        n_sim,
                                        k,
                                        xlim,
                                        ylim,
                                        sequencesdepeche,
                                        zonespersequence,
                                        taillezone,
                                        i0,
                                        i1,
                                        n_fact,
                                        cluster_nb,
                                        n_nodes)
      
      Results <- res[[1]]
      List_param <- res[[2]]
      counter <- res[[3]]
      
      
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next } 
    
        
  }
  
}

