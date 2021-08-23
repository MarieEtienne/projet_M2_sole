
test <- Obj$report()

list_test <- list(y_com_i = Data$y_com_i,
                  E_com=test$E_com,
                  log_pi_j=test$log_pi_j,
                  pi_j=test$pi_j,
                  sum_mu=test$sum_mu,
                  Var_D_com=test$Var_D_com,
                  p_i=test$p_i,
                  
                  Var_Dsup0_com=test$Var_Dsup0_com,
                  E_D_com=test$E_D_com,
                  encounterprob_com=test$encounterprob_com,
                  Sigma_D=test$Sigma_D,
                  log_Sigma_D=test$log_Sigma_D)

data.frame(index_com_i = index_com_i,
           E_com = list_test$E_com,
           p_i = list_test$p_i,
           sum_mu = list_test$sum_mu,
           Var_D_com = list_test$Var_D_com,
           
           E_D_com = list_test$E_D_com,
           pi_j = list_test$pi_j,
           log_pi_j = list_test$log_pi_j,
           
           y_com_i = list_test$y_com_i
           
           )

# Var_D_com
test$p_i / (1 - test$p_i)^2 * test$E_com^2 * ((exp(exp(logSigma_com)^2)) - test$p_i)


## Diagnose model
plot(Strue_x,Report$S_x)

loc_x_2 <- loc_x %>%
        mutate(Strue_x = Strue_x,
               obs_ref = fit_IM_res_ref$Report$S_x,
               obs_non_aggreg = fit_IM_res_non_aggreg$Report$S_x,
               obs_aggreg = fit_IM_res_aggreg$Report$S_x
               )

plot_1 <- ggplot(loc_x_2)+
        geom_point(aes(x=x,y=y,col=Strue_x),size=4,shape=15)+
        scale_color_distiller(palette = "Spectral")+
        theme(aspect.ratio = 1,
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5,size = 12,face = "bold"))+
        ggtitle("Simulated latent field")+
        xlab("")+ylab("")

plot_2 <- ggplot(loc_x_2)+
        geom_point(aes(x=x,y=y,col=obs_ref),size=4,shape=15)+
        scale_color_distiller(palette = "Spectral")+
        theme(aspect.ratio = 1,
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5,size = 12,face = "bold"))+
        ggtitle("Observations: exact values")+
        xlab("")+ylab("")

plot_3 <- ggplot(loc_x_2)+
        geom_point(aes(x=x,y=y,col=obs_non_aggreg),size=4,shape=15)+
        scale_color_distiller(palette = "Spectral")+
        theme(aspect.ratio = 1,
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5,size = 12,face = "bold"))+
        ggtitle("Observations: reallocated values \n \n Model: does not account \n for reallocation")+
        xlab("")+ylab("")


plot_4 <- ggplot(loc_x_2)+
        geom_point(aes(x=x,y=y,col=obs_aggreg),size=4,shape=15)+
        scale_color_distiller(palette = "Spectral")+
        theme(aspect.ratio = 1,
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5,size = 12,face = "bold")
              )+
        ggtitle("Observations: reallocated values \n \n Model: accounts \n for reallocation")+
        xlab("")+ylab("")

plot_grid(plot_1,plot_2,plot_3,plot_4,nrow = 2,align = "hv")
