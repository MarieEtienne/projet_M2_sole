

## Scientific model
load("C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/case_study/sci_df.RData")
prof_sci <- tmbprofile(obj = fit_IM_res$Obj,lincomb = c(0,1,0,0,0,0))
# x11();plot(prof_sci)

prof_sci_2 <- prof_sci %>%
  mutate(Model = "Scientific model") %>%
  mutate(value_2 = (value - mean(value)) / sd(value) ) %>%
  mutate(lkl = "Yi")

prof_sci_2$value_2 <- prof_sci_2$value - min(prof_sci_2$value)

## No realloc
load("C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/case_study/no.realloc_int_df.RData")
prof_no.realloc <- tmbprofile(obj = fit_IM_res$Obj,lincomb = c(0,1,0,0,0,0,0,0,0))
# x11();plot(prof_no.realloc)

prof_no.realloc_2 <- prof_no.realloc %>%
  mutate(Model = "Integrated model") %>%
  mutate(value_2 = (value - mean(value)) / sd(value) ) %>%
  mutate(lkl = "Yi")

prof_no.realloc_2$value_2 <- prof_no.realloc_2$value - min(prof_no.realloc_2$value)

## Reallocation and re-estimation
load("C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/case_study/realloc_int_df_rest.RData")
prof_realloc <- tmbprofile(obj = fit_IM_res$Obj,lincomb = c(0,1,0,0,0,0,0,0))
# x11();plot(prof_no.realloc)

prof_realloc_2 <- prof_realloc %>%
  mutate(Model = "Integrated model") %>%
  mutate(value_2 = (value - min(value))) %>%
  mutate(lkl = "Dj")

prof_realloc_2$value_2 <- prof_realloc_2$value - min(prof_realloc_2$value)

## Test
prof_test <- tmbprofile(obj = fit_IM_res$Obj,lincomb = c(0,1,0,0,0,0,0,0,0))
# x11();plot(prof_no.realloc)

prof_test_2 <- prof_test %>%
  mutate(Model = "Test") %>%
  mutate(value_2 = (value - min(value))) %>%
  mutate(lkl = "Yi")

prof_test_2$value_2 <- prof_test_2$value - min(prof_test_2$value)


lkl_prof_df <- rbind(prof_sci_2,
                     prof_test_2,
                     prof_no.realloc_2,
                     prof_realloc_2) %>%
  mutate(Model_lkl = paste0(Model," - ",lkl))


x11()
ggplot(lkl_prof_df)+
  geom_line(aes(x=parameter,y=-value_2,col=Model_lkl),size=1)+
  theme()+
  theme_bw()+
  ylab("Rescaled log-likelihood values")+
  xlab("Parameter value (beta)")+
  geom_vline(xintercept = 2,col="darkgrey",size=1,linetype="dashed")

##----------------------------------------------------------------------------------------

test <- data.frame(loc_x,S_x=fit_IM_res$Report$S_x)
x11()
ggplot(test)+
  geom_point(aes(x=x,y=y,col=S_x))+
  scale_color_distiller(palette="Spectral")

