
ggplot(loc_x)+
  geom_histogram(aes(x=bathy,fill=factor(substr_Mud_sediment)))


fit_IM_res_no_realloc$Report$jnll_comp
fit_IM_res_realloc$Report$jnll_comp
fit_IM_res_sci$Report$jnll_comp


## Integrated model
obj_int <- fit_IM_res_realloc$Obj
SD_int <- fit_IM_res_realloc$SD
obj_int$fn()
pl_int <- as.list(SD_int,"Est")

## Scientific model
obj_sci <- fit_IM_res_sci$Obj
SD_sci <- fit_IM_res_sci$SD
obj_sci$fn()
pl_sci <- as.list(SD_sci,"Est")
par_sci <- SD_sci$par.fixed
f_sci <- as.numeric(obj_sci$fn(par_sci))

## Fixed
par_int <- SD_int$par.fixed # extractBaseParameters(obj_sci, pl_sci, pl_int)
par_int <- par_int[which(names(par_int) %in% names(par_sci))]

if(F %in% ("logtau" %in% names(par_int))){
  
  par_int <- c(par_int,-2.424969)
  names(par_int)[length(par_int)] <- "logtau"
  par_order <- which(names(par_sci) == "logtau")
  par_int_2 <- c(par_int[1:(par_order-1)],par_int[length(par_int)],par_int[(par_order):(length(par_int)-1)])
  par_int <- par_int_2
  
}

par_int[1] <- par_sci[1]

f_int <- as.numeric(obj_sci$fn(par_int))
def_free <- sum(par_sci != 0)
fixed <- 1 - pchisq( 2 * (f_int-f_sci), df=def_free )

## Random
# par_int.all <- extractBaseParameters(obj_sci, pl_sci, pl_int, all=TRUE)
par_sci.all <- obj_sci$env$last.par.best
par_int.all <- obj_int$env$last.par.best
par_int.all <- par_int.all[which(names(par_int.all) %in% names(par_sci.all))]

if(F %in% ("logtau" %in% names(par_int.all))){
  
  par_int.all <- c(par_int.all,-2.424969)
  names(par_int.all)[length(par_int.all)] <- "logtau"
  par_order <- which(names(par_sci.all) == "logtau")
  par_int.all_2 <- c(par_int.all[1:(par_order-1)],
                     par_int.all[length(par_int.all)],
                     par_int.all[(par_order):(length(par_int.all)-1)])
  par_int.all <- par_int.all_2
  
}

f_sci.all <- obj_sci$env$f(par_sci.all) # Best evaluated parameters
f_int.all <- obj_sci$env$f(par_int.all)
def_free.all <- def_free + length(obj_sci$env$random)
random <- 1 - pchisq( 2 * (f_int.all-f_sci.all), df=def_free.all )
