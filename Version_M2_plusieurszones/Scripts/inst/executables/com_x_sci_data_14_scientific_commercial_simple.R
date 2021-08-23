
      library( TMB )
      dyn.load( dynlib('Scripts/inst/executables/com_x_sci_data_14_scientific_commercial_simple') )
      setwd('C:/R_projects/projet_M2_sole/Version_M2_plusieurszones')
      load( 'All_inputs.RData' )
      Obj = MakeADFun(data=All_inputs[['data']], parameters=All_inputs[['parameters']], random=All_inputs[['random']], All_inputs[['Other_inputs']])
      save( Obj, file='Obj.RData')
    