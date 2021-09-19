// Space time 
#include <TMB.hpp>
//#include <atomic_math.hpp>  //for D_lgamma


template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

template<class Type>
Type dbern(Type x, Type prob, int give_log=1){
  Type logres;
  if( x==0 ) logres = log( 1-prob );
  if( x==1 ) logres = log( prob );
  if(give_log) return logres; else return exp(logres);
}


// dlognorm

template<class Type>

Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=false){
  Type Return;
  if(give_log==false) Return = dnorm( log(x), meanlog, sdlog, false) / x;
  if(give_log==true) Return = dnorm( log(x), meanlog, sdlog, true) - log(x);
  return Return;
}

// square
template<class Type>
Type square(Type x){
  return pow(x,2.0);
}

// square-root
template<class Type>
Type sqrt(Type x){
  return pow(x,0.5);
}

// dzinfgamma, shape = 1/CV^2, scale = mean*CV^2
template<class Type>
Type dzinfgamma(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dgamma( x, pow(cv,-2), posmean*pow(cv,2), false );
    if(give_log==true) Return = log(encounter_prob) + dgamma( x, pow(cv,-2), posmean*pow(cv,2), true );
  } 
  return Return;
}


// dzinflognorm

template<class Type>

Type dzinflognorm(Type x, Type meanlog, Type encounter_prob, Type log_notencounter_prob, Type sdlog, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dlognorm( x, meanlog - square(sdlog)/2, sdlog, false ); // 
    if(give_log==true) Return = log(encounter_prob) + dlognorm( x, meanlog - square(sdlog)/2, sdlog, true );
  } 
  return Return;
}



template<class Type>
Type objective_function<Type>::operator() (){
  using namespace Eigen;
  using namespace density;
  
  ///////////////////////////////////////////////////////////////////
  // Inputs definition 
  /////////////////////
  
//// A changer 
  
  // options vec
  DATA_FACTOR( Options_vec );
  // Slot 0: prior on random effects (0=SPDE_GMRF; 1=ICAR_GMRF)
  // Slot 4: Output S_x in ADREPORT (0=no, 1=yes); needed for unbiased map estimates using bias.correct
  // Slot 5: Estimation model (1 : scientific_commercial, 1 : scientific_only, 3 : commercial_only)
  
  
  //////////////////////////////////// Global /////////////////////////////////////////
  
  ////////
  // Data
  ////////
  
  DATA_MATRIX( Cov_xj );  

  //////////////
  // Parameters 
  //////////////

  PARAMETER_VECTOR(beta_j);

  /////////////////
  // derived values
  /////////////////
  
  int n_j = Cov_xj.row(0).size();
  // int n_k = Cov_xk.row(0).size();
  int n_x = Cov_xj.col(0).size();

  // global stuff
  vector<Type> jnll_comp(5);
  jnll_comp.setZero();


  ////////////////
  // Latent field
  ///////////////

  // Predicted densities
  vector<Type> S_x(n_x);S_x.setZero();
  vector<Type> debug_j(n_j);debug_j.setZero();
  vector<Type> linpredS_x = Cov_xj * beta_j;
  for(int s=0; s<n_x; s++){
    S_x(s) = exp(linpredS_x(s) );
  }

  
  //////////////////////////////////// Commercial /////////////////////////////////////////
  
  
  if(Options_vec(5) == 3 | Options_vec(5) == 1){

    ////////
    // Data
    ////////

    DATA_VECTOR( c_com_x );       	// Response (count) for each observation i (commercial data)
    DATA_VECTOR( y_com_i );       	// Response (0:not surveyed, 1:surveyed) for each site (commercial data)
    DATA_FACTOR( seq_com_i );       	// Response (0:not surveyed, 1:surveyed) for each site (commercial data)
    
    DATA_FACTOR( index_com_i );
   	DATA_VECTOR( q2_com );
    DATA_MATRIX( Cov_xk );  //design matrix for sampling of commercial data
    

    //////////////
    // Parameters
    //////////////

    PARAMETER_VECTOR(q1_com);
    PARAMETER_VECTOR(logSigma_com);
	  PARAMETER_VECTOR( k_com );

	  /////////////////
    // derived values
    /////////////////

    int n_com_i = y_com_i.size();
    int n_ping_i = seq_com_i.size();
    

    /////////////////////
    // Observation model
    ////////////////////
    
    if( Options_vec(9) == 1){
      
      vector<Type> Sigma_com = exp(logSigma_com);
      vector<Type> E_com(n_ping_i);E_com.setZero();
      
      vector<Type> encounterprob_com(n_com_i);encounterprob_com.setZero();
      vector<Type> log_notencounterprob_com(n_com_i);log_notencounterprob_com.setZero();
      
      ///////////////////////////////////
      // Non-aggregated observation model
      //////////////////////////////////
      if( Options_vec(12) == 0 ){
        for(int i=0; i<n_com_i; i++){
          if( !isNA(y_com_i(i)) ){
            
            E_com(i) = q2_com(0) * 1 * S_x(index_com_i(i));
            
            // proba of encountering
            encounterprob_com(i) = ( 1.0 - exp(-1 * E_com(i) * exp(q1_com(0))) );
            log_notencounterprob_com(i) = -1 * E_com(i) * exp(q1_com(0));
            
            jnll_comp(0) -= dzinflognorm(y_com_i(i), log(E_com(i))-log(encounterprob_com(i)), encounterprob_com(i), log_notencounterprob_com(i), Sigma_com(0), true);
            // jnll_comp(0) -= dzinfgamma(y_com_i(i), E_com(i)/encounterprob_com(i), encounterprob_com(i), log_notencounterprob_com(i), Sigma_com(0), true);
            
          }
        }
      }
      
      ///////////////////////////////
      // Aggregated observation model
      ///////////////////////////////
      vector<Type> E_D_com(n_com_i);E_D_com.setZero();
      vector<Type> log_pi_j(n_com_i);log_pi_j.setZero();
      vector<Type> pi_j(n_com_i);pi_j.setZero();
      vector<Type> sum_mu(n_com_i);sum_mu.setZero();
      vector<Type> Var_D_com(n_com_i);Var_D_com.setZero();
      vector<Type> Var_Dsup0_com(n_com_i);Var_Dsup0_com.setZero();
      vector<Type> p_i(n_ping_i);p_i.setZero();
      Type Sigma_D;
      Type log_Sigma_D;
      
      if( Options_vec(12) == 1 ){
        
        for(int i=0; i<n_ping_i; i++){
          
          E_com(i) = q2_com(0) * 1 * S_x(index_com_i(i));
          log_pi_j(seq_com_i(i)) += -1 * E_com(i) * exp(q1_com(0));
          
          if(y_com_i(seq_com_i(i)) > 0){
            
            p_i(i) = exp(-1 * E_com(i) * exp(q1_com(0)));
            sum_mu(seq_com_i(i)) += q2_com(0) * 1 * S_x(index_com_i(i));
            Var_D_com(seq_com_i(i)) += square(E_com(i)) / (1 - p_i(i)) * (exp(square(Sigma_com(0))) - (1 - p_i(i)));
          
          }
          
        }
        
        pi_j = exp(log_pi_j);
        
        for(int j=0; j<n_com_i; j++){
          
          if( y_com_i(j) == 0 ){
            
            jnll_comp(0) -= log_pi_j(j);
            
          }
          
          if( y_com_i(j) > 0 ){
            
            E_D_com(j) = 1 / (1 - pi_j(j)) * sum_mu(j);
            
            Var_Dsup0_com(j) =  Var_D_com(j) / (1 -  pi_j(j)) - pi_j(j) * square(E_D_com(j)) / (1 - pi_j(j)) ;
            
            Sigma_D = sqrt(log(Var_Dsup0_com(j) / square(E_D_com(j)) + 1));
            
            // log_Sigma_D = log(Sigma_D);
            
            jnll_comp(0) -= dlognorm(y_com_i(j),
                      log(E_D_com(j)),
                      Sigma_D,
                      true);
          
          }
          
          
        }
        
      }
      
      REPORT(E_com);
      REPORT(E_D_com);
      REPORT(log_pi_j);
      REPORT(pi_j);
      REPORT(sum_mu);
      REPORT(Var_D_com);
      REPORT(Var_Dsup0_com);
      REPORT(p_i);
      REPORT(Sigma_D);
      REPORT(log_Sigma_D);
      REPORT( encounterprob_com);
      
      REPORT( q1_com );
      REPORT( Sigma_com );
      
    }
    
    //////////////////////////////////////
    // Sampling process of commercial data
    /////////////////////////////////////
    
    if( Options_vec(7) == 1){
      
      // Parameter
      PARAMETER_VECTOR(beta_k);
      PARAMETER_VECTOR( par_b );

      
      // Derived values
      Type b;
      vector<Type> linpredR_x(n_x);linpredR_x.setZero();
      vector<Type> lambda_x(n_x);lambda_x.setZero();
      // matrix<Type> lambda_ref_x(n_x , n_eta );
      vector<Type> fact_S(n_x);fact_S.setZero();
      
      
      // Poisson point process
      for(int s=0; s<n_x; s++){
        linpredR_x = Cov_xk * beta_k;
        fact_S = par_b(0)*(log(S_x(s)));
        lambda_x(s) = exp( linpredR_x(s) + fact_S(s));
        if(!isNA(c_com_x(s)))   jnll_comp(1)-= (Type(1)-lambda_x(s) + c_com_x(s)*log(lambda_x(s))); // log density of poisson process (see Diggle (2013))
      }

      // Outputs
      REPORT( lambda_x );
      REPORT( beta_k );
      REPORT( par_b );
      REPORT( linpredR_x );
      
      ADREPORT(lambda_x);
    }
  }

  //////////////////////////////////// Scientific /////////////////////////////////////////
  
  if(Options_vec(5) != 3){
    
    ////////
    // Data
    ////////
    
    DATA_VECTOR( y_sci_i );       	// Response (0:not surveyed, 1:surveyed) for each site (scientific data)
    DATA_FACTOR( index_sci_i );
    DATA_VECTOR( q2_sci );
    
    //////////////
    // Parameters 
    //////////////
    
    PARAMETER(logSigma_sci);
    PARAMETER_VECTOR(q1_sci);
   	PARAMETER_VECTOR( k_sci );

    /////////////////
    // derived values
    /////////////////
    
    Type Sigma_sci = exp(logSigma_sci);
    int n_sci_i = y_sci_i.size();

    /////////////////////
    // Observation model
    ////////////////////
    
    vector<Type> E_sci(n_sci_i);E_sci.setZero();
    vector<Type> encounterprob_sci(n_sci_i);encounterprob_sci.setZero();
    vector<Type> log_notencounterprob_sci(n_sci_i);log_notencounterprob_sci.setZero();
    
    // Zero-inflated  distribution
    for(int i=0; i<n_sci_i; i++){
      if( !isNA(y_sci_i(i)) ){
        E_sci(i) = k_sci(0) * q2_sci(0) * S_x(index_sci_i(i));
        encounterprob_sci(i) = ( 1.0 - exp(-1 * E_sci(i) * exp(q1_sci(0)) ));
        log_notencounterprob_sci(i) = -1 * E_sci(i) * exp(q1_sci(0));

        jnll_comp(0) -= dzinflognorm(y_sci_i(i), log(E_sci(i))-log(encounterprob_sci(i)), encounterprob_sci(i), log_notencounterprob_sci(i), Sigma_sci, true);
        // if(y_sci_i(i) > 0) jnll_comp(0) -= dlognorm(y_sci_i(i), log(E_sci(i)) - square(Sigma_sci)/2, Sigma_sci, true);
        // jnll_comp(0) -= dzinfgamma(y_sci_i(i), E_sci(i)/encounterprob_sci(i), encounterprob_sci(i), log_notencounterprob_sci(i), Sigma_sci, true);
          
      }
    }

    //////////
    // Outputs
    //////////
    
    REPORT( Sigma_sci );
    REPORT( q1_sci );
    
  }

  Type total_abundance = S_x.sum();

  // Total objective
  // parallel_accumulator<Type> jnll(this);
  Type jnll;
  
  int n_jnll = jnll_comp.size();
  for(int i=0; i<n_jnll; i++){
    jnll += jnll_comp(i);
  }
  
  
  // Reporting
  REPORT( S_x );

  REPORT( total_abundance );
  REPORT( beta_j );
  REPORT( linpredS_x );
  
  // REPORT( Unsampled_x );
  REPORT( jnll_comp );
  // REPORT( jnll );

  // Bias correction output
  //ADREPORT( beta_abs_j );
  //ADREPORT( beta_pos_j );
  ADREPORT( total_abundance );

  if(Options_vec(4)==1){
    ADREPORT(S_x);
  }
  
  
  return jnll;
}

