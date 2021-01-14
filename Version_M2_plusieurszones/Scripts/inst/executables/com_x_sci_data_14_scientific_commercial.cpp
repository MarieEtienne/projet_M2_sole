// Space time 
#include <TMB.hpp>
//#include <atomic_math.hpp>  //for D_lgamma


/** Precision matrix for the anisotropic case, eqn (20) in Lindgren et al. (2011) */    
namespace R_inla_generalized {
using namespace Eigen;
using namespace tmbutils;
using namespace R_inla;

template<class Type>
SparseMatrix<Type> Q_spde_generalized(spde_t<Type> spde, Type kappa, int alpha=2){
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  
  if( alpha==1 ) return kappa_pow2*spde.M0 + spde.M1;
  if( alpha==2 ) return kappa_pow4*spde.M0 + Type(2.0)*kappa_pow2*spde.M1 + spde.M2;
  else EXIT_SUCCESS;
}

template<class Type>
SparseMatrix<Type> Q_spde_generalized(spde_aniso_t<Type> spde, Type kappa, matrix<Type> H, int alpha=2){
  
  int i;
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  
  int n_s = spde.n_s;
  int n_tri = spde.n_tri;
  vector<Type> Tri_Area = spde.Tri_Area;
  matrix<Type> E0 = spde.E0;
  matrix<Type> E1 = spde.E1;
  matrix<Type> E2 = spde.E2;
  matrix<int> TV = spde.TV;
  SparseMatrix<Type> G0 = spde.G0;
  SparseMatrix<Type> G0_inv = spde.G0_inv;
  
  //Type H_trace = H(0,0)+H(1,1);
  //Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);
  SparseMatrix<Type> G1_aniso(n_s,n_s); 
  SparseMatrix<Type> G2_aniso(n_s,n_s); 
  // Calculate adjugate of H
  matrix<Type> adj_H(2,2);
  adj_H(0,0) = H(1,1);
  adj_H(0,1) = -1 * H(0,1);
  adj_H(1,0) = -1 * H(1,0);
  adj_H(1,1) = H(0,0);
  // Calculate new SPDE matrices
  
  // Calculate G1 - pt. 1
  array<Type> Gtmp(n_tri,3,3);
  for(i=0; i<n_tri; i++){    
    // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.    
    Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
  }
  // Calculate G1 - pt. 2
  for(i=0; i<n_tri; i++){
    G1_aniso.coeffRef(TV(i,1),TV(i,0)) = G1_aniso.coeffRef(TV(i,1),TV(i,0)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,0),TV(i,1)) = G1_aniso.coeffRef(TV(i,0),TV(i,1)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,1)) = G1_aniso.coeffRef(TV(i,2),TV(i,1)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,1),TV(i,2)) = G1_aniso.coeffRef(TV(i,1),TV(i,2)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,2),TV(i,0)) = G1_aniso.coeffRef(TV(i,2),TV(i,0)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,2)) = G1_aniso.coeffRef(TV(i,0),TV(i,2)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,0)) = G1_aniso.coeffRef(TV(i,0),TV(i,0)) + (Gtmp(i,0,0));  
    G1_aniso.coeffRef(TV(i,1),TV(i,1)) = G1_aniso.coeffRef(TV(i,1),TV(i,1)) + (Gtmp(i,1,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,2)) = G1_aniso.coeffRef(TV(i,2),TV(i,2)) + (Gtmp(i,2,2));  
  }
  G2_aniso = G1_aniso * G0_inv * G1_aniso; 
  
  if( alpha==1 ) return kappa_pow2*G0 + G1_aniso;
  if( alpha==2 ) return kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1_aniso + G2_aniso;
}
} // end namespace R_inla


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
    if(give_log==false) Return = encounter_prob * dlognorm( x, meanlog - square(sdlog)/2, sdlog, false );
    if(give_log==true) Return = log(encounter_prob) + dlognorm( x, meanlog - square(sdlog)/2, sdlog, true );
  } 
  return Return;
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



// dzinfnorm

template<class Type>

Type dzinfnorm(Type x, Type posmean, Type encounter_prob, Type log_notencounter_prob, Type cv, int give_log=false){
  Type Return;
  if(x==0){
    if(give_log==false) Return = 1.0 - encounter_prob;
    if(give_log==true){
      if( isNA(log_notencounter_prob) ) Return = log(1.0 - encounter_prob);
      if( !isNA(log_notencounter_prob) ) Return = log_notencounter_prob;
    }
  }else{
    if(give_log==false) Return = encounter_prob * dnorm( x, posmean, posmean*cv, false );
    if(give_log==true) Return = log(encounter_prob) + dnorm( x, posmean, posmean*cv, true );
  } 
  return Return;
}

template<class Type>
Type objective_function<Type>::operator() (){
  using namespace R_inla;
  using namespace R_inla_generalized;
  using namespace Eigen;
  using namespace density;
  
  ///////////////////////////////////////////////////////////////////
  // Inputs definition 
  /////////////////////
  
  
  // options vec
  DATA_FACTOR( Options_vec );
  // Slot 0: prior on random effects (0=SPDE_GMRF; 1=ICAR_GMRF)
  // Slot 1: Alpha
  // Slot 2: Include hyperdistribution for delta
  // Slot 3: Include hyperdistribution for eta
  // Slot 4: Output S_x in ADREPORT (0=no, 1=yes); needed for unbiased map estimates using bias.correct
  // Slot 5: Estimation model (1 : scientific_commercial, 1 : scientific_only, 3 : commercial_only)
  
  
  //////////////////////////////////// Global /////////////////////////////////////////
  
  ////////
  // Data
  ////////
  
  DATA_MATRIX( Cov_xj );  

  // Aniso objects
  DATA_STRUCT(spde, spde_t);

  //////////////
  // Parameters 
  //////////////

  PARAMETER_VECTOR( deltainput_x );
  PARAMETER_VECTOR(logtau_S);
  PARAMETER_VECTOR(logkappa_S);
  PARAMETER_VECTOR(beta_j);
  
  /////////////////
  // derived values
  /////////////////
  
  int n_j = Cov_xj.row(0).size();
  // int n_k = Cov_xk.row(0).size();
  int n_x = Cov_xj.col(0).size();
  int n_S = logkappa_S.size();
  
  // global stuff
  vector<Type> jnll_comp(5);
  jnll_comp.setZero();
  vector<Type> MargSD_S(n_S);
  vector<Type> Range_S(n_S);
  for( int s=0; s<n_S; s++){
    MargSD_S(s) = 1 / sqrt(4*M_PI) / exp(logtau_S(s)) / exp(logkappa_S(s));
    Range_S(s) = sqrt(8) / exp( logkappa_S(s) );
  }
  
  // Transform random effects
  vector<Type> delta_x( deltainput_x.size() );
  delta_x = deltainput_x / exp(logtau_S(0));
  
  ////////////////
  // Latent field
  ///////////////
  SparseMatrix<Type> Q;
  //random effects priors
  if( Options_vec(0)==0 ){
    Q = Q_spde_generalized(spde, exp(logkappa_S(0)), Options_vec(1));
    if(Options_vec(2)==1) jnll_comp(2) = GMRF(Q)(deltainput_x); // latent field for presence/absence
  }
  
  // Predicted densities
  vector<Type> S_x(n_x);
  vector<Type> debug_j(n_j);
  vector<Type> linpredS_x = Cov_xj * beta_j;
  for(int s=0; s<n_x; s++){
    S_x(s) = exp(linpredS_x(s) + delta_x(s));
  }
  
  //////////////////////////////////// Commercial /////////////////////////////////////////
  
  
  if(Options_vec(5) == 3 | Options_vec(5) == 1){
    
    ////////
    // Data
    ////////
    
    DATA_MATRIX( c_com_x );       	// Response (count) for each observation i (commercial data)
    // DATA_VECTOR( P_i );      // Proportion of survey unit that is surveyed for each observation i
    // DATA_VECTOR( A_x );      // Relative area of each survey unit s (can set = 1.0 if all the same size)
    //DATA_VECTOR( c_sci_x );       	// Response (count) for each observation i (scientific data)
    DATA_VECTOR( y_com_i );       	// Response (0:not surveyed, 1:surveyed) for each site (commercial data)
    DATA_FACTOR( index_com_i );
   	DATA_FACTOR( b_com_i );
	  DATA_FACTOR( VE_i );
	  
    DATA_VECTOR( q2_com ); 
    DATA_VECTOR( weights_com );

    if(Options_vec(5) == 3) weights_com(0) = 1;
    DATA_MATRIX( Cov_xk );  //design matrix for sampling of commercial data
    
    
    //////////////
    // Parameters 
    //////////////
	
    PARAMETER_MATRIX(q1_com);
    PARAMETER_VECTOR(logSigma_com);
	  PARAMETER_VECTOR( k_com );
	  PARAMETER_VECTOR( logSigma_catch );
	  PARAMETER_VECTOR( logMean_catch );
    /////////////////
    // derived values
    /////////////////
    
    int n_com_i = y_com_i.size();


	/////////////////////
	// Observation model
	////////////////////
	    
    if( Options_vec(9) == 1){
		
		vector<Type> Sigma_com = exp(logSigma_com);
    vector<Type> Sigma_catch = exp(logSigma_catch);

		vector<Type> E_com(n_com_i);
		vector<Type> encounterprob_com(n_com_i);
		vector<Type> log_notencounterprob_com(n_com_i);
		
		// Zero-inflated gamma distribution
		for(int i=0; i<n_com_i; i++){
		  if( !isNA(y_com_i(i)) ){
		    
			if(VE_i(i) == 0) E_com(i) = q2_com(0) * 1 * S_x(index_com_i(i));
			if(VE_i(i) != 0) E_com(i) = q2_com(VE_i(i)) * k_com(VE_i(i)) * S_x(index_com_i(i));
			
			// proba of encountering
			if( Options_vec(8)==0 ){
			  encounterprob_com(i) = plogis( q1_com(VE_i(i),0) + q1_com(VE_i(i),1)*log(E_com(i)) );
			  log_notencounterprob_com(i) = NA_REAL;
			}
			
			if( Options_vec(8)==1 ){
			  encounterprob_com(i) = plogis(q1_com(VE_i(i),1) * ( 1.0 - exp(-1 * E_com(i) * exp(q1_com(VE_i(i),0)))));
			  log_notencounterprob_com(i) = NA_REAL;
			}
			
			if( Options_vec(8)==2 ){
			  encounterprob_com(i) = ( 1.0 - exp(-1 * E_com(i) * exp(q1_com(VE_i(i),0))) );
			  log_notencounterprob_com(i) = -1 * E_com(i) * exp(q1_com(VE_i(i),0));
			}
			//if(y_com_i(i) > 0) jnll_comp(0) -= dgamma( y_com_i(i), 1 / pow(Sigma_com,2) , E_com(i) * pow(Sigma_com,2), true ); // prob_encount_sci(i) * 
			
			if( Options_vec(6)==1 ) jnll_comp(0) -= weights_com(0) * dzinfgamma(y_com_i(i), E_com(i)/encounterprob_com(i), encounterprob_com(i), log_notencounterprob_com(i), Sigma_com(VE_i(i)), true);
			if( Options_vec(6)==2 ) jnll_comp(0) -= weights_com(0) * dzinflognorm(y_com_i(i), log(E_com(i))-log(encounterprob_com(i)), encounterprob_com(i), log_notencounterprob_com(i), Sigma_com(VE_i(i)), true);
			if( Options_vec(6)==3 & y_com_i(i) > 0 ) jnll_comp(0) -= weights_com(0) * dlognorm(y_com_i(i), log(E_com(i)), Sigma_com(VE_i(i)), true);
		  }
		}
		
		if(Options_vec(11) == 1) jnll_comp(0) -= dnorm(log(k_com),logMean_catch(0),Sigma_catch(0),true).sum();;
		
		REPORT( q1_com );
		REPORT( Sigma_com );

	}
    
    //////////////////////////////////////
    // Sampling process of commercial data
    /////////////////////////////////////
    if( Options_vec(7) == 1){
      
      // Parameter
      PARAMETER_MATRIX(beta_k);
      PARAMETER_VECTOR( par_b );
      PARAMETER_MATRIX( etainput_x );
      int n_eta = etainput_x.row(0).size();
      
      // Derived values
      vector<Type> b( n_eta );
      matrix<Type> eta_x( etainput_x.col(0).size(), n_eta );
      matrix<Type> linpredR_x(n_x , n_eta );
      matrix<Type> lambda_x(n_x , n_eta );
      // matrix<Type> lambda_ref_x(n_x , n_eta );
      vector<Type> fact_S(n_x);
      
      // Latent field
      if( Options_vec(0)==0 ){
        for(int l=0; l<n_eta; l++){
          eta_x.col(l) = etainput_x.col(l) / exp(logtau_S(1+l));
          Q = Q_spde_generalized(spde, exp(logkappa_S(1+l)), Options_vec(1));
          if(Options_vec(3)==1) jnll_comp(3+l) = GMRF(Q)(etainput_x.col(l)); // latent field for sampling
        }

      }

      
      // Poisson point process
     for(int l=0; l<n_eta; l++){
		 for(int s=0; s<n_x; s++){
			linpredR_x = Cov_xk * beta_k;
		  if(Options_vec(10)==1) fact_S = exp(par_b(l))*(log(S_x(s)));
		  if(Options_vec(10)==2) fact_S = par_b(l)*(log(S_x(s)));
		   
			lambda_x(s,l) = exp( linpredR_x(s,l) + fact_S(s) + eta_x(s,l) );
			if(!isNA(c_com_x(s,l)) & lambda_x(s,l) > 0)   jnll_comp(1)-= weights_com(0) * (Type(1)-lambda_x(s,l) + c_com_x(s,l)*log(lambda_x(s,l))); // log density of poisson process (see Diggle (2013))
		 }
		 if(Options_vec(10)==1) b(l)=exp(par_b(l));
		 if(Options_vec(10)==2) b(l)=par_b(l);
	 }
      // Outputs
      REPORT( lambda_x );
      REPORT( beta_k );
      REPORT( b );
      REPORT( eta_x );
      REPORT( linpredR_x );
      
      //ADREPORT(eta_x);
      ADREPORT(lambda_x);
      //ADREPORT(lambda_ref_x);
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
    
    vector<Type> E_sci(n_sci_i);
    vector<Type> encounterprob_sci(n_sci_i);
    vector<Type> log_notencounterprob_sci(n_sci_i);
    
    // Zero-inflated gamme distribution
    for(int i=0; i<n_sci_i; i++){
      if( !isNA(y_sci_i(i)) ){
        E_sci(i) = k_sci(0) * q2_sci(0) * S_x(index_sci_i(i));

        // proba of encountering
        if( Options_vec(8)==0 ){
          encounterprob_sci(i) = plogis( q1_sci(0) + q1_sci(1)*log(E_sci(i)) );
          log_notencounterprob_sci(i) = NA_REAL;
        }
        
        if( Options_vec(8)==1 ){
          encounterprob_sci(i) = plogis(q1_sci(1)) * ( 1.0 - exp(-1 * E_sci(i) * exp(q1_sci(0))) );
          log_notencounterprob_sci(i) = NA_REAL;
        }
        
        if( Options_vec(8)==2 ){
          encounterprob_sci(i) = ( 1.0 - exp(-1 * E_sci(i) * exp(q1_sci(0)) ));
          log_notencounterprob_sci(i) = -1 * E_sci(i) * exp(q1_sci(0));
        }
        
        if( Options_vec(6)==1 ) jnll_comp(0) -= dzinfgamma(y_sci_i(i), E_sci(i)/encounterprob_sci(i), encounterprob_sci(i), log_notencounterprob_sci(i), Sigma_sci, true);
        if( Options_vec(6)==2 ) jnll_comp(0) -= dzinflognorm(y_sci_i(i), log(E_sci(i))-log(encounterprob_sci(i)), encounterprob_sci(i), log_notencounterprob_sci(i), Sigma_sci, true);
        if( Options_vec(6)==3 & y_sci_i(i) > 0 ) jnll_comp(0) -= dlognorm(y_sci_i(i), log(E_sci(i)), Sigma_sci, true);
          
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
  REPORT( Range_S );
  REPORT( MargSD_S );
  REPORT( beta_j );
  REPORT( delta_x );
  REPORT( linpredS_x );
  
  // REPORT( Unsampled_x );
  REPORT( jnll_comp );
  // REPORT( jnll );

  // Bias correction output
  //ADREPORT( beta_abs_j );
  //ADREPORT( beta_pos_j );
  ADREPORT( total_abundance );
  ADREPORT( logkappa_S );
  
  if(Options_vec(4)==1){
    ADREPORT(S_x);
    // ADREPORT(delta_x);
  }
  
  
  return jnll;
}
