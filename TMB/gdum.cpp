#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(Y);         // Vector of unshared species
  DATA_VECTOR(Y_den);    // FIXED: missing for binomial family
  DATA_MATRIX(D);         // Matrix mapping Z to sample pairs (n_pairs x 2)
  DATA_MATRIX(X);         // Pair-level variables (n_samples x n_env)
  DATA_MATRIX(W);         // sample-level variables (n_samples x n_env)
  DATA_VECTOR(weights)
  
  DATA_SPARSE_MATRIX(Z);    // Random effects design matrix
  DATA_INTEGER(has_random); // Flag for re
  DATA_INTEGER(mono);       // Flag for monotonic effects
  DATA_INTEGER(link);       // link function
  DATA_INTEGER(family);     // family
  DATA_IVECTOR(map_re);     // Mapping from RE columns
  
  // Parameters
  PARAMETER(intercept);                   // Intercept
  PARAMETER_VECTOR(beta);                 // Slopes
  PARAMETER_VECTOR(lambda);               // Slopes
  PARAMETER_VECTOR(u);                    // random effects
  PARAMETER_VECTOR(log_sigma_re);         // sigma of random effects
  PARAMETER(log_scale);   // scale parameter

  // Derived quantities
  int n_pairs = Y.size();
  int n_X = X.cols();
  int n_W = W.cols();
  
  // Contributions if site-level effects, pair-level effects, and re
  vector<Type> pair_comp_contrib(n_pairs);
  vector<Type> site_comp_contrib(n_pairs);
  vector<Type> re_comp(n_pairs);

  pair_comp_contrib.setZero();
  site_comp_contrib.setZero();
  re_comp.setZero();
  
  // Initiate log likelihood
  Type nll = 0.0;
  
  // Handle random effects if present
  if (has_random) {

    vector<Type> sigma_re = exp(log_sigma_re);
    //Center the random effects
    Type u_mean = u.sum() / u.size();
    vector<Type> u_0 = u - u_mean;
    //Calculate re contribution to nll
    vector<Type> Z_cont = Z * u_0;
    for (int i = 0; i < u.size(); i++) {
      nll -= dnorm(u(i), Type(0), sigma_re(map_re[i]), true);
    }
    //Calculate random effects contribution
    for (int i = 0; i < n_pairs; i++) {
      int s1 = CppAD::Integer(D(i, 0));  
      int s2 = CppAD::Integer(D(i, 1));
      re_comp(i) = Z_cont(s1) + Z_cont(s2);
    }
    ADREPORT(u_0);
    REPORT(u_0);
  }
  
  // Calculate fixed effects components
  if (n_X > 0) {
    matrix<Type> pair_comp(n_pairs, n_X);
    for (int i = 0; i < n_pairs; i++) {
      int s1 = CppAD::Integer(D(i, 0));  
      int s2 = CppAD::Integer(D(i, 1));  
      pair_comp.row(i) = (X.row(s1) - X.row(s2)).array().abs();
    }
    pair_comp_contrib = pair_comp * beta;
  }
  
  if (n_W > 0) {
    matrix<Type> site_comp(n_pairs, n_W);
    for (int i = 0; i < n_pairs; i++) {
      int s1 = CppAD::Integer(D(i, 0));  
      int s2 = CppAD::Integer(D(i, 1));  
      site_comp.row(i) = (W.row(s1) + W.row(s2)).array();
    }
    site_comp_contrib = site_comp * lambda;
  }
  
  // full linear predictor
  vector<Type> eta = intercept + pair_comp_contrib + site_comp_contrib + re_comp;
  vector<Type> mu(n_pairs); 
  
  // LINK 
  if (link == 0) {
    mu = eta;
  }
  
  if (link == 1) {
    mu = exp(eta)/(1+exp(eta));
  }
  
  if (link == 2) {
    mu = Type(1) - exp(-eta);
  }
  
  if (link == 3) {
    mu = Type(1) - exp(-(eta*eta));
  }

  // compute family-dependent neg. log-likelihood 
  if (family == 0) {
    Type scale = exp(log_scale);
    for (int i = 0; i < n_pairs; i++) {
      nll -= weights(i)*dnorm(Y(i), mu(i), scale, true);
    }
    ADREPORT(scale);
  }
  
  else if (family == 1) {
    for (int i = 0; i < n_pairs; i++) {
     nll -= weights(i)*dbinom(Y(i), Y_den(i), mu(i), true);
    }
  }
  
  else if (family == 2) {
    Type scale = exp(log_scale);
    // reparam. alpha beta into mu phi (scale)
    vector<Type> a = mu * scale;
    vector<Type> b = (Type(1) - mu) * scale;
    // compute likelihood
    for (int i = 0; i < n_pairs; i++) {
      nll -= weights(i)*dbeta(Y(i), a(i), b(i), true);
    }
    ADREPORT(scale);
  }
  REPORT(beta);
  REPORT(lambda);
  REPORT(intercept);
  return nll;
}



