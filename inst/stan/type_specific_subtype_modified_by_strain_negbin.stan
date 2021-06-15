// The input data
data {

  int<lower=0> i_max;
  int<lower=0> j_max;
  int<lower=0> k_max;
  int<lower=0> n_obs;
  int<lower=0> i_values[n_obs];
  int<lower=0> j_values[n_obs];
  int<lower=0> k_values[n_obs];
  int<lower=0> c_ij[n_obs];
  int<lower=0> d_ij[n_obs];
  int<lower=0> n_i[n_obs];
  int<lower=0> N_i[n_obs];
  vector<lower=0>[n_obs] t_i;

}

// The parameters accepted by the model
parameters {

  // carriage prevalence ~ U(0,1)
  vector<lower=0.0,upper=1.0>[n_obs] rho_ij;

  // log serotype invasiveness ~ U(-9,0)
  vector<lower=-6.0,upper=0.0>[j_max] log_nu_j;

  // log GPSC invasiveness ~ U(3,-3)
  vector<lower=-3.0,upper=3.0>[k_max] log_nu_k;

  //vector<lower=0.0,upper=1.0>[j_max] nu_j;

  // negative binomial overdispersions
  real<lower=0.0,upper=10.0> phi_nb;

}

transformed parameters {

  // declare transformed parameters
  vector<lower=0,upper=1.0>[j_max] nu_j;
  vector<lower=0.001,upper=1000.0>[k_max] nu_k;

  // calculate serotype invasiveness on a real scale
  for (j in 1:j_max) {
    nu_j[j] = pow(10, log_nu_j[j]);
  }

  // calculate serotype invasiveness on a real scale
  for (k in 1:k_max) {
    nu_k[k] = pow(10, log_nu_k[k]);
  }

}

// The model to be estimated
model {

  // iterate over datasets
  for (index in 1:n_obs) {

    // Get serotype
    int j = j_values[index];

    // Get GPSC
    int k = k_values[index];

    // Get location adjustment
    int i = i_values[index];

    // calculate prior probability
    target += uniform_lpdf( log_nu_j[j] | -9, 0);
    target += uniform_lpdf( log_nu_k[k] | -3, 3);
    target += uniform_lpdf(rho_ij[index] | 0,1);
    target += uniform_lpdf(phi_nb | 0,10);

    // calculate likelihood given data
    target += binomial_lpmf( c_ij[index] | n_i[index], rho_ij[index] );
    target += neg_binomial_2_lpmf( d_ij[index] | nu_j[j]*nu_k[k]*rho_ij[index]*N_i[index]*t_i[index], phi_nb );

  }
}

generated quantities {

  // Calculate and store log likelihood for loo
  vector[2*n_obs] log_lik;

  // Calculate and store predictions for carriage
  vector[n_obs] c_ij_pred;

  // Calculate and store predictions for disease
  vector[n_obs] d_ij_pred;

  // iterate over datasets
  for (index in 1:n_obs) {

    // Get serotype
    int j = j_values[index];

    // Get GPSC
    int k = k_values[index];

    // Get location adjustment
    int i = i_values[index];

    // Store predictions
    c_ij_pred[index] = n_i[index]*rho_ij[index];
    d_ij_pred[index] = nu_j[j]*nu_k[k]*rho_ij[index]*N_i[index]*t_i[index];

    // Calculate likelihood given data
    log_lik[2*(index-1)+1] = binomial_lpmf( c_ij[index] | n_i[index], rho_ij[index] );
    log_lik[2*(index-1)+2] = neg_binomial_2_lpmf(  d_ij[index] | d_ij_pred[index], phi_nb );

  }



}
