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

  // log serotype invasiveness ~ U(-6,1)
  vector<lower=-6.0,upper=1.0>[j_max] log_nu_j;

  // log GPSC invasiveness ~ Cauchy
  vector<lower=-3, upper=3>[k_max-1-1] log_nu_k;

  // negative binomial overdispersions
  real<lower=-3,upper=3> log_phi_nb;

}

transformed parameters {

  // declare transformed parameters
  vector<lower=0,upper=10.0>[j_max] nu_j;
  vector[k_max] nu_k;
  real phi_nb;
  real mu_mod = 0; // position parameter of Cauchy for strain invasiveness
  real tau_mod = 1; // scale parameter of Cauchy for strain invasiveness

  // calculate serotype invasiveness on a real scale
  for (j in 1:j_max) {
    nu_j[j] = pow(10, log_nu_j[j]);
  }

  // calculate serotype invasiveness on a real scale
  nu_k[1] = 1;
  for (k in 2:k_max) {
    nu_k[k] = pow(10, log_nu_k[k-1]);
  }

  // calculate negative binomial overdispersion
  phi_nb = pow(10, log_phi_nb);

}

// The model to be estimated
model {

  // Calculate prior probability for types
  for (j in 1:j_max) {
    target += uniform_lpdf(log_nu_j[j] | -6, 1);
  }

  // Calculate prior probability for strains
  for (k in 1:k_max) {
    target += uniform_lpdf(log_nu_k[k] | -1.25, 1.25);
  }

  // Calculate prior probability for precision parameter
  target += uniform_lpdf(log_phi_nb | -3, 3);

  // iterate over datasets
  for (index in 1:n_obs) {

    // Get serotype
    int j = j_values[index];

    // Get GPSC
    int k = k_values[index];

    // Get location adjustment
    int i = i_values[index];

    // Calculate prior probability for carriage frequency
    target += beta_lpdf(rho_ij[index] | 1, 1);

    // calculate likelihood given data
    target += binomial_lpmf(c_ij[index] | n_i[index], rho_ij[index]);
    target += neg_binomial_2_lpmf(d_ij[index] | nu_j[j]*nu_k[k]*rho_ij[index]*N_i[index]*t_i[index], phi_nb);

  }
}

generated quantities {

  // Calculate and store log likelihood for loo
  vector[n_obs] carriage_log_lik;
  vector[n_obs] disease_log_lik;
  vector[n_obs] log_lik;

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
    carriage_log_lik[index] = binomial_lpmf( c_ij[index] | n_i[index], rho_ij[index] );
    disease_log_lik[index] = poisson_lpmf(  d_ij[index] | d_ij_pred[index] );
    log_lik[index] = carriage_log_lik[index] + disease_log_lik[index];

  }

}
