// The input data
data {

  int<lower=0> i_max;
  int<lower=0> j_max;
  int<lower=0> n_obs;
  int<lower=0> i_values[n_obs];
  int<lower=0> j_values[n_obs];
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

  // log invasiveness ~ U(-6,1)
  real<lower=-6,upper=1.0> log_nu;

  // dataset adjustment
  vector<lower=-pi()/2, upper=pi()/2>[i_max-1] gamma_varying;

  // negative binomial overdispersions
  real phi_nb;

}

transformed parameters {

  // declare transformed parameters
  vector[i_max] gamma_i;
  real mu = 0; // position parameter of Cauchy for gamma
  real tau = 2; // scale parameter of Cauchy for gamma
  real midpoint_inv = log(2.0)/pow(10, -2.5); // midpoint of invasiveness range

  // calculate invasiveness on a real scale
  real<lower=0,upper=10.0> nu;
  nu = pow(10, log_nu);

  // add constant to gamma vector
  gamma_i[1] = 1;
  for (i in 2:i_max) {
    gamma_i[i] = pow(10, mu + tau * tan(gamma_varying[i-1]));
  }

}

// The model to be estimated
model {

  // Calculate prior probability across all types
  target += uniform_lpdf(log_nu | -6, 1);

  // Calculate prior probability for precision parameter
  target += exponential_lpdf(phi_nb | midpoint_inv);

  // Calculate prior probability for study adjustment
  for (i in 2:i_max) {
    target += uniform_lpdf(gamma_varying[i-1] | -pi()/2, pi()/2);
  }

  // iterate over datasets
  for (index in 1:n_obs) {

    // Get location adjustment
    int i = i_values[index];

    // Calculate prior probability for carriage frequency
    target += beta_lpdf(rho_ij[index] | 1, 1);

    // calculate likelihood given data
    target += binomial_lpmf(c_ij[index] | n_i[index], rho_ij[index]);
    target += neg_binomial_2_lpmf(d_ij[index] | gamma_i[i]*nu*rho_ij[index]*N_i[index]*t_i[index], phi_nb);

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

    // Get location adjustment
    int i = i_values[index];

    // Store predictions
    c_ij_pred[index] = n_i[index]*rho_ij[index];
    d_ij_pred[index] = gamma_i[i]*nu*rho_ij[index]*N_i[index]*t_i[index];

    // Calculate likelihood given data
    carriage_log_lik[index] = binomial_lpmf( c_ij[index] | n_i[index], rho_ij[index] );
    disease_log_lik[index] = neg_binomial_2_lpmf(  d_ij[index] | d_ij_pred[index], phi_nb );
    log_lik[index] = carriage_log_lik[index] + disease_log_lik[index];

  }

}
