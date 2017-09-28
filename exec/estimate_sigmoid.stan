functions {

   real hurdle_normal_lpdf(real y, real mu, real sigma, real hu) {
     if (y == 0) {
       return bernoulli_lpmf(1 | hu);
     } else {
       return bernoulli_lpmf(0 | hu) +
              normal_lpdf(y | mu, sigma);
     }
   }

   real hurdle_normal_logit_lpdf(real y, real mu, real sigma, real hu) {
     if (y == 0) {
       return bernoulli_logit_lpmf(1 | hu);
     } else {
       return bernoulli_logit_lpmf(0 | hu) +
              normal_lpdf(y | mu, sigma);
     }
   }

   real prob_fnc(real mu, real location, real scale) {
     return inv_logit((mu - location)/scale);
   }
}
data {
  int<lower=0> N_genes;
  int<lower=0> N_samp;
  int<lower=0> N_cond;
  matrix[N_genes, N_samp] Y;
  int<lower=1> design[N_samp];
  real<lower=0> f_nu;
  real<lower=0> f_sigma;

}

parameters {

  // The Sigmoid Parameters
  real<lower=0> location;
  real<upper=0> scale_inv;

  // The full distribution parameters
  vector[N_genes] sample_mean;
  matrix[N_genes, N_cond-1] beta;

  real f_mean;
  real<lower=0> f_hyper_sd;
  vector<lower=0>[N_genes] f_var;
  real<lower=0> lfc_sd;
}
transformed parameters {
  matrix[N_genes, N_cond] f_x;
  for(g in 1:N_genes){
    for(c in 1:(N_cond-1)){
      f_x[g, c] = sample_mean[g] + beta[g, c];
    }
    f_x[g, N_cond] = sample_mean[g] + (0 - sum(beta[g]));
  }
}

model {

  // If there are multiple groups for a gene with all zero there are multiple solution
  // for beta. I will fix them to the same value. I just don't know how...


  // Estimate the true means with the hurdle model
  f_hyper_sd ~ cauchy(0, 5);
  // f_nu ~ cauchy(0,5);
  // f_sigma ~ cauchy(0,5);
  // f_nu ~ normal(0, N_genes/4.0);
  // f_sigma ~ normal(0, 1);
  lfc_sd ~ normal(0, 1);
  f_var ~ scaled_inv_chi_square(f_nu, sqrt(f_sigma));
  sample_mean ~ normal(f_mean, f_hyper_sd);
  for(g in 1:N_genes){
    for(c in 1:(N_cond-1)){
      beta[g,c]  ~ normal(0, lfc_sd);
    }
  }

  for(g in 1:N_genes){
    for(s in 1:N_samp) {
      Y[g, s] ~ hurdle_normal(f_x[g, design[s]], sqrt(f_var[g]), prob_fnc(f_x[g, design[s]], location, 1/scale_inv));
    }
  }


  // Estimate the location and scale by logstic regression
  location ~ normal(0, 10);
  scale_inv ~ normal(0, 10);

  for(g in 1:N_genes){
    if(sum(Y[g]) != 0){
      for(s in 1:N_samp){
        (Y[g,s] != 0) ~ bernoulli_logit(prob_fnc(f_x[g, design[s]], location, 1/scale_inv));
      }
    }
  }


}
