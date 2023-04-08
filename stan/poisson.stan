data {
  
  int<lower=0> N_obs;
  int<lower=0> P;
  matrix[N_obs, P] X;
  int<lower=0> Y_pos[N_obs];
  int<lower=0> Y_neg[N_obs];
  vector[N_obs] N_samp;
  
}
transformed data {
  
  vector[N_obs] offset = log(N_samp);

}
parameters {
  
  vector[P] beta_pos;
  vector[P] beta_neg;
  
}
model {
  
  beta_pos ~ std_normal();
  beta_neg ~ std_normal();
  
  Y_pos ~ poisson_log(offset + X * beta_pos);
  
  Y_neg ~ poisson_log(offset + X * beta_neg);
  
}
