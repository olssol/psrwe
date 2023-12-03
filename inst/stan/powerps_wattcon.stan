//
//  WATT POWER PRIOR FOR CONTINUOUS DATA
//  Power parameter As follows
//
data {
  //target borrowing
  real<lower = 0>  A;

  //existing data
  real             Y0Tilde;
  real<lower = 0>  SD0;

  //current data
  int<lower = 1> N1;
  array[N1] real Y1;
}

parameters {
  real           thetas;
  real<lower=0>  taus;
}


model {
  //prior
  taus ~ cauchy(0, 2.5);

  if (A > 0) {
    target += normal_lpdf(Y0Tilde | thetas, SD0 / sqrt(A));
  }

  //likelihood
  Y1 ~ normal(thetas, taus);
}
