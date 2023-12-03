//
//  WATT POWER PRIOR FOR CONTINUOUS DATA
//  Power parameter As follows
//
data {
  //target borrowing
  real<lower = 0>   A;
  real<lower = 0> SD0;

  //existing data
  int<lower = 0>            N0;
  array[N0] real            Y0;
  array[N0] real<lower = 0> A_WATT_DI;

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
    for (i in 1:N0) {
      target += A_WATT_DI[i] * normal_lpdf(Y0[i] | thetas, SD0);
    }
  }

  //likelihood
  Y1 ~ normal(thetas, taus);
}
