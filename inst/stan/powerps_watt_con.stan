//
//  WATT POWER PRIOR FOR CONTINUOUS DATA
//  Power parameter As follows
//
data {
  real<lower = 0>   A;

  //existing data
  int<lower = 0>    N0;
  array[N0] real    Y0;
  array[N0] real    A_WATT_DI;

  //current data
  int<lower = 1>    N1;
  array[N1] real    Y1;
}

parameters {
  real              thetas;
  real<lower=0>     taus;
}

transformed parameters {
  array[N0] real<lower = 0>  sd0_watt_di;

  for(i in 1:N0) {
    sd0_watt_di[i] = taus / A_WATT_DI[i];
  }
}

model {
  //prior
  if (A > 0) {
    target += normal_lpdf(Y0 | thetas, sd0_watt_di);
  }

  thetas ~ normal(0, 1000);
  taus ~ cauchy(0, 2.5);

  //likelihood
  Y1 ~ normal(thetas, taus);
}

