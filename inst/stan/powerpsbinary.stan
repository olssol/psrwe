//
//  STRATIFIED POWER PRIOR FOR BINARY DATA
//  Power parameter As follows dirichlet
//
data {
  int<lower = 1> S;

  //existing data
  array[S] int<lower = 0> N0;
  array[S] real<lower = 0, upper = 1> YBAR0;

  //current data
  array[S] int<lower = 1> N1;
  array[S] int<lower = 0> YSUM1;

  //prior of vs
  vector<lower=0>[S] RS;

  //fix vs
  int<lower = 0, upper = 1> FIXVS;

  //target borrowing
  real<lower = 0> A;
}

transformed data {
  row_vector<lower = 0, upper = 1>[S] WS1;
  int<lower = 0> sn1;

  sn1 = sum(N1);
  for (i in 1:S) {
    WS1[i] = N1[i];
    WS1[i] = WS1[i]/sn1;
  }
}

parameters {
  simplex[S] vs;
  vector<lower=0, upper=1>[S]  thetas;
}

transformed parameters {
  array[S] real<lower = 0, upper = 1> as;
  array[S] real<lower = 0> alphas;
  array[S] real<lower = 0> betas;

  for (i in 1:S) {
    if (0 == N0[i]) {
      as[i] = 0;
    } else {
      if (0 == FIXVS) {
        as[i]  = 1 < A*vs[i]/N0[i] ? 1:A*vs[i]/N0[i];
      } else {
        as[i]  = 1 < A*RS[i]/N0[i] ? 1:A*RS[i]/N0[i];
      }
    }

    alphas[i] = as[i] * N0[i] * YBAR0[i]  + 1;
    betas[i]  = as[i] * N0[i] * (1-YBAR0[i]) + 1;
  }
}

model {
  //prior
  if (A > 0) {
    target += beta_lpdf(thetas | alphas, betas);
  } else {
    thetas ~ uniform(0,1);
  }
  vs ~ dirichlet(RS);

  //likelihood
  YSUM1 ~ binomial(N1, thetas);
}

generated quantities {
  real theta;
  theta = WS1 * thetas;
}
