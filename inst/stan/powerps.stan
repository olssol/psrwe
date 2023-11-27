//
//  STRATIFIED POWER PRIOR FOR CONTINUOUS DATA
//  Power parameter As follows
//
data {
  int<lower = 1> S;

  //existing data
  array[S] int<lower = 0>  N0;
  array[S] real            YBAR0;
  array[S] real<lower = 0> SD0;

  //current data
  array[S] int<lower = 1>   N1;
  int<lower = 1>            TN1;
  array[TN1] real           Y1;
  array[TN1] int<lower = 1> INX1;

  //prior of vs
  vector<lower=0>[S] RS;

  //fix vs
  int<lower = 0, upper = 1> FIXVS;

  //target borrowing
  real<lower = 0> A;
}

transformed data {
  row_vector<lower = 0, upper = 1>[S] WS1;
  for (i in 1:S) {
    WS1[i] = N1[i];
    WS1[i] = WS1[i]/TN1;
  }
}

parameters {
  simplex[S]             vs;
  vector[S]              thetas;
  array[S] real<lower=0> taus;
}

transformed parameters {
  array[S] real<lower = 0, upper = 1> as;
  array[S] real<lower = 0> sds;

  for (i in 1:S) {
    if (0 == N0[i]) {
      as[i] = 0;
    } else {
      if (0 == FIXVS) {
        as[i]  = 1 < A * vs[i] / N0[i] ? 1 : A * vs[i] / N0[i];
      } else {
        as[i]  = 1 < A * RS[i] / N0[i] ? 1 : A * RS[i] / N0[i];
      }
    }
    sds[i] = 0 == as[i] ? 0 : SD0[i] / sqrt(as[i] * N0[i]);
  }
}

model {
  //prior
  if (A > 0) {
    target += normal_lpdf(YBAR0 | thetas, sds);
  } else {
    thetas ~ normal(0, 1000);
  }

  vs      ~ dirichlet(RS);
  taus    ~ cauchy(0, 2.5);

  //likelihood
  Y1 ~ normal(thetas[INX1], taus[INX1]);
}

generated quantities {
  real theta;
  theta = WS1 * thetas;
}
