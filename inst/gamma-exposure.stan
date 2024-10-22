data {
  int N;
  vector[N] ncomp;
  vector[N] conc;
}
parameters {
  real<lower=0> shape;
  real<lower=0> scale;
}
model {
  // shape ~ std_normal();
  // scale ~ std_normal();
  
  // conc ~ gamma(ncomp * shape, 1 / (scale * ncomp));
  for (i in 1 : N) {
    conc[i] ~ gamma(ncomp[i] * shape, ncomp[i] / scale);
  }
}
generated quantities {
  array[N] real conc_post;
  // conc_post = gamma_rng(ncomp * shape, 1 / (scale * ncomp));
  for (i in 1 : N) {
    conc_post[i] = gamma_rng(ncomp[i] * shape, ncomp[i] / scale);
  }
}
// TODO Look at including a covariate(s) in the model to see if that can account
// for the excess low values in the current model. That or we should go back to
// the log-normal model.
