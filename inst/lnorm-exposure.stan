functions {
  real lognormal_mean_sdlog(int n, real pop_meanlog, real pop_sdlog) {
    real mean_varlog = log((exp(pop_sdlog ^ 2) - 1) / n + 1);
    return sqrt(mean_varlog);
  }
  real lognormal_mean_meanlog(real pop_meanlog, real pop_sdlog,
                              real mean_sdlog) {
    real mean_meanlog = pop_meanlog + pop_sdlog / 2 - mean_sdlog / 2;
    return mean_meanlog;
  }
}
data {
  int C; // Number of different compositions
  array[C] int ncomp;
  
  int N;
  array[N] int ncomp_idx;
  vector[N] conc;
}
parameters {
  real pop_meanlog;
  real<lower=0> pop_sdlog;
}
transformed parameters {
  vector[C] mean_meanlog;
  vector[C] mean_sdlog;
  
  for (i in 1 : C) {
    mean_sdlog[i] = lognormal_mean_sdlog(ncomp[i], pop_meanlog, pop_sdlog);
    mean_meanlog[i] = lognormal_mean_meanlog(pop_meanlog, pop_sdlog,
                                             mean_sdlog[i]);
  }
}
model {
  target += -log(pop_sdlog);
  conc ~ lognormal(mean_meanlog[ncomp_idx], mean_sdlog[ncomp_idx]);
}
generated quantities {
  array[N] real conc_gen;
  conc_gen = lognormal_rng(mean_meanlog[ncomp_idx], mean_sdlog[ncomp_idx]);
  
  array[N] real conc_gen2 = rep_array(0.0, N);
  for (n in 1 : N) {
    for (c in 1 : ncomp[ncomp_idx[n]]) {
      conc_gen2[n] += lognormal_rng(pop_meanlog, pop_sdlog);
    }
    conc_gen2[n] /= ncomp[ncomp_idx[n]];
  }
}
