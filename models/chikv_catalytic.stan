
data {
  int<lower=1> R;           // number of municipalities
  int<lower=1> T;           // number of years (2014-2025 = 12)
  int<lower=1> T_survey;    // number of survey time points (= 2)
  int<lower=1> A_sero;      // number of age groups in serology (1-yr intervals)
  int<lower=1> A_case;      // number of age groups in cases (5-yr intervals)

  array[T_survey] int<lower=1>   survey_year_idx;  // 1-based index (2023=10, 2024=11)
  array[A_sero]   int<lower=0>   ages_sero;        // actual ages (1-yr), e.g. 16..69
  array[A_case]   int<lower=0>   ages_case;        // lower bound of 5-yr groups (0,5,...,80)

  array[R, A_sero, T_survey] int<lower=0>    pos;     // seropositive counts
  array[R, A_sero, T_survey] int<lower=0>    tested;  // total tested
  array[R, A_case, T]        int<lower=0>    cases;   // notified case counts
  array[R, A_case, T]        real<lower=0>   N;       // population
}

parameters {
  array[R, T] real<lower=0> lambda;   // force of infection (region x year)
  vector<lower=0, upper=1>[R] rho;    // infection-to-reporting ratio (region-specific, time-invariant)
  real<lower=0> phi;                  // NegBin overdispersion
}

transformed parameters {

  // -- Seroprevalence (cumulative exposure probability) -----------
  array[R, A_sero, T_survey] real<lower=0, upper=1> pi_sero;

  for (r in 1:R) {
    for (s in 1:T_survey) {
      int ts = survey_year_idx[s];
      for (a in 1:A_sero) {
        // birth year index (1-based, 1=2014)
        // e.g. age=20, ts=10(2023) -> t_birth = 10-20+1 = -9 -> max(1,.) = 1
        int t_birth = ts - ages_sero[a] + 1;
        int t_start = max(1, t_birth);
        real cum_foi = 0.0;
        for (t in t_start:ts) {
          cum_foi += lambda[r, t];
        }
        pi_sero[r, a, s] = 1.0 - exp(-cum_foi);
      }
    }
  }

  // -- Expected case counts --------------------------------------
  array[R, A_case, T] real<lower=0> mu_cases;

  for (r in 1:R) {
    for (t in 1:T) {
      for (a in 1:A_case) {
        // representative age: midpoint of 5-yr group
        int age_rep = ages_case[a] + 2;
        // cumulative FOI up to (but not including) year t
        real cum_foi_prev = 0.0;
        if (t > 1) {
          int t_birth = t - age_rep;
          int t_start = max(1, t_birth);
          for (t2 in t_start:(t - 1)) {
            cum_foi_prev += lambda[r, t2];
          }
        }
        real s_frac = exp(-cum_foi_prev);  // susceptible fraction
        mu_cases[r, a, t] = rho[r] * N[r, a, t] * s_frac * lambda[r, t];
      }
    }
  }
}

model {
  // -- Priors ----------------------------------------------------
  // TODO: adjust priors based on domain knowledge
  for (r in 1:R)
    for (t in 1:T)
      lambda[r, t] ~ gamma(1, 10);   // mean 0.1, skewed toward low FOI

  rho ~ beta(2, 18);                 // mean ~0.1, low reporting rate prior
  phi ~ exponential(1);

  // -- Likelihood: serology --------------------------------------
  for (r in 1:R)
    for (s in 1:T_survey)
      for (a in 1:A_sero)
        if (tested[r, a, s] > 0)
          pos[r, a, s] ~ binomial(tested[r, a, s], pi_sero[r, a, s]);

  // -- Likelihood: notified cases --------------------------------
  for (r in 1:R)
    for (t in 1:T)
      for (a in 1:A_case)
        cases[r, a, t] ~ neg_binomial_2(mu_cases[r, a, t] + 1e-10, phi);
}

generated quantities {
  array[R, A_sero, T_survey] int pos_pred;
  array[R, A_case, T]        int cases_pred;

  for (r in 1:R) {
    for (s in 1:T_survey)
      for (a in 1:A_sero)
        pos_pred[r, a, s] = binomial_rng(tested[r, a, s], pi_sero[r, a, s]);

    for (t in 1:T)
      for (a in 1:A_case)
        cases_pred[r, a, t] = neg_binomial_2_rng(mu_cases[r, a, t] + 1e-10, phi);
  }
}

