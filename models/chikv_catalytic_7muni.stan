
data {
  int<lower=1> R;        // number of municipalities (= 7)
  int<lower=1> T;        // number of years
  int<lower=1> T_survey; // number of survey time points
  int<lower=1> A_sero;   // number of 1-yr sero age groups
  int          year_start; // first calendar year in the model (e.g. 2014)

  array[T_survey] int<lower=1>               survey_year_idx;
  array[A_sero]   int<lower=0>               ages_sero;        // 1-yr ages
  array[R, A_sero, T_survey] int<lower=0>    pos;              // sero positives
  array[R, A_sero, T_survey] int<lower=0>    tested;           // sero tested
  array[R, T]                int<lower=0>    cases;            // notified cases
  vector<lower=0>[R]                         N_init;           // 2014 population

  // -- Reference literature sero data (age-band aggregated) -------
  int<lower=0>                       N_ref;
  array[N_ref] int<lower=1, upper=R> r_ref;
  array[N_ref] int<lower=0>          age_lo_ref;
  array[N_ref] int<lower=0>          age_hi_ref;
  array[N_ref] int<lower=1, upper=T> survey_yr_ref_idx;
  array[N_ref] int<lower=0>          pos_ref;
  array[N_ref] int<lower=0>          tested_ref;
}

parameters {
  array[R, T] real<lower=0>    lambda;  // FOI (municipality x year)
  vector<lower=0, upper=1>[R]  rho;     // underreporting factor (municipality)
}

transformed parameters {

  // -- Susceptible population --------------------------------------
  array[R, T] real<lower=0> S;

  for (r in 1:R) {
    S[r, 1] = N_init[r];
    for (t in 2:T) {
      S[r, t] = fmax(S[r, t - 1] * exp(-lambda[r, t - 1]), 1.0);
    }
  }

  // -- Expected cases ----------------------------------------------
  array[R, T] real<lower=0> mu_cases;

  for (r in 1:R)
    for (t in 1:T)
      mu_cases[r, t] = fmax(rho[r] * S[r, t] * lambda[r, t], 1e-10);

  // -- Seroprevalence: age-specific cumulative FOI -----------------
  array[R, A_sero, T_survey] real<lower=0, upper=1> pi_sero;

  for (r in 1:R) {
    for (s in 1:T_survey) {
      int ts = survey_year_idx[s];                     // time index of survey year
      int survey_yr = year_start + ts - 1;             // calendar year of survey
      for (a in 1:A_sero) {
        int birth_yr = survey_yr - ages_sero[a];       // calendar birth year
        int t_birth  = birth_yr - year_start + 1;      // time index of birth year
        int t_start  = max(1, t_birth);                // clamp to model start
        real cum_foi = 0.0;
        for (t in t_start:ts)
          cum_foi += lambda[r, t];
        pi_sero[r, a, s] = 1.0 - exp(-cum_foi);
      }
    }
  }
}

model {
  // -- Priors -------------------------------------------------------
  //for (r in 1:R)
  //  for (t in 1:T)
  //    lambda[r, t] ~ exponential(2);
  //rho ~ beta(2, 18);
  // case + sero model is free from identifiability issue. prior is not needed

  // -- Serology likelihood: age-specific binomial (primary) --------
  for (r in 1:R)
    for (s in 1:T_survey)
      for (a in 1:A_sero)
        if (tested[r, a, s] > 0)
          pos[r, a, s] ~ binomial(tested[r, a, s], pi_sero[r, a, s]);

  // -- Case likelihood: Poisson (temporal pattern) -----------------
  for (r in 1:R)
    for (t in 1:T)
      cases[r, t] ~ poisson(mu_cases[r, t]);
}

generated quantities {
  array[R, A_sero, T_survey] int  pos_pred;   // predicted sero positives
  array[R, T]                int  cases_pred; // predicted case counts

  // full age range (0-80) seroprevalence for visualization
  array[R, 81, T_survey] real<lower=0, upper=1> pi_sero_full;

  for (r in 1:R)
    for (s in 1:T_survey)
      for (a in 1:A_sero)
        pos_pred[r, a, s] = binomial_rng(tested[r, a, s], pi_sero[r, a, s]);

  for (r in 1:R)
    for (t in 1:T)
      cases_pred[r, t] = poisson_rng(mu_cases[r, t]);

  // compute seroprevalence for ages 0-80 (index 1 = age 0, index 81 = age 80)
  for (r in 1:R) {
    for (s in 1:T_survey) {
      int ts        = survey_year_idx[s];
      int survey_yr = year_start + ts - 1;
      for (age in 0:80) {
        int birth_yr = survey_yr - age;
        int t_birth  = birth_yr - year_start + 1;
        int t_start  = max(1, t_birth);
        real cum_foi = 0.0;
        for (t in t_start:ts)
          cum_foi += lambda[r, t];
        pi_sero_full[r, age + 1, s] = 1.0 - exp(-cum_foi);
      }
    }
  }
}

