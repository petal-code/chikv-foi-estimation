// ══════════════════════════════════════════════════════════════════
// Model B: Literature review sero data + SINAN case data
// Compare to chikv_catalytic_7muni.stan (Model A) which uses
// primary blood donor sero data from 7 fixed municipalities.
//
// Key difference from Model A:
//   - Sero data is aggregated over age bands (not 1-yr individual bins)
//   - Municipalities come from the literature, not fixed to 7
//   - No 3D sero array — each observation is a flat record
// ══════════════════════════════════════════════════════════════════

data {

  int<lower=1> R;          // number of municipalities (from lit review)
  int<lower=1> T;          // number of years in the model (e.g. 2014–2025 = 12)
  int<lower=0> N_ref;      // number of literature sero observations (rows in dataset)
  int          year_start; // first calendar year of the model (e.g. 2014)

  // -- Literature sero observations --------------------------------
  // Each row i is one published study observation:
  // a municipality, a survey year, an age band, and pos/tested counts.
  array[N_ref] int<lower=1, upper=R> r_ref;             // municipality index for obs i
  array[N_ref] int<lower=0>          age_lo_ref;        // minimum age in study (years)
  array[N_ref] int<lower=0>          age_hi_ref;        // maximum age in study (years)
  array[N_ref] int<lower=1, upper=T> survey_yr_ref_idx; // time index of survey year
  array[N_ref] int<lower=0>          pos_ref;           // observed seropositives
  array[N_ref] int<lower=0>          tested_ref;        // total tested

  // -- Case data (SINAN) -------------------------------------------
  array[R, T] int<lower=0> cases;   // notified cases per municipality per year

  // -- Population --------------------------------------------------
  vector<lower=0>[R] N_init;        // 2014 population per municipality
}

parameters {

  // Annual force of infection per municipality per year
  array[R, T] real<lower=0> lambda;

  // Underreporting factor per municipality (fraction of true cases notified)
  vector<lower=0, upper=1>[R] rho;
}

transformed parameters {

  // -- Susceptible population dynamics -----------------------------
  // S[r, 1] = N_init[r] (full population susceptible at model start)
  // S[r, t] decays each year as people get infected
  array[R, T] real<lower=0> S;

  for (r in 1:R) {
    S[r, 1] = N_init[r];
    for (t in 2:T)
      S[r, t] = fmax(S[r, t - 1] * exp(-lambda[r, t - 1]), 1.0);
  }

  // -- Expected case counts ----------------------------------------
  // Observed cases = rho (reporting fraction) x new infections
  // New infections in year t = S[r,t] * lambda[r,t]
  array[R, T] real<lower=0> mu_cases;

  for (r in 1:R)
    for (t in 1:T)
      mu_cases[r, t] = fmax(rho[r] * S[r, t] * lambda[r, t], 1e-10);

  // -- Seroprevalence: reference-aggregated data fitting -----------
  // For each literature observation i, we don't know individual ages —
  // only that participants were aged [age_lo_ref, age_hi_ref].
  // We compute the catalytic model seroprevalence for every integer
  // age in that band, then average across the band.
  array[N_ref] real<lower=0, upper=1> pi_ref;

  for (i in 1:N_ref) {

    int ts        = survey_yr_ref_idx[i];  // time index of the survey year
    int survey_yr = year_start + ts - 1;   // convert back to calendar year

    // Number of integer ages in the band (used to compute the average)
    int  n_ages = age_hi_ref[i] - age_lo_ref[i] + 1;
    real cum_pi = 0.0;  // running sum of pi(a) across the band

    for (a in age_lo_ref[i]:age_hi_ref[i]) {

      // A person aged a at survey time was born in:
      int birth_yr = survey_yr - a;

      // Convert birth year to model time index
      int t_birth = birth_yr - year_start + 1;

      // Clamp to model start — we assume no CHIKV exposure before year_start
      int t_start = max(1, t_birth);

      // Sum annual FOI from birth (or model start) to survey year
      real cum_foi = 0.0;
      for (t in t_start:ts)
        cum_foi += lambda[r_ref[i], t];

      // Catalytic model: P(seropositive) = 1 - P(never infected)
      cum_pi += 1.0 - exp(-cum_foi);
    }

    // Average seroprevalence across the age band
    pi_ref[i] = cum_pi / n_ages;
  }
}

model {

  // -- Seroprevalence likelihood (literature age-band binomial) ----
  // For each literature observation: observed positives ~ Binomial(tested, pi_ref)
  for (i in 1:N_ref)
    pos_ref[i] ~ binomial(tested_ref[i], pi_ref[i]);

  // -- Case likelihood (Poisson) -----------------------------------
  // Notified cases in each municipality and year ~ Poisson(expected cases)
  for (r in 1:R)
    for (t in 1:T)
      cases[r, t] ~ poisson(mu_cases[r, t]);
}

generated quantities {

  // Posterior predictive: simulated seropositives for each lit observation
  array[N_ref] int pos_pred;

  // Posterior predictive: simulated case counts
  array[R, T] int cases_pred;

  for (i in 1:N_ref)
    pos_pred[i] = binomial_rng(tested_ref[i], pi_ref[i]);

  for (r in 1:R)
    for (t in 1:T)
      cases_pred[r, t] = poisson_rng(mu_cases[r, t]);
}