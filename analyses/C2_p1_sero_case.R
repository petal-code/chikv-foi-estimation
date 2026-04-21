# ══════════════════════════════════════════════════════════════════
# Brazil CHIKV Catalytic Model  (7 sero municipalities only)
# FOI      : municipality x year  (R = 7)
# Reporting: municipality-level underreporting factor (rho)
# Cases    : Poisson(mu_cases)
# Serology : age-specific binomial likelihood (1-yr age bins)
# S[r,1]   = N_init[r]  (2014 population)
# S[r,t+1] = S[r,t] * exp(-lambda[r,t])
# Data     : Serology (2023, 2024, 2025) + Notified cases (2014-2025)
# ══════════════════════════════════════════════════════════════════

# ── Libraries ─────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(cmdstanr)
library(brpop)

# ══════════════════════════════════════════════════════════════════
# 1. Stan Model
# ══════════════════════════════════════════════════════════════════

stan_code <- "
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
"

dir.create("models", showWarnings = FALSE)
writeLines(stan_code, "models/chikv_catalytic_7muni.stan")
mod <- cmdstan_model("models/chikv_catalytic_7muni.stan", force_recompile = TRUE)

# ══════════════════════════════════════════════════════════════════
# 2. Helper Functions
# ══════════════════════════════════════════════════════════════════

parse_age_years <- function(x) {
  unit  <- as.integer(substr(as.character(x), 1, 1))
  value <- as.integer(substr(as.character(x), 2, 4))
  case_when(
    unit == 4 ~ as.numeric(value),
    unit == 3 ~ as.numeric(value) / 12,
    unit == 2 ~ as.numeric(value) / 365,
    TRUE      ~ NA_real_
  )
}

# ══════════════════════════════════════════════════════════════════
# 3. Reference Tables
# ══════════════════════════════════════════════════════════════════

estado_map <- tibble(
  sg_uf      = c("13", "23", "26", "31", "33", "35", "41"),
  co_estado  = c("AM", "CE", "PE", "MG", "RJ", "SP", "PR"),
  estado_idx = 1:7
)

sero_muni_info <- tibble(
  r            = 1:7,
  region       = c("BeloHorizonte", "Curitiba", "Fortaleza",
                   "Manaus", "Recife", "RiodeJaneiro", "SãoPaulo"),
  co_municipio = c(310620, 410690, 230440, 130260, 261160, 330455, 355030),
  co_estado    = c("MG",   "PR",   "CE",   "AM",   "PE",   "RJ",   "SP"),
  sg_uf        = c("31",   "41",   "23",   "13",   "26",   "33",   "35")
)

# ══════════════════════════════════════════════════════════════════
# 4. Serology Data (2023, 2024, 2025)
# ══════════════════════════════════════════════════════════════════

load("data/cohort/Brazil_CHIKV serology survey_Nov 2023.RData")
load("data/cohort/Brazil_CHIKV serology survey_Jun 2024.RData")
load("data/cohort/Brazil_CHIKV serology survey_Jun 2025.RData")

parse_sero_df <- function(df, municipality, survey_year) {
  df <- df %>% rename_with(~ gsub("CHKG[-.]1_", "CHKG1_", .x))
  df %>%
    transmute(
      region      = municipality,
      survey_year = survey_year,
      age_int     = as.integer(floor(
        as.numeric(difftime(as.Date(DONATION_DATE),
                            as.Date(BIRTH_DATE),
                            units = "days")) / 365.25
      )),
      sero_pos = as.integer(CHKG1_INTERPRETATION == "Positivo")
    ) %>%
    filter(!is.na(age_int), !is.na(sero_pos), age_int >= 0)
}

sero_individual <- bind_rows(
  map_dfr(sero_muni_info$region,
          ~ parse_sero_df(get(paste0("df_", .x, "_Nov23")), .x, 2023)),
  map_dfr(sero_muni_info$region,
          ~ parse_sero_df(get(paste0("df_", .x, "_Jun24")), .x, 2024)),
  map_dfr(sero_muni_info$region,
          ~ parse_sero_df(get(paste0("df_", .x, "_Jun25")), .x, 2025))
)

# Age-aggregated for PPC plotting
sero_agg_age <- sero_individual %>%
  group_by(region, survey_year, age_int) %>%
  summarise(pos = sum(sero_pos), tested = n(), .groups = "drop")

# Municipality x survey total (for reference)
sero_agg_total <- sero_individual %>%
  group_by(region, survey_year) %>%
  summarise(
    pos_total    = sum(sero_pos),
    tested_total = n(),
    obs_prev     = pos_total / tested_total,
    .groups      = "drop"
  )

cat("Sero summary (total per municipality x survey):\n")
print(sero_agg_total)

# ══════════════════════════════════════════════════════════════════
# 5. Case Data  (7 municipalities only)
# ══════════════════════════════════════════════════════════════════

df_raw_slim <- readRDS("data/case/df_raw_SINAN_CHIKV_2013_2025_slim.rds")

cases_agg <- df_raw_slim %>%
  filter(CLASSI_FIN == "13", CRITERIO %in% c("1", "2")) %>%
  mutate(
    co_municipio = as.integer(ID_MUNICIP),
    year         = year(DT_NOTIFIC)
  ) %>%
  filter(
    co_municipio %in% sero_muni_info$co_municipio,
    !is.na(year), year >= 2014, year <= 2025
  ) %>%
  group_by(co_municipio, year) %>%
  summarise(cases = n(), .groups = "drop")

cat("Cases total:", sum(cases_agg$cases), "\n")

# ══════════════════════════════════════════════════════════════════
# 6. Population Data  (2014 total population, 7 municipalities)
# ══════════════════════════════════════════════════════════════════

popdata <- readRDS("data/pop/pop_all_raw_datasus.rds")

N_init_df <- popdata %>%
  mutate(code_muni = as.integer(code_muni)) %>%
  filter(
    code_muni %in% sero_muni_info$co_municipio,
    age_group != "Total",
    year == 2014
  ) %>%
  group_by(co_municipio = code_muni) %>%
  summarise(N_init = sum(pop), .groups = "drop")

stopifnot(all(sero_muni_info$co_municipio %in% N_init_df$co_municipio))

N_init_vec <- sero_muni_info %>%
  left_join(N_init_df, by = "co_municipio") %>%
  arrange(r) %>%
  pull(N_init) %>%
  as.numeric()

cat("N_init (2014 population):\n")
print(sero_muni_info %>%
        mutate(N_init = N_init_vec) %>%
        select(region, co_estado, N_init))

# ══════════════════════════════════════════════════════════════════
# 7. Build Stan Data Arrays
# ══════════════════════════════════════════════════════════════════

years        <- 2014:2025
survey_years <- c(2023, 2024, 2025)
ages_sero    <- sort(unique(sero_agg_age$age_int))

R        <- 7L
T        <- length(years)
T_survey <- length(survey_years)
A_sero   <- length(ages_sero)

cat("R:", R, "| T:", T, "| T_survey:", T_survey, "| A_sero:", A_sero, "\n")

# Serology arrays [R, A_sero, T_survey]
pos_arr    <- array(0L, dim = c(R, A_sero, T_survey))
tested_arr <- array(0L, dim = c(R, A_sero, T_survey))

for (i in seq_len(nrow(sero_agg_age))) {
  r <- match(sero_agg_age$region[i],      sero_muni_info$region)
  a <- match(sero_agg_age$age_int[i],     ages_sero)
  s <- match(sero_agg_age$survey_year[i], survey_years)
  if (!anyNA(c(r, a, s))) {
    pos_arr[r, a, s]    <- sero_agg_age$pos[i]
    tested_arr[r, a, s] <- sero_agg_age$tested[i]
  }
}

# Case array [R, T]
cases_arr <- array(0L, dim = c(R, T))

for (i in seq_len(nrow(cases_agg))) {
  r <- match(cases_agg$co_municipio[i], sero_muni_info$co_municipio)
  t <- match(cases_agg$year[i],         years)
  if (!anyNA(c(r, t))) {
    cases_arr[r, t] <- cases_agg$cases[i]
  }
}

stan_data <- list(
  R               = R,
  T               = T,
  T_survey        = T_survey,
  A_sero          = A_sero,
  year_start      = min(years),  
  survey_year_idx = as.array(match(survey_years, years)),
  ages_sero       = as.array(ages_sero),
  pos             = pos_arr,
  tested          = tested_arr,
  cases           = cases_arr,
  N_init          = N_init_vec
)

cat("survey_year_idx:", stan_data$survey_year_idx, "\n")
cat("Total tested:", sum(tested_arr), "| Positive:", sum(pos_arr), "\n")
cat("Total cases:", sum(cases_arr), "\n")
cat("N_init:", round(N_init_vec), "\n")

# ══════════════════════════════════════════════════════════════════
# 8. MCMC Sampling
# ══════════════════════════════════════════════════════════════════

init_fn <- function() {
  list(
    lambda = matrix(0.05, nrow = R, ncol = T),
    rho    = rep(0.05, R)
  )
}

fit <- mod$sample(
  data            = stan_data,
  seed            = 42,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 100,
  init            = init_fn
)

fit$diagnostic_summary()

# ══════════════════════════════════════════════════════════════════
# 9. Parameter Summaries
# ══════════════════════════════════════════════════════════════════

# N_init lookup for PPC
n_init_tbl <- tibble(r = 1:R, N_init = N_init_vec)

# -- 9a. Underreporting factor -------------------------------------------
rho_summary <- fit$summary("rho") %>%
  mutate(r = row_number(), region = sero_muni_info$region[r]) %>%
  select(region, mean, sd, q5, q95, rhat, ess_bulk)

cat("=== Underreporting factor (rho) ===\n")
print(rho_summary)

# -- 9b. FOI time series ------------------------------------------
lambda_summary <- fit$summary("lambda") %>%
  mutate(
    idx    = str_match(variable, "\\[(\\d+),(\\d+)\\]"),
    r      = as.integer(idx[, 2]),
    t      = as.integer(idx[, 3]),
    region = sero_muni_info$region[r],
    year   = years[t]
  ) %>%
  select(region, year, mean, sd, q5, q95, rhat, ess_bulk)

cat("High-Rhat lambda (rhat > 1.05):", sum(lambda_summary$rhat > 1.05), "\n")

lambda_summary %>%
  ggplot(aes(x = year, y = mean, ymin = q5, ymax = q95)) +
  geom_ribbon(alpha = 0.3) +
  geom_line() +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Force of Infection by municipality",
       y = "FOI (lambda)", x = NULL) +
  theme_bw()

# ══════════════════════════════════════════════════════════════════
# 10. Posterior Predictive Checks
# ══════════════════════════════════════════════════════════════════

parse_idx2 <- function(summary_df) {
  idx <- str_match(summary_df$variable, "\\[(\\d+),(\\d+)\\]")
  summary_df %>%
    mutate(i1 = as.integer(idx[, 2]),
           i2 = as.integer(idx[, 3]))
}

parse_idx3 <- function(summary_df) {
  idx <- str_match(summary_df$variable, "\\[(\\d+),(\\d+),(\\d+)\\]")
  summary_df %>%
    mutate(i1 = as.integer(idx[, 2]),
           i2 = as.integer(idx[, 3]),
           i3 = as.integer(idx[, 4]))
}

# -- 10a. Seroprevalence PPC: full age range 0-80 (primary) -------
# pi_sero_full[r, age+1, s]: index 1 = age 0, index 81 = age 80
sero_pred_full_summary <- fit$summary("pi_sero_full") %>%
  parse_idx3() %>%
  mutate(
    region      = sero_muni_info$region[i1],
    age_int     = i2 - 1L,        
    survey_year = survey_years[i3]
  )

sero_obs_age <- sero_agg_age %>%
  mutate(obs_prev = pos / tested)

ggplot(sero_pred_full_summary, aes(x = age_int)) +
  geom_ribbon(aes(ymin = q5, ymax = q95),
              alpha = 0.3, fill = "tomato") +
  geom_line(aes(y = mean), color = "tomato") +
  geom_point(data = sero_obs_age,
             aes(y = obs_prev), size = 0.8, alpha = 0.6) +
  facet_grid(survey_year ~ region) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  labs(title = "Seroprevalence PPC: observed (black) vs predicted (red)",
       x = "Age", y = "Seroprevalence") +
  theme_bw()

# -- 10b. Incidence PPC (per 100k) --------------------------------
incidence_pred_summary <- fit$summary("cases_pred") %>%
  parse_idx2() %>%
  mutate(region = sero_muni_info$region[i1], year = years[i2]) %>%
  left_join(n_init_tbl, by = c("i1" = "r")) %>%
  mutate(
    pred_mean = mean / N_init * 100000,
    pred_q5   = q5   / N_init * 100000,
    pred_q95  = q95  / N_init * 100000
  )

incidence_obs <- cases_agg %>%
  left_join(sero_muni_info %>% select(co_municipio, region), by = "co_municipio") %>%
  left_join(tibble(co_municipio = sero_muni_info$co_municipio,
                   N_init = N_init_vec), by = "co_municipio") %>%
  mutate(incidence = cases / N_init * 100000)

ggplot(incidence_obs, aes(x = year)) +
  geom_ribbon(data = incidence_pred_summary,
              aes(ymin = pred_q5, ymax = pred_q95),
              alpha = 0.3, fill = "steelblue") +
  geom_line(data = incidence_pred_summary,
            aes(y = pred_mean), color = "steelblue") +
  geom_point(aes(y = incidence), color = "black", size = 1.5) +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Incidence PPC: observed (black) vs predicted (blue)",
       y = "Cases per 100,000 population", x = NULL) +
  theme_bw()

# -- 10c. Underreporting factor (rho) ----------------------------
rho_summary %>%
  mutate(region = fct_reorder(region, mean)) %>%
  ggplot(aes(x = region, y = mean, ymin = q5, ymax = q95)) +
  geom_pointrange(color = "darkgreen", size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     limits = c(0, NA)) +
  coord_flip() +
  labs(title = "Underreporting factor by municipality",
       x = NULL, y = "Underreporting factor") +
  theme_bw()