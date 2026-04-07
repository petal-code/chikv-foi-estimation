# ══════════════════════════════════════════════════════════════════
# Brazil CHIKV Catalytic Model
# Estimating Force of Infection (FOI) and Infection-to-Reporting Ratio
# Data: Serology (2023, 2024) + Notified cases (2014-2025)
# ══════════════════════════════════════════════════════════════════

# ── Libraries ─────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(cmdstanr)
# remotes::install_github("rfsaldanha/brpop")
library(brpop)

# ══════════════════════════════════════════════════════════════════
# 1. Stan Model (embedded as string)
# ══════════════════════════════════════════════════════════════════

stan_code <- "
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
      lambda[r, t] ~ exponential(2);   //

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
"

# Write Stan model to file and compile
writeLines(stan_code, "models/chikv_catalytic.stan")
mod <- cmdstan_model("models/chikv_catalytic.stan")

# ══════════════════════════════════════════════════════════════════
# 2. Serology Data
# ══════════════════════════════════════════════════════════════════

load("data/cohort/Brazil_CHIKV serology survey_Jun 2024.RData")
load("data/cohort/Brazil_CHIKV serology survey_Nov 2023.RData")

municipalities <- c("BeloHorizonte", "Curitiba", "Fortaleza",
                    "Manaus", "Recife", "RiodeJaneiro", "SãoPaulo")

parse_sero_df <- function(df, municipality, survey_year) {
  # Normalize column names: CHKG.1_ or CHKG-1_ -> CHKG1_
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
  map_dfr(municipalities, ~ parse_sero_df(get(paste0("df_", .x, "_Jun24")), .x, 2024)),
  map_dfr(municipalities, ~ parse_sero_df(get(paste0("df_", .x, "_Nov23")), .x, 2023))
)

# Aggregate to region x survey_year x age_int level
sero_agg <- sero_individual %>%
  group_by(region, survey_year, age_int) %>%
  summarise(pos = sum(sero_pos), tested = n(), .groups = "drop")

# ══════════════════════════════════════════════════════════════════
# 3. Case Data
# ══════════════════════════════════════════════════════════════════

load("data/case/bd_CHIKV_cases_2014-2025 assigned week.RData")

municipality_map <- tibble(
  co_municipio = c(310620, 410690, 230440, 130260, 261160, 330455, 355030),
  co_estado    = c("MG",   "PR",   "CE",   "AM",   "PE",   "RJ",   "SP"),
  region       = c("BeloHorizonte", "Curitiba", "Fortaleza",
                   "Manaus", "Recife", "RiodeJaneiro", "SãoPaulo")
)

# 5-year age group boundaries and labels
age_breaks <- c(0, 5, 10, 15, 20, 25, 30, 35, 40,
                45, 50, 55, 60, 65, 70, 75, 80, Inf)
age_labels <- c(0, 5, 10, 15, 20, 25, 30, 35, 40,
                45, 50, 55, 60, 65, 70, 75, 80)

cases_agg <- df_CHIKV_full %>%
  inner_join(municipality_map, by = c("co_municipio", "co_estado")) %>%
  mutate(
    year      = year(dt_coleta),
    age_int   = as.integer(floor(age)),
    age_lower = age_labels[findInterval(age_int, age_breaks,
                                        rightmost.closed = TRUE)]
  ) %>%
  filter(!is.na(age_int), !is.na(year), year >= 2014, year <= 2025) %>%
  group_by(region, year, age_lower) %>%
  summarise(cases = n(), .groups = "drop")

# ══════════════════════════════════════════════════════════════════
# 4. Population Data
# ══════════════════════════════════════════════════════════════════

popdata <- mun_pop_age(source = "datasus", sex = "all")

parse_age_group <- function(ag) {
  case_when(
    ag == "Total"                 ~ NA_integer_,
    ag == "From 80 years or more" ~ 80L,
    TRUE ~ as.integer(str_extract(ag, "^From (\\d+)", group = 1))
  )
}

pop_processed <- popdata %>%
  filter(code_muni %in% municipality_map$co_municipio) %>%
  left_join(municipality_map %>% select(co_municipio, region),
            by = c("code_muni" = "co_municipio")) %>%
  mutate(age_lower = parse_age_group(age_group)) %>%
  filter(!is.na(age_lower)) %>%
  select(region, year, age_lower, pop)

# Linear extrapolation for 2022-2025 (data only available up to 2021)
pop_extrap <- pop_processed %>%
  filter(year %in% 2019:2021) %>%
  group_by(region, age_lower) %>%
  arrange(year) %>%
  summarise(
    slope    = (pop[year == 2021] - pop[year == 2019]) / 2,
    pop_2021 = pop[year == 2021],
    .groups  = "drop"
  ) %>%
  crossing(year = 2022:2025) %>%
  mutate(pop = pmax(0L, as.integer(pop_2021 + slope * (year - 2021)))) %>%
  select(region, year, age_lower, pop)

pop_full <- bind_rows(
  pop_processed %>% filter(year >= 2014),
  pop_extrap
) %>%
  arrange(region, year, age_lower)

# ══════════════════════════════════════════════════════════════════
# 5. Validation
# ══════════════════════════════════════════════════════════════════

stopifnot(all(unique(cases_agg$age_lower) %in% unique(pop_full$age_lower)))

cat("Serology:", nrow(sero_agg), "rows |",
    n_distinct(sero_agg$region), "municipalities |",
    n_distinct(sero_agg$survey_year), "time points\n")
cat("Cases:", nrow(cases_agg), "rows |",
    n_distinct(cases_agg$region), "municipalities |",
    n_distinct(cases_agg$year), "years\n")
cat("Population:", nrow(pop_full), "rows |",
    n_distinct(pop_full$year), "years\n")

# ══════════════════════════════════════════════════════════════════
# 6. Build Stan Data List
# ══════════════════════════════════════════════════════════════════

regions      <- sort(unique(cases_agg$region))
years        <- 2014:2025
survey_years <- c(2023, 2024)
ages_sero    <- sort(unique(sero_agg$age_int))    # 1-yr intervals
ages_case    <- sort(unique(cases_agg$age_lower)) # 5-yr intervals

R        <- length(regions)
T        <- length(years)
T_survey <- length(survey_years)
A_sero   <- length(ages_sero)
A_case   <- length(ages_case)

cat("R:", R, "| T:", T, "| T_survey:", T_survey,
    "| A_sero:", A_sero, "| A_case:", A_case, "\n")

# Serology arrays [R, A_sero, T_survey]
pos_arr    <- array(0L, dim = c(R, A_sero, T_survey))
tested_arr <- array(0L, dim = c(R, A_sero, T_survey))

for (i in seq_len(nrow(sero_agg))) {
  r <- match(sero_agg$region[i],      regions)
  a <- match(sero_agg$age_int[i],     ages_sero)
  s <- match(sero_agg$survey_year[i], survey_years)
  if (!anyNA(c(r, a, s))) {
    pos_arr[r, a, s]    <- sero_agg$pos[i]
    tested_arr[r, a, s] <- sero_agg$tested[i]
  }
}

# Case array [R, A_case, T]
cases_arr <- array(0L, dim = c(R, A_case, T))

for (i in seq_len(nrow(cases_agg))) {
  r <- match(cases_agg$region[i],    regions)
  a <- match(cases_agg$age_lower[i], ages_case)
  t <- match(cases_agg$year[i],      years)
  if (!anyNA(c(r, a, t))) {
    cases_arr[r, a, t] <- cases_agg$cases[i]
  }
}

# Population array [R, A_case, T]
N_arr <- array(0.0, dim = c(R, A_case, T))

for (i in seq_len(nrow(pop_full))) {
  r <- match(pop_full$region[i],    regions)
  a <- match(pop_full$age_lower[i], ages_case)
  t <- match(pop_full$year[i],      years)
  if (!anyNA(c(r, a, t))) {
    N_arr[r, a, t] <- pop_full$pop[i]
  }
}

stan_data <- list(
  R               = R,
  T               = T,
  T_survey        = T_survey,
  A_sero          = A_sero,
  A_case          = A_case,
  ages_sero       = as.array(ages_sero),
  ages_case       = as.array(ages_case),
  survey_year_idx = as.array(match(survey_years, years)),  # {10, 11}
  pos             = pos_arr,
  tested          = tested_arr,
  cases           = cases_arr,
  N               = N_arr
)

# Quick sanity checks
cat("survey_year_idx:", stan_data$survey_year_idx, "\n")
cat("Total tested:", sum(tested_arr), "| Positive:", sum(pos_arr), "\n")
cat("Total cases:", sum(cases_arr), "\n")
cat("Population NA:", sum(is.na(N_arr)), "| Zero:", sum(N_arr == 0), "\n")

# ══════════════════════════════════════════════════════════════════
# 7. MCMC Sampling
# ══════════════════════════════════════════════════════════════════

fit <- mod$sample(
  data            = stan_data,
  seed            = 42,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 500,
  iter_sampling   = 500,
  refresh         = 100
)

fit$diagnostic_summary()

# ══════════════════════════════════════════════════════════════════
# 8. Parameter Summaries
# ══════════════════════════════════════════════════════════════════

# Reporting ratio by region
fit$summary("rho") %>%
  mutate(region = regions) %>%
  select(region, mean, sd, q5, q95, rhat, ess_bulk)

# NegBin overdispersion
fit$summary("phi")

# FOI summary
lambda_summary <- fit$summary("lambda") %>%
  mutate(
    r      = rep(1:R, times = T),
    t      = rep(1:T, each  = R),
    region = regions[r],
    year   = years[t]
  ) %>%
  select(region, year, mean, sd, q5, q95, rhat, ess_bulk)

# Check convergence
lambda_summary %>% filter(rhat > 1.05)

# ── FOI time series plot ──────────────────────────────────────────
lambda_summary %>%
  ggplot(aes(x = year, y = mean, ymin = q5, ymax = q95)) +
  geom_ribbon(alpha = 0.3) +
  geom_line() +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Force of Infection by municipality and year",
       y = "FOI (lambda)", x = NULL) +
  theme_bw()

# ══════════════════════════════════════════════════════════════════
# 9. Posterior Predictive Checks
# ══════════════════════════════════════════════════════════════════

# Helper: parse 3-index Stan variable names [i1, i2, i3]
parse_idx3 <- function(summary_df) {
  idx <- str_match(summary_df$variable, "\\[(\\d+),(\\d+),(\\d+)\\]")
  summary_df %>%
    mutate(
      i1 = as.integer(idx[, 2]),
      i2 = as.integer(idx[, 3]),
      i3 = as.integer(idx[, 4])
    )
}

# ── 9a. Cases PPC ─────────────────────────────────────────────────
cases_pred_summary <- fit$summary("cases_pred") %>%
  parse_idx3() %>%
  mutate(
    region    = regions[i1],
    age_lower = ages_case[i2],
    year      = years[i3]
  )

# Total cases per region x year
cases_fit <- cases_agg %>%
  group_by(region, year) %>%
  summarise(obs = sum(cases), .groups = "drop") %>%
  left_join(
    cases_pred_summary %>%
      group_by(region, year) %>%
      summarise(pred_mean = sum(mean),
                pred_q5   = sum(q5),
                pred_q95  = sum(q95),
                .groups   = "drop"),
    by = c("region", "year")
  )

ggplot(cases_fit, aes(x = year)) +
  geom_ribbon(aes(ymin = pred_q5, ymax = pred_q95),
              alpha = 0.3, fill = "steelblue") +
  geom_line(aes(y = pred_mean), color = "steelblue") +
  geom_point(aes(y = obs), color = "black", size = 1.5) +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Cases: observed (black) vs predicted (blue)",
       y = "Cases", x = NULL) +
  theme_bw()

# Cases by age group (per municipality)
cases_fit_age <- cases_agg %>%
  rename(obs = cases) %>%
  left_join(
    cases_pred_summary %>%
      group_by(region, year, age_lower) %>%
      summarise(pred_mean = sum(mean),
                pred_q5   = sum(q5),
                pred_q95  = sum(q95),
                .groups   = "drop"),
    by = c("region", "year", "age_lower")
  )

plot_cases_age <- function(target_region) {
  
  full_grid <- expand_grid(
    region    = target_region,
    year      = years,
    age_lower = ages_case
  )
  
  pred <- cases_pred_summary %>%
    filter(region == target_region) %>%
    group_by(region, year, age_lower) %>%
    summarise(pred_mean = sum(mean),
              pred_q5   = sum(q5),
              pred_q95  = sum(q95),
              .groups   = "drop")
  
  obs <- cases_agg %>%
    filter(region == target_region) %>%
    rename(obs = cases)
  
  plot_df <- full_grid %>%
    left_join(pred, by = c("region", "year", "age_lower")) %>%
    left_join(obs,  by = c("region", "year", "age_lower"))
  
  ggplot(plot_df, aes(x = year)) +
    geom_ribbon(aes(ymin = pred_q5, ymax = pred_q95),
                alpha = 0.3, fill = "steelblue") +
    geom_line(aes(y = pred_mean), color = "steelblue") +
    geom_point(data = plot_df %>% filter(!is.na(obs)),
               aes(y = obs), color = "black", size = 1.5) +
    facet_wrap(~ age_lower, scales = "free_y",
               labeller = label_bquote(.(age_lower) ~ "yr+")) +
    labs(title = paste(target_region, "- Cases by age group"),
         y = "Cases", x = NULL) +
    theme_bw()
}


for (r in regions) print(plot_cases_age(r))

# ── 9b. Seroprevalence PPC (ages 0-80, computed from posterior draws) ──
lambda_draws <- fit$draws("lambda", format = "matrix")

lambda_cols <- fit$summary("lambda") %>%
  mutate(
    idx = str_match(variable, "\\[(\\d+),(\\d+)\\]"),
    r   = as.integer(idx[, 2]),
    t   = as.integer(idx[, 3])
  ) %>%
  select(variable, r, t)

ages_full         <- 0:80
survey_years_calc <- c(2023, 2024)
n_draws           <- nrow(lambda_draws)

# Compute posterior pi_sero with CI for all ages 0-80
pi_full_ci <- expand_grid(
  region      = regions,
  age_int     = ages_full,
  survey_year = survey_years_calc
) %>%
  mutate(
    r  = match(region, regions),
    ts = match(survey_year, years),   # 1-based index
    draws_cum = pmap(list(r, age_int, ts), function(r_i, a_i, ts_i) {
      t_start <- max(1L, ts_i - a_i)
      cols <- lambda_cols %>%
        filter(r == r_i, t >= t_start, t <= ts_i) %>%
        pull(variable)
      if (length(cols) == 0) return(rep(0, n_draws))
      rowSums(lambda_draws[, cols, drop = FALSE])
    }),
    pi_draws = map(draws_cum, ~ 1 - exp(-.x)),
    pi_mean  = map_dbl(pi_draws, mean),
    pi_q5    = map_dbl(pi_draws, ~ quantile(.x, 0.05)),
    pi_q95   = map_dbl(pi_draws, ~ quantile(.x, 0.95))
  ) %>%
  select(-draws_cum, -pi_draws)

sero_obs <- sero_agg %>% mutate(obs_prev = pos / tested)

ggplot(pi_full_ci, aes(x = age_int)) +
  geom_ribbon(aes(ymin = pi_q5, ymax = pi_q95),
              alpha = 0.3, fill = "tomato") +
  geom_line(aes(y = pi_mean), color = "tomato") +
  geom_point(data = sero_obs, aes(y = obs_prev),
             size = 0.8, alpha = 0.6) +
  facet_grid(survey_year ~ region) +
  labs(title = "Seroprevalence: observed (black) vs predicted (red)",
       y = "Seroprevalence", x = "Age") +
  theme_bw()

# ══════════════════════════════════════════════════════════════════
# 10. Population Structure (2014-2025, observed + extrapolated)
# ══════════════════════════════════════════════════════════════════

# Flag which years are observed vs extrapolated
pop_plot <- pop_full %>%
  mutate(
    source     = if_else(year <= 2021, "Observed", "Extrapolated"),
    age_label  = paste0(age_lower, "-", age_lower + 4),
    age_label  = if_else(age_lower == 80, "80+", age_label)
  )

# ── 10a. Total population trajectory by municipality ─────────────────────
pop_plot %>%
  group_by(region, year, source) %>%
  summarise(total_pop = sum(pop), .groups = "drop") %>%
  ggplot(aes(x = year, y = total_pop / 1e6, color = source, group = region)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 2021.5, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Observed" = "steelblue",
                                "Extrapolated" = "tomato")) +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Total population by municipality (2014-2025)",
       subtitle = "Dashed line: observed / extrapolated boundary",
       y = "Population (millions)", x = NULL, color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

# ── 10b. Age structure by municipality (selected years) ───────────────────
pop_plot %>%
  filter(year %in% c(2014, 2017, 2020, 2023, 2025)) %>%
  ggplot(aes(x = age_lower, y = pop / 1e3,
             color = factor(year), linetype = source)) +
  geom_line(linewidth = 0.7) +
  scale_linetype_manual(values = c("Observed" = "solid",
                                   "Extrapolated" = "dashed")) +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Age structure by municipality (selected years)",
       y = "Population (thousands)", x = "Age group (lower bound)",
       color = "Year", linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

# ── 10c. Population heatmap: age x year per municipality ─────────────────
for (target_region in regions) {
  p <- pop_plot %>%
    filter(region == target_region) %>%
    ggplot(aes(x = year, y = factor(age_lower), fill = pop / 1e3)) +
    geom_tile() +
    geom_vline(xintercept = 2021.5, linetype = "dashed",
               color = "white", linewidth = 0.8) +
    scale_fill_viridis_c(option = "mako", direction = -1) +
    scale_y_discrete(labels = function(x) {
      x_num <- as.integer(x)
      if_else(x_num == 80, "80+", paste0(x_num, "-", x_num + 4))
    }) +
    labs(title = paste(target_region, "- Population by age group (2014-2025)"),
         subtitle = "Dashed line: observed / extrapolated boundary",
         x = NULL, y = "Age group", fill = "Population\n(thousands)") +
    theme_bw()
  print(p)
}

# ══════════════════════════════════════════════════════════════════
# 11. Cumulative Cases by Age Group and Municipality
# ══════════════════════════════════════════════════════════════════

# ── 11a. Cumulative cases over time by age group (per municipality) ───────
cases_cumulative <- cases_agg %>%
  arrange(region, age_lower, year) %>%
  group_by(region, age_lower) %>%
  mutate(cum_cases = cumsum(cases)) %>%
  ungroup() %>%
  mutate(
    age_label = if_else(age_lower == 80, "80+",
                        paste0(age_lower, "-", age_lower + 4))
  )

plot_cum_cases <- function(target_region) {
  cases_cumulative %>%
    filter(region == target_region) %>%
    ggplot(aes(x = year, y = cum_cases,
               color = factor(age_lower), group = factor(age_lower))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.2) +
    scale_color_viridis_d(option = "turbo",
                          labels = function(x) {
                            x_num <- as.integer(x)
                            if_else(x_num == 80, "80+",
                                    paste0(x_num, "-", x_num + 4))
                          }) +
    labs(title = paste(target_region, "- Cumulative cases by age group"),
         y = "Cumulative cases", x = NULL, color = "Age group") +
    theme_bw() +
    theme(legend.position = "right")
}

for (r in regions) print(plot_cum_cases(r))

# ── 11b. Total cumulative cases by municipality ───────────────────────────
cases_agg %>%
  group_by(region, year) %>%
  summarise(cases = sum(cases), .groups = "drop") %>%
  arrange(region, year) %>%
  group_by(region) %>%
  mutate(cum_cases = cumsum(cases)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = cum_cases / 1e3, color = region)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(title = "Cumulative notified cases by municipality (2014-2025)",
       y = "Cumulative cases (thousands)", x = NULL, color = "Municipality") +
  theme_bw() +
  theme(legend.position = "bottom")

# ── 11c. Age distribution of cumulative cases per municipality (bar) ──────
cases_agg %>%
  group_by(region, age_lower) %>%
  summarise(cum_cases = sum(cases), .groups = "drop") %>%
  mutate(
    age_label = if_else(age_lower == 80, "80+",
                        paste0(age_lower, "-", age_lower + 4)),
    age_label = fct_reorder(age_label, age_lower)
  ) %>%
  ggplot(aes(x = age_label, y = cum_cases, fill = region)) +
  geom_col() +
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Total notified cases by age group (2014-2025)",
       y = "Total cases", x = "Age group") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ── 11d. Cumulative cases by age group: observed + model range ────
# Sum predicted cases across all years per region x age group
cases_pred_cum <- cases_pred_summary %>%
  group_by(region, age_lower) %>%
  summarise(
    pred_mean = sum(mean),
    pred_q5   = sum(q5),
    pred_q95  = sum(q95),
    .groups   = "drop"
  )

# Observed cumulative cases across all years
cases_obs_cum <- cases_agg %>%
  group_by(region, age_lower) %>%
  summarise(obs = sum(cases), .groups = "drop")

# Join observed and predicted
cases_cum_combined <- full_join(cases_obs_cum, cases_pred_cum,
                                by = c("region", "age_lower")) %>%
  mutate(
    age_label = if_else(age_lower == 80, "80+",
                        paste0(age_lower, "-", age_lower + 4)),
    age_label = fct_reorder(age_label, age_lower)
  )

ggplot(cases_cum_combined, aes(x = age_label)) +
  geom_col(aes(y = pred_mean), fill = "steelblue", alpha = 0.5) +
  geom_errorbar(aes(ymin = pred_q5, ymax = pred_q95),
                color = "steelblue", width = 0.4, linewidth = 0.6) +
  geom_point(aes(y = obs), color = "black", size = 1.8) + 
  facet_wrap(~ region, scales = "free_y") +
  labs(title = "Cumulative cases by age group: observed (black) vs predicted (blue)",
       subtitle = "Bar = posterior mean, error bar = 90% CI, dot = observed",
       y = "Cumulative cases", x = "Age group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ══════════════════════════════════════════════════════════════════
# 12. Underreporting factor by Municipality
# ══════════════════════════════════════════════════════════════════

fit$summary("rho") %>%
  mutate(municipality = regions) %>%
  ggplot(aes(x = fct_reorder(municipality, mean), y = mean,
             ymin = q5, ymax = q95)) +
  geom_pointrange(color = "steelblue", linewidth = 0.8, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Underreporting factor by municipality",
       subtitle = "Point = posterior mean, range = 90% CI",
       x = NULL, y = "Underreporting factor") +
  theme_bw()



