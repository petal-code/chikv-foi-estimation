# Brazil CHIKV Catalytic Model  — Model B (Literature Sero Data)

# ── Libraries ─────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(cmdstanr)
library(brpop)

cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.38.0")

# ══════════════════════════════════════════════════════════════════
# 1. Compile Stan Model
# ══════════════════════════════════════════════════════════════════

mod <- cmdstan_model("models/chikv_catalytic_lit.stan",
                     force_recompile = TRUE)

# ══════════════════════════════════════════════════════════════════
# 2. Load Literature Sero Data
# ══════════════════════════════════════════════════════════════════

# CSV has 7 columns: Municipality, co_municipio, age_lo, age_hi,
# year_model, pos, tested — one row per study observation.
ref_serodata <- read_csv("data/ref_serodata.csv") %>%
  mutate(co_municipio = as.integer(co_municipio)) %>%
  filter(
    !is.na(co_municipio),
    !is.na(pos),
    !is.na(tested),
    !is.na(age_lo),
    !is.na(age_hi),
    tested > 0
  )

cat("Loaded", nrow(ref_serodata), "literature sero observations\n")
cat("Unique municipalities:", n_distinct(ref_serodata$co_municipio), "\n")

# ══════════════════════════════════════════════════════════════════
# 3. Build Municipality Index
# ══════════════════════════════════════════════════════════════════

# Build Municipality lookup table from whatever is in ref_serodata.
# r = 1, 2, ... R is the Stan index used throughout.
lit_muni_info <- ref_serodata %>%
  distinct(Municipality, co_municipio) %>%
  arrange(co_municipio) %>%
  mutate(r = row_number())

R <- nrow(lit_muni_info)
cat("Number of municipalities (R):", R, "\n")
print(lit_muni_info)

# ══════════════════════════════════════════════════════════════════
# 4. Prepare Sero Observations for Stan
# ══════════════════════════════════════════════════════════════════

# Model years — same range as Model A for comparability
years      <- 2014:2025
year_start <- min(years)
T          <- length(years)

lit_sero <- ref_serodata %>%
  left_join(lit_muni_info %>% select(r, co_municipio), by = "co_municipio") %>%
  mutate(
    # Round fractional years (e.g. 2018.5) to nearest integer
    survey_yr     = as.integer(round(year_model)),
    # Convert calendar year to Stan time index (1 = year_start = 2014)
    survey_yr_idx = survey_yr - year_start + 1L
  ) %>%
  filter(survey_yr_idx >= 1L, survey_yr_idx <= T)

N_ref <- nrow(lit_sero)
cat("Sero observations within model time window (N_ref):", N_ref, "\n")

# ══════════════════════════════════════════════════════════════════
# 5. Case Data (SINAN — all lit municipalities)
# ══════════════════════════════════════════════════════════════════

df_raw_slim <- readRDS("data/case/df_raw_SINAN_CHIKV_2013_2025_slim.rds")

cases_agg <- df_raw_slim %>%
  filter(CLASSI_FIN == "13", CRITERIO %in% c("1", "2")) %>%
  mutate(
    co_municipio = as.integer(ID_MUNICIP),
    year         = year(DT_NOTIFIC)
  ) %>%
  filter(
    co_municipio %in% lit_muni_info$co_municipio,
    !is.na(year),
    year >= year_start,
    year <= max(years)
  ) %>%
  group_by(co_municipio, year) %>%
  summarise(cases = n(), .groups = "drop")

cat("Total notified cases across lit municipalities:", sum(cases_agg$cases), "\n")

# ══════════════════════════════════════════════════════════════════
# 6. Population Data (2014, all lit municipalities)
# ══════════════════════════════════════════════════════════════════

popdata <- readRDS("data/pop/pop_all_raw_datasus.rds")

N_init_df <- popdata %>%
  mutate(code_muni = as.integer(code_muni)) %>%
  filter(
    code_muni %in% lit_muni_info$co_municipio,
    age_group != "Total",
    year == 2014
  ) %>%
  group_by(co_municipio = code_muni) %>%
  summarise(N_init = sum(pop), .groups = "drop")

# Warn if any Municipality is missing population data
missing_pop <- setdiff(lit_muni_info$co_municipio, N_init_df$co_municipio)
if (length(missing_pop) > 0) {
  warning("Missing 2014 population for co_municipio: ",
          paste(missing_pop, collapse = ", "))
}

# Build N_init vector ordered by r
N_init_vec <- lit_muni_info %>%
  left_join(N_init_df, by = "co_municipio") %>%
  arrange(r) %>%
  pull(N_init) %>%
  as.numeric()

cat("N_init (2014 population):\n")
print(lit_muni_info %>% mutate(N_init = N_init_vec) %>%
        select(r, Municipality, N_init))

# ══════════════════════════════════════════════════════════════════
# 7. Build Stan Data Arrays
# ══════════════════════════════════════════════════════════════════

# Case array [R, T] — zero-filled, then populated from cases_agg
cases_arr <- array(0L, dim = c(R, T))

for (i in seq_len(nrow(cases_agg))) {
  r <- match(cases_agg$co_municipio[i], lit_muni_info$co_municipio)
  t <- match(cases_agg$year[i],         years)
  if (!anyNA(c(r, t)))
    cases_arr[r, t] <- cases_agg$cases[i]
}

stan_data <- list(
  R          = R,
  T          = T,
  N_ref      = N_ref,
  year_start = year_start,
  
  r_ref             = as.array(lit_sero$r),
  age_lo_ref        = as.array(as.integer(lit_sero$age_lo)),
  age_hi_ref        = as.array(as.integer(lit_sero$age_hi)),
  survey_yr_ref_idx = as.array(lit_sero$survey_yr_idx),
  pos_ref           = as.array(as.integer(lit_sero$pos)),
  tested_ref        = as.array(as.integer(lit_sero$tested)),
  
  cases  = cases_arr,
  N_init = N_init_vec
)

cat("Stan data summary:\n")
cat("  R:", R, "| T:", T, "| N_ref:", N_ref, "\n")
cat("  Total tested (sero):", sum(lit_sero$tested),
    "| Positive:", sum(lit_sero$pos), "\n")
cat("  Total cases:", sum(cases_arr), "\n")

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

# -- 9a. Underreporting factor (rho) -------------------------------
rho_summary <- fit$summary("rho") %>%
  mutate(r = row_number()) %>%
  left_join(lit_muni_info, by = "r") %>%
  select(Municipality, mean, sd, q5, q95, rhat, ess_bulk)

cat("=== Underreporting factor (rho) ===\n")
print(rho_summary)

# -- 9b. FOI time series (lambda) ----------------------------------
lambda_summary <- fit$summary("lambda") %>%
  mutate(
    idx  = str_match(variable, "\\[(\\d+),(\\d+)\\]"),
    r    = as.integer(idx[, 2]),
    t    = as.integer(idx[, 3]),
    year = years[t]
  ) %>%
  left_join(lit_muni_info, by = "r") %>%
  select(Municipality, year, mean, sd, q5, q95, rhat, ess_bulk)

cat("High-Rhat lambda (rhat > 1.05):", sum(lambda_summary$rhat > 1.05), "\n")

lambda_summary %>%
  ggplot(aes(x = year, y = mean, ymin = q5, ymax = q95)) +
  geom_ribbon(alpha = 0.3) +
  geom_line() +
  facet_wrap(~ Municipality, scales = "free_y") +
  labs(title = "Model B — Force of Infection by Municipality (literature sero)",
       y = "FOI (lambda)", x = NULL) +
  theme_bw()

# ══════════════════════════════════════════════════════════════════
# 10. Posterior Predictive Checks
# ══════════════════════════════════════════════════════════════════

# -- 10a. Seroprevalence PPC ---------------------------------------
sero_pred_summary <- fit$summary("pos_pred") %>%
  mutate(i = row_number()) %>%
  left_join(
    lit_sero %>%
      mutate(i = row_number()) %>%
      left_join(lit_muni_info, by = "r"),
    by = "i"
  ) %>%
  mutate(
    obs_prev  = pos / tested,
    pred_mean = mean / tested,
    pred_q5   = q5  / tested,
    pred_q95  = q95 / tested
  )

ggplot(sero_pred_summary, aes(x = survey_yr, y = obs_prev)) +
  geom_pointrange(aes(y = pred_mean, ymin = pred_q5, ymax = pred_q95),
                  color = "tomato", alpha = 0.7) +
  geom_point(color = "black", size = 2) +
  facet_wrap(~ Municipality, scales = "free_y") +
  labs(title = "Model B — Seroprevalence PPC: observed (black) vs predicted (red)",
       x = "Survey year", y = "Seroprevalence") +
  theme_bw()

# -- 10b. Incidence PPC (per 100k) ---------------------------------
n_init_tbl <- tibble(r = 1:R, N_init = N_init_vec)

incidence_pred_summary <- fit$summary("cases_pred") %>%
  mutate(
    idx  = str_match(variable, "\\[(\\d+),(\\d+)\\]"),
    r    = as.integer(idx[, 2]),
    t    = as.integer(idx[, 3]),
    year = years[t]
  ) %>%
  left_join(n_init_tbl, by = "r") %>%
  left_join(lit_muni_info, by = "r") %>%
  mutate(
    pred_mean = mean / N_init * 1e5,
    pred_q5   = q5   / N_init * 1e5,
    pred_q95  = q95  / N_init * 1e5
  )

incidence_obs <- cases_agg %>%
  left_join(lit_muni_info %>% select(co_municipio, Municipality),
            by = "co_municipio") %>%
  left_join(n_init_tbl %>%
              left_join(lit_muni_info %>% select(r, co_municipio), by = "r"),
            by = "co_municipio") %>%
  mutate(incidence = cases / N_init * 1e5)

ggplot(incidence_obs, aes(x = year)) +
  geom_ribbon(data = incidence_pred_summary,
              aes(ymin = pred_q5, ymax = pred_q95),
              alpha = 0.3, fill = "steelblue") +
  geom_line(data = incidence_pred_summary,
            aes(y = pred_mean), color = "steelblue") +
  geom_point(aes(y = incidence), color = "black", size = 1.5) +
  facet_wrap(~ Municipality, scales = "free_y") +
  labs(title = "Model B — Incidence PPC: observed (black) vs predicted (blue)",
       y = "Cases per 100,000", x = NULL) +
  theme_bw()

# -- 10c. Underreporting factor forest plot ------------------------
rho_summary %>%
  mutate(Municipality = fct_reorder(Municipality, mean)) %>%
  ggplot(aes(x = Municipality, y = mean, ymin = q5, ymax = q95)) +
  geom_pointrange(color = "darkgreen", size = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     limits = c(0, NA)) +
  coord_flip() +
  labs(title = "Model B — Underreporting factor by Municipality",
       x = NULL, y = "Underreporting factor (rho)") +
  theme_bw()