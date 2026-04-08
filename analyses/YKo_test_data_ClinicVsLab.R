library(microdatasus)
remotes::install_github("rfsaldanha/microdatasus")
df_raw <- fetch_datasus(
  year_start=  2013,
  year_end = 2025,
  uf = "all",
  information_system = "SINAN-CHIKUNGUNYA"
)

# --- Step 1: Classify cases ---
df_classified <- df_raw |>
  filter(CLASSI_FIN == "13") |>
  mutate(
    case_type = case_when(
      CRITERIO == "1" ~ "laboratory_confirmed",
      CRITERIO == "2" ~ "clinical_confirmed",
      TRUE            ~ NA_character_
    )
  ) |>
  filter(!is.na(case_type))

df_lab      <- df_classified |> filter(case_type == "laboratory_confirmed")
df_clin     <- df_classified |> filter(case_type == "clinical_confirmed")
df_combined <- df_classified  # both groups

# --- Step 2: Parse age from NU_IDADE_N (unit-prefix encoding) ---
parse_age_years <- function(df) {
  df |>
    mutate(
      idade_unit  = as.integer(substr(NU_IDADE_N, 1, 1)),
      idade_value = as.integer(substr(NU_IDADE_N, 2, 4)),
      age_years   = case_when(
        idade_unit == 4 ~ as.numeric(idade_value),          # years
        idade_unit == 3 ~ as.numeric(idade_value) / 12,     # months -> years
        idade_unit == 2 ~ as.numeric(idade_value) / 365,    # days -> years
        TRUE            ~ NA_real_
      )
    )
}

df_lab      <- parse_age_years(df_lab)
df_clin     <- parse_age_years(df_clin)
df_combined <- parse_age_years(df_combined)

# --- Step 3: Assign age groups ---
age_breaks  <- c(0, 5, 10, 20, 30, 40, 50, 60, 70, Inf)
age_labels  <- c("<5", "5-9", "10-19", "20-29", "30-39",
                 "40-49", "50-59", "60-69", "70+")

assign_age_group <- function(df) {
  df |>
    mutate(
      age_group = cut(age_years,
                      breaks = age_breaks,
                      labels = age_labels,
                      right  = FALSE,
                      include.lowest = TRUE)
    )
}

df_lab      <- assign_age_group(df_lab)
df_clin     <- assign_age_group(df_clin)
df_combined <- assign_age_group(df_combined)

# --- Step 4: Summarise age distribution ---
summarise_age <- function(df, label) {
  df |>
    filter(!is.na(age_group)) |>
    count(age_group) |>
    mutate(
      pct       = n / sum(n) * 100,
      case_type = label
    )
}

age_summary <- bind_rows(
  summarise_age(df_lab,      "Laboratory confirmed"),
  summarise_age(df_clin,     "Clinical confirmed"),
  summarise_age(df_combined, "Combined")
)

# --- Step 5: Plot ---
ggplot(age_summary, aes(x = age_group, y = pct, fill = case_type)) +
  geom_col(position = "dodge") +
  facet_wrap(~case_type, ncol = 1) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Age distribution of chikungunya cases by confirmation type",
    x     = "Age group (years)",
    y     = "Proportion (%)",
    fill  = "Confirmation type"
  ) +
  theme_minimal() +
  theme(legend.position = "none")