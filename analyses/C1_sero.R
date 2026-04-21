# ══════════════════════════════════════════════════════════════════
# Brazil CHIKV - Seroprevalence EDA: 7 regions x 3 survey years x sex
# Color = sex (M/F), facet = survey_year (rows) × region (cols)
# ══════════════════════════════════════════════════════════════════

library(tidyverse)
library(broom)


# ══════════════════════════════════════════════════════════════════
# 1. Region lookup
# ══════════════════════════════════════════════════════════════════

region_info <- tibble(
  sero_region  = c("Manaus", "Fortaleza", "Recife",
                   "BeloHorizonte", "RiodeJaneiro", "SãoPaulo", "Curitiba"),
  region_label = c("Manaus (AM)", "Fortaleza (CE)", "Recife (PE)",
                   "Belo Horizonte (MG)", "Rio de Janeiro (RJ)",
                   "São Paulo (SP)", "Curitiba (PR)")
)

# ══════════════════════════════════════════════════════════════════
# 2. Load & parse sero data
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
      sex      = SEX,
      sero_pos = as.integer(CHKG1_INTERPRETATION == "Positivo")
    ) %>%
    filter(!is.na(age_int), !is.na(sero_pos), !is.na(sex), age_int >= 0)
}

sero_all <- bind_rows(
  lapply(region_info$sero_region, function(reg) {
    bind_rows(
      parse_sero_df(get(paste0("df_", reg, "_Nov23")), reg, 2023),
      parse_sero_df(get(paste0("df_", reg, "_Jun24")), reg, 2024),
      parse_sero_df(get(paste0("df_", reg, "_Jun25")), reg, 2025)
    )
  })
) %>%
  mutate(sex = case_when(
    toupper(sex) %in% c("M", "MALE",    "MASCULINO") ~ "M",
    toupper(sex) %in% c("F", "FEMALE",  "FEMININO")  ~ "F",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(sex))

# ══════════════════════════════════════════════════════════════════
# 3. Aggregate for plot
# ══════════════════════════════════════════════════════════════════

year_levels  <- c(2023, 2024, 2025)
year_labels  <- c("Nov 2023", "Jun 2024", "Jun 2025")

age_breaks_10 <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, Inf)

sero_by_age <- sero_all %>%
  group_by(region, survey_year, sex, age_int) %>%
  summarise(prev = mean(sero_pos), .groups = "drop") %>%
  left_join(region_info, by = c("region" = "sero_region")) %>%
  mutate(
    survey_year  = factor(survey_year, levels = year_levels, labels = year_labels),
    region_label = factor(region_label, levels = region_info$region_label)
  )

sero_bin <- sero_all %>%
  mutate(
    bin_idx  = findInterval(age_int, age_breaks_10, rightmost.closed = TRUE),
    bin_low  = age_breaks_10[bin_idx],
    bin_high = pmin(bin_low + 10, 85)
  ) %>%
  group_by(region, survey_year, sex, bin_low, bin_high) %>%
  summarise(bin_prev = mean(sero_pos), .groups = "drop") %>%
  left_join(region_info, by = c("region" = "sero_region")) %>%
  mutate(
    survey_year  = factor(survey_year, levels = year_levels, labels = year_labels),
    region_label = factor(region_label, levels = region_info$region_label)
  )

# ══════════════════════════════════════════════════════════════════
# 4. Plot: rows = survey_year, cols = region
# ══════════════════════════════════════════════════════════════════

sex_colors <- c("F" = "#D6604D", "M" = "#2166AC")

p_sero_sex <- ggplot() +
  geom_point(
    data  = sero_by_age,
    aes(x = age_int, y = prev, color = sex),
    size  = 0.8, alpha = 0.45
  ) +
  geom_segment(
    data = sero_bin,
    aes(x = bin_low, xend = bin_high,
        y = bin_prev, yend = bin_prev,
        color = sex),
    linewidth = 1.2, alpha = 0.9
  ) +
  scale_color_manual(
    values = sex_colors,
    name   = "Sex",
    labels = c("F" = "Female", "M" = "Male")
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_x_continuous(
    breaks = seq(0, 80, by = 20),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_grid(survey_year ~ region_label) +
  labs(x = "Age (years)", y = "Seroprevalence") +
  theme_bw(base_size = 20) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(2, "cm"),
    panel.grid.minor = element_blank()
  )

p_sero_sex

# ══════════════════════════════════════════════════════════════════
# 5. Sex heterogeneity: chi-square, region-stratified
#    2×2 table (sex × sero_pos) per region × survey_year
#    BH correction across 21 cells
# ══════════════════════════════════════════════════════════════════

chisq_results <- sero_all %>%
  group_by(region, survey_year) %>%
  summarise(
    n_F_pos = sum(sero_pos == 1 & sex == "F"),
    n_F_neg = sum(sero_pos == 0 & sex == "F"),
    n_M_pos = sum(sero_pos == 1 & sex == "M"),
    n_M_neg = sum(sero_pos == 0 & sex == "M"),
    prev_F  = round(n_F_pos / (n_F_pos + n_F_neg), 3),
    prev_M  = round(n_M_pos / (n_M_pos + n_M_neg), 3),
    p_raw   = tryCatch(
      chisq.test(matrix(c(n_F_pos, n_F_neg,
                          n_M_pos, n_M_neg), nrow = 2))$p.value,
      warning = function(w) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj   = p.adjust(p_raw, method = "BH"),
    sig_adj = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    survey_year = factor(survey_year, levels = year_levels, labels = year_labels)
  ) %>%
  left_join(region_info, by = c("region" = "sero_region")) %>%
  select(region_label, survey_year, prev_F, prev_M, p_raw, p_adj, sig_adj)

cat("\n── Sex heterogeneity: chi-square per region × year (BH-corrected) ──\n")
print(chisq_results, n = Inf)

# ══════════════════════════════════════════════════════════════════
# 6. Age heterogeneity: logistic regression
#    6a. Pooled (all regions)
#    6b. Region-stratified (separate model per region)
#    Reference: age group 25-29
# ══════════════════════════════════════════════════════════════════

age_breaks_5 <- c(0, 5, 10, 15, 20, 25, 30, 35, 40,
                  45, 50, 55, 60, 65, 70, 75, 80, Inf)
age_labels_5 <- c("0-4","5-9","10-14","15-19","20-24",
                  "25-29","30-34","35-39","40-44",
                  "45-49","50-54","55-59","60-64",
                  "65-69","70-74","75-79","80+")

sero_age <- sero_all %>%
  mutate(
    age_group   = cut(age_int, breaks = age_breaks_5,
                      labels = age_labels_5, right = FALSE),
    age_group   = relevel(factor(age_group), ref = "25-29"),
    survey_year = factor(survey_year)
  ) %>%
  filter(!is.na(age_group))

format_age_results <- function(fit) {
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("age_group", term)) %>%
    mutate(
      term    = gsub("age_group", "", term),
      across(c(estimate, conf.low, conf.high), ~ round(.x, 3)),
      p_value = signif(p.value, 3),
      sig     = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    ) %>%
    select(age_group = term, OR = estimate,
           CI_low = conf.low, CI_high = conf.high, p_value, sig)
}

# ── 6a. Pooled ────────────────────────────────────────────────────
fit_pooled <- glm(
  sero_pos ~ age_group + region + survey_year,
  data   = sero_age,
  family = binomial(link = "logit")
)

cat("\n── Age heterogeneity: pooled (ref = 25-29) ──\n")
print(format_age_results(fit_pooled), n = Inf)

# ── 6b. Region-stratified ─────────────────────────────────────────
cat("\n── Age heterogeneity: region-stratified (ref = 25-29) ──\n")
region_info$sero_region %>%
  set_names() %>%
  lapply(function(reg) {
    fit <- glm(
      sero_pos ~ age_group + survey_year,
      data   = sero_age %>% filter(region == reg),
      family = binomial(link = "logit")
    )
    format_age_results(fit) %>% mutate(region = reg, .before = 1)
  }) %>%
  bind_rows() %>%
  left_join(region_info, by = c("region" = "sero_region")) %>%
  select(region_label, age_group, OR, CI_low, CI_high, p_value, sig) %>%
  print(n = Inf)