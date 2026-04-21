# ══════════════════════════════════════════════════════════════════
# Brazil CHIKV - Exploratory Analysis: All Estados
# ══════════════════════════════════════════════════════════════════

library(tidyverse)
library(lubridate)
library(brpop)

# ══════════════════════════════════════════════════════════════════
# 0. Prerequisites
# ══════════════════════════════════════════════════════════════════

df_raw_slim <- readRDS("data/case/df_raw_SINAN_CHIKV_2013_2025_slim.rds")

# Population data (save after first run)
# pop_all_raw <- mun_pop_age(source = "datasus", sex = "all")
# saveRDS(pop_all_raw, "data/pop/pop_all_raw_datasus.rds")
pop_all_raw <- readRDS("data/pop/pop_all_raw_datasus.rds")


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

parse_age_group <- function(ag) {
  case_when(
    ag == "Total"                 ~ NA_integer_,
    ag == "From 80 years or more" ~ 80L,
    TRUE ~ as.integer(str_extract(ag, "^From (\\d+)", group = 1))
  )
}

# ══════════════════════════════════════════════════════════════════
# 1. Mappings & Shared Objects
# ══════════════════════════════════════════════════════════════════

estado_name_map <- tibble(
  co_estado = c("11","12","13","14","15","16","17",
                "21","22","23","24","25","26","27","28","29",
                "31","32","33","35",
                "41","42","43",
                "50","51","52","53"),
  uf        = c("RO","AC","AM","RR","PA","AP","TO",
                "MA","PI","CE","RN","PB","PE","AL","SE","BA",
                "MG","ES","RJ","SP",
                "PR","SC","RS",
                "MS","MT","GO","DF")
)

regiao_map <- tibble(
  uf     = c("RO","AC","AM","RR","PA","AP","TO",
             "MA","PI","CE","RN","PB","PE","AL","SE","BA",
             "MG","ES","RJ","SP",
             "PR","SC","RS",
             "MS","MT","GO","DF"),
  regiao = c(rep("North",        7),
             rep("Northeast",    9),
             rep("Southeast",    4),
             rep("South",        3),
             rep("Central-West", 4))
)

regiao_order  <- c("North", "Northeast", "Central-West", "Southeast", "South")
regiao_colors <- c(
  "North"        = "#1565C0",
  "Northeast"    = "#E65100",
  "Central-West" = "#6A1B9A",
  "Southeast"    = "#AD1457",
  "South"        = "#2E7D32"
)

uf_order <- regiao_map %>%
  mutate(regiao = factor(regiao, levels = regiao_order)) %>%
  arrange(regiao, uf) %>%
  pull(uf)

age_breaks_model <- c(0, 5, 10, 15, 20, 25, 30, 35, 40,
                      45, 50, 55, 60, 65, 70, 75, 80, Inf)
age_labels_model <- c(0, 5, 10, 15, 20, 25, 30, 35, 40,
                      45, 50, 55, 60, 65, 70, 75, 80)

age_axis <- scale_x_continuous(
  breaks = age_labels_model,
  labels = if_else(age_labels_model == 80, "80+",
                   paste0(age_labels_model, "-", age_labels_model + 4)),
  expand = expansion(mult = c(0.02, 0.05))
)

# ══════════════════════════════════════════════════════════════════
# 2. Data Preparation
# ══════════════════════════════════════════════════════════════════

# ── 2a. Case data ─────────────────────────────────────────────────
cases_all <- df_raw_slim %>%
  filter(
    CLASSI_FIN == "13",
    CRITERIO   %in% c("1", "2")
  ) %>%
  mutate(
    co_municipio = as.integer(ID_MUNICIP),
    co_estado    = SG_UF_NOT,
    year         = year(DT_NOTIFIC),
    age_years    = parse_age_years(NU_IDADE_N),
    age_lower    = age_labels_model[
      findInterval(floor(age_years), age_breaks_model,
                   rightmost.closed = TRUE)
    ]
  ) %>%
  filter(!is.na(age_years), !is.na(year), year >= 2014, year <= 2025) %>%
  left_join(estado_name_map, by = "co_estado")

# ── 2b. Population data ───────────────────────────────────────────
all_munis <- cases_all %>% distinct(co_municipio) %>% pull()

pop_all <- pop_all_raw %>%
  filter(code_muni %in% all_munis) %>%
  mutate(
    co_municipio = code_muni,
    age_lower    = parse_age_group(age_group)
  ) %>%
  filter(!is.na(age_lower)) %>%
  left_join(cases_all %>% distinct(co_municipio, co_estado, uf),
            by = "co_municipio") %>%
  select(co_municipio, co_estado, uf, year, age_lower, pop)

pop_all_extrap <- pop_all %>%
  filter(year %in% 2019:2021) %>%
  group_by(co_municipio, co_estado, uf, age_lower) %>%
  arrange(year) %>%
  summarise(
    slope    = (pop[year == 2021] - pop[year == 2019]) / 2,
    pop_2021 = pop[year == 2021],
    .groups  = "drop"
  ) %>%
  crossing(year = 2022:2025) %>%
  mutate(pop = pmax(0L, as.integer(pop_2021 + slope * (year - 2021)))) %>%
  select(co_municipio, co_estado, uf, year, age_lower, pop)

pop_all_full <- bind_rows(
  pop_all %>% filter(year >= 2014),
  pop_all_extrap
) %>%
  arrange(uf, co_municipio, year, age_lower)

pop_2019_muni <- pop_all_full %>%
  filter(year == 2019) %>%
  select(co_municipio, uf, age_lower, pop)

# ── 2c. Aggregations ──────────────────────────────────────────────

# Municipality x year: incidence (plot 1)
cases_muni_year <- cases_all %>%
  group_by(co_municipio, uf, year) %>%
  summarise(cases = n(), .groups = "drop") %>%
  left_join(
    pop_all_full %>%
      group_by(co_municipio, uf, year) %>%
      summarise(pop = sum(pop), .groups = "drop"),
    by = c("co_municipio", "uf", "year")
  ) %>%
  mutate(incidence = cases / pop * 1e5)

# Municipality x age: cumulative cases and incidence
cases_muni_age <- cases_all %>%
  filter(!is.na(age_lower)) %>%
  group_by(co_municipio, uf, age_lower) %>%
  summarise(cum_cases = n(), .groups = "drop") %>%
  left_join(pop_2019_muni, by = c("co_municipio", "uf", "age_lower")) %>%
  mutate(cum_incidence = cum_cases / pop * 1e5)

# State x age: sum across municipalities, then normalize within state
cases_estado_norm <- cases_muni_age %>%
  group_by(uf, age_lower) %>%
  summarise(
    cum_cases = sum(cum_cases, na.rm = TRUE),
    pop       = sum(pop,       na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(cum_incidence = cum_cases / pop * 1e5) %>%
  group_by(uf) %>%
  mutate(cum_inc_norm = cum_incidence / sum(cum_incidence, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(regiao_map, by = "uf") %>%
  mutate(
    regiao = factor(regiao, levels = regiao_order),
    uf     = factor(uf,     levels = uf_order)
  )

# Region x age: mean of state-level values, then re-normalize within region
cases_regiao_norm <- cases_estado_norm %>%
  group_by(regiao, age_lower) %>%
  summarise(
    mean_cum_cases     = mean(cum_cases,      na.rm = TRUE),
    mean_cum_incidence = mean(cum_incidence,  na.rm = TRUE),
    mean_norm          = mean(cum_inc_norm,   na.rm = TRUE),
    .groups            = "drop"
  ) %>%
  group_by(regiao) %>%
  mutate(mean_norm = mean_norm / sum(mean_norm, na.rm = TRUE)) %>%
  ungroup()

# State total incidence (all ages combined, for bar plot)
cases_estado_total <- cases_estado_norm %>%
  group_by(uf, regiao) %>%
  summarise(
    cum_cases = sum(cum_cases, na.rm = TRUE),
    pop       = sum(pop,       na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(cum_incidence = cum_cases / pop * 1e5)

# ══════════════════════════════════════════════════════════════════
# 3. Plots
# ══════════════════════════════════════════════════════════════════

# ── Plot 1: Incidence time series by state (municipality-level lines) ──
p1 <- cases_muni_year %>%
  ggplot(aes(x = year, y = incidence, group = co_municipio)) +
  geom_line(alpha = 0.2, linewidth = 0.4, color = "steelblue") +
  facet_wrap(~ uf, scales = "free_y", ncol = 5) +
  scale_x_continuous(
    breaks = seq(2015, 2025, by = 5),
    limits = c(2015, 2025)
  ) +
  labs(y = "Incidence (per 100,000)", x = NULL) +
  theme_bw()

print(p1)

# ── Plot 2: Cumulative cases by age group ─────────────────────────
p2 <- ggplot() +
  geom_line(data = cases_estado_norm,
            aes(x = age_lower, y = cum_cases,
                group = uf, color = regiao),
            alpha = 0.4, linewidth = 0.6) +
  geom_line(data = cases_regiao_norm,
            aes(x = age_lower, y = mean_cum_cases,
                group = regiao, color = regiao),
            linewidth = 1.5) +
  geom_point(data = cases_regiao_norm,
             aes(x = age_lower, y = mean_cum_cases, color = regiao),
             size = 2.0) +
  scale_color_manual(values = regiao_colors) +
  facet_wrap(~ regiao, scales = "free_y", ncol = 5) +
  age_axis +
  labs(y = "Cumulative cases", x = "Age group") +
  theme_bw(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p2)

# ── Plot 3: Cumulative incidence by age group ─────────────────────
p3 <- ggplot() +
  geom_line(data = cases_estado_norm,
            aes(x = age_lower, y = cum_incidence,
                group = uf, color = regiao),
            alpha = 0.4, linewidth = 0.6) +
  geom_line(data = cases_regiao_norm,
            aes(x = age_lower, y = mean_cum_incidence,
                group = regiao, color = regiao),
            linewidth = 1.5) +
  geom_point(data = cases_regiao_norm,
             aes(x = age_lower, y = mean_cum_incidence, color = regiao),
             size = 2.0) +
  scale_color_manual(values = regiao_colors) +
  facet_wrap(~ regiao, scales = "free_y", ncol = 5) +
  age_axis +
  labs(y = "Cumulative incidence (per 100,000)", x = "Age group") +
  theme_bw(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p3)

# ── Plot 4: Normalized incidence by age group, faceted by region ──
p4 <- ggplot() +
  geom_line(data = cases_estado_norm,
            aes(x = age_lower, y = cum_inc_norm,
                group = uf, color = regiao),
            alpha = 0.4, linewidth = 0.6) +
  geom_line(data = cases_regiao_norm,
            aes(x = age_lower, y = mean_norm,
                group = regiao, color = regiao),
            linewidth = 1.5) +
  geom_point(data = cases_regiao_norm,
             aes(x = age_lower, y = mean_norm, color = regiao),
             size = 2.0) +
  scale_color_manual(values = regiao_colors) +
  facet_wrap(~ regiao, scales = "free_y", ncol = 5) +
  age_axis +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(y = "Proportion of total incidence", x = "Age group") +
  theme_bw(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p4)

# ── Plot 5: Normalized incidence by age group, single panel ───────
p5 <- ggplot() +
  geom_line(data = cases_estado_norm,
            aes(x = age_lower, y = cum_inc_norm,
                color = regiao, group = uf),
            alpha = 0.3, linewidth = 0.6) +
  geom_line(data = cases_regiao_norm,
            aes(x = age_lower, y = mean_norm,
                color = regiao, group = regiao),
            linewidth = 1.8) +
  geom_point(data = cases_regiao_norm,
             aes(x = age_lower, y = mean_norm, color = regiao),
             size = 2.5) +
  scale_color_manual(values = regiao_colors) +
  age_axis +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    limits = c(0, 0.25)
  ) +
  labs(x = "Age group", y = "Proportion of total incidence",
       color = "Region") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p5)

# ── Plot 5 (no South): Normalized incidence, South region excluded ─
p5_nosouth <- ggplot() +
  geom_line(data = cases_estado_norm %>% filter(regiao != "South"),
            aes(x = age_lower, y = cum_inc_norm,
                color = regiao, group = uf),
            alpha = 0.3, linewidth = 0.6) +
  geom_line(data = cases_regiao_norm %>% filter(regiao != "South"),
            aes(x = age_lower, y = mean_norm,
                color = regiao, group = regiao),
            linewidth = 1.8) +
  geom_point(data = cases_regiao_norm %>% filter(regiao != "South"),
             aes(x = age_lower, y = mean_norm, color = regiao),
             size = 2.5) +
  scale_color_manual(values = regiao_colors) +
  age_axis +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    limits = c(0, 0.13)
  ) +
  labs(x = "Age group", y = "Proportion of total incidence",
       color = "Region") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p5_nosouth)

# ── Plot 6: Cumulative incidence bar chart by state ───────────────
p6 <- ggplot(cases_estado_total,
             aes(x = uf, y = cum_incidence, fill = regiao)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = regiao_colors) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "State", y = "Cumulative incidence (per 100,000)",
       fill = "Region") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p6)