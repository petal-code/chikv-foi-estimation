# ══════════════════════════════════════════════════════════════════
# Brazil CHIKV - Map Visualization
# ══════════════════════════════════════════════════════════════════

library(tidyverse)
library(sf)
library(geobr)
library(dplyr)

# ══════════════════════════════════════════════════════════════════
# 1. Municipality geometry
# ══════════════════════════════════════════════════════════════════

# Download once and cache locally (takes ~30s first time)
# muni_geo <- read_municipality(year = 2020, simplified = TRUE)
# muni_geo <- read_municipality(year = 2020, simplified = TRUE)
# saveRDS(muni_geo, "data/geo/muni_geo_2020.rds")
muni_geo <- readRDS("data/geo/muni_geo_2020.rds")

# geobr uses 7-digit code_muni; trim to 6 digits to match co_municipio
muni_geo <- muni_geo %>%
  mutate(co_municipio = as.integer(code_muni) %/% 10L)

all_muni_geo <- muni_geo %>% st_drop_geometry() %>% pull(co_municipio)

# ══════════════════════════════════════════════════════════════════
# 2. Sero study municipalities (7 state capitals)
# ══════════════════════════════════════════════════════════════════

sero_info <- tibble(
  co_municipio = c(130260, 230440, 261160, 310620, 330455, 355030, 410690),
  muni_name    = c("Manaus", "Fortaleza", "Recife",
                   "Belo Horizonte", "Rio de Janeiro", "São Paulo", "Curitiba")
)

sero_geo <- muni_geo %>%
  filter(co_municipio %in% sero_info$co_municipio)

sero_centroids <- sero_geo %>%
  st_centroid() %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  left_join(sero_info, by = "co_municipio") %>%
  mutate(
    label_lon = case_when(
      muni_name == "Manaus"          ~ lon + 4,
      muni_name == "Fortaleza"       ~ lon + 5,
      muni_name == "Recife"          ~ lon + 5,
      muni_name == "Belo Horizonte"  ~ lon + 9,
      muni_name == "Rio de Janeiro"  ~ lon + 5,
      muni_name == "São Paulo"       ~ lon + 5,
      muni_name == "Curitiba"        ~ lon + 5
    ),
    label_lat = case_when(
      muni_name == "Manaus"          ~ lat + 1,
      muni_name == "Fortaleza"       ~ lat + 3,
      muni_name == "Recife"          ~ lat + 0,
      muni_name == "Belo Horizonte"  ~ lat + 2,
      muni_name == "Rio de Janeiro"  ~ lat - 2,
      muni_name == "São Paulo"       ~ lat - 4,
      muni_name == "Curitiba"        ~ lat - 3
    )
  )

# ══════════════════════════════════════════════════════════════════
# 3. Cumulative incidence per municipality (all years)
# ══════════════════════════════════════════════════════════════════

cum_inc_muni <- cases_all %>%
  group_by(co_municipio, uf) %>%
  summarise(cum_cases = n(), .groups = "drop") %>%
  left_join(
    pop_all_full %>%
      group_by(co_municipio) %>%
      summarise(pop = sum(pop[year == 2019], na.rm = TRUE), .groups = "drop"),
    by = "co_municipio"
  ) %>%
  mutate(cum_incidence = cum_cases / pop * 1e5)

map_cum <- muni_geo %>%
  left_join(cum_inc_muni, by = "co_municipio") %>%
  mutate(cum_incidence = replace_na(cum_incidence, 0))

# ══════════════════════════════════════════════════════════════════
# 4. Period-specific incidence (2-year bands)
# ══════════════════════════════════════════════════════════════════

period_breaks <- list(
  "2014-2015" = 2014:2015,
  "2016-2017" = 2016:2017,
  "2018-2019" = 2018:2019,
  "2020-2021" = 2020:2021,
  "2022-2023" = 2022:2023,
  "2024-2025" = 2024:2025
)

inc_period <- map_dfr(names(period_breaks), function(pname) {
  yrs <- period_breaks[[pname]]
  
  cases_p <- cases_all %>%
    filter(year %in% yrs) %>%
    group_by(co_municipio) %>%
    summarise(cases = n(), .groups = "drop")
  
  pop_p <- pop_all_full %>%
    filter(year %in% yrs) %>%
    group_by(co_municipio) %>%
    summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")
  
  tibble(co_municipio = all_muni_geo) %>%
    left_join(cases_p, by = "co_municipio") %>%
    left_join(pop_p,   by = "co_municipio") %>%
    mutate(
      cases     = replace_na(cases, 0L),
      pop       = replace_na(pop,   0L),
      incidence = if_else(pop > 0, cases / pop * 1e5, 0),
      period    = pname
    )
}) %>%
  mutate(period = factor(period, levels = names(period_breaks)))

map_period <- muni_geo %>%
  left_join(inc_period, by = "co_municipio")

inc_max <- quantile(inc_period$incidence, 0.995, na.rm = TRUE)

# ══════════════════════════════════════════════════════════════════
# 5. Shared plot elements
# ══════════════════════════════════════════════════════════════════

theme_map <- theme_void(base_size = 11) +
  theme(
    legend.position   = "right",
    legend.key.width  = unit(0.4, "cm"),
    legend.key.height = unit(1.5, "cm"),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(size = 9, color = "grey40"),
    strip.text        = element_text(face = "bold", size = 10)
  )

sero_color <- "#FF3333"

# Reusable sero site overlay layers
sero_layers <- list(
  # Hollow circle at centroid
  geom_point(
    data   = sero_centroids,
    aes(x = lon, y = lat),
    shape  = 21,
    size   = 3.0,
    color  = sero_color,
    fill   = NA,
    stroke = 1.1
  ),
  # Line segment from label to centroid
  geom_segment(
    data = sero_centroids,
    aes(x    = label_lon - 0.8, y    = label_lat,
        xend = lon,             yend = lat),
    color     = sero_color,
    linewidth = 0.4
  ),
  # White label box
  geom_label(
    data       = sero_centroids,
    aes(x = label_lon, y = label_lat, label = muni_name),
    fill       = "white",
    color      = "black",
    size       = 2.2,
    label.size = 0.15,
    label.r    = unit(0.08, "cm"),
    fontface   = "bold"
  )
)

# ══════════════════════════════════════════════════════════════════
# 6. Map 1: Cumulative incidence — plain
# ══════════════════════════════════════════════════════════════════

p_map1_plain <- ggplot(map_cum) +
  geom_sf(aes(fill = cum_incidence), color = NA) +
  scale_fill_viridis_c(
    option   = "plasma",
    name     = "Cumulative\nincidence\n(per 100k)",
    na.value = viridis::viridis(1, option = "plasma", begin = 0),
    trans    = "log1p",
    breaks   = c(0, 10, 50, 100, 250, 500, 1000, 3000, 10000, 20000),
    labels   = scales::comma_format(accuracy = 1)
  ) +
  theme_map

print(p_map1_plain)

# ══════════════════════════════════════════════════════════════════
# 7. Map 1: Cumulative incidence — with sero sites
# ══════════════════════════════════════════════════════════════════

p_map1_sero <- ggplot(map_cum) +
  geom_sf(aes(fill = cum_incidence), color = NA) +
  sero_layers +
  scale_fill_viridis_c(
    option   = "plasma",
    name     = "Cumulative\nincidence\n(per 100k)",
    na.value = viridis::viridis(1, option = "plasma", begin = 0),
    trans    = "log1p",
    breaks   = c(0, 10, 50, 100, 250, 500, 1000, 3000, 10000, 30000),
    labels   = scales::comma_format(accuracy = 1)
  ) +
  theme_map

print(p_map1_sero)

# ══════════════════════════════════════════════════════════════════
# 8. Map 2: Period incidence — plain
# ══════════════════════════════════════════════════════════════════

p_map2_plain <- ggplot(map_period) +
  geom_sf(aes(fill = incidence), color = NA) +
  scale_fill_viridis_c(
    option   = "plasma",
    name     = "Incidence\n(per 100k)",
    na.value = viridis::viridis(1, option = "plasma", begin = 0),
    trans    = "log1p",
    limits   = c(0, inc_max),
    oob      = scales::squish,
    breaks   = c(0, 10, 50, 100, 250, 500, 1000),
    labels   = scales::comma_format(accuracy = 1)
  ) +
  facet_wrap(~ period, ncol = 3) +
  theme_map +
  theme(panel.spacing = unit(0.3, "cm"))

print(p_map2_plain)

# ══════════════════════════════════════════════════════════════════
# 9. Map 2: Period incidence — with sero sites
# ══════════════════════════════════════════════════════════════════

p_map2_sero <- ggplot(map_period) +
  geom_sf(aes(fill = incidence), color = NA) +
  sero_layers +
  scale_fill_viridis_c(
    option   = "plasma",
    name     = "Incidence\n(per 100k)",
    na.value = viridis::viridis(1, option = "plasma", begin = 0),
    trans    = "log1p",
    limits   = c(0, inc_max),
    oob      = scales::squish,
    breaks   = c(0, 10, 50, 100, 250, 500, 1000),
    labels   = scales::comma_format(accuracy = 1)
  ) +
  facet_wrap(~ period, ncol = 3) +
  theme_map +
  theme(panel.spacing = unit(0.3, "cm"))

print(p_map2_sero)