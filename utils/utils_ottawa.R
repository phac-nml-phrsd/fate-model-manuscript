


# Reproducible helpers for the Ottawa neighbourhood-WWTP benchmark ------------

library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(janitor)
library(stringr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(gridExtra)


#--- helper function
ML_PER_L <- 1000
wwtp_pop <- 910000#29580#910000
site_pop <- 8000#1020#10000

get_site_population <- function(site){
  dplyr::case_when(
    site == "WWTP" ~ wwtp_pop,
    site == "D" ~ site_pop,
    TRUE ~ NA_real_
  )
}

# =============================================================================
# 1) Import wastewater data and format + add derived targets
# =============================================================================
load_data <- function(file_path){
  
  sheets <- excel_sheets(file_path)
  sheets <- sheets[!str_detect(tolower(sheets), "residence|hrt|travel|retention")]
  
  # --- robust date fixer (Excel serial + R numeric + many character formats) ---
  fix_mixed_date <- function(x) {
    if (inherits(x, "Date")) return(x)
    
    xn <- suppressWarnings(as.numeric(x))
    
    out <- dplyr::case_when(
      !is.na(xn) & xn > 30000 ~ as.Date(xn, origin = "1899-12-30"),  # Excel serial
      !is.na(xn) & xn > 10000 ~ as.Date(xn, origin = "1970-01-01"),  # R numeric
      TRUE ~ as.Date(NA)
    )
    
    needs <- is.na(out) & !is.na(x) & as.character(x) != ""
    if (any(needs)) {
      out[needs] <- as.Date(
        parse_date_time(
          as.character(x)[needs],
          orders = c("Ymd", "Y-m-d", "mdY", "m/d/Y", "d/m/Y", "dmY", "d-b-Y", "d-b-y")
        )
      )
    }
    
    out
  }
  
  # ---- read all sheets ----
  df_all <- purrr::map_dfr(
    sheets,
    function(sh) {
      readxl::read_excel(file_path, sheet = sh, col_types = "text") %>%
        janitor::clean_names() %>%
        dplyr::select(dplyr::where(~ !all(is.na(.x)))) %>%
        dplyr::mutate(sheet = sh)
    }
  )
  
  # ---- tidy dataset + de-duplicate ----
  df_long <- df_all %>%
    dplyr::mutate(
      date = fix_mixed_date(date),
      note = if ("note" %in% names(.)) note else NA_character_
    ) %>%
    tidyr::pivot_longer(

      cols = dplyr::any_of(c("wwtp", "d")),
      names_to = "site",
      values_to = "value",
      values_transform = list(value = ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    dplyr::mutate(
      site = dplyr::recode(site, wwtp = "WWTP", d = "D"),
      target = as.character(target)
    ) %>%
    tidyr::separate(
      col = target,
      into = c("target", "type"),
      sep = "_",
      extra = "merge",
      fill = "right"
    ) %>%
    dplyr::mutate(
      target = stringr::str_to_upper(target),
      type   = stringr::str_to_upper(type)
    ) %>%
    dplyr::select(date, site, value, target, type, note) %>%
    dplyr::filter(!is.na(date), !is.na(site), !is.na(target)) %>%
    dplyr::group_by(date, site, target, type) %>%
    dplyr::summarise(
      value = {
        v <- value[!is.na(value)]
        if (length(v) == 0) NA_real_
        else if (length(v) == 1) v
        else stats::median(v)
      },
      note = {
        n <- note[!is.na(note)]
        n <- stringr::str_trim(n)
        n <- n[n != ""]
        if (length(n) == 0) NA_character_ else paste(unique(n), collapse = "; ")
      },
      .groups = "drop"
    ) %>%
    dplyr::arrange(target, type, site, date)
  
  # ---- add derived targets: N1/N2 over TS, VS, FLOW, PMMOV (within each date/site/type) ----
  safe_ratio <- function(num, den) {
    dplyr::if_else(!is.na(num) & !is.na(den) & num > 0 & den > 0, num / den, NA_real_)
  }
  
  # ---- add derived targets: per (date, site, type) ----
  #   - TS/VS normalization: divide concentration by TS/VS
  #   - PMMOV normalization: divide by PMMOV (same date/site/type)
  #   - FLOW: flow-weighted load = conc * flow * 1000  (gc/s)
  #   - FLOW per-capita: (gc/s) / population           (gc/person/s)
  
  df_long_aug <- df_long %>%
    dplyr::mutate(
      date   = as.Date(date),
      site   = trimws(as.character(site)),
      type   = toupper(trimws(as.character(type))),
      target = toupper(trimws(as.character(target))),
      value  = as.numeric(value)
    ) %>%
    {
      df_orig <- .
      
      
      
      safe_ratio <- function(num, den) {
        dplyr::if_else(!is.na(num) & !is.na(den) & num > 0 & den > 0, num / den, NA_real_)
      }
      
      # Denominators by (date, site) only (TS/VS/FLOW usually have unit-type in "type")
      df_denoms <- df_orig %>%
        dplyr::filter(target %in% c("TS","VS","FLOW"), !is.na(value)) %>%
        dplyr::group_by(date, site, target) %>%
        dplyr::summarise(value = stats::median(value, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = target, values_from = value) %>%
        dplyr::rename(TS_d = TS, VS_d = VS, FLOW_d = FLOW)
      
      # Biomarkers by (date, site, type)  âœ… keeps CF and UF separate
      df_bio <- df_orig %>%
        dplyr::filter(target %in% c("N1","N2","PMMOV"), !is.na(value)) %>%
        dplyr::group_by(date, site, type, target) %>%

        dplyr::summarise(
          value = stats::median(value, na.rm = TRUE),
          note  = dplyr::first(note),
          .groups = "drop"
        ) %>%
        tidyr::pivot_wider(names_from = target, values_from = value) %>%
        dplyr::left_join(df_denoms, by = c("date","site")) %>%
        dplyr::mutate(
          population = get_site_population(site),
          
          # ---- TS/VS normalizations (gc/mL per (mg/L) etc) ----
          N1_TS    = safe_ratio(N1, TS_d),
          N2_TS    = safe_ratio(N2, TS_d),
          PMMOV_TS = safe_ratio(PMMOV, TS_d),
          
          N1_VS    = safe_ratio(N1, VS_d),
          N2_VS    = safe_ratio(N2, VS_d),
          PMMOV_VS = safe_ratio(PMMOV, VS_d),
          
          # ---- PMMOV normalization (within CF/UF type) ----
          N1_PMMOV = safe_ratio(N1, PMMOV),
          N2_PMMOV = safe_ratio(N2, PMMOV),
          
          # ---- FLOW-weighted load (gc/s) ----
          N1_FLOW    = dplyr::if_else(!is.na(N1)    & !is.na(FLOW_d) & N1    > 0 & FLOW_d > 0,
                                      N1    * FLOW_d * ML_PER_L, NA_real_),
          N2_FLOW    = dplyr::if_else(!is.na(N2)    & !is.na(FLOW_d) & N2    > 0 & FLOW_d > 0,
                                      N2    * FLOW_d * ML_PER_L, NA_real_),
          PMMOV_FLOW = dplyr::if_else(!is.na(PMMOV) & !is.na(FLOW_d) & PMMOV > 0 & FLOW_d > 0,
                                      PMMOV * FLOW_d * ML_PER_L, NA_real_),
          
          # ---- FLOW per-capita ONLY (gc/person/s) ----
          N1_FLOW_PC    = dplyr::if_else(!is.na(N1_FLOW)    & !is.na(population) & population > 0,
                                         N1_FLOW / population, NA_real_),
          N2_FLOW_PC    = dplyr::if_else(!is.na(N2_FLOW)    & !is.na(population) & population > 0,
                                         N2_FLOW / population, NA_real_),
          PMMOV_FLOW_PC = dplyr::if_else(!is.na(PMMOV_FLOW) & !is.na(population) & population > 0,
                                         PMMOV_FLOW / population, NA_real_)
        ) %>%
        tidyr::pivot_longer(
          cols = dplyr::any_of(c(
            "N1_TS","N2_TS","PMMOV_TS",
            "N1_VS","N2_VS","PMMOV_VS",
            "N1_PMMOV","N2_PMMOV",
            "N1_FLOW","N2_FLOW","PMMOV_FLOW",
            "N1_FLOW_PC","N2_FLOW_PC","PMMOV_FLOW_PC"
          )),
          names_to = "target",
          values_to = "value"
        ) %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::transmute(
          date, site, type,
          target,
          value,
          note = dplyr::coalesce(note, NA_character_)
        )
      
      dplyr::bind_rows(df_orig, df_bio) %>%
        dplyr::arrange(date, site, type, target)
    }
  
  
  return(df_long_aug)
}


#----- ADD TOTAL (liquid + Solid) to the data frame
# =============================================================================
# Add TOTAL (CF+UF) for biomarkers, then add TOTAL-normalized ratios
#   - TOTAL biomarkers: N1_TOTAL = N1_UF + N1_CF (same date/site)
#                      N2_TOTAL = N2_UF + N2_CF
#                      PMMOV_TOTAL = PMMOV_UF + PMMOV_CF
#   - Then compute (within TOTAL):
#       N1_TS, N2_TS, N1_VS, N2_VS, N1_FLOW, N2_FLOW, N1_PMMOV, N2_PMMOV
# =============================================================================
add_total_and_total_normalized <- function(df_long,
                                           biomarkers = c("N1","N2","PMMOV"),
                                           type_liquid = "UF",
                                           type_solid  = "CF",

                                           type_total  = "TOTAL",
                                           denom_targets = c("TS","VS","FLOW")) {
  
  stopifnot(all(c("date","site","target","type","value") %in% names(df_long)))
  
  safe_ratio <- function(num, den) {
    dplyr::if_else(!is.na(num) & !is.na(den) & num > 0 & den > 0, num / den, NA_real_)
  }
  
  df0 <- df_long %>%
    dplyr::mutate(
      date   = as.Date(date),
      site   = trimws(as.character(site)),
      target = toupper(trimws(as.character(target))),
      type   = toupper(trimws(as.character(type))),
      value  = as.numeric(value)
    )
  
  # ------------------------------------------------------------
  # 1) Build TOTAL biomarker rows by summing CF + UF (date/site)
  # ------------------------------------------------------------
  df_bio_wide <- df0 %>%
    dplyr::filter(target %in% biomarkers, type %in% c(type_liquid, type_solid)) %>%
    dplyr::group_by(date, site, target, type) %>%
    dplyr::summarise(
      value = stats::median(value, na.rm = TRUE),
      note  = dplyr::first(dplyr::coalesce(note, NA_character_)),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = type, values_from = value)
  
  # ensure both columns exist even if missing entirely
  if (!(type_liquid %in% names(df_bio_wide))) df_bio_wide[[type_liquid]] <- NA_real_
  if (!(type_solid  %in% names(df_bio_wide))) df_bio_wide[[type_solid ]] <- NA_real_
  
  df_total_bio <- df_bio_wide %>%
    dplyr::mutate(
      value = dplyr::coalesce(.data[[type_liquid]], 0) + dplyr::coalesce(.data[[type_solid]], 0),
      type  = type_total
    ) %>%
    dplyr::select(date, site, target, type, value) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::mutate(
      note = paste0("TOTAL = ", type_liquid, " + ", type_solid)
    )
  
  # ------------------------------------------------------------
  # 2) Attach denominators (TS/VS/FLOW) by (date, site) only
  #    because their `type` is units (MG_PER_L, RATE_L_PER_S)
  # ------------------------------------------------------------
  df_denoms <- df0 %>%
    dplyr::filter(target %in% denom_targets, !is.na(value)) %>%
    dplyr::group_by(date, site, target) %>%
    dplyr::summarise(value = stats::median(value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = target, values_from = value) %>%
    dplyr::rename_with(~ paste0(.x, "_d"), dplyr::any_of(denom_targets))
  
  # ------------------------------------------------------------
  # 3) Compute TOTAL-normalized ratios using TOTAL biomarker rows
  #    - N1/PMMOV and N2/PMMOV use TOTAL PMMOV as denominator
  # ------------------------------------------------------------
  df_total_wide <- df_total_bio %>%
    dplyr::select(date, site, target, value) %>%
    tidyr::pivot_wider(names_from = target, values_from = value) %>%
    dplyr::left_join(df_denoms, by = c("date","site"))
  
  # make sure biomarker columns exist
  for (b in biomarkers) if (!(b %in% names(df_total_wide))) df_total_wide[[b]] <- NA_real_
  
  # denom columns (may be missing if not in data)
  TS_d   <- if ("TS_d"   %in% names(df_total_wide)) df_total_wide$TS_d   else NA_real_
  VS_d   <- if ("VS_d"   %in% names(df_total_wide)) df_total_wide$VS_d   else NA_real_
  FLOW_d <- if ("FLOW_d" %in% names(df_total_wide)) df_total_wide$FLOW_d else NA_real_
  
  ML_PER_L <- 1000
  
  df_total_ratios <- df_total_wide %>%
    dplyr::mutate(
      population = get_site_population(site)
    ) %>%

    dplyr::transmute(
      date, site,
      type = type_total,
      
      N1_TS    = safe_ratio(.data$N1, TS_d),
      N2_TS    = safe_ratio(.data$N2, TS_d),
      N1_VS    = safe_ratio(.data$N1, VS_d),
      N2_VS    = safe_ratio(.data$N2, VS_d),
      
      # FLOW-weighted TOTAL loads (gc/s)
      N1_FLOW = dplyr::if_else(!is.na(.data$N1) & !is.na(FLOW_d) & .data$N1 > 0 & FLOW_d > 0,
                               .data$N1 * FLOW_d * ML_PER_L, NA_real_),
      N2_FLOW = dplyr::if_else(!is.na(.data$N2) & !is.na(FLOW_d) & .data$N2 > 0 & FLOW_d > 0,
                               .data$N2 * FLOW_d * ML_PER_L, NA_real_),
      PMMOV_FLOW = dplyr::if_else(!is.na(.data$PMMOV) & !is.na(FLOW_d) & .data$PMMOV > 0 & FLOW_d > 0,
                                  .data$PMMOV * FLOW_d * ML_PER_L, NA_real_),
      
      # per-capita ONLY for FLOW (gc/person/s)
      N1_FLOW_PC = dplyr::if_else(!is.na(.data$N1) & !is.na(FLOW_d) & !is.na(population) &
                                    .data$N1 > 0 & FLOW_d > 0 & population > 0,
                                  .data$N1 * FLOW_d * ML_PER_L / population, NA_real_),
      N2_FLOW_PC = dplyr::if_else(!is.na(.data$N2) & !is.na(FLOW_d) & !is.na(population) &
                                    .data$N2 > 0 & FLOW_d > 0 & population > 0,
                                  .data$N2 * FLOW_d * ML_PER_L / population, NA_real_),
      PMMOV_FLOW_PC = dplyr::if_else(!is.na(.data$PMMOV) & !is.na(FLOW_d) & !is.na(population) &
                                       .data$PMMOV > 0 & FLOW_d > 0 & population > 0,
                                     .data$PMMOV * FLOW_d * ML_PER_L / population, NA_real_),
      
      N1_PMMOV = safe_ratio(.data$N1, .data$PMMOV),
      N2_PMMOV = safe_ratio(.data$N2, .data$PMMOV),
      PMMOV_TS = safe_ratio(.data$PMMOV, TS_d),
      PMMOV_VS = safe_ratio(.data$PMMOV, VS_d)
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        N1_TS, N2_TS, N1_VS, N2_VS,
        N1_FLOW, N2_FLOW, PMMOV_FLOW,
        N1_FLOW_PC, N2_FLOW_PC, PMMOV_FLOW_PC,
        N1_PMMOV, N2_PMMOV,
        PMMOV_TS, PMMOV_VS
      ),
      names_to = "target",
      values_to = "value"
    ) %>%
    dplyr::filter(!is.na(value))%>%
    dplyr::mutate(note = "Computed from TOTAL biomarkers; TS/VS/FLOW joined by (date, site). FLOW-derived values are flow-weighted loads (gc/s).")
  
  # ------------------------------------------------------------
  # 4) Bind back to original df
  # ------------------------------------------------------------
  dplyr::bind_rows(
    df0,
    df_total_bio %>% dplyr::select(dplyr::any_of(names(df0))),  # keep same columns where possible
    df_total_ratios %>% dplyr::select(dplyr::any_of(names(df0)))
  ) %>%
    dplyr::arrange(date, site, type, target)
}




#========= Plotting biomarkers by methods az a box pot
# =============================================================================
# Publication-ready styling + updated plotting functions
# Applies:
# - Panel tags (A, B, â€¦) on the final combined figure
# - Consistent visual grammar (TS/VS/FLOW colors; CF/UF/TOTAL fill)
# - Reduced non-data ink (clean theme)
# - Dual-axis note "VS scaled to TS axis"
# - Harmonized x-limits across time-series panels
# - Cleaner boxplots (no outliers, stronger median, jittered points)
# - Robust â€œdrop last date per siteâ€
# =============================================================================

# ---- global helper theme (use everywhere) ----
theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),

      axis.line = ggplot2::element_line(color = "black", linewidth = 0.35),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank()
    )
}

# ---- consistent palettes ----
pal_hydro <- list(
  TS   = "#ff7f00",  
  VS   = "gray",  
  FLOW = "#1f78b4"   
)

pal_type_fill <- c(
  CF    = "#08519c",  # solid
  UF    = "#6baed6",  # liquid
  TOTAL = "#8c510a"   # total
)

# ---- display labels for sample fractions ----
fraction_labels <- c(
  CF    = "Solids",
  UF    = "Liquid",
  TOTAL = "Total"
)

# =============================================================================
# 1) Boxplots (with points) for N1/N2/PMMOV by type (CF/UF/TOTAL)
#    - returns two ggplots: p_wwtp and p_d
# =============================================================================
plot_biomarker_boxplots_WWTP_and_D_pub <- function(df_long,
                                                   biomarkers = c("N1","N2","PMMOV"),
                                                   sites = c("WWTP","D"),
                                                   types = c("CF","UF","TOTAL"),
                                                   log10_y = TRUE) {
  stopifnot(all(c("date","site","target","type","value") %in% names(df_long)))
  
  library(dplyr)
  library(ggplot2)
  
  dfp <- df_long %>%
    mutate(
      date   = as.Date(date),
      site   = toupper(trimws(as.character(site))),
      target = toupper(trimws(as.character(target))),
      type   = toupper(trimws(as.character(type))),
      value  = as.numeric(value)
    ) %>%
    filter(
      site %in% toupper(sites),
      target %in% toupper(biomarkers),
      type %in% toupper(types),
      !is.na(value)
    ) %>%
    mutate(
      site   = factor(site, levels = toupper(sites)),
      target = factor(target, levels = toupper(biomarkers)),
      type   = factor(type, levels = toupper(types)),
      value_plot = if (log10_y) ifelse(value > 0, log10(value), NA_real_) else value
    ) %>%
    filter(!is.na(value_plot))
  
  # Ensure same y-scale for WWTP vs D within each biomarker by using identical facet structure,
  # but separate plots per site. (Comparison happens visually because same transform/labels/theme.)
  make_one <- function(site_id) {
    
    site_label <- dplyr::case_when(
      site_id == "WWTP" ~ "WWTP",
      site_id == "D"    ~ "Neighborhood",
      TRUE              ~ site_id
    )
    
    ggplot(dfp %>% filter(site == site_id),
           aes(x = type, y = value_plot, fill = type)) +
      geom_boxplot(
        outlier.shape = NA,
        alpha = 0.35,

        linewidth = 0.7,
        fatten = 1.2
      ) +
      geom_point(
        position = position_jitter(width = 0.14, height = 0),
        alpha = 0.55,
        size = 1.9
      ) +
      facet_wrap(~ target, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = pal_type_fill, drop = FALSE) +
      scale_x_discrete(labels = fraction_labels) +
      labs(
        x = NULL,
        y = if (log10_y) "log10(gc / mL)" else "concentration",
        title = paste0("Sampling Location: ", site_label)
      ) +
      theme_pub() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(face = "bold")
      )
  }
  
  list(
    p_wwtp = make_one("WWTP"),
    p_d    = make_one("D")
  )
}

# =============================================================================
# 2) 2x2 time-series: TS & VS on same plot (dual y-axis), FLOW separate
#    - WWTP and D separately; stacked as 2x2 (top TS/VS, bottom FLOW)
#    - removes last date per site
#    - harmonizes x-limits across all four panels
# =============================================================================
plot_2x2_ts_vs_dual_and_flow_pub <- function(df_long,
                                             sites = c("WWTP", "D"),
                                             ts_target = "TS",
                                             vs_target = "VS",
                                             flow_target = "FLOW",
                                             drop_last_point = TRUE) {
  
  stopifnot(all(c("date","site","target","type","value") %in% names(df_long)))
  
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  df0 <- df_long %>%
    mutate(
      date   = as.Date(date),
      site   = toupper(trimws(as.character(site))),
      target = toupper(trimws(as.character(target))),
      value  = as.numeric(value)
    ) %>%
    filter(
      site %in% toupper(sites),
      target %in% c(ts_target, vs_target, flow_target),
      !is.na(date),
      !is.na(value)
    )
  
  # Drop last date per site (applies to TS/VS/FLOW together)
  if (isTRUE(drop_last_point)) {
    last_dates <- df0 %>%
      group_by(site) %>%
      summarise(last_date = max(date), .groups = "drop")
    
    df0 <- df0 %>%
      left_join(last_dates, by = "site") %>%
      filter(date < last_date) %>%
      dplyr::select(-last_date)
  }
  
  # Harmonize x-limits across panels (after dropping)
  x_rng <- range(df0$date, na.rm = TRUE)
  
  make_site_plots <- function(site_id) {
    
    dsite <- df0 %>% filter(site == site_id)

    
    d_ts <- dsite %>% filter(target == ts_target) %>%
      group_by(date) %>%
      summarise(TS = median(value, na.rm = TRUE), .groups = "drop")
    
    d_vs <- dsite %>% filter(target == vs_target) %>%
      group_by(date) %>%
      summarise(VS = median(value, na.rm = TRUE), .groups = "drop")
    
    d_flow <- dsite %>% filter(target == flow_target) %>%
      group_by(date) %>%
      summarise(FLOW = median(value, na.rm = TRUE), .groups = "drop")
    
    d_ts_vs <- full_join(d_ts, d_vs, by = "date") %>% arrange(date)
    
    # Scale factor for dual axis
    max_ts <- suppressWarnings(max(d_ts_vs$TS, na.rm = TRUE))
    max_vs <- suppressWarnings(max(d_ts_vs$VS, na.rm = TRUE))
    k <- if (is.finite(max_ts) && is.finite(max_vs) && max_vs > 0) max_ts / max_vs else 1
    
    # Annotation position
    ann_x <- x_rng[1]
    ann_y <- suppressWarnings(max(d_ts_vs$TS, na.rm = TRUE))
    
    p_ts_vs <- ggplot(d_ts_vs, aes(x = date)) +
      geom_line(aes(y = TS, color = "TS"), linewidth = 0.85, na.rm = TRUE) +
      geom_point(aes(y = TS, color = "TS"), size = 2.0, na.rm = TRUE) +
      geom_line(aes(y = VS * k, color = "VS"), linetype = "dashed",
                linewidth = 0.85, na.rm = TRUE) +
      geom_point(aes(y = VS * k, color = "VS"), shape = 1, size = 2.0, na.rm = TRUE) +
      scale_color_manual(values = c(TS = pal_hydro$TS, VS = pal_hydro$VS), name = NULL) +
      scale_y_continuous(
        name = "TS (mg / L)",
        sec.axis = sec_axis(~ . / k, name = "VS (mg / L)")
      ) +
      scale_x_date(limits = x_rng) +
      labs(
        #title = paste0(site_id, ": TS & VS"),
        x = NULL
      ) +
      annotate(
        "text",
        x = ann_x,
        y = ann_y,
        label = "VS scaled to TS axis",
        hjust = 0,
        vjust = -0.7,
        size = 3.0
      ) +
      theme_pub() +
      theme(
        legend.position = "top",
        legend.justification = "left",
        legend.direction = "horizontal"
      )
    
    p_flow <- ggplot(d_flow, aes(x = date, y = FLOW)) +
      geom_line(color = pal_hydro$FLOW, linewidth = 0.85, na.rm = TRUE) +
      geom_point(color = pal_hydro$FLOW, size = 2.0, na.rm = TRUE) +
      scale_x_date(limits = x_rng) +
      labs(
        #title = paste0(site_id, ": Flow rate"),
        x = NULL,
        y = "Flow rate (L / s)"
      ) +
      theme_pub()
    
    list(p_ts_vs = p_ts_vs, p_flow = p_flow)
  }
  
  pW <- make_site_plots("WWTP")
  pD <- make_site_plots("D")
  
  # 2x2 layout
  (pW$p_ts_vs | pD$p_ts_vs) /
    (pW$p_flow | pD$p_flow)
}

# =============================================================================
# 3) Combined figure:

#    TOP: boxplots (WWTP | D)
#    BOT: 2x2 TS/VS + Flow
#    + panel tags A, B, ...
# =============================================================================
plot_combined_biomarker_and_hydraulics_pub <- function(df_ww,
                                                       top_height = 1,
                                                       bottom_height = 0.8) {
  
  bx <- plot_biomarker_boxplots_WWTP_and_D_pub(df_ww)
  top_panel <- bx$p_wwtp | bx$p_d
  
  bottom_panel <- plot_2x2_ts_and_flow_pub(df_ww)
  
  (top_panel / bottom_panel) +
    plot_layout(heights = c(top_height, bottom_height)) +
    plot_annotation(
      tag_levels = "A",
      tag_suffix = ")",
      theme = theme(
        plot.margin = margin(6, 6, 6, 6),
        plot.tag = element_text(face = "bold")
      )
    )
}


#=======================================================
# =============================================================================
# Boxplots (with points) of log10(WWTP / D) ratios, filtered to a date window,
# arranged as:
#   - TOP ROW facets: N1 (raw) + N1 normalized (PMMOV, TS, VS, FLOW)
#   - BOTTOM ROW facets: N2 (raw) + N2 normalized (PMMOV, TS, VS, FLOW)
# By fraction/type: CF (solid), UF (liquid), TOTAL
#
# Requires df_long to already contain derived targets:
#   N1_PMMOV, N2_PMMOV, N1_TS, N2_TS, N1_VS, N2_VS, N1_FLOW, N2_FLOW
# =============================================================================
plot_log10_ratio_WWTP_over_D_boxplots <- function(
    df_long,
    box_start = as.Date("2021-03-22"),
    box_end   = as.Date("2021-05-18"),
    types = c("CF","UF","TOTAL"),
    show_points = TRUE,
    point_alpha = 0.55,
    box_alpha = 0.35
) {
  
  stopifnot(all(c("date","site","target","type","value") %in% names(df_long)))
  
  
  ratio_targets <- c(
    "N1", "N2", "PMMOV",
    "N1_PMMOV", "N2_PMMOV",
    "N1_TS", "N2_TS", "PMMOV_TS",
    "N1_FLOW_PC", "N2_FLOW_PC", "PMMOV_FLOW_PC"
  )
  
  df0 <- df_long %>%
    mutate(
      date   = as.Date(date),
      site   = toupper(trimws(as.character(site))),
      target = toupper(trimws(as.character(target))),
      type   = toupper(trimws(as.character(type))),
      value  = as.numeric(value)
    ) %>%
    filter(
      date >= box_start, date <= box_end,
      site %in% c("WWTP","D"),
      type %in% toupper(types),
      target %in% ratio_targets,
      !is.na(value)
    )
  
  # Pair WWTP and D per date/target/type then compute log10 ratio
  df_ratio <- df0 %>%
    group_by(date, site, target, type) %>%
    summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%

    pivot_wider(names_from = site, values_from = value) %>%
    mutate(
      log10_ratio = dplyr::if_else(
        !is.na(WWTP) & !is.na(D) & WWTP > 0 & D > 0,
        log10(WWTP / D),
        NA_real_
      )
    ) %>%
    filter(!is.na(log10_ratio)) %>%
    mutate(
      type = factor(type, levels = toupper(types)),
      virus = dplyr::case_when(
        grepl("^N1", target) ~ "N1",
        grepl("^N2", target) ~ "N2",
        grepl("^PMMOV", target) ~ "PMMoV",
        TRUE ~ NA_character_
      ),
      norm = dplyr::case_when(
        target %in% c("N1","N2","PMMOV") ~ "raw",
        grepl("_PMMOV$", target) ~ "PMMoV",
        grepl("_TS$", target)    ~ "TS",
        grepl("_VS$", target)    ~ "VS",
        grepl("_FLOW_PC$", target)  ~ "Flow (per cap)",
        TRUE ~ "other"
      ),
      norm = factor(norm, levels = c("raw","PMMoV","TS","VS","Flow (per cap)")),
      virus = factor(virus, levels = c("N1","N2","PMMoV"))
    ) %>%
    filter(!is.na(virus), norm %in% levels(norm))
  
  # Revised manuscript layout: normalization on x, sample fraction in columns.
  p <- ggplot(df_ratio, aes(x = norm, y = log10_ratio, fill = type)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.6) +
    geom_boxplot(outlier.shape = NA, alpha = box_alpha, linewidth = 0.7, fatten = 1.2) +
    { if (show_points)
      geom_point(
        aes(color = type),
        position = position_jitter(width = 0.14, height = 0),
        alpha = point_alpha,
        size = 1.9
      ) else NULL } +
    facet_grid(
      rows = vars(virus), cols = vars(type), scales = "free_y",
      labeller = labeller(type = c(CF = "Solids", UF = "Liquid", TOTAL = "Total"))
    ) +
    scale_fill_manual(values = c(CF = "#08519c", UF = "#6baed6", TOTAL = "#8c510a"), drop = FALSE) +
    scale_color_manual(values = c(CF = "#08519c", UF = "#6baed6", TOTAL = "#8c510a"), drop = FALSE) +
    labs(
      x = NULL,
      y = expression(log[10]*"(WWTP / neighborhood)")
    ) +
    theme_pub() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "bold", angle = 0, vjust = 0.5),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.9, "lines")
    )
  
  list(data = df_ratio, plot = p)
}



#=============== Estimate Loss =============================
estimate_IQR_loss_from_log10_ratio <- function(df_ratio) {
  
  stopifnot(all(c("virus","norm","type","log10_ratio") %in% names(df_ratio)))
  
  df_ratio %>%
    dplyr::group_by(virus, norm, type) %>%
    dplyr::summarise(
      n = sum(!is.na(log10_ratio)),
      
      q25_log10 = quantile(log10_ratio, 0.25, na.rm = TRUE),
      q50_log10 = quantile(log10_ratio, 0.50, na.rm = TRUE),
      q75_log10 = quantile(log10_ratio, 0.75, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    dplyr::mutate(

      # Convert back to linear ratios
      ratio_q25 = 10^q25_log10,
      ratio_q50 = 10^q50_log10,
      ratio_q75 = 10^q75_log10,
      
      # Estimated signal loss (%)
      loss_q25_pct = (1 - ratio_q25) * 100,
      loss_q50_pct = (1 - ratio_q50) * 100,
      loss_q75_pct = (1 - ratio_q75) * 100
    ) %>%
    dplyr::arrange(virus, norm, type)
}



#============= Heatmap of loss
plot_loss_heatmap <- function(df_ratio, facet_var = "type") {
  
  df_sum <- estimate_IQR_loss_from_log10_ratio(df_ratio) %>%
    dplyr::mutate(
      virus_chr = as.character(virus),
      norm_chr  = as.character(norm),
      facet_chr = as.character(.data[[facet_var]])
    )
  
  # ---- Detect PMMoV spelling (case-insensitive)
  all_levels  <- unique(c(df_sum$virus_chr, df_sum$norm_chr))
  pmmov_match <- all_levels[toupper(trimws(all_levels)) == "PMMOV"][1]
  
  # ---- Biomarker order
  desired <- c("N1", "N2", pmmov_match)
  virus_levels <- c(desired, setdiff(unique(df_sum$virus_chr), desired))
  
  # ---- Clean normalization labels for display/order
  df_sum <- df_sum %>%
    dplyr::mutate(
      norm_clean = dplyr::case_when(
        toupper(trimws(norm_chr)) %in% c("rAW", "UNNORMALIZED") ~ "raw",
        toupper(trimws(norm_chr)) == "PMMOV" ~ "PMMoV",
        toupper(trimws(norm_chr)) == "TS" ~ "TS",
        TRUE ~ norm_chr
      )
    )
  
  norm_levels <- c(
    "raw", "PMMoV", "TS",
    setdiff(unique(df_sum$norm_clean), c("raw", "PMMoV", "TS"))
  )
  
  # ---- Facet order
  df_sum <- df_sum %>%
    dplyr::mutate(
      virus_plot = factor(virus_chr, levels = virus_levels),
      norm_plot  = factor(norm_clean, levels = norm_levels),
      facet_key  = dplyr::case_when(
        facet_chr == "CF"    ~ "1_CF",
        facet_chr == "UF"    ~ "2_UF",
        facet_chr == "TOTAL" ~ "3_TOTAL",
        TRUE                 ~ paste0("9_", facet_chr)
      ),
      facet_key = factor(facet_key, levels = c("1_CF", "2_UF", "3_TOTAL")),
      y_dummy = 1
    )
  
  # ---- Grey PMMoV/PMMoV tile in every facet column
  overlay_df <- tibble::tibble(
    virus_plot = factor(pmmov_match, levels = virus_levels),
    norm_plot  = factor("PMMoV", levels = levels(df_sum$norm_plot)),
    y_dummy    = 1
  ) %>%
    tidyr::crossing(facet_key = levels(df_sum$facet_key)) %>%
    dplyr::filter(!is.na(facet_key))
  
  ggplot(df_sum, aes(x = norm_plot, y = y_dummy, fill = loss_q50_pct)) +
    geom_tile(color = "white", linewidth = 0.4) +
    
    geom_tile(
      data = overlay_df,
      inherit.aes = FALSE,
      aes(x = norm_plot, y = y_dummy),

      fill = "grey80",
      color = "white",
      linewidth = 0.4
    ) +
    
    geom_text(
      data = df_sum %>%
        dplyr::filter(!(toupper(trimws(virus_chr)) == "PMMOV" &
                          toupper(trimws(norm_chr)) == "PMMOV")),
      aes(label = sprintf("%.0f%%\n[%.0f, %.0f]",
                          loss_q50_pct, loss_q25_pct, loss_q75_pct)),
      size = 3
    ) +
    
    facet_grid(
      virus_plot ~ facet_key,
      labeller = ggplot2::labeller(
        facet_key = c(
          "1_CF"    = "Solids",
          "2_UF"    = "Liquid",
          "3_TOTAL" = "Total"
        )
      )
    ) +
    
    scale_fill_gradientn(
      colours = c("grey95", "red3", "red4"),
      values = scales::rescale(c(0, 50, 100)),
      limits = c(0, 100),
      name = "Median apparent loss (%)"
    ) +
    
    scale_y_continuous(NULL, breaks = NULL) +
    labs(x = NULL, y = NULL) +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      strip.background = element_rect(fill = "grey96", colour = NA),
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title.align = 0.5
    ) +
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    ))
}
#===== combined figure for results
combine_boxplot_and_loss_heatmap <- function(df_raw,
                                             heatmap_fun = plot_loss_heatmap,
                                             heights = c(1.2, 1),
                                             add_tags = TRUE) {
  
  res <- plot_log10_ratio_WWTP_over_D_boxplots(df_raw)
  
  p_top    <- res$plot
  p_bottom <- heatmap_fun(res$data)
  
  combined <- p_top / p_bottom +
    patchwork::plot_layout(heights = heights)
  
  if (add_tags) {
    combined <- combined +
      patchwork::plot_annotation(tag_levels = "A")
  }
  
  return(combined)
  
  
}

#------------------------ just flow time series
plot_2x2_ts_and_flow_pub <- function(df_long,
                                     sites = c("WWTP", "D"),
                                     ts_target = "TS",
                                     flow_target = "FLOW",
                                     drop_last_point = TRUE) {
  
  stopifnot(all(c("date","site","target","type","value") %in% names(df_long)))
  
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  df0 <- df_long %>%
    mutate(
      date   = as.Date(date),
      site   = toupper(trimws(as.character(site))),
      target = toupper(trimws(as.character(target))),
      value  = as.numeric(value)
    ) %>%
    filter(
      site %in% toupper(sites),
      target %in% c(ts_target, flow_target),
      !is.na(date),
      !is.na(value)
    )
  
  # Drop last date per site (applies to TS/FLOW together)
  if (isTRUE(drop_last_point)) {
    last_dates <- df0 %>%
      group_by(site) %>%
      summarise(last_date = max(date), .groups = "drop")
    
    df0 <- df0 %>%
      left_join(last_dates, by = "site") %>%
      filter(date < last_date) %>%
      dplyr::select(-last_date)
  }
  
  # Harmonize x-limits across panels (after dropping)
  x_rng <- range(df0$date, na.rm = TRUE)
  
  make_site_plots <- function(site_id) {
    
    dsite <- df0 %>% filter(site == site_id)
    
    d_ts <- dsite %>%
      filter(target == ts_target) %>%
      group_by(date) %>%
      summarise(TS = median(value, na.rm = TRUE), .groups = "drop")
    
    d_flow <- dsite %>%
      filter(target == flow_target) %>%
      group_by(date) %>%
      summarise(FLOW = median(value, na.rm = TRUE), .groups = "drop")
    
    p_ts <- ggplot(d_ts, aes(x = date, y = TS)) +
      geom_line(color = pal_hydro$TS, linewidth = 0.85, na.rm = TRUE) +
      geom_point(color = pal_hydro$TS, size = 2.0, na.rm = TRUE) +
      scale_x_date(limits = x_rng) +
      labs(x = NULL, y = "TS (mg / L)") +
      theme_pub() +
      theme(
        legend.position = "none"
      )
    
    p_flow <- ggplot(d_flow, aes(x = date, y = FLOW)) +
      geom_line(color = pal_hydro$FLOW, linewidth = 0.85, na.rm = TRUE) +
      geom_point(color = pal_hydro$FLOW, size = 2.0, na.rm = TRUE) +
      scale_x_date(limits = x_rng) +
      labs(x = NULL, y = "Flow rate (L / s)") +
      theme_pub()
    
    list(p_ts = p_ts, p_flow = p_flow)

  }
  
  pW <- make_site_plots("WWTP")
  pD <- make_site_plots("D")
  
  # 2x2 layout
  (pW$p_ts | pD$p_ts) /
    (pW$p_flow | pD$p_flow)
}
