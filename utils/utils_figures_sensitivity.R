# Revised stochastic surveillance-coverage and minimum-detection figures.
# Extracted from the manuscript source and kept separate from core plotting helpers.

plot_particle_loss_distribution_summary <- function(summary_file, labels_file, level, filename) {
    dat <- readr::read_csv(summary_file, show_col_types = FALSE) %>%
        dplyr::filter(hydraulic_level == level)
    labels <- readRDS(labels_file)
    velocity_labels <- stats::setNames(paste0("Class ", labels$part.class, "\n(V=", labels$vel.set, " m/d)"), labels$part.class)
    g <- ggplot2::ggplot(dat, ggplot2::aes(bin_mid_proportion, observation_count, fill = wwtp)) +
        ggplot2::geom_col(width = 0.05) +
        ggplot2::facet_grid(wwtp ~ particle_class, scales = "free_y", labeller = ggplot2::labeller(particle_class = ggplot2::as_labeller(velocity_labels))) +
        ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
        ggplot2::labs(x = paste("Proportion of solid particle lost per", tolower(level)), y = "Frequency") +
        ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "none")
    ggplot2::ggsave(filename, g, width = 9, height = 4)
    g
}

calculate_simulation_boxplot_metrics <- function(sim_loss, df_polygons, demographic, fsa_polygons, stochastic_parameters) {
    required_loss_columns <- c("sim", "node_id")
    missing_loss_columns <- setdiff(required_loss_columns, names(sim_loss))
    if (length(missing_loss_columns) > 0) {
        stop("sim_loss is missing required column(s): ", paste(missing_loss_columns, collapse = ", "))
    }
    if (!"sim" %in% names(stochastic_parameters)) {
        stop("stochastic_parameters must contain a sim column.")
    }
    polygon_lookup <- df_polygons %>% dplyr::select(node_id, geometry) %>% dplyr::distinct(node_id, .keep_all = TRUE)
    simulation_ids <- sort(unique(sim_loss$sim))
    results <- lapply(simulation_ids, function(sim_id) {
        message("Calculating boxplot metrics for simulation ", sim_id, "...")
        loss_one_sim <- sim_loss %>% dplyr::filter(sim == sim_id)
        parameters_one_sim <- stochastic_parameters %>% dplyr::filter(sim == sim_id)
        if (nrow(parameters_one_sim) == 0) {
            stop("No stochastic parameters found for simulation ", sim_id, ".")
        }
        loss_one_sim_mean <- get_mean_sim(loss_one_sim) %>% dplyr::left_join(polygon_lookup, by = "node_id")
        loss_fsa_one_sim <- calcu_loss_fsa(loss = loss_one_sim_mean, demo = demographic, wpg = fsa_polygons, total.parms = parameters_one_sim)
        population_one_sim <- calculate_pop_loss(loss_fsa_one_sim, parameters_one_sim)
        calcu_inf_rate(population_one_sim, parameters_one_sim) %>% dplyr::mutate(sim = sim_id)
    })
    dplyr::bind_rows(results)
} 

boxplot_eff_pop_coverage <- function(df) {
    df_wwtp <- df %>% sf::st_drop_geometry() %>% dplyr::transmute(sim, wwtp = stringr::str_to_title(wwtp), Liquid = eff.pop.bio.deg.path.wwtp, Solids = eff.pop.set.slow.wwtp, Total = eff.pop.total.slow.wwtp) %>% dplyr::distinct(sim, wwtp, .keep_all = TRUE)
    df_city <- df %>% sf::st_drop_geometry() %>% dplyr::transmute(sim, wwtp = "City", Liquid = eff.pop.bio.deg.path.city, Solids = eff.pop.set.slow.city, Total = eff.pop.total.slow.city) %>% dplyr::distinct(sim, .keep_all = TRUE)
    df_long <- dplyr::bind_rows(df_wwtp, df_city) %>% tidyr::pivot_longer(cols = c(Solids, Liquid, Total), names_to = "fraction", values_to = "value") %>% dplyr::mutate(fraction = factor(fraction, levels = c("Solids", "Liquid", "Total")))
    ggplot2::ggplot(df_long, ggplot2::aes(x = wwtp, y = value, fill = wwtp)) + ggplot2::geom_boxplot(width = 0.6, outlier.alpha = 0.25, outlier.size = 1) + ggplot2::facet_grid(~fraction) + ggplot2::scale_fill_manual(values = c(North = "#DEEBF7", South = "#9ECAE1", West = "#3182BD", City = "#33a02c")) + ggplot2::scale_y_continuous(limits = c(0, 100)) + ggplot2::labs(x = "City and catchments of wastewater treatment plants", y = "Effective population coverage (%)") + ggplot2::guides(fill = "none") + 
        ggplot2::theme(legend.position = "none", panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", color = "white"), strip.background = ggplot2::element_rect(fill = "steelblue4"), axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, face = "bold"), axis.text.y = ggplot2::element_text(face = "bold"), axis.title = ggplot2::element_text(face = "bold"))
} 

boxplot_detection_fractiontype <- function(df) {
    df_wwtp <- df %>% sf::st_drop_geometry() %>% dplyr::transmute(sim, wwtp = stringr::str_to_title(wwtp), Liquid = wwtp.inf.rate.bio.deg.highshed, Solids = wwtp.inf.rate.set.slow.highshed, Total = wwtp.inf.rate.slow.highshed) %>% dplyr::distinct(sim, wwtp, .keep_all = TRUE)
    df_city <- df %>% sf::st_drop_geometry() %>% dplyr::transmute(sim, wwtp = "City", Liquid = city.inf.rate.bio.deg.highshed, Solids = city.inf.rate.set.slow.highshed, Total = city.inf.rate.slow.highshed) %>% dplyr::distinct(sim, .keep_all = TRUE)
    df_long <- dplyr::bind_rows(df_wwtp, df_city) %>% tidyr::pivot_longer(cols = c(Solids, Liquid, Total), names_to = "fraction", values_to = "value") %>% dplyr::mutate(fraction = factor(fraction, levels = c("Solids", "Liquid", "Total")))
    ggplot2::ggplot(df_long, ggplot2::aes(x = wwtp, y = value, fill = wwtp)) + ggplot2::geom_boxplot(width = 0.6, outlier.alpha = 0.25, outlier.size = 1) + ggplot2::facet_wrap(~fraction, scales = "free_y", nrow = 1) + ggplot2::scale_fill_manual(values = c(North = "#DEEBF7", South = "#9ECAE1", West = "#3182BD", City = "goldenrod3")) + ggplot2::labs(x = "City and catchments of wastewater treatment plants", y = "Minimum infection rate\nfor detection (%)") + ggplot2::guides(fill = "none") + ggplot2::theme(legend.position = "none", 
        panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", color = "white"), strip.background = ggplot2::element_rect(fill = "steelblue4"), axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, face = "bold"), axis.text.y = ggplot2::element_text(face = "bold"), axis.title = ggplot2::element_text(face = "bold"))
} 

summarise_stochastic_boxplot <- function(df, type = c("population_coverage", "infection_rate")) {
    type <- match.arg(type)
    if (type == "population_coverage") {
        df_wwtp <- df %>% sf::st_drop_geometry() %>% dplyr::transmute(sim, location = stringr::str_to_title(wwtp), Solids = eff.pop.set.slow.wwtp, Liquid = eff.pop.bio.deg.path.wwtp, Total = eff.pop.total.slow.wwtp) %>% dplyr::distinct(sim, location, .keep_all = TRUE)
        df_city <- df %>% sf::st_drop_geometry() %>% dplyr::transmute(sim, location = "City", Solids = eff.pop.set.slow.city, Liquid = eff.pop.bio.deg.path.city, Total = eff.pop.total.slow.city) %>% dplyr::distinct(sim, .keep_all = TRUE)
    }
    else {
        df_wwtp <- df %>% sf::st_drop_geometry() %>% dplyr::transmute(sim, location = stringr::str_to_title(wwtp), Solids = wwtp.inf.rate.set.slow.highshed, Liquid = wwtp.inf.rate.bio.deg.highshed, Total = wwtp.inf.rate.slow.highshed) %>% dplyr::distinct(sim, location, .keep_all = TRUE)
        df_city <- df %>% sf::st_drop_geometry() %>% dplyr::transmute(sim, location = "City", Solids = city.inf.rate.set.slow.highshed, Liquid = city.inf.rate.bio.deg.highshed, Total = city.inf.rate.slow.highshed) %>% dplyr::distinct(sim, .keep_all = TRUE)
    }
    boxplot_long <- dplyr::bind_rows(df_wwtp, df_city) %>% tidyr::pivot_longer(cols = c(Solids, Liquid, Total), names_to = "fraction", values_to = "value") %>% dplyr::mutate(fraction = factor(fraction, levels = c("Solids", "Liquid", "Total")))
    boxplot_long %>% dplyr::filter(!is.na(value)) %>% dplyr::group_by(fraction, location) %>% dplyr::summarise(n = dplyr::n(), minimum = min(value), q1 = as.numeric(stats::quantile(value, 0.25, names = FALSE)), median = stats::median(value), q3 = as.numeric(stats::quantile(value, 0.75, names = FALSE)), maximum = max(value), iqr = stats::IQR(value), lower_whisker = min(value[value >= q1 - 1.5 * iqr]), upper_whisker = max(value[value <= q3 + 1.5 * iqr]), .groups = "drop") %>% dplyr::select(fraction, 
        location, n, minimum, lower_whisker, q1, median, q3, upper_whisker, maximum, iqr) %>% dplyr::arrange(fraction, location)
} 

summarise_eff_pop_coverage_boxplot <- function(df) {
    summarise_stochastic_boxplot(df, type = "population_coverage")
} 

summarise_detection_boxplot <- function(df) {
    summarise_stochastic_boxplot(df, type = "infection_rate")
} 

figure_combo_map_loss_fsa_boxplot <- function(df_map, df_boxplot, filename) {
    size.text <- 2
    col.fsa <- "white"
    g.pipe <- figure_map_loss_fsa(df_map, size.text, col.fsa)
    g.pop <- figure_map_pop_loss_fsa(df_map, size.text, col.fsa)
    g.box <- boxplot_eff_pop_coverage(df_boxplot)
    top_row <- g.pipe + g.pop + patchwork::plot_layout(ncol = 2)
    g.combo <- top_row/g.box + patchwork::plot_layout(heights = c(2, 1)) + patchwork::plot_annotation(tag_levels = "A")
    ggplot2::ggsave(filename, plot = g.combo, width = 9, height = 10, dpi = 400, units = "in", bg = "white", device = "png")
    return(g.combo)
} 

figure_combo_detection_boxplot <- function(df_map, df_boxplot, filename) {
    g.map <- figure_map_inf_rate(df_map, shed.type = "high") + ggplot2::theme(legend.position = "right", legend.justification = "center", legend.box.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5), plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5))
    g.box <- boxplot_detection_fractiontype(df_boxplot)
    top_row <- g.map + patchwork::guide_area() + patchwork::plot_layout(ncol = 2, widths = c(3, 1.25), guides = "collect")
    g.combo <- top_row/g.box + patchwork::plot_layout(heights = c(1.5, 1)) + patchwork::plot_annotation(tag_levels = "A")
    ggplot2::ggsave(filename, plot = g.combo, width = 7, height = 6, dpi = 250, units = "in", bg = "white", device = "png")
    return(g.combo)
}
