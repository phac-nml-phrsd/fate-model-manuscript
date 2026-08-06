# Revised core manuscript figures. Source this after utils_figures.R so the
# revised histogram definition overrides the legacy public implementation.

figure_loss_histo_summary <- function(summary_data, filename = NULL) {
    dat <- if (is.character(summary_data)) readr::read_csv(summary_data, show_col_types = FALSE) else summary_data
    dat <- dat %>% dplyr::mutate(
        wwtp = stringr::str_to_title(wwtp),
        fate_process = factor(fate_process, levels = c("Settling / Resuspension", "Biofilm Adsorption", "Biodegradation"))
    )
    stats <- dat %>% dplyr::distinct(wwtp, fate_process, mean_loss_proportion, median_loss_proportion)
    g <- ggplot2::ggplot(dat, ggplot2::aes(bin_mid_proportion, flow_path_count)) +
        ggplot2::geom_col(width = 0.05, fill = "#7da0b1", color = "white") +
        ggplot2::geom_vline(data = stats, ggplot2::aes(xintercept = mean_loss_proportion), color = "#1f78b4", linetype = "dashed") +
        ggplot2::geom_vline(data = stats, ggplot2::aes(xintercept = median_loss_proportion), color = "#ff7f00") +
        ggplot2::facet_grid(wwtp ~ fate_process, scales = "free_y") +
        ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
        ggplot2::labs(x = "Proportion of in-sewer viral biomarker loss", y = "Count of flow paths (500 simulations)") +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "#2b4c7e"), strip.text = ggplot2::element_text(color = "white", face = "bold"))
    if (!is.null(filename)) ggplot2::ggsave(filename, g, width = 9, height = 7, dpi = 250, bg = "white")
    g
}

figure_loss_histo_fate_processes_summary <- function(fig1, summary_data, filename) {
    fate_img <- magick::image_read_pdf(fig1, density = 300)[1]
    fate_img <- magick::image_background(magick::image_trim(fate_img), "white", flatten = TRUE)
    info <- magick::image_info(fate_img)
    g_fate <- ggplot2::ggplot() +
        ggplot2::annotation_raster(as.raster(fate_img), 0, info$width[1], 0, info$height[1]) +
        ggplot2::coord_fixed(xlim = c(0, info$width[1]), ylim = c(0, info$height[1]), expand = FALSE) +
        ggplot2::theme_void()
    g <- (g_fate / figure_loss_histo_summary(summary_data)) +
        patchwork::plot_layout(heights = c(1, 0.8)) + patchwork::plot_annotation(tag_levels = "A")
    ggplot2::ggsave(filename, g, width = 8, height = 9, dpi = 250, bg = "white")
    g
}

figure_sensitivity_loss_histo_summary <- function(rna_summary, conduit_summary, filename) {
    dat <- readr::read_csv(rna_summary, show_col_types = FALSE) %>%
        dplyr::mutate(wwtp = stringr::str_to_title(wwtp), scenario = factor(scenario, levels = c("Skewed (Base)", "Highly Skewed", "Homogeneous")))
    g_set <- ggplot2::ggplot(dat, ggplot2::aes(bin_mid, count)) +
        ggplot2::geom_col(width = 0.05, fill = "#7da0b1", color = "white") +
        ggplot2::facet_grid(wwtp ~ scenario, scales = "free_y") +
        ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
        ggplot2::labs(x = "Proportion of in-sewer settling/resuspension viral loss", y = "Count of flow paths (500 simulations)") +
        ggplot2::theme_bw()
    g_total <- get_histo_plot(conduit_summary, type = "total")
    g <- (g_set / g_total) + patchwork::plot_annotation(tag_levels = "A")
    ggplot2::ggsave(filename, g, width = 7.5, height = 8, dpi = 250, bg = "white")
    g
}

figure_loss_histo <- function(sim_df_loss_total, filename = NULL) {
    dat <- helper_sim_df(sim_df_loss_total) %>% mutate(nameplot = case_when(name == "loss.set.path.slow" ~ "Settling / Resuspension", name == "loss.bio.path" ~ "Biofilm Adsorption", name == "loss.deg.path" ~ "Biodegradation"), wwtp = case_when(wwtp == "north" ~ "North", wwtp == "south" ~ "South", wwtp == "west" ~ "West", TRUE ~ wwtp))
    stats_df <- dat %>% group_by(wwtp, nameplot) %>% summarise(mean_val = mean(value, na.rm = TRUE), median_val = median(value, na.rm = TRUE), .groups = "drop")
    g <- dat %>% ggplot(aes(x = value)) + geom_histogram(binwidth = 0.05, fill = "#7da0b1", color = "white", alpha = 0.9) + geom_vline(data = stats_df, aes(xintercept = mean_val), color = "#1f78b4", linetype = "dashed", size = 1) + geom_vline(data = stats_df, aes(xintercept = median_val), color = "#ff7f00", linetype = "solid", size = 1) + facet_grid(wwtp ~ nameplot, scales = "free_y") + theme_bw(base_size = 14) + theme(panel.grid.major.x = element_line(color = "gray90"), panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_rect(fill = "#2b4c7e", color = NA), strip.text = element_text(color = "white", face = "bold", size = 13), axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)), axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)), axis.text = element_text(size = 12), plot.margin = margin(12, 18, 10, 10)) + scale_x_continuous(breaks = seq(0, 1, 0.25), labels = scales::label_number(accuracy = 0.01)) + 
        labs(x = "Proportion of in-sewer viral biomarker loss", y = "Count of flow paths (500 simulations)")
    return(g)
    stats_df
    if (!is.null(filename) && length(filename) == 1 && !is.na(filename) && nzchar(filename)) {
        ggsave(filename, plot = g, width = 9, height = 7, dpi = 250, units = "in", bg = "white", device = "png")
    }
} 

figure_biofilm_pipe <- function(df, filename) {
    df = df %>% mutate(wwtp = str_to_title(wwtp))
    stats <- df %>% group_by(wwtp) %>% summarise(mean_val = mean(conduit_av, na.rm = TRUE), median_val = median(conduit_av, na.rm = TRUE))
    g = ggplot(df, aes(x = conduit_av, fill = wwtp)) + geom_histogram(color = "black") + geom_vline(data = stats, aes(xintercept = mean_val), color = "red", linetype = "dashed", size = 1) + geom_vline(data = stats, aes(xintercept = median_val), color = "purple", linetype = "dotted", size = 1) + scale_fill_manual(values = c(North = "#DEEBF7", South = "#9ECAE1", West = "#3182BD")) + facet_wrap(~wwtp, scales = "free_y") + labs(x = expression("Biofilm Area/Volume ratio (1/" * m * ")"), y = "Counts of Conduits") + 
        geom_text(data = stats, aes(x = mean_val, y = Inf, label = paste("Mean:", round(mean_val, 2))), color = "red", vjust = 3.5, hjust = -0.25) + geom_text(data = stats, aes(x = median_val, y = Inf, label = paste("Median:", round(median_val, 2))), color = "purple", vjust = 5, hjust = -0.3) + theme_minimal() + theme(strip.text = element_text(face = "bold"))
    plot(g)
    ggsave(filename, plot = g, width = 8, height = 4, dpi = 400, units = "in", bg = "white", device = "png")
} 

figure_loss_histo_fate_processes <- function(fig1 = "doc/figs/figure_fate_processes.pdf", sim_df_loss_total, filename = "doc/figs/figure_loss_histo_fate.png", panel_widths = c(1, 0.8), width = 8, height = 9) {
    if (!file.exists(fig1)) {
        stop("Figure file does not exist: ", fig1)
    }
    extension <- tolower(tools::file_ext(fig1))
    if (extension == "pdf") {
        fate_img <- magick::image_read_pdf(path = fig1, density = 300)
        fate_img <- fate_img[1]
    }
    else if (extension %in% c("png", "jpg", "jpeg", "tif", "tiff")) {
        fate_img <- magick::image_read(fig1)
        fate_img <- fate_img[1]
    }
    else {
        stop("Unsupported file format: ", extension, ". Use PDF, PNG, JPG, JPEG, TIF, or TIFF.")
    }
    if (!inherits(fate_img, "magick-image")) {
        stop("The fate-process figure could not be read as a Magick image.")
    }
    fate_img <- magick::image_trim(fate_img)
    fate_img <- magick::image_background(image = fate_img, color = "white", flatten = TRUE)
    fate_info <- magick::image_info(fate_img)
    fate_width <- fate_info$width[1]
    fate_height <- fate_info$height[1]
    fate_raster <- as.raster(fate_img)
    g.fate <- ggplot2::ggplot() + ggplot2::annotation_raster(raster = fate_raster, xmin = 0, xmax = fate_width, ymin = 0, ymax = fate_height) + ggplot2::coord_fixed(ratio = 1, xlim = c(0, fate_width), ylim = c(0, fate_height), expand = FALSE, clip = "off") + ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))
    g.loss <- figure_loss_histo(sim_df_loss_total = sim_df_loss_total, filename = NULL)
    if (!inherits(g.loss, c("ggplot", "patchwork"))) {
        stop("figure_loss_histo() did not return a ggplot or patchwork object. ", "Add return(g.histo) at the end of figure_loss_histo().")
    }
    g.loss <- g.loss + ggplot2::theme(plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))
    g.combo <- (patchwork::wrap_plots(g.fate, g.loss, ncol = 1, heights = panel_widths) + patchwork::plot_annotation(tag_levels = "A")) & ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 18, colour = "black", hjust = 0, vjust = 1), plot.tag.position = c(0.005, 0.995))
    output_directory <- dirname(filename)
    if (!dir.exists(output_directory)) {
        dir.create(output_directory, recursive = TRUE)
    }
    ggsave(filename, plot = g.combo, width = width, height = height, dpi = 250, units = "in", bg = "white", device = "png")
    return(g.combo)
}
