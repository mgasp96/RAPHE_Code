# Plotting helpers used by Sequencing/scripts (001_QC, 003_Annotation_analysis).
# Minimal version: only the functions these scripts call.

QC_VlnPlot <- function(seurat_obj, hline_values = c(1000, 2000, 10)) {
  features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  
  # Generate violin plots without combining them into one plot
  vln_plots <- VlnPlot(seurat_obj, features = features, ncol = 3, 
                       cols = "#025e8d", pt.size = 0.1, combine = FALSE)
  
  # Add jitter, horizontal lines, and annotations to each plot
  vln_plots <- lapply(seq_along(vln_plots), function(i) {
    vln_plots[[i]] + 
      geom_jitter(width = 0.2, size = 0.5, alpha = 0.2, color = "#363633") +
      geom_hline(yintercept = hline_values[i], linetype = "dashed", 
                 color = "#be1818", size = 0.7) +
      annotate("text", x = 0.5, 
               y = hline_values[i] + (0.25 * hline_values[i]), 
               label = paste0(hline_values[i]), color = "#be1818", size = 4)
  })
  
  # Arrange plots in a single row using patchwork and return the result
  wrap_plots(vln_plots, ncol = 3)
}


QC_HistPlot <- function(seurat_obj) {
  # List of features to plot
  features <- c("nFeature_RNA", "nCount_RNA")
  
  # Function to generate the histogram plot for a single feature
  make_hist_plot <- function(feature) {
    # Extract the data from the metadata
    dat <- seurat_obj@meta.data[[feature]]
    
    # Define breaks in steps of 100
    brks <- seq(0, max(dat) + 100, by = 100)
    
    # Compute the histogram (without plotting)
    hist_data <- hist(dat, breaks = brks, plot = FALSE)
    
    # Log-transform the counts (adding 1 to avoid log10(0))
    log_counts <- log10(hist_data$counts + 1)
    
    # Smooth the log-transformed histogram using a spline
    fit <- smooth.spline(x = hist_data$mids, y = log_counts, spar = 0.5)
    smoothed <- fit$y
    
    # Identify local minima in the smoothed curve
    # Local minima are identified as positions where the sign of the second difference is 2.
    local_min_idx <- which(diff(sign(diff(smoothed))) == 2) + 1
    # Convert indices to x-values using histogram midpoints
    local_min_vals <- hist_data$mids[local_min_idx]
    
    # Create a data frame for plotting
    df_hist <- data.frame(mid = hist_data$mids,
                          log_counts = log_counts)
    
    # Determine a y position for annotations (90% of maximum log-count)
    label_y <- max(df_hist$log_counts) * 0.9
    
    # Build the plot
    p <- ggplot(df_hist, aes(x = mid, y = log_counts)) +
      geom_col(fill = "#025e8d") +
      labs(x = feature, y = "Log10(Frequency + 1)", 
           title = paste("Log-Transformed Histogram of", feature)) +
      theme_minimal()
    
    # Add vertical lines and annotations for each local minimum
    for(i in seq_along(local_min_vals)){
      p <- p +
        geom_vline(xintercept = local_min_vals[i], linetype = "dashed", 
                   color = "#be1818", size = 1) +
        annotate("text", x = local_min_vals[i], y = label_y, 
                 label = paste("min:", round(local_min_vals[i], 1)), 
                 color = "#be1818", vjust = -0.5, angle = 90, size = 3)
    }
    return(p)
  }
  
  # Generate a list of plots for the two features
  plot_list <- lapply(features, make_hist_plot)
  
  # Combine the plots into a patchwork object (side-by-side)
  combined_plot <- wrap_plots(plot_list, ncol = 2)
  
  return(combined_plot)
}



QC_ScatterPlots <- function(seurat_obj, 
                            low_nFeature = 1500, 
                            high_nFeature = 9000, 
                            mt_threshold = 10) {
  # Extract metadata for convenience
  md <- seurat_obj@meta.data
  
  # Determine extremes for annotations
  max_nFeature <- max(md$nFeature_RNA, na.rm = TRUE)
  max_percent  <- max(md$percent.mt, na.rm = TRUE)
  min_percent  <- min(md$percent.mt, na.rm = TRUE)
  
  # Calculate counts for each region based on thresholds
  # Lower row: percent.mt < mt_threshold
  reg1 <- sum(md$nFeature_RNA < low_nFeature & md$percent.mt < mt_threshold, na.rm = TRUE)
  reg2 <- sum(md$nFeature_RNA >= low_nFeature & md$nFeature_RNA < high_nFeature & md$percent.mt < mt_threshold, na.rm = TRUE)
  reg3 <- sum(md$nFeature_RNA >= high_nFeature & md$percent.mt < mt_threshold, na.rm = TRUE)
  # Upper row: percent.mt >= mt_threshold
  reg4 <- sum(md$nFeature_RNA < low_nFeature & md$percent.mt >= mt_threshold, na.rm = TRUE)
  reg5 <- sum(md$nFeature_RNA >= low_nFeature & md$nFeature_RNA < high_nFeature & md$percent.mt >= mt_threshold, na.rm = TRUE)
  reg6 <- sum(md$nFeature_RNA >= high_nFeature & md$percent.mt >= mt_threshold, na.rm = TRUE)
  
  # Create p1: FeatureScatter of nFeature_RNA vs. percent.mt
  p1 <- FeatureScatter(seurat_obj, 
                       feature1 = "nFeature_RNA", 
                       feature2 = "percent.mt", 
                       cols = "#025e8d") +
    # Add horizontal and vertical dashed lines
    geom_hline(yintercept = mt_threshold, linetype = "dashed", color = "#be1818", size = 1) +
    geom_vline(xintercept = low_nFeature, linetype = "dashed", color = "#be1818", size = 1) +
    geom_vline(xintercept = high_nFeature, linetype = "dashed", color = "#be1818", size = 1) +
    # Annotations for threshold lines
    annotate("text", 
             x = max_nFeature * 0.95, 
             y = mt_threshold + 3, 
             label = paste("y =", mt_threshold), 
             color = "#be1818", size = 4, hjust = 1) +
    annotate("text", 
             x = low_nFeature + low_nFeature, 
             y = min_percent + 50, 
             label = paste("x =", low_nFeature), 
             color = "#be1818", size = 4, vjust = -1) +
    annotate("text", 
             x = high_nFeature - low_nFeature, 
             y = min_percent + 50, 
             label = paste("x =", high_nFeature), 
             color = "#be1818", size = 4, vjust = -1) +
    # Annotations for cell counts in each region
    # Lower row (percent.mt < mt_threshold)
    annotate("text", 
             x = low_nFeature / 2, 
             y = (min_percent + mt_threshold) / 2, 
             label = reg1, 
             color = "black", size = 4) +
    annotate("text", 
             x = (low_nFeature + high_nFeature) / 2, 
             y = (min_percent + mt_threshold) / 2, 
             label = reg2, 
             color = "#ffcc00", size = 4) +
    annotate("text", 
             x = (high_nFeature + max_nFeature) / 2, 
             y = (min_percent + mt_threshold) / 2, 
             label = reg3, 
             color = "black", size = 4) +
    # Upper row (percent.mt >= mt_threshold)
    annotate("text", 
             x = low_nFeature / 2, 
             y = (mt_threshold + max_percent) / 2, 
             label = reg4, 
             color = "black", size = 4) +
    annotate("text", 
             x = (low_nFeature + high_nFeature) / 2, 
             y = (mt_threshold + max_percent) / 2, 
             label = reg5, 
             color = "black", size = 4) +
    annotate("text", 
             x = (high_nFeature + max_nFeature) / 2, 
             y = (mt_threshold + max_percent) / 2, 
             label = reg6, 
             color = "black", size = 4)
  
  # Create p2: Scatter plot of log10(nCount_RNA) vs. log10(nFeature_RNA)
  # Extract metadata and compute log10-transformed values (adding 1 to avoid log10(0))
  md <- seurat_obj@meta.data
  md$log_nCount_RNA <- log10(md$nCount_RNA + 1)
  md$log_nFeature_RNA <- log10(md$nFeature_RNA + 1)
  
  p2 <- ggplot(md, aes(x = log_nCount_RNA, y = log_nFeature_RNA)) +
    geom_point(color = "#025e8d", alpha = 0.7) +
    theme_minimal() +
    labs(x = "log10(nCount_RNA + 1)", y = "log10(nFeature_RNA + 1)",
         title = "Scatter: log10 Transformed") +
    geom_smooth(method = "lm", se = FALSE, color = "#be1818")
  
  # Combine the two plots side-by-side using patchwork and return the patchwork object
  return(p1 + p2)
}




plot_GO_barplot_v2 <- function(
    go_results,
    comparison,
    ontology,
    color_palette,
    n_terms = 5,
    xlim_range = NULL,
    sort.by = "FoldEnrichment",
    scale.by = "-log10(p.adjust)"
) {
  go_df <- as.data.frame(go_results)
  
  if (nrow(go_df) == 0) {
    print(paste("No enriched GO terms for", ontology, "in", comparison))
    return(NULL)
  }
  
  go_df <- go_df %>%
    filter(p.adjust < 0.05)
  
  if (nrow(go_df) == 0) {
    print(paste("No significant enriched GO terms for", ontology, "in", comparison))
    return(NULL)
  }
  
  allowed_sort <- c("FoldEnrichment", "zScore", "p.adjust")
  if (!sort.by %in% allowed_sort) {
    stop("sort.by must be one of: ", paste(allowed_sort, collapse = ", "))
  }
  
  allowed_scale <- c("FoldEnrichment", "zScore", "p.adjust", "-log10(p.adjust)")
  if (!scale.by %in% allowed_scale) {
    stop("scale.by must be one of: ", paste(allowed_scale, collapse = ", "))
  }
  
  if (sort.by %in% c("FoldEnrichment", "zScore") && !sort.by %in% colnames(go_df)) {
    stop(sort.by, " column not found in go_results")
  }
  
  if (scale.by %in% c("FoldEnrichment", "zScore", "p.adjust") && !scale.by %in% colnames(go_df)) {
    stop(scale.by, " column not found in go_results")
  }
  
  go_df <- go_df %>%
    mutate(`-log10(p.adjust)` = -log10(p.adjust))
  
  if (sort.by == "p.adjust") {
    go_top <- go_df %>%
      arrange(p.adjust) %>%
      slice_head(n = n_terms) %>%
      mutate(x_value = -log10(p.adjust))
    
    x_lab <- "-log10(p.adjust)"
    
  } else {
    go_top <- go_df %>%
      arrange(desc(.data[[sort.by]])) %>%
      slice_head(n = n_terms) %>%
      mutate(x_value = .data[[sort.by]])
    
    x_lab <- sort.by
  }
  
  go_top$Description <- factor(go_top$Description, levels = rev(go_top$Description))
  
  color_range <- range(go_top[[scale.by]], na.rm = TRUE)
  colors <- color_palette(10)
  
  p <- ggplot(go_top, aes(x = x_value, y = Description, fill = .data[[scale.by]])) +
    geom_col() +
    scale_fill_gradientn(
      colors = colors,
      name = scale.by,
      limits = color_range
    ) +
    theme_minimal(base_size = 10) +
    xlab(x_lab) +
    ylab("") +
    ggtitle(paste(ontology, "-", comparison))
  
  if (!is.null(xlim_range)) {
    p <- p + xlim(xlim_range)
  }
  
  return(p)
}
