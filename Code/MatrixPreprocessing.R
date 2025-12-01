# Load required libraries
library(tidyverse)
library(zoo)
library(mgcv)

# ============================================================================
# CONFIG
# ============================================================================
DATA_PATH <- "/Users/benjaminmakhlouf/Research_repos/PNNL/PinkSalmonArcticRangeExpansion/Data/FromMegan/OtoDataNotCleaned.csv"
OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/PNNL/PinkSalmonArcticRangeExpansion/Figures/Trimmed"
CSV_OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/PNNL/PinkSalmonArcticRangeExpansion/Data/Cleaned"
MA_WINDOW <- 60
VALUE_MIN <- 0.7040
VALUE_MAX <- 0.7110
TRIM_MODE <- 2  # Set to 1 for interactive trimming, 2 for no trimming
GAMMA_VALUE <- 1.4  # Smoothing parameter for GAM

# Create output directories if they don't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(CSV_OUTPUT_DIR, showWarnings = FALSE)

# ============================================================================
# LOAD AND RESHAPE DATA
# ============================================================================
alldat <- read.csv(DATA_PATH)

oto_long <- alldat %>%
  select(Measurement.number, HART.14:Elson.26.6) %>%
  pivot_longer(
    cols = -Measurement.number,
    names_to = "Individual",
    values_to = "Value"
  ) %>%
  mutate(
    Value = as.numeric(Value),
    Group = str_extract(Individual, "^[A-Za-z]+")
  ) %>%
  drop_na(Value) %>%
  arrange(Individual, Measurement.number)

# ============================================================================
# FILTER AND CALCULATE MOVING AVERAGE
# ============================================================================
oto_filtered <- oto_long %>%
  filter(Value >= VALUE_MIN & Value <= VALUE_MAX)

oto_with_gam <- oto_filtered %>%
  group_by(Individual) %>%
  mutate(Index = row_number()) %>%
  ungroup()

# ============================================================================
# INTERACTIVE PLOTTING AND POINT SELECTION
# ============================================================================
individuals_list <- unique(oto_with_gam$Individual)
cleaned_data_list <- list()
gam_data_list <- list()  # Store GAM data for later interpolation

for (i in seq_along(individuals_list)) {
  
  individual <- individuals_list[i]
  individual_data <- oto_with_gam %>% filter(Individual == individual)
  
  # Print instructions
  cat("\n========== ", i, "of", length(individuals_list), "==========\n")
  cat("Individual:", individual, "\n")
  cat("Measurement range:", min(individual_data$Measurement.number), 
      "to", max(individual_data$Measurement.number), "\n")
  
  if (TRIM_MODE == 1) {
    cat("Click on the plot to select the trim point (keep data AFTER this point)\n\n")
    
    # Create base R plot for interactive selection
    dev.new(width = 10, height = 6)
    with(individual_data, {
      plot(Measurement.number, Value, 
           main = paste("Otolith Data -", individual),
           xlab = "Measurement Number", 
           ylab = "Value",
           pch = 16, col = "steelblue", cex = 1.5)
      abline(h = 0.7091, col = "gray50", lty = 2, lwd = 1.5)
    })
    
    # Get user click location
    trim_loc <- locator(1)
    dev.off()
    
    # Trim data: keep measurements >= selected measurement number
    trim_measurement <- round(trim_loc$x)
    trimmed_data <- individual_data %>%
      filter(Measurement.number >= trim_measurement) %>%
      select(Measurement.number, Individual, Value, Group, Index)
    
    # Fit GAM
    k <- min(30, floor(15 * (nrow(trimmed_data)^(2/9))))
    k <- max(k, 3)
    gam_model <- gam(Value ~ s(Index, bs = "tp", k = k), gamma = GAMMA_VALUE, data = trimmed_data)
    gam_pred <- predict(gam_model, se.fit = TRUE)
    trimmed_data$GAM <- gam_pred$fit
    trimmed_data$GAM_SE <- gam_pred$se.fit
    trimmed_data$GAM_upper <- trimmed_data$GAM + 1.96 * trimmed_data$GAM_SE
    trimmed_data$GAM_lower <- trimmed_data$GAM - 1.96 * trimmed_data$GAM_SE
    
  } else if (TRIM_MODE == 2) {
    # No trimming: keep all filtered data
    trimmed_data <- individual_data %>%
      select(Measurement.number, Individual, Value, Group, Index)
    
    # Fit GAM
    k <- min(30, floor(15 * (nrow(trimmed_data)^(2/9))))
    k <- max(k, 3)
    gam_model <- gam(Value ~ s(Index, bs = "tp", k = k), gamma = GAMMA_VALUE, data = trimmed_data)
    gam_pred <- predict(gam_model, se.fit = TRUE)
    trimmed_data$GAM <- gam_pred$fit
    trimmed_data$GAM_SE <- gam_pred$se.fit
    trimmed_data$GAM_upper <- trimmed_data$GAM + 1.96 * trimmed_data$GAM_SE
    trimmed_data$GAM_lower <- trimmed_data$GAM - 1.96 * trimmed_data$GAM_SE
    
  } else {
    stop("Invalid TRIM_MODE. Please set TRIM_MODE to 1 or 2 in CONFIG section.")
  }
  
  # Create ggplot of data with GAM and confidence intervals
  p_cleaned <- trimmed_data %>%
    ggplot(aes(x = Measurement.number, y = Value)) +
    geom_ribbon(aes(ymin = GAM_lower, ymax = GAM_upper), fill = "darkred", alpha = 0.2) +
    geom_point(size = 2.5, alpha = 0.5, color = "steelblue") +
    geom_line(aes(y = GAM), size = 1, color = "darkred", alpha = 0.8) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste("Otolith Data -", individual),
      x = "Measurement Number",
      y = "Value"
    )
  
  # Store cleaned data (without GAM-related columns)
  cleaned_data_list[[individual]] <- trimmed_data %>% select(-GAM, -GAM_SE, -GAM_upper, -GAM_lower)
  
  # Store GAM data for interpolation
  gam_data_list[[individual]] <- trimmed_data %>% select(Index, GAM)
  
  cat("✓ GAM fitted with k =", k, "\n")
  cat("  GAM range:", round(min(trimmed_data$GAM, na.rm = TRUE), 4), "to", round(max(trimmed_data$GAM, na.rm = TRUE), 4), "\n")
  
  # Save plot to appropriate folder based on trim mode
  if (TRIM_MODE == 1) {
    plot_folder <- OUTPUT_DIR
  } else {
    plot_folder <- gsub("Trimmed$", "Untrimmed", OUTPUT_DIR)
  }
  
  dir.create(plot_folder, showWarnings = FALSE, recursive = TRUE)
  
  ggsave(
    filename = file.path(plot_folder, paste0(individual, "_cleaned.png")),
    plot = p_cleaned,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  cat("✓ Saved plot:", individual, "\n")
}

# ============================================================================
# COMBINE AND EXPORT
# ============================================================================
cleaned_data <- bind_rows(cleaned_data_list)

# Save combined cleaned data
write.csv(
  cleaned_data,
  file = file.path(CSV_OUTPUT_DIR, "OtoDataCleaned_Combined.csv"),
  row.names = FALSE
)

# Save individual files
for (individual in individuals_list) {
  if (individual %in% names(cleaned_data_list)) {
    write.csv(
      cleaned_data_list[[individual]],
      file = file.path(CSV_OUTPUT_DIR, paste0(individual, "_cleaned.csv")),
      row.names = FALSE
    )
  }
}

cat("\n========== SUMMARY ==========\n")
cat("Trimming mode:", if (TRIM_MODE == 1) "Interactive" else "No trimming", "\n")
cat("Total individuals processed:", length(individuals_list), "\n")
cat("Total rows in cleaned data:", nrow(cleaned_data), "\n")
cat("Cleaned files saved to:", CSV_OUTPUT_DIR, "\n")
cat("Plots saved to:", if (TRIM_MODE == 1) OUTPUT_DIR else gsub("Trimmed$", "Untrimmed", OUTPUT_DIR), "\n")

# ============================================================================
# CREATE INTERPOLATED MATRIX FROM GAM DATA
# ============================================================================
# Calculate average length of all GAM timeseries
gam_timeseries_lengths <- sapply(gam_data_list, nrow)
avg_length <- round(mean(gam_timeseries_lengths))

cat("\n========== INTERPOLATION (FROM GAM DATA) ==========\n")
cat("Individual GAM timeseries lengths:\n")
for (i in seq_along(gam_timeseries_lengths)) {
  cat("  ", names(gam_timeseries_lengths)[i], ": ", gam_timeseries_lengths[i], " observations\n", sep = "")
}
cat("Average length:", avg_length, "\n")

# Create interpolated matrix from GAM data
interpolated_matrix <- matrix(NA, nrow = length(individuals_list), ncol = avg_length)
rownames(interpolated_matrix) <- individuals_list

for (i in seq_along(individuals_list)) {
  individual <- individuals_list[i]
  gam_data <- gam_data_list[[individual]]
  
  # Get the GAM column
  original_gam_values <- gam_data$GAM
  original_length <- length(original_gam_values)
  
  # Create new sequence for interpolation
  original_seq <- seq(1, original_length, length.out = original_length)
  new_seq <- seq(1, original_length, length.out = avg_length)
  
  # Interpolate GAM values using linear interpolation
  interpolated_gam_values <- approx(original_seq, original_gam_values, xout = new_seq)$y
  
  # Store in matrix
  interpolated_matrix[i, ] <- interpolated_gam_values
}

# Save interpolated GAM matrix
write.csv(
  interpolated_matrix,
  file = file.path(CSV_OUTPUT_DIR, "OtoDataInterpolated_GAM_Matrix.csv"),
  row.names = TRUE
)

cat("✓ Interpolated GAM matrix saved with dimensions:", nrow(interpolated_matrix), "x", ncol(interpolated_matrix), "\n")
cat("  File: OtoDataInterpolated_GAM_Matrix.csv\n")