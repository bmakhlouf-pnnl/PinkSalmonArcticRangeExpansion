# Load required libraries
library(tidyverse)
library(zoo)

# ============================================================================
# CONFIG
# ============================================================================
DATA_PATH <- "/Users/benjaminmakhlouf/Research_repos/PNNL/PinkSalmonArcticRangeExpansion/Data/FromMegan/OtoDataNotCleaned.csv"
OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/PNNL/PinkSalmonArcticRangeExpansion/Figures/Trimmed"
CSV_OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/PNNL/PinkSalmonArcticRangeExpansion/Data/Cleaned"
MA_WINDOW <- 7
VALUE_MIN <- 0.7040
VALUE_MAX <- 0.7110

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

oto_with_ma <- oto_filtered %>%
  group_by(Individual) %>%
  mutate(MA = rollmean(Value, k = MA_WINDOW, fill = NA)) %>%
  ungroup()

# ============================================================================
# INTERACTIVE PLOTTING AND POINT SELECTION
# ============================================================================
individuals_list <- unique(oto_with_ma$Individual)
cleaned_data_list <- list()

for (i in seq_along(individuals_list)) {
  
  individual <- individuals_list[i]
  individual_data <- oto_with_ma %>% filter(Individual == individual)
  
  # Print instructions
  cat("\n========== ", i, "of", length(individuals_list), "==========\n")
  cat("Individual:", individual, "\n")
  cat("Measurement range:", min(individual_data$Measurement.number), 
      "to", max(individual_data$Measurement.number), "\n")
  cat("Click on the plot to select the trim point (keep data AFTER this point)\n\n")
  
  # Create base R plot for interactive selection
  dev.new(width = 10, height = 6)
  with(individual_data, {
    plot(Measurement.number, Value, 
         main = paste("Otolith Data -", individual),
         xlab = "Measurement Number", 
         ylab = "Value",
         pch = 16, col = "steelblue", cex = 1.5)
    lines(Measurement.number, MA, col = "darkred", lwd = 2)
  })
  
  # Get user click location
  trim_loc <- locator(1)
  dev.off()
  
  # Trim data: keep measurements >= selected measurement number
  trim_measurement <- round(trim_loc$x)
  trimmed_data <- individual_data %>%
    filter(Measurement.number >= trim_measurement) %>%
    select(Measurement.number, Individual, Value, Group) %>%
    mutate(MA = rollmean(Value, k = MA_WINDOW, fill = NA))
  
  # Store in list (without MA column)
  cleaned_data_list[[individual]] <- trimmed_data %>% select(-MA)
  
  # Create ggplot of trimmed data with moving average
  p_trimmed <- trimmed_data %>%
    ggplot(aes(x = Measurement.number, y = Value)) +
    geom_point(size = 2.5, alpha = 0.5, color = "steelblue") +
    geom_line(aes(y = MA), size = 1, color = "darkred", alpha = 0.8) +
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
  
  # Save plot
  ggsave(
    filename = file.path(OUTPUT_DIR, paste0(individual, "_cleaned.png")),
    plot = p_trimmed,
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
cat("Total individuals processed:", length(individuals_list), "\n")
cat("Total rows in cleaned data:", nrow(cleaned_data), "\n")
cat("Cleaned files saved to:", CSV_OUTPUT_DIR, "\n")
cat("Plots saved to:", OUTPUT_DIR, "\n")