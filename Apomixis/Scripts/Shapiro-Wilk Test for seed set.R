# Shapiro-Wilk Test for seed set and save results as CSV

# Clean slate
cat("\014"); graphics.off(); rm(list = ls())

# Load data
df <- read.csv("C:/Users/Maarten/Documents/Plantstage/Seed set/Seed setGIZ001.csv", row.names=1)

# Select only M1–M10 columns
seed_data <- df[, grepl("^M[0-9]+$", names(df))]

# Function to safely perform the test (returns NA if not possible)
safe_shapiro <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if(length(unique(x)) < 3) return(NA)  # Shapiro needs at least 3 unique values
  shapiro.test(x)$p.value
}

# Apply the safe function per row (sample)
shapiro_pvals <- apply(seed_data, 1, safe_shapiro)

# Create a results data frame
shapiro_results <- data.frame(
  sample = rownames(seed_data),
  p_value = shapiro_pvals,
  stringsAsFactors = FALSE
)

# Write results to CSV
output_file <- "C:/Users/Maarten/Documents/Plantstage/Seed set/shapiroSS001.csv"
write.csv(shapiro_results, output_file, row.names = FALSE)
cat(sprintf("Results saved to: %s\n", output_file))

# Print all p-values > 0.05
over_05 <- shapiro_results[!is.na(shapiro_results$p_value) & shapiro_results$p_value > 0.05, ]
if(nrow(over_05) > 0) {
  cat("Samples with p-values > 0.05:\n")
  print(over_05)
} else {
  cat("All p-values are under 0.05\n")
}

