# Shapiro-Wilk Test for pollen viability (R version)
# Clean slate
cat("\014"); graphics.off(); rm(list = ls())

# Load required libraries
library(readr)
library(dplyr)
library(tidyr)

# Load data
df <- read.csv("C:/Users/Maarten/Documents/Plantstage/Amphasys/Datasheet001/CrossGatingFreq2Summary_with_second_sheet.csv", 
               na.strings = c("", "NA"))

# Combine Construct and Plant into one identifier
df <- df %>%
  mutate(CP = paste0(Callus, Plant))

# Gather measurement columns into long format
df_long <- df %>%
  pivot_longer(
    cols = starts_with("Meas_"),
    names_to = "Measurement",
    values_to = "VPercent"
  ) %>%
  filter(!is.na(VPercent))

# Perform Shapiro-Wilk test for each CP combination
shapiro_results <- df_long %>%
  group_by(CP) %>%
  summarise(
    n = n(),
    W = if(n >= 3) shapiro.test(VPercent)$statistic else NA,
    p.value = if(n >= 3) shapiro.test(VPercent)$p.value else NA,
    .groups = "drop"
  )

# Save results to CSV
write.csv(shapiro_results, 
          file = "C:/Users/Maarten/Documents/Plantstage/Amphasys/shapiro/shapiroGIZ001_shapiro_results.csv", 
          row.names = FALSE)

# Optional: Print the results
print(n = 99, shapiro_results)

