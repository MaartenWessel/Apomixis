library(readr)
library(dplyr)
library(stringr)
library(openxlsx)

set.seed(42) # For reproducibility

# Load the CSV
df <- read_csv(
  "C:/Users/Maarten/Documents/Plantstage/Pollen_staining/Totalstain/total_pollen_results.csv",
  name_repair = "unique"
)

# Prepare data
df_clean <- df %>%
  select(Img, Feret) %>%
  mutate(
    Img = as.character(Img),
    Feret = as.numeric(Feret),
    group = str_extract(Img, "^[^-]+")
  ) %>%
  filter(group == "GIZ39") %>%
  mutate(callusplant = str_extract(Img, "(?<=-)[^_]+"))

# Shapiro-Wilk per callusplant (subsample if n > 5000)
unique_callusplants <- unique(df_clean$callusplant)

shapiro_results <- data.frame(
  callusplant = character(),
  n = integer(),
  used_n = integer(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (cp in unique_callusplants) {
  cp_data <- df_clean %>% filter(callusplant == cp) %>% pull(Feret)
  n <- length(cp_data)
  used_n <- n
  if (n >= 3 && n <= 5000) {
    sw <- shapiro.test(cp_data)
    p_val <- sw$p.value
  } else if (n > 5000) {
    # Random sample of 5000
    cp_sample <- sample(cp_data, 5000)
    sw <- shapiro.test(cp_sample)
    p_val <- sw$p.value
    used_n <- 5000
  } else {
    p_val <- NA
  }
  shapiro_results <- rbind(
    shapiro_results,
    data.frame(callusplant = cp, n = n, used_n = used_n, p_value = p_val, stringsAsFactors = FALSE)
  )
}

print(shapiro_results)

# Optionally save to CSV
# write.csv(shapiro_results, "shapiro_results.csv", row.names = FALSE)
cat("Callusplants with p-value > 0.05 (normality not rejected):\n")
print(shapiro_results$callusplant[shapiro_results$p_value > 0.05])

# Write results to CSV
output_dir <- "C:/Users/Maarten/Documents/Plantstage/Pollen_staining/ShapiroWilk"
output_file <- file.path(output_dir, "GIZ39.csv")
write.csv(shapiro_results, output_file, row.names = FALSE)
cat(sprintf("Results saved to: %s\n", output_file))
