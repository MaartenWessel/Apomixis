# Clean up environment
cat("\014"); graphics.off(); rm(list = ls())

# Load libraries
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# Set your folder path
folder <- "C:/Users/Maarten/Documents/Plantstage/Pollen_staining/Totalstain"

# Define file paths
csv_main   <- file.path(folder, "total_pollen_results_28042025.csv")
csv_v2     <- file.path(folder, "total_pollen_results_v2.csv")
csv_orig   <- file.path(folder, "total_pollen_results.csv")

# Read CSVs
df         <- read_csv(csv_main, name_repair = "unique")
df_v2      <- read_csv(csv_v2, name_repair = "unique")
df_orig    <- read_csv(csv_orig, name_repair = "unique")

# Find and add measurements
imgs_df      <- unique(df$Img)
imgs_v2      <- unique(df_v2$Img)
imgs_orig    <- unique(df_orig$Img)

# Find images in v2 and orig that are NOT in the main df
imgs_v2_only   <- setdiff(imgs_v2, imgs_df)
imgs_orig_only <- setdiff(imgs_orig, imgs_df)

# Combine unique new images from v2 and orig 
imgs_unique_new <- union(imgs_v2_only, imgs_orig_only)

# Get the rows for these new images from v2 and orig
df_v2_new   <- df_v2 %>% filter(Img %in% imgs_v2_only)
df_orig_new <- df_orig %>% filter(Img %in% imgs_orig_only)

# Combine new measurements into one data frame
df_new_measurements <- bind_rows(df_v2_new, df_orig_new)

# Combine with the main data frame for downstream analysis
df_combined <- bind_rows(df, df_new_measurements)

# Define columns
img_col <- "Img"
feret_col <- "Feret"

# Prepare data, label Wildtype and MiMe
df_clean <- df_combined %>%
  select(Img = all_of(img_col), Feret = all_of(feret_col)) %>%
  mutate(
    Img = as.character(Img),
    Feret = as.numeric(Feret),
    group = str_extract(Img, "^[^-]+")
  ) %>%
  filter(
    group == "GIZ39" |  
      Img %in% c("WT_MT_merged_RAW_ch00.tif", "MiMe 182 3-1_RAW_ch00.tif", "WT2_RAW_ch00.tif", "WT3_RAW_ch00.tif", "WT_RAW_ch00.tif"),
    is.finite(Feret)
  ) %>%
  mutate(
    callusplant = case_when(
      Img %in% c("WT_MT_merged_RAW_ch00.tif", "WT2_RAW_ch00.tif", "WT3_RAW_ch00.tif", "WT_RAW_ch00.tif") ~ "Wild-type",
      Img == "MiMe 182 3-1_RAW_ch00.tif" ~ "MiMe",
      TRUE ~ str_extract(Img, "(?<=-)[^_]+")
    )
  )

# Get all unique callusplants (controls + samples)
all_callusplants <- unique(df_clean$callusplant)
callusplants_no_ctrl <- setdiff(all_callusplants, c("Wild-type", "MiMe"))

# Get mean pollen diameter for each callusplant
callusplant_means <- df_clean %>%
  filter(!(callusplant %in% c("Wild-type", "MiMe"))) %>%
  group_by(callusplant) %>%
  summarise(mean_feret = mean(Feret, na.rm = TRUE)) %>%
  arrange(desc(mean_feret))

# Create new order
callusplant_levels <- c(
  "Wild-type",
  "MiMe",
  callusplant_means$callusplant
)

# Mann-Whitney U test all: vs Wildtype
wildtype_data <- df_clean$Feret[df_clean$callusplant == "Wild-type"]
pvals <- setNames(rep(NA, length(callusplant_levels)), callusplant_levels)
for (cp in callusplant_levels) {
  if (cp == "Wild-type") next
  cp_data <- df_clean$Feret[df_clean$callusplant == cp]
  if (length(cp_data) > 0 && length(wildtype_data) > 0) {
    pvals[cp] <- wilcox.test(cp_data, wildtype_data)$p.value
  }
}

# Mann-Whitney U test: all vs MiMe
mime_data <- df_clean$Feret[df_clean$callusplant == "MiMe"]
pvals_mime <- setNames(rep(NA, length(callusplant_levels)), callusplant_levels)
for (cp in callusplant_levels) {
  if (cp == "MiMe") next
  cp_data <- df_clean$Feret[df_clean$callusplant == cp]
  if (length(cp_data) > 0 && length(mime_data) > 0) {
    pvals_mime[cp] <- wilcox.test(cp_data, mime_data)$p.value
  }
}

# Prepare plotting data 
plot_data <- df_clean %>%
  filter(callusplant %in% callusplant_levels)
plot_data$callusplant <- factor(plot_data$callusplant, levels = callusplant_levels)

# Define color palette 
color_palette <- setNames(
  rep("lightblue", length(callusplant_levels)),
  callusplant_levels
)
if ("Wild-type" %in% names(color_palette)) color_palette["Wild-type"] <- "tomato"
if ("MiMe" %in% names(color_palette)) color_palette["MiMe"] <- "palegreen"

# 6. A/B/C label logic 
WT_labels <- names(pvals)[pvals > 0.05 & !is.na(pvals) & names(pvals) != "Wild-type"]
MM_labels <- names(pvals_mime)[pvals_mime > 0.05 & !is.na(pvals_mime) & names(pvals_mime) != "MiMe"]

label_df <- tibble(
  callusplant = callusplant_levels,
  WT = callusplant %in% WT_labels,
  MM = callusplant %in% MM_labels
) %>%
  mutate(label = case_when(
    WT & MM ~ "AB",
    WT ~ "A",
    MM ~ "B",
    !WT & !MM & !(callusplant %in% c("Wild-type", "MiMe")) ~ "C",
    TRUE ~ ""
  )) %>%
  filter(label != "" | callusplant %in% c("Wild-type", "MiMe")) %>%
  mutate(
    label = case_when(
      callusplant == "Wild-type" ~ "A",
      callusplant == "MiMe" ~ "B",
      TRUE ~ label
    ),
    ypos = 65
  )

# Plot 
p <- ggplot(plot_data, aes(x = callusplant, y = Feret, fill = callusplant)) +
  geom_violin(color = "black", trim = FALSE) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(limits = c(0, 68), breaks = seq(0, 60, by = 10)) +
  theme_minimal() +
  labs(
    title = "STDM1-P19L - Pollen diameter",
    x = NULL,
    y = "Pollen diameter (µm)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

if (nrow(label_df) > 0) {
  p <- p + geom_text(
    data = label_df,
    aes(x = callusplant, y = ypos, label = label, color = label),
    size = 6,
    vjust = 0,
    show.legend = FALSE
  ) +
    scale_color_manual(
      values = c("A" = "red", "B" = "green", "AB" = "black", "C" = "blue")
    )
}
print(p)

## Print summary stats
cat("Wild-type average pollen diameter:", mean(df_clean$Feret[df_clean$callusplant == "Wild-type"], na.rm=TRUE), "\n")
cat("MiMe average pollen diameter:", mean(df_clean$Feret[df_clean$callusplant == "MiMe"], na.rm=TRUE), "\n")
cat("Highest pollen diameter:", max(df_clean$Feret, na.rm=TRUE), "\n")
cat("Lowest pollen diameter:", min(df_clean$Feret, na.rm=TRUE), "\n")
avg_n_measurements <- df_clean %>%
  group_by(callusplant) %>%
  summarise(n = n(), .groups = "drop") %>%
  summarise(avg = mean(n)) %>%
  pull(avg)
cat("Average number of measurements per plant:", avg_n_measurements, "\n")

# Add number of unique plants measured
n_plants <- length(unique(df_clean$callusplant))
cat("Number of unique plants measured:", n_plants, "\n")