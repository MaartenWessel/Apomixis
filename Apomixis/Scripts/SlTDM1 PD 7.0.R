# Clean up environment
cat("\014"); graphics.off(); rm(list = ls())

# Load libraries
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)

# Set folder path
folder <- "C:/Users/Maarten/Documents/Plantstage/Pollen_staining/Totalstain"

# Define file paths
csv_main   <- file.path(folder, "total_pollen_results_28042025.csv")
csv_v2     <- file.path(folder, "total_pollen_results_v2.csv")
csv_orig   <- file.path(folder, "total_pollen_results.csv")

# Read CSVs
df         <- read_csv(csv_main, name_repair = "unique")
df_v2      <- read_csv(csv_v2, name_repair = "unique")
df_orig    <- read_csv(csv_orig, name_repair = "unique")

# add measurements from v2 and orig
imgs_df      <- unique(df$Img)
imgs_v2      <- unique(df_v2$Img)
imgs_orig    <- unique(df_orig$Img)

imgs_v2_only   <- setdiff(imgs_v2, imgs_df)
imgs_orig_only <- setdiff(imgs_orig, imgs_df)

df_v2_new   <- df_v2 %>% filter(Img %in% imgs_v2_only)
df_orig_new <- df_orig %>% filter(Img %in% imgs_orig_only)

df_new_measurements <- bind_rows(df_v2_new, df_orig_new)
df_combined <- bind_rows(df, df_new_measurements)

# Prepare data and assign callusplant (original names) 
img_col <- "Img"
feret_col <- "Feret"

df_clean <- df_combined %>%
  select(Img = all_of(img_col), Feret = all_of(feret_col)) %>%
  mutate(
    Img = as.character(Img),
    Feret = as.numeric(Feret),
    group = str_extract(Img, "^[^-]+")
  ) %>%
  filter(
    group == "GIZ32" | group == "32" | 
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

# All stats and ordering with original callusplant names
callusplant_means <- df_clean %>%
  filter(!(callusplant %in% c("Wild-type", "MiMe"))) %>%
  group_by(callusplant) %>%
  summarise(mean_feret = mean(Feret, na.rm = TRUE)) %>%
  arrange(desc(mean_feret))

callusplant_levels <- c(
  "Wild-type",
  "MiMe",
  callusplant_means$callusplant
)

# Mann-Whitney U test: all vs Wild-type 
wildtype_data <- df_clean$Feret[df_clean$callusplant == "Wild-type"]
pvals <- setNames(rep(NA, length(callusplant_levels)), callusplant_levels)
for (cp in callusplant_levels) {
  if (cp == "Wild-type") next
  if (!(cp %in% df_clean$callusplant)) next
  cp_data <- df_clean$Feret[df_clean$callusplant == cp]
  if (length(cp_data) > 0 && length(wildtype_data) > 0 &&
      length(unique(cp_data)) > 1 && length(unique(wildtype_data)) > 1) {
    pvals[cp] <- wilcox.test(cp_data, wildtype_data)$p.value
  }
}

# Mann-Whitney U test: all vs MiMe 
mime_data <- df_clean$Feret[df_clean$callusplant == "MiMe"]
pvals_mime <- setNames(rep(NA, length(callusplant_levels)), callusplant_levels)
for (cp in callusplant_levels) {
  if (cp == "MiMe") next
  if (!(cp %in% df_clean$callusplant)) next
  cp_data <- df_clean$Feret[df_clean$callusplant == cp]
  if (length(cp_data) > 0 && length(mime_data) > 0 &&
      length(unique(cp_data)) > 1 && length(unique(mime_data)) > 1) {
    pvals_mime[cp] <- wilcox.test(cp_data, mime_data)$p.value
  }
}

########################################################
# apply callusplant renaming

callusplants_original <- unique(df_clean$callusplant)
callusplants_renamed <- ifelse(
  grepl("^\\d+-\\d+$", callusplants_original),
  paste0("C", sub("^(\\d+)-(\\d+)$", "\\1P\\2", callusplants_original)),
  callusplants_original
)

# Ensure uniqueness
# If duplicate, add -2, 
for (cp in unique(callusplants_renamed)) {
  idx <- which(callusplants_renamed == cp)
  if (length(idx) > 1) {
    callusplants_renamed[idx] <- paste0(cp, "-", seq_along(idx))
    callusplants_renamed[idx[1]] <- cp  # Keep the first as is
  }
}
rename_lookup <- setNames(callusplants_renamed, callusplants_original)
df_clean$callusplant_plot <- rename_lookup[df_clean$callusplant]

# Redefine callusplant_levels for plotting (renamed)
callusplant_levels_renamed <- rename_lookup[callusplant_levels]

# Prepare plotting data 
plot_data <- df_clean %>%
  filter(callusplant %in% callusplant_levels) %>%
  mutate(callusplant_plot = factor(callusplant_plot, levels = callusplant_levels_renamed))

# Define color palette
color_palette <- setNames(
  rep("lightblue", length(callusplant_levels_renamed)),
  callusplant_levels_renamed
)
if ("Wild-type" %in% names(rename_lookup)) color_palette[rename_lookup["Wild-type"]] <- "tomato"
if ("MiMe" %in% names(rename_lookup)) color_palette[rename_lookup["MiMe"]] <- "palegreen"

# A/B/C labels
WT_labels <- names(pvals)[pvals > 0.05 & !is.na(pvals) & names(pvals) != "Wild-type"]
MM_labels <- names(pvals_mime)[pvals_mime > 0.05 & !is.na(pvals_mime) & names(pvals_mime) != "MiMe"]

label_df <- tibble(
  callusplant = callusplant_levels,
  callusplant_plot = callusplant_levels_renamed,
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

# Plot figure
p <- ggplot(plot_data, aes(x = callusplant_plot, y = Feret, fill = callusplant_plot)) +
  geom_violin(color = "black", trim = FALSE) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(limits = c(0, 68), breaks = seq(0, 60, by = 10)) +
  theme_minimal() +
  labs(
    title = "SlTDM1-P19L - Pollen diameter",
    x = NULL,
    y = "Pollen diameter (µm)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

if (nrow(label_df) > 0) {
  p <- p + geom_text(
    data = label_df,
    aes(x = callusplant_plot, y = ypos, label = label, color = label),
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

n_plants <- length(unique(df_clean$callusplant))
cat("Number of unique plants measured:", n_plants, "\n")

target_groups <- c("C4P4", "C18P1", "C24P1", "C31P1")

for (grp in target_groups) {
  avg <- mean(df_clean$Feret[df_clean$callusplant_plot == grp], na.rm = TRUE)
  cat(sprintf("Average pollen diameter for %s: %.2f\n", grp, avg))
}