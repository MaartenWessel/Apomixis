# Clean up environment
cat("\014"); graphics.off(); rm(list = ls())

# Load libraries
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

#load data

df <- read_csv(
  "C:/Users/Maarten/Documents/Plantstage/Pollen_staining/Totalstain/total_pollen_results_v2.csv",
  name_repair = "unique"
)

# Define columns
img_col <- "Img"
feret_col <- "Feret"

# Clean and group original data
df_clean <- df %>%
  select(Img = all_of(img_col), Feret = all_of(feret_col)) %>%
  mutate(
    Img = as.character(Img),
    Feret = as.numeric(Feret),
    group = str_extract(Img, "^[^-]+"),
    callusplant = str_extract(Img, "(?<=-)[^_]+")
  )

# Filter for the three plants of interest
plants_of_interest <- tribble(
  ~group, ~callusplant, ~full_name,
  "GIZ39", "C2P1",    "GIZ39-C2P1",
  "GIZ32", "C13P1",   "GIZ32-C13P1",
  "GIZ32", "C16P2",   "GIZ32-C16P2"
)

df_original <- df_clean %>%
  inner_join(plants_of_interest, by = c("group", "callusplant")) %>%
  mutate(source = "automatic")


# 2. Load manual data 

load_manual <- function(file, full_name) {
  dat <- read_csv(file)
  # Select column 3 (C) for Feret/length
  dat <- dat %>%
    mutate(
      Feret = .[[3]],  # column 3 is "length"
      full_name = full_name,
      source = "manual"
    ) %>%
    select(Feret, full_name, source)
  return(dat)
}

manual_dir <- "C:/Users/Maarten/Documents/Plantstage/Pollen_staining"
manual_files <- list(
  "GIZ32-C13P1" = file.path(manual_dir, "Manual_Meas_GIZ32-C13P1.csv"),
  "GIZ32-C16P2" = file.path(manual_dir, "Manual_Meas_GIZ32-C16P2.csv"),
  "GIZ39-C2P1"  = file.path(manual_dir, "Manual_Meas_GIZ39-C2P1.csv")
)

manual_data <- bind_rows(
  load_manual(manual_files[["GIZ32-C13P1"]], "GIZ32-C13P1"),
  load_manual(manual_files[["GIZ32-C16P2"]], "GIZ32-C16P2"),
  load_manual(manual_files[["GIZ39-C2P1"]], "GIZ39-C2P1")
)

# Make sure original data has same columns
df_original_plot <- df_original %>%
  select(Feret, full_name, source)

# Combine manual and automatic data
plot_data <- bind_rows(df_original_plot, manual_data)

# Set factor order for plotting
plot_data$full_name <- factor(
  plot_data$full_name,
  levels = c("GIZ39-C2P1","GIZ32-C13P1","GIZ32-C16P2")
)
plot_data$source <- factor(plot_data$source, levels = c("automatic", "manual"))


# Run Mann-Whitney U tests and prepare annotation for stars

annot_df <- plot_data %>%
  group_by(full_name) %>%
  summarize(
    y_pos = max(Feret[source == "manual"], na.rm = TRUE) + 3, # adjust offset as needed
    pval = {
      auto <- Feret[source == "automatic"]
      man  <- Feret[source == "manual"]
      if (length(auto) > 0 && length(man) > 0) {
        wilcox.test(auto, man)$p.value
      } else {
        NA
      }
    }
  ) %>%
  mutate(
    stars = case_when(
      is.na(pval)     ~ "",
      pval < 0.001    ~ "***",
      pval < 0.01     ~ "**",
      pval < 0.05     ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  filter(stars != "") %>%
  mutate(
    source = "manual"
  )


# 4. Plot violinplots with significance stars

color_map <- c("automatic" = "lightblue", "manual" = "orange")

p <- ggplot(plot_data, aes(x = full_name, y = Feret, fill = source)) +
  geom_violin(width=1, color="black", trim=TRUE, alpha=0.9, position=position_dodge(width=0.9)) + # trim=TRUE!
  scale_fill_manual(values = color_map) +
  coord_cartesian(ylim = c(0, 30)) + # restrict y axis to measurement range
  theme_minimal() +
  labs(
    title = "Pollen diameter validation, Automatic vs Manual",
    x = "",
    y = "Pollen diameter (µm)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add significance stars above manual violins where significant
if (nrow(annot_df) > 0) {
  p <- p + geom_text(
    data = annot_df,
    aes(x = full_name, y = y_pos, label = stars),
    color = "black",
    size = 8,
    fontface = "bold",
    inherit.aes = FALSE
  )
}
print(p)


# 5. Print Mann-Whitney U results and summary stats

cat("\nMann-Whitney U test results (automatic vs manual):\n")
for (plant in levels(plot_data$full_name)) {
  auto <- plot_data$Feret[plot_data$full_name == plant & plot_data$source == "automatic"]
  man  <- plot_data$Feret[plot_data$full_name == plant & plot_data$source == "manual"]
  if(length(auto) > 0 && length(man) > 0){
    test <- wilcox.test(auto, man)
    cat(sprintf("%s: p = %.4g (n_auto = %d, n_manual = %d)\n",
                plant, test$p.value, length(auto), length(man)))
  }
}

cat("\nSummary statistics:\n")
summary_stats <- plot_data %>%
  group_by(full_name, source) %>%
  summarise(
    n = n(),
    mean = mean(Feret, na.rm=TRUE),
    sd = sd(Feret, na.rm=TRUE),
    min = min(Feret, na.rm=TRUE),
    max = max(Feret, na.rm=TRUE),
    .groups = "drop"
  )
print(summary_stats)