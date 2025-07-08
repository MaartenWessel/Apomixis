# Clean slate
rm(list = ls())
graphics.off()
cat("\014") 

# Read the CSV file (skip first line)
file_path <- "C:/Users/Maarten/Documents/Plantstage/Seed set/Seed setGIZ32.csv"
df_raw <- read.csv(file_path, skip = 1, header = FALSE, stringsAsFactors = FALSE)

# Variable names
varnames <- c("Construct_Plant", paste0("M", 1:10), "Average", "Stdev", "Total")
colnames(df_raw) <- varnames

# Convert measurement columns to numeric
for (k in 1:10) {
  col <- paste0("M", k)
  if (!is.numeric(df_raw[[col]])) {
    df_raw[[col]] <- as.numeric(df_raw[[col]])
  }
}

# Exclude rows with less than 10 measurements
num_measurements <- rowSums(!is.na(df_raw[, paste0("M", 1:10)]))
df10 <- df_raw[num_measurements == 10, ]

# Original Figure
plants_all <- as.character(df10$Construct_Plant)
means_all <- rowMeans(df10[, paste0("M", 1:10)], na.rm = TRUE)
stdevs_all <- apply(df10[, paste0("M", 1:10)], 1, sd, na.rm = TRUE)

# Find wildtype row
wildtypeGroupName <- "wildtype"
wildtype_idx <- which(plants_all == wildtypeGroupName)
non_wildtype_idx <- setdiff(seq_along(plants_all), wildtype_idx)

# Sort non-wildtype indices by mean (descending)
sorted_non_wt_idx <- non_wildtype_idx[order(means_all[non_wildtype_idx], decreasing = TRUE)]

# Combine for plotting order: wildtype first, then sorted by mean
plot_idx <- c(wildtype_idx, sorted_non_wt_idx)

# Build sorted vectors for plotting
plants_all_sorted <- plants_all[plot_idx]
means_all_sorted <- means_all[plot_idx]
stdevs_all_sorted <- stdevs_all[plot_idx]

# Axis labels
x_labels_all <- plants_all_sorted
x_labels_all[x_labels_all == "wildtype"] <- "Wild-type"

# Define bar colors
bar_cols <- rep(rgb(0.6784, 0.8510, 0.9020), length(x_labels_all))
bar_cols[x_labels_all == "Wild-type"] <- "tomato"

# Plot
bar_centers <- barplot(
  means_all_sorted,
  col = bar_cols,
  names.arg = x_labels_all,
  las = 2, # x labels vertical
  ylim = c(0, 100),
  ylab = "Seeds/fruit",
  main = "", # No main, add title manually to control font
)

title("SlTDM1-P19L - Seed set", font.main = 1) # font.main = 1 for plain (not bold) title

arrows(
  x0 = bar_centers,
  y0 = means_all_sorted - stdevs_all_sorted,
  x1 = bar_centers,
  y1 = means_all_sorted + stdevs_all_sorted,
  code = 3, angle = 90, length = 0.05, col = "black", lwd = 1.2
)
box()

# Mann-Whitney U test vs wildtype for original data
uniqueGroupsForPlot_all <- plants_all_sorted
wildtypeIdx_all <- which(plants_all_sorted == wildtypeGroupName)
wildtypeData_all <- as.numeric(unlist(df10[plants_all == wildtypeGroupName, paste0("M", 1:10)]))
wildtypeData_all <- wildtypeData_all[!is.na(wildtypeData_all)]
pvals_all <- rep(NA, length(uniqueGroupsForPlot_all))

for (i in seq_along(uniqueGroupsForPlot_all)) {
  if (uniqueGroupsForPlot_all[i] == wildtypeGroupName) next
  thisGroupData <- as.numeric(unlist(df10[plants_all == uniqueGroupsForPlot_all[i], paste0("M", 1:10)]))
  thisGroupData <- thisGroupData[!is.na(thisGroupData)]
  if (length(thisGroupData) == 0 | length(wildtypeData_all) == 0) {
    pvals_all[i] <- NA
  } else {
    pvals_all[i] <- wilcox.test(thisGroupData, wildtypeData_all)$p.value
  }
}

sigStars_all <- rep("", length(uniqueGroupsForPlot_all))
for (i in seq_along(sigStars_all)) {
  if (uniqueGroupsForPlot_all[i] != wildtypeGroupName && !is.na(pvals_all[i])) {
    if (pvals_all[i] < 0.001) {
      sigStars_all[i] <- "***"
    } else if (pvals_all[i] < 0.01) {
      sigStars_all[i] <- "**"
    } else if (pvals_all[i] < 0.05) {
      sigStars_all[i] <- "*"
    }
  }
}

# Add significance stars directly above whiskers
for (i in seq_along(means_all_sorted)) {
  if (sigStars_all[i] != "" && uniqueGroupsForPlot_all[i] != wildtypeGroupName) {
    text(bar_centers[i], means_all_sorted[i] + stdevs_all_sorted[i], sigStars_all[i],
         cex = 1.3, col = "black", xpd = TRUE, pos = 3, font = 2)
  }
}

cat("\nMann-Whitney U test p-values (vs wildtype) - All groups:\n")
for (i in seq_along(uniqueGroupsForPlot_all)) {
  if (uniqueGroupsForPlot_all[i] == wildtypeGroupName) next
  cat(sprintf("%s vs wildtype: p = %.4g\n", uniqueGroupsForPlot_all[i], pvals_all[i]))
}

# Exemplary Figure
exemplaryGroups <- c("wildtype", "C6P1", "C3P3", "C12P1", "C3P1", "C15P1", "C16P1", "C7P3", "C17P1", "C3P4", "C2P2", "C2P3", "C2P6")

plantsAll <- as.character(df10$Construct_Plant)
idxOrder <- sapply(exemplaryGroups, function(g) {
  idx <- match(g, plantsAll)
  if (is.na(idx)) {
    warning(sprintf("Group %s not found in data!", g))
    return(NA)
  } else {
    idx
  }
})
valid <- !is.na(idxOrder)
dfEx <- df10[idxOrder[valid], ]
exemplaryGroups <- exemplaryGroups[valid]
plants <- as.character(dfEx$Construct_Plant)
means <- rowMeans(dfEx[, paste0("M", 1:10)], na.rm = TRUE)
stdevs <- apply(dfEx[, paste0("M", 1:10)], 1, sd, na.rm = TRUE)

# Axis labels
x_labels_ex <- plants
x_labels_ex[x_labels_ex == "wildtype"] <- "Wild-type"

# Bar colors
bar_cols_ex <- rep(rgb(0.6784, 0.8510, 0.9020), length(x_labels_ex))
bar_cols_ex[x_labels_ex == "Wild-type"] <- "tomato"

bar_centers_ex <- barplot(
  means,
  col = bar_cols_ex,
  names.arg = x_labels_ex,
  las = 2,
  ylim = c(0, 100),
  ylab = "Seeds/fruit",
  main = "" # No main
)
title("SlTDM1 - Seed set - overlap", font.main = 1) 
arrows(
  x0 = bar_centers_ex,
  y0 = means - stdevs,
  x1 = bar_centers_ex,
  y1 = means + stdevs,
  code = 3, angle = 90, length = 0.05, col = "black", lwd = 1.2
)
box()

uniqueGroupsForPlot <- plants
wildtypeIdx <- which(plants == wildtypeGroupName)
wildtypeData <- as.numeric(unlist(dfEx[wildtypeIdx, paste0("M", 1:10)]))
wildtypeData <- wildtypeData[!is.na(wildtypeData)]

pvals <- rep(NA, length(uniqueGroupsForPlot))
for (i in seq_along(uniqueGroupsForPlot)) {
  if (uniqueGroupsForPlot[i] == wildtypeGroupName) next
  thisGroupIdx <- which(plants == uniqueGroupsForPlot[i])
  thisGroupData <- as.numeric(unlist(dfEx[thisGroupIdx, paste0("M", 1:10)]))
  thisGroupData <- thisGroupData[!is.na(thisGroupData)]
  if (length(thisGroupData) == 0 | length(wildtypeData) == 0) {
    pvals[i] <- NA
  } else {
    pvals[i] <- wilcox.test(thisGroupData, wildtypeData)$p.value
  }
}

sigStars <- rep("", length(uniqueGroupsForPlot))
for (i in seq_along(sigStars)) {
  if (uniqueGroupsForPlot[i] != wildtypeGroupName && !is.na(pvals[i])) {
    if (pvals[i] < 0.001) {
      sigStars[i] <- "***"
    } else if (pvals[i] < 0.01) {
      sigStars[i] <- "**"
    } else if (pvals[i] < 0.05) {
      sigStars[i] <- "*"
    }
  }
}

# Add significance stars directly above whiskers
for (i in seq_along(means)) {
  if (sigStars[i] != "" && uniqueGroupsForPlot[i] != wildtypeGroupName) {
    text(bar_centers_ex[i], means[i] + stdevs[i], sigStars[i],
         cex = 1.3, col = "black", xpd = TRUE, pos = 3, font = 2)
  }
}

cat("\nMann-Whitney U test p-values (vs wildtype) - Exemplary groups:\n")
for (i in seq_along(uniqueGroupsForPlot)) {
  if (uniqueGroupsForPlot[i] == wildtypeGroupName) next
  cat(sprintf("%s vs wildtype: p = %.4g\n", uniqueGroupsForPlot[i], pvals[i]))
}