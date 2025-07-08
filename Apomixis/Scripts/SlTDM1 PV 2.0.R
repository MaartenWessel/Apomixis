# Clean slate
rm(list = ls())
graphics.off()

# File paths 
dataFile <- "C:/Users/Maarten/Documents/Plantstage/Amphasys/Datasheet32/CrossGatingFreq2Summary_with_second_sheet.xlsx"
wildtypeFile <- "C:/Users/Maarten/Documents/Plantstage/Amphasys/Wildtype.xlsx"
wildtypeGroupName <- "Wildtype"   # define wildtype name 
wildtypeDisplayName <- "Wild-type" # Displayed on x-axis

# Load data 
library(readxl)
library(ggplot2)

data <- read_excel(dataFile, sheet = 2)
wildtypeData <- read_excel(wildtypeFile)

# Preprocess group column
data$Group <- paste0(data$Construct, data$Plant)
meas_cols <- grep("^Meas_", names(data))
measurementColumns <- data[, meas_cols]

# Prepare data containers 
boxplotData <- numeric()
groupLabels <- character()
averagePercentages <- numeric()
measurementCounts <- numeric()

# Process main groups
uniqueGroups <- unique(data$Group)
for (currentGroup in uniqueGroups) {
  groupRows <- which(data$Group == currentGroup)
  groupMeasurements <- as.numeric(unlist(measurementColumns[groupRows, ]))
  reshapedMeasurements <- groupMeasurements[!is.na(groupMeasurements)]
  
  if (length(reshapedMeasurements) >= 3) {
    boxplotData <- c(boxplotData, reshapedMeasurements)
    groupLabels <- c(groupLabels, rep(currentGroup, length(reshapedMeasurements)))
    averagePercentages <- c(averagePercentages, mean(reshapedMeasurements))
    measurementCounts <- c(measurementCounts, length(reshapedMeasurements))
  }
}

# Add Wildtype data
wildtypeMeasurements <- as.numeric(unlist(wildtypeData[, 2]))
if (length(wildtypeMeasurements) >= 3) {
  boxplotData <- c(boxplotData, wildtypeMeasurements)
  groupLabels <- c(groupLabels, rep(wildtypeGroupName, length(wildtypeMeasurements)))
  averagePercentages <- c(averagePercentages, mean(wildtypeMeasurements))
  measurementCounts <- c(measurementCounts, length(wildtypeMeasurements))
}

# Build display labels 
uniqueGroupsForPlot <- unique(groupLabels)
countsForPlot <- sapply(uniqueGroupsForPlot, function(x) sum(groupLabels == x))
uniqueGroupsForPlot_display <- uniqueGroupsForPlot
uniqueGroupsForPlot_display[uniqueGroupsForPlot_display == wildtypeGroupName] <- wildtypeDisplayName

# Remove count from display label for x axis
displayLabels <- uniqueGroupsForPlot_display

# Mann-Whitney U vs Wildtype 
pvals <- rep(NA_real_, length(uniqueGroupsForPlot))
wildtypeDataForTest <- boxplotData[groupLabels == wildtypeGroupName]

for (i in seq_along(uniqueGroupsForPlot)) {
  if (uniqueGroupsForPlot[i] == wildtypeGroupName) next
  thisGroupData <- boxplotData[groupLabels == uniqueGroupsForPlot[i]]
  if (length(thisGroupData) > 0) {
    pvals[i] <- wilcox.test(thisGroupData, wildtypeDataForTest, exact = FALSE)$p.value
  }
}

# Assign significance stars
sigStars <- rep("", length(pvals))
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

# Sorting 
wildtypeIdx <- which(uniqueGroupsForPlot == wildtypeGroupName)
nonWildtypeIdx <- setdiff(seq_along(uniqueGroupsForPlot), wildtypeIdx)
nonWTsort <- order(averagePercentages[nonWildtypeIdx], decreasing = TRUE)
finalOrder <- c(wildtypeIdx, nonWildtypeIdx[nonWTsort])

# Sort everything consistently
uniqueGroupsForPlot <- uniqueGroupsForPlot[finalOrder]
displayLabels <- displayLabels[finalOrder]
sigStars <- sigStars[finalOrder]
averagePercentages <- averagePercentages[finalOrder]

# Reorder data points
sortedBoxplotData <- numeric()
sortedGroupLabels <- character()
for (i in seq_along(uniqueGroupsForPlot)) {
  mask <- groupLabels == uniqueGroupsForPlot[i]
  sortedBoxplotData <- c(sortedBoxplotData, boxplotData[mask])
  sortedGroupLabels <- c(sortedGroupLabels, rep(displayLabels[i], sum(mask)))
}

# General trend figure

all_groupLabels <- sortedGroupLabels
all_boxplotData <- sortedBoxplotData
all_displayLabels <- displayLabels

all_groupLabels_factor <- factor(all_groupLabels, levels = displayLabels)

# Colors: tomato for wild-type, lightblue for others
all_fill_vals <- setNames(
  rep("lightblue", length(displayLabels)),
  displayLabels
)
wildtype_lab_all <- displayLabels[displayLabels == wildtypeDisplayName]
all_fill_vals[wildtype_lab_all] <- "tomato"

# Data frame for plot
all_plot_df <- data.frame(
  value = all_boxplotData,
  group = all_groupLabels_factor
)

# Compute y-position for significance stars: top of box + 3, max 100
all_max_per_group <- aggregate(value ~ group, data = all_plot_df, FUN = max, na.rm = TRUE)
all_max_per_group$star <- sigStars
all_max_per_group$labelName <- as.character(all_max_per_group$group)
all_max_per_group$ypos <- pmin(all_max_per_group$value + 3, 100)
all_star_df <- all_max_per_group[all_max_per_group$star != "" & all_max_per_group$labelName != wildtypeDisplayName, ]

# Plot with x labels, with significance stars
p2 <- ggplot(all_plot_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(color = "black", outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_fill_manual(values = all_fill_vals) +
  guides(fill = "none") +
  theme_minimal() +
  labs(
    title = "SlTDM1-P19L - Pollen viability",
    x = NULL,
    y = "Pollen viability%"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 8)),
    plot.title = element_text(size = 12)
  )

if (nrow(all_star_df) > 0) {
  p2 <- p2 + geom_text(
    data = all_star_df,
    aes(x = group, y = ypos, label = star),
    color = "black",
    size = 6,
    vjust = 0
  )
}

print(p2)
ggsave("AtTDM1-P17L_viability_general_trend.png", p2, width = 14, height = 6, dpi = 300)