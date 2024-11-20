library(data.table)
library(readxl)
library(flowCore)
library(Biobase)
library(flowStats)
library(ggplot2)
library(umap) # For visualization
library(uwot)
library(Spectre)
library(sva)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
library(dplyr)
library(ggsignif)

############################################ Global Variables #######################################################
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/FLORAH")
in.path <-"C:/Users/axi313/Documents/NK_Manuscript/R_dat/"

all.flow <- readRDS(paste0(in.path,"FLORAH_batch_corrected_filtered_FlowSOM.rds"))



######
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/FLORAH/Results/")
dir.create("Comparison_Plots")
setwd("Comparison_Plots")

### Check cells
sample.col <- 'Sample_Name'
treatment.col <- "Condition"
timepoint.col <- "Timepoint"
HIV.col <- "HIV_Status"
batch.col <- "Batch"
all.flow$Condition <- gsub("HUT\\+IL-15", "HUT78+IL-15", all.flow$Condition)
levels(as.factor(all.flow$Condition))
data.frame(table(all.flow[[treatment.col]])) # Check number of cells per sample.
cluster.cols <- names(all.flow)[c(44:48,50:59,61,63:67)]

### Dimensionality reduction
sub.targets <- c(50000,50000,50000,50000,50000,50000,50000, 50000,
                 50000,50000,50000,50000,50000,50000,50000, 50000)
flow.sub <- do.subsample(all.flow, sub.targets, treatment.col)
flow.sub <- run.umap(flow.sub, cluster.cols)


make.multi.plot(flow.sub, 'UMAP_X', 'UMAP_Y', 'Ordered_Clusters',col.type = 'factor'
                ,divide.by = 'Condition', figure.title= "Split By Treatment")

make.multi.plot(flow.sub, 'UMAP_X', 'UMAP_Y', 'Ordered_Clusters',col.type = 'factor'
                ,divide.by = 'HIV_Status', figure.title= "Split By HIV Status")


### Summary Stats
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/FLORAH/Results")

dir.create("Output 4 - summary data")
setwd("Output 4 - summary data")

### Setup

variance.test <- 'kruskal.test'
pairwise.test <- "wilcox.test"

comparisons.1 <- list(c("HEI", "HEU"))
comparisons.2 <- list(c("CEM", "CEM+IL-15", "HUT78","HUT78+IL-15", "K562", "K562+IL-15","Untreated","IL-15"))

grp.order.1 <- c("HEI", "HEU")
grp.order.2 <- c("CEM", "CEM+IL-15", "HUT78","HUT78+IL-15", "K562", "K562+IL-15","Untreated","IL-15")

### Select columns to measure MFI

sum.dat <- create.sumtable(dat = all.flow, 
                           sample.col = sample.col,
                           pop.col = "Ordered_Clusters",
                           use.cols = cluster.cols, 
                           annot.cols = c(treatment.col, HIV.col,timepoint.col)
)
as.matrix(names(sum.dat))
annot.cols <- c(treatment.col, HIV.col,timepoint.col)

plot.cols <- names(sum.dat)[c(5:124)]
plot.cols

######### Correlation Analysis ##############

hei_data <- subset(all.flow, HIV_Status == "HEI" & Condition == 'K562')
markers <- cluster.cols
clinical_vars <-  "Specific_Killing"  # Replace with your clinical data columns

cor_results <- sapply(markers, function(marker) {
  sapply(clinical_vars, function(clinical_var) {
    cor.test(hei_data[[marker]], hei_data[[clinical_var]], method = "spearman")$estimate
  })
})

# Convert the correlation results to a matrix
cor_matrix <- matrix(unlist(cor_results), nrow = length(markers), byrow = TRUE)
rownames(cor_matrix) <- markers
colnames(cor_matrix) <- clinical_vars

# Step 6: Adjust p-values for multiple comparisons (Benjamini-Hochberg correction)
p_values <- sapply(markers, function(marker) {
  sapply(clinical_vars, function(clinical_var) {
    cor.test(hei_data[[marker]], hei_data[[clinical_var]], method = "spearman")$p.value
  })
})
adjusted_p_values <- p.adjust(p_values, method = "BH")
adjusted_p_matrix <- matrix(unlist(adjusted_p_values), nrow = length(markers), byrow = TRUE)
rownames(adjusted_p_matrix) <- markers
colnames(adjusted_p_matrix) <- clinical_vars

# Step 8: Create a correlation heatmap with significance annotations

# Round the correlation matrix and p-values for display
cor_matrix_rounded <- round(cor_matrix, 2)
adjusted_p_matrix_rounded <- round(adjusted_p_matrix, 3)

# Create custom annotations for significance: display p-values for each cell
annotations <- matrix("", nrow = nrow(adjusted_p_matrix), ncol = ncol(adjusted_p_matrix))
annotations[adjusted_p_matrix < 0.05] <- "*"
annotations[adjusted_p_matrix < 0.01] <- "**"
annotations[adjusted_p_matrix < 0.001] <- "***"

# Use pheatmap to plot the heatmap with significance annotations
# Save heatmap as PNG
png("Marker_Clinical_Variable_Correlation_Heatmap.png", width = 1200, height = 800)

# Generate the heatmap
pheatmap(cor_matrix_rounded, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         display_numbers = annotations,  # Use significance annotations
         main = "Marker-Clinical Variable Correlation Heatmap")

# Close the PNG device
dev.off()

# Step 7: Visualize results


# Scatter plot for specific killing vs Marker1
ggplot(hei_data, aes(x = `FASL (CD178)`, y = Specific_Killing)) +
  geom_point(color = "blue", alpha = 0.6) +  # Points on the plot
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear regression line
  theme_minimal() +  # Use a clean theme
  labs(title = "Scatter Plot: Marker1 vs Specific Killing", 
       x = "Marker1 Expression", 
       y = "Specific Killing")
ggsave("FLORAH_SK_Pheatmap.png", width = 7, height = 9, bg='White')

# Step 8: Create a correlation heatmap
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, 
         display_numbers = TRUE, main = "Marker-Clinical Variable Correlation Heatmap")


# Step 10: Log-transform viral load for skewed data (optional)
datatable$log_viral_load <- log(datatable$viral_load + 1)

# You can also repeat the correlation analysis using the log-transformed viral load
cor_results_log <- sapply(markers, function(marker) {
  cor.test(hei_data[[marker]], hei_data$log_viral_load, method = "spearman")$estimate
})

# Optional: scatter plot for log-transformed viral load
ggplot(hei_data, aes(x = Marker1, y = log_viral_load)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Marker1 vs Log Viral Load in HEI group", x = "Marker1 Expression", y = "Log Viral Load")

####
# Heatmap/Table of all the significant differences between HEI/HEU and between treatment 
#
#
#

library(tidyr)

# Step 1: Summarize the cluster.cols by Sample_Name and Ordered_Clusters
summarized_data <- all.flow %>%
  group_by(Sample_Name, Ordered_Clusters) %>%  # Group by Sample_Name and Cluster
  summarize(across(all_of(cluster.cols), ~ median(.x, na.rm = TRUE)),  # Calculate median for each cluster and marker
            SampleID = first(SampleID),       # Keep metadata columns as before
            Batch = first(Batch),
            PID = first(PID),
            Cohort = first(Cohort),
            Condition=first(Condition),
            HIV_Status = first(HIV_Status),
            Specific_Killing = first(Specific_Killing)) %>%
  ungroup()

# Step 2: Pivot the data so that each cluster has its own column for each marker
summarized_wide <- summarized_data %>%
  pivot_wider(names_from = Ordered_Clusters,  # Create new columns based on clusters
              values_from = all_of(cluster.cols),  # The marker columns to split by cluster
              names_glue = "{.value}_Cluster{Ordered_Clusters}")  # Create new column names (e.g., Marker1_Cluster1)

#### Statistical Testing: Differences in Marker Expression by HIV_Status ###

# Step 1: Define a function to run the Wilcoxon rank-sum test (or use t-test if appropriate) for each cluster
test_marker_by_hiv <- function(marker_cluster_col) {
  # Run the Wilcoxon test on the full dataset without grouping by HIV_Status
  result <- wilcox.test(summarized_wide[[marker_cluster_col]] ~ summarized_wide$HIV_Status)
  return(result$p.value)
}
# Step 2: Apply the function to each cluster for Marker1 (e.g., Marker1_Cluster1, Marker1_Cluster2, etc.)
marker1_cluster1_test <- test_marker_by_hiv("NKG2A_Cluster1")
marker1_cluster2_test <- test_marker_by_hiv("NKG2A_Cluster2")

wilcox.test(summarized_wide$NKG2A_Cluster1 ~ summarized_wide$HIV_Status)

# You can repeat this for other markers and clusters
# Step 1: Get all cluster columns (assuming "_Cluster" is in the name of the cluster columns)
cluster_cols <- grep("_Cluster", colnames(summarized_wide), value = TRUE)

# Step 2: Define a function to run the Wilcoxon test for each cluster column
run_wilcox_tests <- function(cluster_cols) {
  results <- lapply(cluster_cols, function(col) {
    # Run the Wilcoxon test for the current marker-cluster column
    test_result <- wilcox.test(summarized_wide[[col]] ~ summarized_wide$HIV_Status)
    
    # Store the p-value and column name
    data.frame(marker_cluster = col, p_value = test_result$p.value)
  })
  
  # Combine the results into a single data frame
  results_df <- bind_rows(results)
  
  return(results_df)
}

# Step 3: Run the Wilcoxon tests for all cluster columns
wilcox_results <- run_wilcox_tests(cluster_cols)
# Adjust p-values for multiple testing using Benjamini-Hochberg (FDR control)
wilcox_results$adjusted_p_value <- p.adjust(wilcox_results$p_value, method = "BH")

# Filter by adjusted p-value
significant_results <- wilcox_results %>%
  filter(adjusted_p_value < 0.05)


# View the significant results
print(significant_results)

### Plotting
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/FLORAH/Results/Output 4 - summary data/HEIvsHEU")

# Step 1: Loop through significant results
for (i in 1:nrow(significant_results)) {
  # Get marker-cluster and p-value
  marker_cluster <- significant_results$marker_cluster[i]
  p_value <- significant_results$p_value[i]
  
  # Create the plot for each significant marker-cluster
  plot <- ggplot(summarized_wide, aes(x = HIV_Status, y = get(marker_cluster), fill = HIV_Status)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, alpha = 0.7) +  # Boxplot
    geom_jitter(position = position_jitter(width = 0.2), size = 2, color = "black", alpha = 0.6) +  # Jittered points
    scale_fill_manual(values = c("HEU" = "#5bbae3", "HEI" = "#fc913f")) +  # Custom colors for HEU and HEI
    theme_minimal() +
    labs(title = paste("Expression of", marker_cluster, "by HIV Status"),
         x = "HIV Status", y = paste(marker_cluster, "Expression")) +
    annotate("text", x = 1.5, y = max(summarized_wide[[marker_cluster]], na.rm = TRUE), 
             label = paste("p-value =", signif(p_value, digits = 3)), size = 5, color = "red") +  # P-value annotation
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.position = "none"  # Hide legend
    )
  
  # Save the plot
  file_name <- paste0(marker_cluster, "_plot.png")
  ggsave(filename = file_name, plot = plot, width = 8, height = 7,bg='white')
}

# Step 2: Save the p-values table to a CSV file
write.csv(wilcox_results, file = "wilcox__results_p_values.csv", row.names = FALSE)

# View the saved p-values
head(significant_results)
levels(as.factor(summarized_wide$Condition))
?pairwise.wilcox.test

##### Treatment KW and Wilcoxins ###


# Step 1: Filter the data for relevant conditions (CEM, HUT78, K562, Untreated)
baseline_conditions <- summarized_wide %>%
  filter(Condition %in% c("CEM", "HUT78", "K562", "Untreated"))

# Step 2: Define the columns for clusters
cluster_cols <- grep("_Cluster", colnames(baseline_conditions), value = TRUE)

# Step 3: Function to run KW and pairwise Wilcoxon tests
run_kw_and_pairwise_tests <- function(cluster_cols) {
  results <- lapply(cluster_cols, function(col) {
    # Kruskal-Wallis test for overall comparison
    kw_test <- kruskal.test(baseline_conditions[[col]] ~ baseline_conditions$Condition)
    
    # Pairwise Wilcoxon rank-sum tests for each condition pair
    pairwise_tests <- pairwise.wilcox.test(baseline_conditions[[col]], baseline_conditions$Condition, p.adjust.method = "BH")
    
    # Flatten the pairwise p-value matrix into individual columns (ensure NA handling)
    pairwise_p_value_df <- as.data.frame(as.table(pairwise_tests$p.value)) %>%
      pivot_wider(names_from = Var2, values_from = Freq, names_prefix = "vs_")  # Rename pairs for clarity
    
    # Spread the `Var1` values across the rows and align the p-values
    pairwise_p_value_df <- pivot_wider(pairwise_p_value_df, names_from = Var1, values_from = starts_with("vs_"))
    
    # Combine KW p-value and pairwise p-values into one data frame
    combined_results <- data.frame(
      marker_cluster = col,
      kw_p_value = kw_test$p.value
    )
    
    combined_results <- cbind(combined_results, pairwise_p_value_df)
    
    return(combined_results)
  })
  
  bind_rows(results)
}

# Step 4: Run the tests for all cluster columns
baseline_results <- run_kw_and_pairwise_tests(cluster_cols)

# View the result
print(baseline_results)

# Step 5: Remove columns that contain only NA values
baseline_results_clean <- baseline_results %>%
  select(where(~ !all(is.na(.))))  # Removes columns where all values are NA

# View the cleaned result
print(baseline_results_clean)

# Save p-values
write.csv(baseline_results_clean, file = "significant_baseline_treatment_p_values.csv", row.names = FALSE)


#### TEST
# Set the marker-cluster for the test plot
marker_cluster <- "CD107a_Cluster2"

# Extract the Kruskal-Wallis p-value from the cleaned results
kw_p_value <- baseline_results_clean %>%
  filter(marker_cluster == !!marker_cluster) %>%
  pull(kw_p_value)

# Extract the pairwise Wilcoxon p-values for this marker-cluster
pairwise_p_values <- baseline_results_clean %>%
  filter(marker_cluster == !!marker_cluster) %>%
  select(starts_with("vs_")) %>%
  unlist()

# Define all pairwise comparisons between the four conditions
pairwise_conditions <- list(
  "CEM_HUT78" = c("CEM", "HUT78"),
  "CEM_K562" = c("CEM", "K562"),
  "CEM_Untreated" = c("CEM", "Untreated"),
  "HUT78_K562" = c("HUT78", "K562"),
  "HUT78_Untreated" = c("HUT78", "Untreated"),
  "K562_Untreated" = c("K562", "Untreated")
)
pairwise_p_values
# Filter only significant pairwise comparisons (p-value < 0.05)
significant_pairs <- names(pairwise_p_values)[pairwise_p_values < 0.05]
significant_condition_pairs <- lapply(significant_pairs, function(pair) pairwise_conditions[[gsub("vs_", "", pair)]])
significant_pairs
### Create filtered plotting data
# Reorder the Condition factor levels
filtered_data <- summarized_wide %>%
  filter(Condition %in% c("Untreated", "CEM", "HUT78", "K562"))
filtered_data$Condition <- factor(filtered_data$Condition, levels = c("Untreated", "CEM", "HUT78", "K562"))

# Set custom colors for the Condition groups
condition_colors <- c("Untreated" = "#8c510a", "CEM" = "#01665e", "HUT78" = "#d8b365", "K562" = "#abdda4")

# Calculate maximum y value for dynamic placement of significance bars
y_max <- max(filtered_data[[marker_cluster]], na.rm = TRUE)


# Calculate maximum y value for dynamic placement of significance bars
y_max <- max(filtered_data[[marker_cluster]], na.rm = TRUE)
# Create the plot with the filtered data
# Create the plot with the filtered data

# Create the plot with the filtered data

# Create the plot with the filtered data
ggplot(filtered_data, aes(x = Condition, y = get(marker_cluster), fill = Condition)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +  # Hide outliers and remove Condition legend
  geom_jitter(aes(color = HIV_Status), width = 0.2, size = 3, alpha = 0.6) +  # Color jittered dots by HIV_Status
  scale_color_manual(values = c("HEU" = "#5bbae3", "HEI" = "#fc913f")) +  # Custom colors for HIV_Status
  scale_fill_manual(values = condition_colors, guide = "none") +  # Custom colors for Conditions, hide legend
  theme_minimal() +
  labs(
    title = paste("Expression of", marker_cluster, "by Condition"),
    subtitle = paste("KW p-value =", signif(kw_p_value, digits = 3)),  # Move KW p-value to subtitle
    x = "Condition", y = paste(marker_cluster, "Expression")
  ) +
  geom_signif(
    comparisons = significant_condition_pairs, 
    map_signif_level = TRUE, 
    y_position = seq(y_max + 0.2, y_max + 0.5, length.out = length(significant_condition_pairs)),  # Dynamic y positions
    step_increase = 0.1,  # Spacing between bars
    textsize = 5
  ) +
  ylim(NA, y_max + 1.5) +  # Increase plot space dynamically
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "red"),  # Add subtitle for KW p-value
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

###
# Loop through each significant marker-cluster
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/FLORAH/Results/Output 4 - summary data/Treatment/Significant_No_Scale")

for (marker_cluster in baseline_results_clean$marker_cluster[baseline_results_clean$kw_p_value < 0.05]) {
  
  # Extract the Kruskal-Wallis p-value for the current marker-cluster
  kw_p_value <- baseline_results_clean %>%
    filter(marker_cluster == !!marker_cluster) %>%
    pull(kw_p_value)
  
  # Extract the pairwise Wilcoxon p-values for this marker-cluster
  pairwise_p_values <- baseline_results_clean %>%
    filter(marker_cluster == !!marker_cluster) %>%
    select(starts_with("vs_")) %>%
    unlist()
  
  # Filter only significant pairwise comparisons (p-value < 0.05)
  significant_pairs <- names(pairwise_p_values)[pairwise_p_values < 0.05]
  significant_condition_pairs <- lapply(significant_pairs, function(pair) pairwise_conditions[[gsub("vs_", "", pair)]])
  
  # Subset the data to include only the specified conditions
  filtered_data <- summarized_wide %>%
    filter(Condition %in% c("Untreated", "CEM", "HUT78", "K562"))
  
  # Reorder the Condition factor levels
  filtered_data$Condition <- factor(filtered_data$Condition, levels = c("Untreated", "CEM", "HUT78", "K562"))
  
  # Set custom colors for the Condition groups
  condition_colors <- c("Untreated" = "#8c510a", "CEM" = "#01665e", "HUT78" = "#d8b365", "K562" = "#abdda4")
  
  # Calculate maximum y value for dynamic placement of significance bars
  y_max <- max(filtered_data[[marker_cluster]], na.rm = TRUE)
  
  # Create the plot with the filtered data
  plot <- ggplot(filtered_data, aes(x = Condition, y = get(marker_cluster), fill = Condition)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +  # Hide outliers and remove Condition legend
    geom_jitter(aes(color = HIV_Status), width = 0.2, size = 3, alpha = 0.6) +  # Color jittered dots by HIV_Status
    scale_color_manual(values = c("HEU" = "#5bbae3", "HEI" = "#fc913f")) +  # Custom colors for HIV_Status
    scale_fill_manual(values = condition_colors, guide = "none") +  # Custom colors for Conditions, hide legend
    theme_minimal() +
    labs(
      title = paste("Expression of", marker_cluster, "by Condition"),
      subtitle = paste("KW p-value =", signif(kw_p_value, digits = 3)),  # Move KW p-value to subtitle
      x = "Condition", y = paste(marker_cluster, "Expression")
    ) +
    geom_signif(
      comparisons = significant_condition_pairs, 
      map_signif_level = TRUE, 
      y_position = seq(y_max + 0.2, y_max + 0.5, length.out = length(significant_condition_pairs)),  # Dynamic y positions
      step_increase = 0.1,  # Spacing between bars
      textsize = 5
    ) +
    ylim(NA, y_max + 1.5) +  # Increase plot space dynamically
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "red"),  # Add subtitle for KW p-value
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13)
    )
  
  # Save the plot to the working directory
  ggsave(filename = paste0(marker_cluster, "_expression_plot.png"), plot = plot, width = 8, height = 6, bg='white')
}
write.csv(baseline_results_clean , file = "Wilcox_KW_results.csv", row.names = FALSE)

#### Timepoint
# Define the list of conditions
conditions <- c("Untreated", "CEM", "HUT78", "K562")

# Filter the data for HEI first
# Subset the data for HEI and the specified conditions
hei_data <- subset(summarized_wide, HIV_Status == "HEI" & Condition %in% conditions)

# Identify PIDs that have both timepoints for each condition
paired_PIDs <- hei_data %>%
  group_by(PID, Condition) %>%
  filter(n_distinct(Timepoint) == 2) %>%
  pull(PID) %>%
  unique()


# Subset the data to include only paired PIDs
paired_hei_data <- hei_data[hei_data$PID %in% paired_PIDs, ]

# Get the list of marker columns (assuming markers start from the 5th column onwards)
marker_columns <- colnames(paired_hei_data)[11:ncol(paired_hei_data)]

# Create a list to store p-values for each marker
p_value_list <- list()

# Loop through each marker cluster and calculate p-values
for (marker_cluster in marker_columns) {
  
  # Perform paired Wilcoxon test for each condition and extract p-values
  p_values <- paired_hei_data %>%
    group_by(Condition) %>%
    summarise(p_value = wilcox.test(get(marker_cluster)[Timepoint == "Entry"], 
                                    get(marker_cluster)[Timepoint == "12 Months"], 
                                    paired = TRUE)$p.value)
  
  # Store p-values in the list, using marker_cluster as the key
  p_value_list[[marker_cluster]] <- p_values
}

### Plotting
paired_hei_data$Condition <- factor(paired_hei_data$Condition, levels = c("Untreated", "CEM", "HUT78", "K562"))

setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/FLORAH/Results/Output 4 - summary data/Timepoint/HEI")
# Define the colors for Conditions
condition_colors <- c("Untreated" = "#8c510a", "CEM" = "#01665e", "HUT78" = "#d8b365", "K562" = "#abdda4")
paired_hei_data$Condition_Timepoint <- factor(interaction(paired_hei_data$Condition, paired_hei_data$Timepoint, sep = "_"), 
                                              levels = c("Untreated_Entry", "Untreated_12 Months", 
                                                         "CEM_Entry", "CEM_12 Months", 
                                                         "HUT78_Entry", "HUT78_12 Months", 
                                                         "K562_Entry", "K562_12 Months"))
# Loop through each marker cluster, generate the plot, and save it
for (marker_cluster in marker_columns) {
  
  # Set the maximum y-axis value dynamically based on the data
  y_max <- max(paired_hei_data[[marker_cluster]], na.rm = TRUE)
  
  plot <- ggplot(paired_hei_data, aes(x = Condition_Timepoint, y = get(marker_cluster), fill = Condition)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +  # Remove Condition legend
    geom_jitter(aes(color = Timepoint), width = 0.2, size = 2, alpha = 0.6) +  # Color jittered dots by Timepoint
    scale_color_manual(values = c("Entry" = "#2E8B57", "12 Months" = "#fc8d62")) +  # Custom colors for Timepoint
    scale_fill_manual(values = condition_colors) +  # Custom colors for Conditions
    theme_minimal() +
    labs(
      title = paste("Expression of", marker_cluster, "by Condition and Timepoint (HEI)"),
      x = "Condition and Timepoint", y = paste(marker_cluster, "Expression")
    ) +
    ylim(NA, y_max + 1.5) +  # Increase plot space dynamically
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.position = "none"  # Remove legends for Timepoint and Condition
    ) +
    scale_x_discrete(labels = c("Untreated_Entry" = "Untreated (Entry)", 
                                "Untreated_12 Months" = "Untreated (12 Months)", 
                                "CEM_Entry" = "CEM (Entry)", 
                                "CEM_12 Months" = "CEM (12 Months)", 
                                "HUT78_Entry" = "HUT78 (Entry)", 
                                "HUT78_12 Months" = "HUT78 (12 Months)", 
                                "K562_Entry" = "K562 (Entry)", 
                                "K562_12 Months" = "K562 (12 Months)"))  # Custom x-axis labels
  
  # Add significance bars comparing timepoints within each condition
  for (cond in c("Untreated", "CEM", "HUT78", "K562")) {
    p_value <- p_value_list[[marker_cluster]]$p_value[p_value_list[[marker_cluster]]$Condition == cond]
    
    # If p-value is significant, add the geom_signif between timepoints for each condition
    if (p_value < 0.05) {
      plot <- plot + 
        geom_signif(
          comparisons = list(c(paste0(cond, "_Entry"), paste0(cond, "_12 Months"))),  # Compare timepoints within each condition
          annotations = "*",  # Add star for significant p-values
          map_signif_level = FALSE,  # Use custom annotations
          y_position = y_max + 0.5,  # Set position for significance bar above the plot
          step_increase = 0.1,  # Spacing between bars
          textsize = 5 )
    }
  }
  # Save the plot as a PNG file in the working directory
  ggsave(filename = paste0(marker_cluster, "_plot.png"), plot = plot, width = 14, height = 7, bg='White')
}

### Save p-values
# Initialize an empty data frame to store the p-values
p_value_df <- data.frame(Marker = character(), Condition = character(), P_Value = numeric(), stringsAsFactors = FALSE)

# Loop through each marker in the p_value_list
for (marker_cluster in names(p_value_list)) {
  # Extract the data frame for each marker
  marker_p_values <- p_value_list[[marker_cluster]]
  
  # Add the marker name to the data frame
  marker_p_values$Marker <- marker_cluster
  
  # Combine this marker's p-values with the overall p_value_df
  p_value_df <- rbind(p_value_df, marker_p_values)
}

# Save the combined p-value data frame to a CSV file
write.csv(p_value_df, file = "p_values_output.csv", row.names = FALSE)

setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/FLORAH/Results/Output 4 - summary data/Timepoint/HEU")

# Subset the data for HEU
heu_data <- subset(summarized_wide, HIV_Status == "HEU"& Condition %in% conditions)

# Ensure correct order for Condition and Timepoint
heu_data$Condition_Timepoint <- factor(interaction(heu_data$Condition, heu_data$Timepoint, sep = "_"), 
                                       levels = c("Untreated_Entry", "Untreated_12 Months", 
                                                  "CEM_Entry", "CEM_12 Months", 
                                                  "HUT78_Entry", "HUT78_12 Months", 
                                                  "K562_Entry", "K562_12 Months"))

# Loop through each marker cluster, generate the plot, and save it
for (marker_cluster in marker_columns) {
  
  # Set the maximum y-axis value dynamically based on the data
  y_max <- max(heu_data[[marker_cluster]], na.rm = TRUE)
  
  plot <- ggplot(heu_data, aes(x = Condition_Timepoint, y = get(marker_cluster), fill = Condition)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +  # Remove Condition legend
    geom_jitter(aes(color = Timepoint), width = 0.2, size = 2, alpha = 0.6) +  # Color jittered dots by Timepoint
    scale_color_manual(values = c("Entry" = "#2E8B57", "12 Months" = "#fc8d62")) +  # Custom colors for Timepoint
    scale_fill_manual(values = condition_colors) +  # Custom colors for Conditions
    theme_minimal() +
    labs(
      title = paste("Expression of", marker_cluster, "by Condition and Timepoint (HEU)"),
      x = "Condition and Timepoint", y = paste(marker_cluster, "Expression")
    ) +
    ylim(NA, y_max + 1.5) +  # Increase plot space dynamically
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.position = "none"  # Remove legends for Timepoint and Condition
    ) +
    scale_x_discrete(labels = c("Untreated_Entry" = "Untreated (Entry)", 
                                "Untreated_12 Months" = "Untreated (12 Months)", 
                                "CEM_Entry" = "CEM (Entry)", 
                                "CEM_12 Months" = "CEM (12 Months)", 
                                "HUT78_Entry" = "HUT78 (Entry)", 
                                "HUT78_12 Months" = "HUT78 (12 Months)", 
                                "K562_Entry" = "K562 (Entry)", 
                                "K562_12 Months" = "K562 (12 Months)"))  # Custom x-axis labels
  
  # Add significance bars comparing timepoints within each condition
  for (cond in c("Untreated", "CEM", "HUT78", "K562")) {
    p_value <- p_value_list[[marker_cluster]]$p_value[p_value_list[[marker_cluster]]$Condition == cond]
    
    # If p-value is significant, add the geom_signif between timepoints for each condition
    if (p_value < 0.05) {
      plot <- plot + 
        geom_signif(
          comparisons = list(c(paste0(cond, "_Entry"), paste0(cond, "_12 Months"))),  # Compare timepoints within each condition
          annotations = "*",  # Add star for significant p-values
          map_signif_level = FALSE,  # Use custom annotations
          y_position = y_max + 0.5,  # Set position for significance bar above the plot
          step_increase = 0.1,  # Spacing between bars
          textsize = 5 )
    }
  }
  # Save the plot as a PNG file in the working directory
  ggsave(filename = paste0("HEU_", marker_cluster, "_plot.png"), plot = plot, width = 14, height = 7, bg='White')
}

# Create an empty data frame to store the p-values for HEU
p_value_df_heu <- data.frame(Marker = character(), Condition = character(), P_Value = numeric(), stringsAsFactors = FALSE)

# Loop through each marker in the p_value_list for HEU
for (marker_cluster in names(p_value_list)) {
  # Extract the data frame for each marker
  marker_p_values <- p_value_list[[marker_cluster]]
  
  # Add the marker name to the data frame
  marker_p_values$Marker <- marker_cluster
  
  # Combine this marker's p-values with the overall p_value_df for HEU
  p_value_df_heu <- rbind(p_value_df_heu, marker_p_values)
}

# Save the p-values data frame to a CSV file for HEU
write.csv(p_value_df_heu, file = "p_values_output_HEU.csv", row.names = FALSE)


## Corelation to clinical variables
# Subset the data for HUT78 and K562, include both HEU and HEI
subset_data <- summarized_wide %>% 
  filter(Condition %in% c("HUT78", "K562")) %>%
  select(Condition, Timepoint, HIV_Status, Specific_Killing, all_of(marker_columns))  # Include specific Killing and marker columns

# Define the marker cluster you're interested in (for example, "marker_1_Cluster1")
marker_cluster <- "NKG2A_Cluster1"

# Create a scatter plot for both HUT78 and K562
ggplot(subset_data, aes_string(x = marker_cluster, y = "Specific_Killing", color = "Condition")) +
  geom_point(alpha = 0.6) +  # Add points
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
  labs(
    title = paste("Relation between", marker_cluster, "Expression and Specific Killing"),
    x = paste(marker_cluster, "Expression"),
    y = "Specific Killing"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  ) +
  scale_color_manual(values = c("HUT78" = "#d8b365", "K562" = "#abdda4"))  # Custom colors for conditions


#####
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/FLORAH/Results/Output 4 - summary data/Specific_Killing")

# Initialize a data frame to store correlations
correlation_results <- data.frame(Marker = character(), Condition = character(), Correlation = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

# Loop through each marker cluster and calculate correlations for both conditions
for (marker_cluster in marker_columns) {
  for (cond in c("HUT78", "K562")) {
    # Subset data for each condition
    condition_data <- subset(subset_data, Condition == cond)
    
    # Calculate the correlation between marker expression and Specific_Killing
    correlation_test <- cor.test(condition_data[[marker_cluster]], condition_data$Specific_Killing, method = "pearson")  # Use Pearson if you prefer linear correlation
    
    # Store the results in the data frame
    correlation_results <- rbind(correlation_results, 
                                 data.frame(Marker = marker_cluster, 
                                            Condition = cond, 
                                            Correlation = correlation_test$estimate, 
                                            P_Value = correlation_test$p.value))
  }
}

# Add a column for significance stars based on p-value
correlation_results$Significance <- ifelse(correlation_results$P_Value < 0.001, "***",
                                           ifelse(correlation_results$P_Value < 0.01, "**",
                                                  ifelse(correlation_results$P_Value < 0.05, "*", "")))

# Filter the results to keep only significant correlations (p < 0.05)
significant_results <- correlation_results %>% 
  filter(P_Value < 0.05)

# Create a heatmap for significant correlations only with stars for significance
cp <- ggplot(significant_results, aes(x = Condition, y = Marker, fill = Correlation)) +
  geom_tile(color = "white") +  # Create the heatmap tiles
  scale_fill_gradient2(low = "#66c2a5", high = "#fc8d62", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Pearsons Correlation") +
  geom_text(aes(label = Significance), color = "black", size = 5) +  # Add stars for significance
  theme_minimal() +
  labs(
    title = "Significant Correlations of Marker Expression with Specific Killing",
    x = "Condition", 
    y = "Marker Cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
ggsave("Pearsons_Cor_plot.png", plot = cp, width = 7, height = 9, bg='White')


### Individual Plotting
# Filter for significant marker clusters in either HUT78 or K562
significant_markers <- correlation_results %>%
  filter(P_Value < 0.05 & Condition %in% c("HUT78", "K562")) %>%
  distinct(Marker)  # Get unique marker clusters that are significant in either condition

# Reshape the data from wide to long format
long_data <- subset_data %>%
  pivot_longer(cols = all_of(marker_columns), 
               names_to = "Marker", 
               values_to = "Expression")
# Subset the data to include only significant marker clusters
significant_data <- long_data %>%
  filter(Marker %in% significant_markers$Marker)




# Loop through each significant marker cluster to plot the relationship
for (marker_cluster in significant_markers$Marker) {
  
  # Extract p-values for both HUT78 and K562 for the current marker
  p_value_hut78 <- correlation_results %>%
    filter(Marker == marker_cluster & Condition == "HUT78") %>%
    pull(P_Value)
  
  p_value_k562 <- correlation_results %>%
    filter(Marker == marker_cluster & Condition == "K562") %>%
    pull(P_Value)
  
  # Create a scatter plot faceted by condition, showing the p-values
  plot <- ggplot(significant_data %>% filter(Marker == marker_cluster), aes(x = Expression, y = Specific_Killing)) +
    geom_point(alpha = 0.6, aes(color = Condition)) +  # Scatter plot points
    geom_smooth(method = "lm", se = FALSE, aes(color = Condition)) +  # Linear regression line
    facet_wrap(~Condition) +  # Facet side by side for HUT78 and K562
    labs(
      title = paste("Relationship between", marker_cluster, "Expression and Specific Killing"),
      subtitle = paste("HUT78 p-value:", signif(p_value_hut78, 3), " | K562 p-value:", signif(p_value_k562, 3)),
      x = paste(marker_cluster, "Expression"),
      y = "Specific Killing"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13)
    )
  
  # Display the plot
  print(plot)
}
###
# Loop through each significant marker cluster to plot the relationship
# Define custom colors for the conditions
custom_colors <- c("HUT78" = "#d8b365", "K562" = "#abdda4")
for (marker_cluster in significant_markers$Marker) {
  
  # Extract p-values for both HUT78 and K562 for the current marker
  p_value_hut78 <- correlation_results %>%
    filter(Marker == marker_cluster & Condition == "HUT78") %>%
    pull(P_Value)
  
  p_value_k562 <- correlation_results %>%
    filter(Marker == marker_cluster & Condition == "K562") %>%
    pull(P_Value)
  
  # Create a scatter plot faceted by condition, showing the p-values
  plot <- ggplot(significant_data %>% filter(Marker == marker_cluster), aes(x = Expression, y = Specific_Killing, color = Condition)) +
    geom_point(size = 3, alpha = 0.8) +  # Scatter plot points
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +  # Linear regression line
    facet_wrap(~Condition) +  # Facet side by side for HUT78 and K562
    scale_color_manual(values = custom_colors) +  # Apply custom colors
    labs(
      title = paste(marker_cluster, "Expression Vs Specific Killing"),
      subtitle = paste("HUT78 p-value:", signif(p_value_hut78, 3), " | K562 p-value:", signif(p_value_k562, 3)),
      x = paste(marker_cluster, "Expression"),
      y = "Specific Killing"
    ) +
    theme_minimal(base_size = 14) +  # Minimal theme with larger text
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Center and bold title
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "darkblue"),  # Subtitle with color
      axis.text = element_text(size = 12, face = "bold"),  # Bold axis text
      axis.title = element_text(size = 15, face = "bold"),  # Bold axis labels
      legend.position = "none",  # Remove legend (since conditions are faceted)
      panel.grid.major = element_line(color = "gray80", linewidth = 0.5),  # Custom grid lines
      strip.background = element_rect(fill = "#606060"),  # Facet label background
      strip.text = element_text(size = 14, face = "bold", color = "white")  # Customize facet label text
    )
  
  # Display the plot
  ggsave(filename = paste0("Plot_", marker_cluster, ".png"), plot = plot, width = 12, height = 7, dpi = 300, bg='White')
  
}
