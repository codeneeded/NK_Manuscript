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

############################################ Global Variables #######################################################
setwd("C:/Users/ammas/Documents/NK_Manuscript/FLOWSOM/TARA/")
in.path <-"C:/Users/ammas/Documents/NK_Manuscript/R_dat/"

all.flow <- readRDS(paste0(in.path,"batch_corrected_filtered_FlowSOM.rds"))



######
setwd("C:/Users/ammas/Documents/NK_Manuscript/FLOWSOM/TARA/Results/")
dir.create("Comparison_Plots")
setwd("Comparison_Plots")

### Check cells
sample.col <- 'Sample_Name'
treatment.col <- "Condition"
timepoint.col <- "Timepoint"
HIV.col <- "HIV_Status"
batch.col <- "Batch"
all.flow$Condition <- gsub("HUT\\+IL-15", "HUT78+IL-15", all.flow$Condition)

data.frame(table(all.flow[[treatment.col]])) # Check number of cells per sample.
cluster.cols <- names(all.flow)[c(44:50,52:61,63,65:69)]

### Dimensionality reduction
sub.targets <- c(50000,50000,50000,50000,50000,50000,50000, 50000)
flow.sub <- do.subsample(all.flow, sub.targets, treatment.col)
flow.sub <- run.umap(flow.sub, cluster.cols)


make.multi.plot(flow.sub, 'UMAP_X', 'UMAP_Y', 'Ordered_Clusters',col.type = 'factor'
                ,divide.by = 'Condition', figure.title= "Split By Treatment")

make.multi.plot(flow.sub, 'UMAP_X', 'UMAP_Y', 'Ordered_Clusters',col.type = 'factor'
                ,divide.by = 'HIV_Status', figure.title= "Split By Treatment")


### Summary Stats
setwd("C:/Users/ammas/Documents/NK_Manuscript/FLOWSOM/TARA/Results")

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
### Reorder summary data and SAVE

sum.dat <- do.reorder(sum.dat, treatment.col, grp.order.2)
fwrite(sum.dat, 'sum.dat_treatment.csv')

### Autographs

for(i in plot.cols){
  
  measure <- gsub("\\ --.*", "", i)
  measure
  
  pop <- gsub("^[^--]*.-- ", "", i)
  pop
  
  make.autograph(sum.dat,
                 x.axis = HIV.col,
                 y.axis = i,
                 y.axis.label = measure, 
                 violin = FALSE, 
                 colour.by = treatment.col,
                 colours = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F"),
                 
                 grp.order = grp.order.1,
                 my_comparisons = comparisons.1,
                 
                 Variance_test = variance.test,
                 Pairwise_test = pairwise.test,
                 
                 dot.size = 3,
                 width = 7,
                 
                 title = paste0(pop," ", measure),
                 filename = paste0(i, '_by_HIV_Status.pdf'))
  
}
dir.create("Treatment")
setwd("Treatment")
# Define the vector of treatments
treatments <- c("CEM", "HUT78","K562",'Untreated')

# Generate combinations of pairs
pair_list <- combn(treatments, 2, simplify = FALSE)

grp.order.3 <- c("Untreated","CEM","HUT78", "K562","CEM+IL-15" ,"HUT78+IL-15" , "K562+IL-15","IL-15")

for(i in plot.cols){
  
  measure <- gsub("\\ --.*", "", i)
  measure
  
  pop <- gsub("^[^--]*.-- ", "", i)
  pop
  
  make.autograph(sum.dat,
                 x.axis = treatment.col,
                 y.axis = i,
                 y.axis.label = measure, 
                 violin = FALSE, 
                 colour.by = HIV.col,

                 grp.order = grp.order.3,
                 my_comparisons = pair_list,
                 
                 Variance_test = variance.test,
                 Pairwise_test = pairwise.test,
                 
                 dot.size = 3,
                 width = 9,
                 height = 7,
                 max.y = 2,
                 
                 title = paste0(pop," ", measure),
                 filename = paste0(i, '_by_Treatment_Status.pdf'))
  }


# Display the list of pairs
print("K562")

il15_comparison_list <- list(
  c("CEM", "CEM+IL-15"),
  c("HUT78", "HUT78+IL-15"),
  c("K562", "K562+IL-15"),
  c("Untreated", "IL-15")  # This pair compares untreated with pure IL-15 treatment
)

grp.order.2 <- c("CEM", "CEM+IL-15", "HUT78","HUT78+IL-15", "K562", "K562+IL-15","Untreated","IL-15")


setwd("C:/Users/ammas/Documents/NK_Manuscript/FLOWSOM/TARA/Results")
setwd("Output 4 - summary data")

dir.create("IL-15")
setwd("IL-15")
for(i in plot.cols){
  
  measure <- gsub("\\ --.*", "", i)
  measure
  
  pop <- gsub("^[^--]*.-- ", "", i)
  pop
  
  make.autograph(sum.dat,
                 x.axis = treatment.col,
                 y.axis = i,
                 y.axis.label = measure, 
                 violin = FALSE, 
                 colour.by = HIV.col,
                 
                 grp.order = grp.order.2,
                 my_comparisons = il15_comparison_list,
                 
                 Variance_test = variance.test,
                 Pairwise_test = pairwise.test,
                 
                 dot.size = 3,
                 width = 9,
                 height = 7,
                 max.y = 2,
                 
                 title = paste0(pop," ", measure),
                 filename = paste0(i, '_by_IL-15Treatment_Status.pdf'))
}

######### Correlation Analysis
hei_data <- subset(all.flow, HIV_Status == "HEI" & Timepoint == "Entry" & Condition == 'HUT78')
markers <- cluster.cols
clinical_vars <- c("Viral_Load", "Specific_Killing")  # Replace with your clinical data columns

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
pheatmap(cor_matrix_rounded, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = annotations,  # Use significance annotations
         main = "Marker-Clinical Variable Correlation Heatmap")
# Step 7: Visualize results

# Scatter plot for one marker (e.g., Marker1) vs viral load
ggplot(hei_data, aes(x = `FASL (CD178)`, y = Specific_Killing)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Marker1 vs Viral Load in HEI group", x = "Marker1 Expression", y = "Viral Load")
# Scatter plot for specific killing vs Marker1
ggplot(hei_data, aes(x = `FASL (CD178)`, y = Specific_Killing)) +
  geom_point(color = "blue", alpha = 0.6) +  # Points on the plot
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear regression line
  theme_minimal() +  # Use a clean theme
  labs(title = "Scatter Plot: Marker1 vs Specific Killing", 
       x = "Marker1 Expression", 
       y = "Specific Killing")

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
            Timepoint = first(Timepoint),
            HIV_Status = first(HIV_Status),
            Viral_Load = first(Viral_Load),
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
marker1_cluster1_test <- test_marker_by_hiv("NKG2D_Cluster1")
marker1_cluster2_test <- test_marker_by_hiv("NKG2D_Cluster2")

wilcox.test(summarized_wide$NKG2D_Cluster1 ~ summarized_wide$HIV_Status)

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
setwd("C:/Users/ammas/Documents/NK_Manuscript/FLOWSOM/TARA/Results/Output 4 - summary data/HEIvsHEU/Significant_No_Scale")

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
