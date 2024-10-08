library(data.table)
library(readxl)
library(flowCore)
library(Biobase)
library(flowStats)
library(CytoNorm)
library(ggplot2)
library(umap) # For visualization
library(cyCombine)
library(uwot)
library(Spectre)
library(sva)

############################################ Global Variables #######################################################
setwd("C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/FLOWSOM/TARA/FCS/All")
in.path <-"C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/FLOWSOM/TARA/FCS/"
fcs_path <-"C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/FLOWSOM/TARA/FCS/All"

##### Read in flow files ###
T_metadata <- fread(paste0(in.path,"T_Metadata.csv"))
T_metadata[, File_Name := gsub("\\.fcs$", "_Total_NK.fcs", File_Name)]

# Get the list of FCS files that correspond to the metadata
files_to_load <- file.path(fcs_path, T_metadata$File_Name)

# Load only the FCS files specified in the metadata
data_list <- LoadFCS(files_to_load)

# Check the structure of the loaded data to ensure correct files are loaded
str(data_list)
# Ensure the row names of T_metadata match the sample names in the flowSet

T_metadata <- as.data.frame(T_metadata)
rownames(T_metadata) <- T_metadata$File_Name



################################## Quality Control ####################################################
###
#Flow Rate Check (flow_rate_qc): Detects fluctuations in flow rate, which can indicate technical issues.
#Signal Acquisition Check (signal_qc): Identifies abrupt changes in signal intensities.
#Dynamic Range Check (dynamic_range_qc): Identifies issues with fluorescence channels, such as abrupt changes in intensity.

### Flow AI
#library(flowAI)

#flow_auto_qc('H7 CP008 V12 NK+ CEM+IL-15.fcs', output = 0, remove_from = "all", second_fractionFR = 0.02,fcs_QC=FALSE, fcs_highQ="_hQC")



#### Create Flowset out of QCed files
##########

# Function to read and correct the $TOT value
correct_and_read_fcs <- function(file_path) {
  # Read the FCS file without any transformations
  fcs_data <- read.FCS(file_path, transformation = FALSE,truncate_max_range = FALSE )
  
  # Get the actual number of events in the data
  actual_events <- nrow(exprs(fcs_data))
  
  # Update the $TOT keyword to the actual number of events
  fcs_data@description[["$TOT"]] <- as.character(actual_events)
  
  # Return the corrected flowFrame
  return(fcs_data)
}

# Read and correct each FCS file into a list of flowFrames
fcs_files_corrected <- lapply(T_metadata$File_Name, function(file_name) {
  file_path <- file.path(fcs_path, file_name)
  correct_and_read_fcs(file_path)
})

# Name the list elements based on File_Name for easier identification
names(fcs_files_corrected) <- T_metadata$File_Name

###########
#Specter Pipeline
# Convert the list of flowFrames to a list of data frames
data_list <- lapply(fcs_files_corrected, function(flowFrame) {
  # Convert flowFrame to a data frame
  df <- as.data.frame(exprs(flowFrame))
  # Add a SampleID column based on the file name or any metadata of choice
  df$SampleID <- identifier(flowFrame)
  return(df)
})

# Combine all data frames into a single aggregated data frame
combined_data <- do.call(rbind, data_list)

# Add additional metadata from T_metadata if necessary (e.g., Batch, Treatment)
# Ensure that the column 'SampleID' in combined_data matches with 'File_Name' in T_metadata
combined_data <- merge(combined_data, T_metadata, by.x = "SampleID", by.y = "File_Name", all.x = TRUE)

# Check the combined data structure
str(combined_data)

### Data Transformation


setwd("C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/FLOWSOM/TARA/Results/QC")
dir.create("Output 1 - transformed plots")
setwd("Output 1 - transformed plots")

as.matrix(names(combined_data))


to.asinh <- names(combined_data)[c(8:33)]
cofactor <- 500


combined_data <- as.data.table(combined_data)
combined_data_norm <- do.asinh(combined_data, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")




for(i in transformed.cols){
  make.colour.plot(do.subsample(combined_data_norm, 20000), i, plot.against)
}

### Define cellular columns
# Create a named vector for the mapping
col_mapping <- c(
  "FJComp-APC-A_asinh" = "NKG2D",
  "FJComp-APC-Fire 750-A_asinh" = "SIGLEC-7",
  "FJComp-Alexa Fluor 594-A_asinh" = "NKG2A",
  "FJComp-Alexa Fluor 647-A_asinh" = "Granzyme B",
  "FJComp-Alexa Fluor 700-A_asinh" = "TNF-a",
  "FJComp-BB700-A_asinh" = "IFN-g",
  "FJComp-BUV395-A_asinh" = "NKG2C",
  "FJComp-BUV496-A_asinh" = "CD3",
  "FJComp-BUV563-A_asinh" = "CD56",
  "FJComp-BUV615-A_asinh" = "KLRG1",
  "FJComp-BUV737-A_asinh" = "CD25",
  "FJComp-BV421-A_asinh" = "CD57",
  "FJComp-BV480-A_asinh" = "CD16",
  "FJComp-BV510-A_asinh" = "CD122",
  "FJComp-BV605-A_asinh" = "Ki-67",
  "FJComp-BV650-A_asinh" = "NKp46",
  "FJComp-BV711-A_asinh" = "NKp30",
  "FJComp-BV786-A_asinh" = "CD127",
  "FJComp-Cell Trace Violet-A_asinh" = "CT Violet",
  "FJComp-FITC-A_asinh" = "MIP-1b",
  "FJComp-LIVE DEAD Blue-A_asinh" = "LIVE/DEAD",
  "FJComp-PE-A_asinh" = "Perforin",
  "FJComp-PE-Cy5-A_asinh" = "CD107a",
  "FJComp-PE-Cy7-A_asinh" = "FASL (CD178)",
  "FJComp-PerCP-Cy5.5-A_asinh" = "CX3CR1",
  "FJComp-PerCP-eFluor 710-A_asinh" = "TIGIT"
)

# Rename the columns in combined_data_norm using the mapping
names(combined_data_norm)[names(combined_data_norm) %in% names(col_mapping)] <- col_mapping[names(combined_data_norm)[names(combined_data_norm) %in% names(col_mapping)]]

# Verify the new column names
names(combined_data_norm)

as.matrix(names(combined_data_norm))
cellular.cols <- names(combined_data_norm)[c(44:69)]
as.matrix(cellular.cols)
cluster.cols <- names(combined_data_norm)[c(44:50,52:61,63,65:69)]
as.matrix(cluster.cols)
names(combined_data_norm)[c(44:50,52:61,63,65:69)]
sample.col <- 'Sample_Name'
treatment.col <- "Condition"
timepoint.col <- "Timepoint"
HIV.col <- "HIV_Status"
batch.col <- "Batch"

data.frame(table(combined_data_norm[[batch.col]])) # Check number of cells per sample.

### Remove samples with <1000 cells

# Create a table of the number of cells per sample
cell_counts <- data.frame(table(combined_data_norm[[sample.col]]))

# Filter the samples with fewer than 1000 cells
samples_to_remove <- cell_counts$Var1[cell_counts$Freq < 1000]
samples_to_remove
# Step 2: Remove these samples from the combined data
combined_data_filtered <- combined_data_norm[!combined_data_norm[[sample.col]] %in% samples_to_remove, ]

data.frame(table(combined_data_filtered[[sample.col]])) # Check number of cells per sample.

data.frame(table(combined_data_filtered[[batch.col]]))

#### Run Flowsom NO BATCH EFFECT CORRECTION 
combined_data_filtered_f <-run.flowsom(combined_data_filtered, cluster.cols)

setwd("C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/FLOWSOM/TARA/Results/QC")
dir.create("Output 2 - Batch Effect")
setwd("Output 2 - Batch Effect")

### Dimensionality reduction
sub.targets <- c(5000,5000,5000,5000,5000,5000,5000,5000)
batch.sub <- do.subsample(combined_data_filtered_f, sub.targets, batch.col)
batch.sub <- run.umap(batch.sub, cluster.cols)

make.colour.plot(batch.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE,
                 title = 'Pre Batch effect Correction - By Flowsom',dot.size=0.5,plot.width	=13)

make.colour.plot(batch.sub, 'UMAP_X', 'UMAP_Y', 'Batch', 'factor', title = 'Pre Batch effect Correction - By Batch'
                 ,dot.size=0.5,plot.width	=13)

### Remove Batch 3 and HUU
# Remove samples with Batch 3 and HIV_Status HUU
combined_data_filtered_2 <- combined_data_filtered[!(combined_data_filtered$Batch == 'Batch 3' | combined_data_filtered$HIV_Status == 'HUU'), ]

# Check the structure to confirm the removal
data.frame(table(combined_data_filtered_2$Batch))
data.frame(table(combined_data_filtered_2$HIV_Status))

combined_data_filtered_f_2 <-run.flowsom(combined_data_filtered_2, cluster.cols)

### Dimensionality reduction
sub.targets <- c(5000,5000,5000,5000,5000,5000,5000)
batch.sub <- do.subsample(combined_data_filtered_f_2, sub.targets, batch.col)
batch.sub <- run.umap(batch.sub, cluster.cols)

make.colour.plot(batch.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE,
                 title = 'Pre Batch effect Correction - Flowsom (-3)',dot.size=0.5,plot.width	=13)

make.colour.plot(batch.sub, 'UMAP_X', 'UMAP_Y', 'Batch', 'factor', title = 'Pre Batch effect Correction - By Batch (-3)'
                 ,dot.size=0.5,plot.width	=13)


### Correct for Batch Effect

# Scale by batch 
combined_data_filtered_2 <- as.data.frame(combined_data_filtered_2)
# Normalize each marker/column by z-score normalization
combined_data_filtered_2[, cluster.cols] <- scale(combined_data_filtered_2[, cluster.cols])

# Create the matrix of columns to correct (using the vector cluster.cols)
expr_matrix <- as.matrix(combined_data_filtered_2[, cluster.cols])

# Create the batch vector (from the Batch column)
batch_vector <- combined_data_filtered_2$Batch
mod <- model.matrix(~ HIV_Status + Timepoint + Condition, data = combined_data_filtered_f_2)


# Apply ComBat to adjust for batch effects
combat_corrected <- ComBat(dat = t(expr_matrix), batch = batch_vector, mod = mod, par.prior = TRUE, prior.plots = FALSE)

# Transpose the result back to the original format
combat_corrected <- t(combat_corrected)

# Replace the original data in combined_data_filtered_f_2 with the batch-corrected values
batch.corrected <- combined_data_filtered_2
batch.corrected[, cluster.cols] <- combat_corrected
batch.corrected <- as.data.table(batch.corrected)
batch.corrected_f <-run.flowsom(batch.corrected, cluster.cols)

# Calculate the size of each FlowSOM_metacluster
cluster_sizes <- table(batch.corrected_f$FlowSOM_metacluster)

# Plot Cluster sizes
# Create a data frame from the cluster_sizes vector
cluster_sizes_df <- data.frame(
  Cluster = factor(names(cluster_sizes)),  # Cluster names (as factors for discrete x-axis)
  Size = as.numeric(cluster_sizes)         # Cluster sizes
)

# Create the bar plot using ggplot2
cluster_size_plot <- ggplot(cluster_sizes_df, aes(x = Cluster, y = Size)) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Bar plot with blue color
  geom_text(aes(label = Size), vjust = -0.5) +  # Add labels above the bars
  labs(title = "Cluster Sizes", x = "Cluster", y = "Size (Number of Cells)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability


# Save the plot using ggsave
ggsave("cluster_sizes_plot.png", plot = cluster_size_plot, width = 8, height = 6, dpi = 300, bg='white')
# Map the cluster sizes back to the original data frame as a new column
batch.corrected_f$Cluster_Size <- cluster_sizes[batch.corrected_f$FlowSOM_metacluster]

large_clusters <- names(cluster_sizes[cluster_sizes >= 20000])

# Step 3: Filter the data to keep only rows from large clusters
batch.corrected_filtered <- batch.corrected_f[batch.corrected_f$FlowSOM_metacluster %in% large_clusters, ]

# Step 4: Verify the filtered data
table(batch.corrected_filtered$FlowSOM_metacluster)
# Verify the new column
head(batch.corrected_f)

# Order the clusters by size (from largest to smallest)
ordered_clusters <- names(sort(cluster_sizes, decreasing = TRUE))
cluster_ranks <- setNames(seq_along(ordered_clusters), ordered_clusters)
# Create a new column with clusters ordered by size
batch.corrected_filtered$Ordered_Clusters <- cluster_ranks[as.character(batch.corrected_filtered$FlowSOM_metacluster)]


### Dimensionality reduction
sub.targets <- c(5000,5000,5000,5000,5000,5000,5000)
batch.sub <- do.subsample(batch.corrected_filtered, sub.targets, batch.col)
batch.sub <- run.umap(batch.sub, cluster.cols)


make.colour.plot(batch.sub, "UMAP_X", "UMAP_Y", "Ordered_Clusters", col.type = 'factor', add.label = TRUE,
                 title = 'Post-Batch Correction - Flowsom_size',dot.size=0.5,plot.width	=13)
make.colour.plot(batch.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE,
                 title = 'Post-Batch Correction - Flowsom_MST',dot.size=0.5,plot.width	=13)

make.colour.plot(batch.sub, 'UMAP_X', 'UMAP_Y', 'Batch', 'factor', title = 'Post Batch effect Correction - By Batch'
                 ,dot.size=0.5,plot.width	=13)

setwd("C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/FLOWSOM/TARA/Results/QC")
dir.create("Output 3 - Annotation")
setwd("Output 3 - Annotation")
make.multi.plot(batch.sub, 'UMAP_X', 'UMAP_Y', cluster.cols)

  ## Make heatmap
### Expression heatmap

exp <- do.aggregate(batch.corrected_filtered, cluster.cols, by = "Ordered_Clusters")
make.pheatmap(exp, "Ordered_Clusters", cluster.cols,transpose=TRUE)

setwd("C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/FLOWSOM")
saveRDS(batch.corrected_filtered, file = "batch_corrected_filtered_FlowSOM.rds")
saveRDS(batch.corrected, file = "batch_corrected_data.rds")
