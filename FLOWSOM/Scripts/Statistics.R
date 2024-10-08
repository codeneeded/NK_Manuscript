library(data.table)
library(readxl)
library(flowCore)
library(Biobase)
library(flowStats)
library(ggplot2)
library(umap) # For visualization
library(cyCombine)
library(uwot)
library(Spectre)
library(sva)

############################################ Global Variables #######################################################
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/TARA/")
in.path <-"C:/Users/axi313/Documents/NK_Manuscript/R_dat/"

all.flow <- readRDS(paste0(in.path,"batch_corrected_filtered_FlowSOM.rds"))



######
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/TARA/Results/")
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


exp <- do.aggregate(all.flow, cluster.cols, by = "Ordered_Clusters")
make.pheatmap(exp, "Ordered_Clusters", cluster.cols,transpose=TRUE)