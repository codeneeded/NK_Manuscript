library(readxl)
library(dplyr)
library(lme4)
library(ggplot2)
library(ggeffects)
library(lmerTest)
library(ggpubr)
library(gridExtra)
library(tidyr)
library(reshape2)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(plotly)
library(readr)
library(broom.mixed)
library(extrafont)
loadfonts(device = "win")

############################################ Global Variables #######################################################
setwd("C:/Users/ammas/Documents/NK_Manuscript")
in.path <-"C:/Users/ammas/Documents/NK_Manuscript/Saved_R_Data/"

load(paste0(in.path,"tara_freq_clean.RDS"))
load(paste0(in.path,"florah_freq_clean.RDS"))


### Comparison of Frequencies across subsets across conditions ###

### Subset untreated condition only

tara_untreated <- tara_Freq %>%
  filter(Treatment == "untreated",HIV != "HUU")
colnames(tara_untreated)

florah_untreated <- florah_Freq %>%
  filter(Treatment == "untreated",HIV != "HUU")
colnames(tara_untreated)

################ Wilcoxins HEI vs HEU plot ################################

### TARA
#plot wilcoxins test of all columns between HEi and HEU
setwd("C:/Users/ammas/Documents/NK_Manuscript/Boxplots/HEIvsHEU/TARA")

# Step 1: Get the list of columns to plot (all columns after column 20)
columns_to_plot <- colnames(tara_untreated)[21:ncol(tara_untreated)]



# Step 2: Loop through each column and create/save the plot
for (column_name in columns_to_plot) {
  plot <- ggplot(tara_untreated %>% filter(!is.na(!!sym(column_name))), 
                 aes(x = HIV, y = !!sym(column_name), fill = HIV)) +
    geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(color = HIV), width = 0.2, size = 2, alpha = 0.7) +
    facet_wrap(~ Timepoint, scales = "free_x", nrow = 1) +
    stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                       method = "wilcox.test", 
                       label = 'p.format') +
    labs(title = column_name,
         x = "HIV Status",
         y = "Expression Frequency (% of Total NK)") +
    theme_minimal() +
    scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3"),
                      labels = c("HEI" = "pHIV", "HEU" = "pHEU")) +
    scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3"),
                       labels = c("HEI" = "pHIV", "HEU" = "pHEU")) +
    scale_x_discrete(labels = c("HEI" = "pHIV", "HEU" = "pHEU")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", colour = 'black'),
      axis.title.x = element_text(size = 18, margin = margin(t = 15), colour = 'black'),
      axis.title.y = element_text(size = 18, colour = 'black'),
      axis.text.x = element_text(size = 18, colour = 'black'),
      axis.text.y = element_text(size = 16, colour = 'black'),
      strip.text = element_text(size = 17, colour = 'black'),
      strip.placement = "outside",
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
    )
  
  safe_column_name <- gsub("[^[:alnum:]_]", "_", column_name)
  
  ggsave(filename = paste0(safe_column_name, "_pHIVvspHEU.png"), 
         plot = plot, 
         width = 8.5, 
         height = 5, 
         dpi = 300, 
         bg = "white")
}

### FLORAH

#plot wilcoxins test of all columns between HEi and HEU
setwd("C:/Users/ammas/Documents/NK_Manuscript/Boxplots/HEIvsHEU/Florah")

# Step 1: Get the list of columns to plot (all columns after column 20)
columns_to_plot <- colnames(florah_untreated)[21:ncol(florah_untreated)]



# Step 2: Loop through each column and create/save the plot
for (column_name in columns_to_plot) {
  # Step 3: Create the plot for the current column
  plot <- ggplot(florah_untreated %>% filter(!is.na(!!sym(column_name))), 
                 aes(x = HIV, y = !!sym(column_name), fill = HIV)) +
    geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +  # Boxplot with adjusted width
    geom_jitter(aes(color = HIV), width = 0.2, size = 2, alpha = 0.7) +  # Jitter for individual points
    stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                       method = "wilcox.test", 
                       label = "p.format") +  # Adjust y position for significance stars
    labs(title = column_name,
         x = "HIV Status",
         y = paste("Expression Frequency (% of Total NK)")) +
    theme_minimal() +
    scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom colors for HEI and HEU
    scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom colors for jitter points
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold",colour = 'black'),  # Title formatting
      axis.title.x = element_text(size = 18, margin = margin(t = 15),colour = 'black'),  # Shifts x-axis title down
      axis.title.y = element_text(size = 18,colour = 'black'),  # y-axis title formatting
      axis.text.x = element_text(size = 18,colour = 'black'),  # Larger x-axis text for Group labels
      axis.text.y = element_text(size = 16,colour = 'black'),  # y-axis text size
      strip.text = element_text(size = 17,colour = 'black'),  # Facet label size
      strip.placement = "outside",  # Moves facet label outside
      panel.spacing = unit(1, "lines"),  # Adjusts space between facets
      legend.position = "none"  # Remove legend if not needed
    )
  # Create a safe filename by replacing special characters with an underscore
  safe_column_name <- gsub("[^[:alnum:]_]", "_", column_name)
  
  # Step 4: Save the plot using ggsave
  ggsave(filename = paste0( safe_column_name, "_HEIvsHEU.png"), 
         plot = plot, 
         width = 8.5, 
         height = 5, 
         dpi = 300, 
         bg = "white")  # Specify background color as white
}

### NOTE: NS. is displayed for p-values equal to 1, the string ns in all other (non-significant)

#################### Wilcoxins Splt by Timepoint (TARA ONLY) ####################################

#### Read in Files ####
setwd("C:/Users/ammas/Documents/NK_Manuscript")

nk_data <- read_csv("NK_Cell_data_cleaned_normalized_to_CD45_Lesley Reduced 061124.csv")


demographics <- read_csv("Viral_Titres.csv")
out.path.2 <- "C:/Users/ammas/Documents/NK_Manuscript/Boxplots/Timepoint_HIV/"

#### Clean and combine datasets ####

# Get the names of all columns
column_names <- names(demographics)

# Swap the second and third column names
column_names[c(2, 3)] <- column_names[c(3, 2)]

# Reorder the dataframe based on the updated column names
demographics <- demographics[column_names]

# Inner join the demographics and nk_data dataframes based on PID and Age
nkiller <- inner_join(demographics, nk_data, by = c("PID", "Age"))

nkiller <- as.data.frame(nkiller)

# Assuming 'group_variable' is your factor variable for groups in 'data_frame'
nkiller$Group <- as.factor(nkiller$Group)
nkiller$Group <- relevel(nkiller$Group, ref = "HEU")
str(nkiller)

# Create the plot
### Make plot df
nkiller_plot <- nkiller %>%
  mutate(timepoint.2 = case_when(
    Age %in% c(1, 2) ~ "Entry",      # Ages 1 and 2 become 'Entry'
    Age == 5 ~ "5 months",           # Age 5 becomes '5 months'
    Age == 10 ~ "10 months",         # Age 10 becomes '10 months'
    Age == 18 ~ "18 months"
  ))
# Reorder the 'timepoint.2' column as a factor

nkiller_plot$timepoint.2 <- factor(nkiller_plot$timepoint.2, 
                                   levels = c("Entry", "5 months", "10 months", "18 months"))

# Create the plot with -test and significance stars, split by Timepoint
# Create the plot with all facets in one row and custom order
library(ggplot2)
library(ggpubr)

# Define common fill/color mapping
group_colors <- c("HEI" = "#fc913f", "HEU" = "#5bbae3")
group_labels <- c("HEI" = "pHIV", "HEU" = "pHEU")

# List of plots with y-variable, title, y-label, and output name
plot_info <- list(
  list(y = "Total_NK", title = expression("Total NK in CD45"^"+"~"Lymphocytes"),
       ylab = "Total NK (%)", file = "Total_NK_percent"),
  
  list(y = "NKbright", title = expression("CD56"^"bright"*"CD16"^"-"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD56"^"bright"*"CD16"^"-"~"(%)"), file = "CD56bright_percent_CD45L"),
  
  list(y = "NK_CD56-CD16+", title = expression("CD56"^"-"*"CD16"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD56"^"-"*"CD16"^"+"~"(%)"), file = "CD56-CD16+_percent_CD45L"),
  
  list(y = "NKdim_CD16+", title = expression("CD56"^"dim"*"CD16"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD56"^"dim"*"CD16"^"+"~"(%)"), file = "CD56dimCD16+_percent_CD45L"),
  
  list(y = "Total_NK/NKG2A+", title = expression("NKG2A"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("NKG2A"^"+"~"(%)"), file = "NKG2A+_percent_CD45L"),
  
  list(y = "Total_NK/CD38+", title = expression("CD38"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD38"^"+"~"(%)"), file = "CD38+_percent_CD45L"),
  
  list(y = "NK_CD56-CD16+/NKG2A+", title = expression("CD56"^"-"*"CD16"^"+"~"NKG2A"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD56"^"-"*"CD16"^"+"~"NKG2A"^"+"~"(%)"), file = "CD56-CD16+NKG2A+_percent_CD45L"),
  
  list(y = "NK_CD56-CD16+/CD38+", title = expression("CD56"^"-"*"CD16"^"+"~"CD38"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD56"^"-"*"CD16"^"+"~"CD38"^"+"~"(%)"), file = "CD56-CD16+CD38+_percent_CD45L"),
  
  list(y = "NKdim_CD16+/NKG2A+", title = expression("CD56"^"dim"*"CD16"^"+"~"NKG2A"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD56"^"dim"*"CD16"^"+"~"NKG2A"^"+"~"(%)"), file = "CD56dimCD16+NKG2A+_percent_CD45L"),
  
  list(y = "NKdim_CD16+/CD38+", title = expression("CD56"^"dim"*"CD16"^"+"~"CD38"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD56"^"dim"*"CD16"^"+"~"CD38"^"+"~"(%)"), file = "CD56dimCD16+CD38+_percent_CD45L"),
  
  list(y = "NKbright/NKG2A+", title = expression("CD56"^"bright"*"CD16"^"-"*"NKG2A"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD56"^"bright"*"CD16"^"-"*"NKG2A"^"+"~"(%)"), file = "CD56bright_NKG2A+_percent_CD45L"),
  
  list(y = "NKbright/CD38+", title = expression("CD56"^"bright"*"CD16"^"-"*"CD38"^"+"~"in CD45"^"+"~"Lymphocytes"),
       ylab = expression("CD56"^"bright"*"CD16"^"-"*"CD38"^"+"~"(%)"), file = "CD38+_percent_CD45L")
)

# Loop through and generate plots
for (info in plot_info) {
  p <- ggplot(nkiller_plot, aes(x = Group, y = .data[[info$y]], fill = Group)) +
    geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +
    labs(title = info$title,
         x = "HIV Status",
         y = info$ylab) +
    theme_minimal() +
    scale_fill_manual(values = group_colors, labels = group_labels) +
    scale_color_manual(values = group_colors, labels = group_labels) +
    scale_x_discrete(labels = group_labels) +
    stat_compare_means(comparisons = list(c("HEI", "HEU")),
                       method = "wilcox.test",
                       label = "p.format") +
    facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", colour = 'black'),
      axis.title.x = element_text(size = 18, margin = margin(t = 15), colour = 'black'),
      axis.title.y = element_text(size = 18, colour = 'black'),
      axis.text.x = element_text(size = 18, colour = 'black'),
      axis.text.y = element_text(size = 16, colour = 'black'),
      strip.text = element_text(size = 17, colour = 'black'),
      strip.placement = "outside",
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
    )
  
  ggsave(filename = paste0(out.path.2, "HEIvsHEU_", info$file, ".png"),
         plot = p, height = 5, width = 8.5, dpi = 300, bg = "white")
}

########################### % of total NK   ##################################################################
### Adjust columns of interest to percent of Total NK from % of total CD45 Lymphocytes
nkiller_plot.2 <- nkiller_plot
nkiller_plot.2$Total_NK <- (nkiller_plot.2$Total_NK/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts
nkiller_plot.2$`NKdim_CD16+` <- ((nkiller_plot.2$`NKdim_CD16+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`NK_CD56-CD16+` <- ((nkiller_plot.2$`NK_CD56-CD16+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$NKbright<- ((nkiller_plot.2$NKbright/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`Total_NK/CD38+`<- ((nkiller_plot.2$`Total_NK/CD38+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`Total_NK/NKG2A+`<- ((nkiller_plot.2$`Total_NK/NKG2A+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`NKdim_CD16+/CD38+`<- ((nkiller_plot.2$`NKdim_CD16+/CD38+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`NKdim_CD16+/NKG2A+`<- ((nkiller_plot.2$`NKdim_CD16+/NKG2A+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`NK_CD56-CD16+/CD38+`<- ((nkiller_plot.2$`NK_CD56-CD16+/CD38+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`NK_CD56-CD16+/NKG2A+`<- ((nkiller_plot.2$`NK_CD56-CD16+/NKG2A+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`NKbright/CD38+`<- ((nkiller_plot.2$`NKbright/CD38+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`NKbright/NKG2A+`<- ((nkiller_plot.2$`NKbright/NKG2A+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)


group_colors <- c("HEI" = "#fc913f", "HEU" = "#5bbae3")
group_labels <- c("HEI" = "pHIV", "HEU" = "pHEU")

# Define all the plots for Total NK %
plot_totalnk_info <- list(
  list(y = "NKbright",
       title = expression("CD56"^"bright"*"CD16"^"-"~"in Total NK"),
       ylab = expression("CD56"^"bright"*"CD16"^"-"~"(%)"),
       file = "CD56bright_percent_TotalNK"),
  
  list(y = "NK_CD56-CD16+",
       title = expression("CD56"^"-"*"CD16"^"+"~"in Total NK"),
       ylab = expression("CD56"^"-"*"CD16"^"+"~"(%)"),
       file = "CD56-CD16+_percent_TotalNK"),
  
  list(y = "NKdim_CD16+",
       title = expression("CD56"^"dim"*"CD16"^"+"~"in Total NK"),
       ylab = expression("CD56"^"dim"*"CD16"^"+"~"(%)"),
       file = "CD56dimCD16+_percent_TotalNK"),
  
  list(y = "Total_NK/NKG2A+",
       title = expression("NKG2A"^"+"~"in Total NK"),
       ylab = expression("NKG2A"^"+"~"(%)"),
       file = "NKG2A+_percent_TotalNK"),
  
  list(y = "Total_NK/CD38+",
       title = expression("CD38"^"+"~"in Total NK"),
       ylab = expression("CD38"^"+"~"(%)"),
       file = "CD38+_percent_TotalNK"),
  
  list(y = "NK_CD56-CD16+/NKG2A+",
       title = expression("CD56"^"-"*"CD16"^"+"~"NKG2A"^"+"~"in Total NK"),
       ylab = expression("CD56"^"-"*"CD16"^"+"~"NKG2A"^"+"~"(%)"),
       file = "CD56-CD16+NKG2A+_percent_TotalNK"),
  
  list(y = "NK_CD56-CD16+/CD38+",
       title = expression("CD56"^"-"*"CD16"^"+"~"CD38"^"+"~"in Total NK"),
       ylab = expression("CD56"^"-"*"CD16"^"+"~"CD38"^"+"~"(%)"),
       file = "CD56-CD16+CD38+_percent_TotalNK"),
  
  list(y = "NKdim_CD16+/NKG2A+",
       title = expression("CD56"^"dim"*"CD16"^"+"~"NKG2A"^"+"~"in Total NK"),
       ylab = expression("CD56"^"dim"*"CD16"^"+"~"NKG2A"^"+"~"(%)"),
       file = "CD56dimCD16+NKG2A+_percent_TotalNK"),
  
  list(y = "NKdim_CD16+/CD38+",
       title = expression("CD56"^"dim"*"CD16"^"+"~"CD38"^"+"~"in Total NK"),
       ylab = expression("CD56"^"dim"*"CD16"^"+"~"CD38"^"+"~"(%)"),
       file = "CD56dimCD16+CD38+_percent_TotalNK"),
  
  list(y = "NKbright/NKG2A+",
       title = expression("CD56"^"bright"*"CD16"^"-"*"NKG2A"^"+"~"in Total NK"),
       ylab = expression("CD56"^"bright"*"CD16"^"-"*"NKG2A"^"+"~"(%)"),
       file = "CD56bright_NKG2A+_percent_TotalNK"),
  
  list(y = "NKbright/CD38+",
       title = expression("CD56"^"bright"*"CD16"^"-"*"CD38"^"+"~"in Total NK"),
       ylab = expression("CD56"^"bright"*"CD16"^"-"*"CD38"^"+"~"(%)"),
       file = "CD56bright_CD38+_percent_TotalNK")
)

# Loop through and generate plots
for (info in plot_totalnk_info) {
  p <- ggplot(nkiller_plot.2, aes(x = Group, y = .data[[info$y]], fill = Group)) +
    geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +
    geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +
    labs(title = info$title,
         x = "HIV Status",
         y = info$ylab) +
    theme_minimal() +
    scale_fill_manual(values = group_colors, labels = group_labels) +
    scale_color_manual(values = group_colors, labels = group_labels) +
    scale_x_discrete(labels = group_labels) +
    stat_compare_means(comparisons = list(c("HEI", "HEU")),
                       method = "wilcox.test",
                       label = "p.format") +
    facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", colour = 'black'),
      axis.title.x = element_text(size = 18, margin = margin(t = 15), colour = 'black'),
      axis.title.y = element_text(size = 18, colour = 'black'),
      axis.text.x = element_text(size = 18, colour = 'black'),
      axis.text.y = element_text(size = 16, colour = 'black'),
      strip.text = element_text(size = 17, colour = 'black'),
      strip.placement = "outside",
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
    )
  
  ggsave(filename = paste0(out.path.2, "HEIvsHEU_", info$file, ".png"),
         plot = p, height = 5, width = 8.5, dpi = 300, bg = "white")
}




################ NKG2A and CD38 for each subset #####################

out.path.3 <- "C:/Users/ammas/Documents/NK_Manuscript/Boxplots/NKG2A_CD38/"

# Label and color settings
group_colors <- c("HEI" = "#fc913f", "HEU" = "#5bbae3")
group_labels <- c("HEI" = "pHIV", "HEU" = "pHEU")

# Subsets and markers to plot
subsets <- c("NKdim_CD16+", "NK_CD56-CD16+", "NKbright")
markers <- c("NKG2A", "CD38")

# Loop over combinations
for (subset in subsets) {
  for (marker in markers) {
    # Column names
    subset_col <- subset
    marker_col <- paste0(subset, "/", marker, "+")
    percent_col <- paste0(subset, "_", marker, "_percent")
    
    # Calculate % marker within subset
    nkiller_plot[[percent_col]] <- (nkiller_plot[[marker_col]] / nkiller_plot[[subset_col]]) * 100
    
    # Plot
    p <- ggplot(nkiller_plot, aes(x = Group, y = .data[[percent_col]], fill = Group)) +
      geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +
      geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +
      labs(title = bquote(.(marker)^"+"~"in"~.(subset)),
           x = "HIV Status",
           y = bquote(.(marker)^"+"~"(%)")) +
      theme_minimal() +
      scale_fill_manual(values = group_colors, labels = group_labels) +
      scale_color_manual(values = group_colors, labels = group_labels) +
      scale_x_discrete(labels = group_labels) +
      stat_compare_means(comparisons = list(c("HEI", "HEU")),
                         method = "wilcox.test",
                         label = "p.format") +
      facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", colour = 'black'),
        axis.title.x = element_text(size = 18, margin = margin(t = 15), colour = 'black'),
        axis.title.y = element_text(size = 18, colour = 'black'),
        axis.text.x = element_text(size = 18, colour = 'black'),
        axis.text.y = element_text(size = 16, colour = 'black'),
        strip.text = element_text(size = 17, colour = 'black'),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        legend.position = "none"
      )
    
    # Save
    ggsave(
      filename = paste0(out.path.3, "HEIvsHEU_", subset, "_", marker, "_percent.png"),
      plot = p, height = 5, width = 8.5, dpi = 300, bg = "white"
    )
  }
}


################################ Correlate Viral Load with Variables ##########

### Select Viral Load Only

nkiller_HEI <- nkiller_plot %>%
  filter(Group == 'HEI')

#### CORR FUNCTIONS ############
#
# Function to calculate correlations
calculate_correlations <- function(data, target_var, cols_to_check) {
  # Create an empty results data frame
  results <- data.frame(
    Variable = character(),
    Correlation = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through the columns to check correlations
  for (col in cols_to_check) {
    var_name <- colnames(data)[col]
    correlation_test <- cor.test(data[[var_name]], data[[target_var]], method = "spearman")
    
    # Store the results
    results <- rbind(results, data.frame(
      Variable = var_name,
      Correlation = correlation_test$estimate,
      P_value = correlation_test$p.value
    ))
  }
  
  # Return the results
  return(results)
}

# Function to plot significant correlations
plot_significant_correlations <- function(correlations, title) {
  # Filter significant correlations
  significant_correlations <- correlations %>%
    filter(P_value < 0.05) %>%
    mutate(
      Direction = ifelse(Correlation > 0, "Positive", "Negative")
    )
  
  # Plot the significant correlations
  ggplot(significant_correlations, aes(x = reorder(Variable, Correlation), y = Correlation, fill = Direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +  # Flip coordinates for better readability
    labs(
      title = title,
      x = "Variable",
      y = "Spearman Correlation"
    ) +
    scale_fill_manual(values = c("Positive" = "#5ab4ac", "Negative" = "#d8b365"), name = "Direction") +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
}




setwd("C:/Users/ammas/Documents/NK_Manuscript/Correlations/TARA/Viral_Load_HEI_Only/All_TP_Dataset")
data_cols <- 7:(ncol(nkiller_HEI)-1)
correlations_VL_Untreated <- calculate_correlations(nkiller_HEI, "Viral Titre", data_cols)

# Step 2: Plot the significant correlations
plot_VL_Untreated <- plot_significant_correlations(
  correlations_VL_Untreated,
  "Significant Correlations with Viral Load - Untreated"
)

# Step 3: Save the plot
ggsave(
  filename = "HEI_Untreated_Viral_Load_Correlations_Combined_Timepoints.png",
  plot = plot_VL_Untreated,
  width = 10,
  height = 8,
  bg='white',
  dpi = 300
)

###### Split by Timepoint

for (timepoint in unique(nkiller_HEI$timepoint.2)) {
  # Filter data for the current timepoint
  timepoint_data <- nkiller_HEI %>%
    filter(timepoint.2 == timepoint)
  
  # Step 2: Calculate correlations for Viral Load
  correlations_VL_timepoint <- calculate_correlations(timepoint_data, "Viral Titre", data_cols)
  
  # Step 3: Plot the significant correlations
  plot_VL_timepoint <- plot_significant_correlations(
    correlations_VL_timepoint,
    paste("Significant Correlations with Viral Load - Untreated - Timepoint", timepoint)
  )
  
  # Step 4: Save the plot
  ggsave(
    filename = paste0("HEI_Untreated_Viral_Load_Correlations_Timepoint_", timepoint, ".png"),
    plot = plot_VL_timepoint,
    width = 10,
    height = 8,
    bg = 'white',
    dpi = 300
  )
}







#################################
### Compare viral Load with variables

nkiller_plot <- nkiller_plot %>%
  mutate(Scaled_Viral_Load = scale(`Viral Titre`))

nkiller_subset_hei <- nkiller_plot %>%
  filter(Group == "HEI")

# Step 3: Create scatter plots of Scaled Viral Load vs. Total NK, split by Timepoint and Group
ggplot(nkiller_subset_hei, aes(x = Scaled_Viral_Load, y = Total_NK, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot with points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line
  labs(title = "Total NK Frequency vs Scaled Viral Load",
       x = "Scaled Viral Load",
       y = "Total NK Frequency (%)") +
  theme_minimal() +
  facet_wrap(~ Age, scales = "free_x", nrow = 1) +  # Split by Timepoint
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom colors for HEI and HEU
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold",colour = 'black'),  # Title formatting
    axis.title.x = element_text(size = 18, margin = margin(t = 15),colour = 'black'),  # Shifts x-axis title down
    axis.title.y = element_text(size = 18,colour = 'black'),  # y-axis title formatting
    axis.text.x = element_text(size = 18,colour = 'black'),  # Larger x-axis text for Group labels
    axis.text.y = element_text(size = 16,colour = 'black'),  # y-axis text size
    strip.text = element_text(size = 16,colour = 'black'),  # Facet label size
    strip.placement = "outside",  # Moves facet label outside
    panel.spacing = unit(1, "lines"),  # Adjusts space between facets
    legend.position = "none"  # Remove legend if not needed
  )

ggplot(nkiller_subset_hei, aes(x = Scaled_Viral_Load, y = `NKdim_CD16+`, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot with points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line
  labs(title = "`NKdim_CD16+` vs Scaled Viral Load",
       x = "Scaled Viral Load",
       y = "NKdim_CD16+ (% of CD45+ Lymphocytes)") +
  theme_minimal() +
  facet_wrap(~ Age, scales = "free_x", nrow = 1) +  # Split by Timepoint
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom colors for HEI and HEU
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold",colour = 'black'),  # Title formatting
    axis.title.x = element_text(size = 18, margin = margin(t = 15),colour = 'black'),  # Shifts x-axis title down
    axis.title.y = element_text(size = 18,colour = 'black'),  # y-axis title formatting
    axis.text.x = element_text(size = 18,colour = 'black'),  # Larger x-axis text for Group labels
    axis.text.y = element_text(size = 16,colour = 'black'),  # y-axis text size
    strip.text = element_text(size = 16,colour = 'black'),  # Facet label size
    strip.placement = "outside",  # Moves facet label outside
    panel.spacing = unit(1, "lines"),  # Adjusts space between facets
    legend.position = "none"  # Remove legend if not needed
  )

ggplot(nkiller_subset_hei, aes(x = Scaled_Viral_Load, y = `NK_CD56-CD16+`, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot with points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line
  labs(title = "`NK_CD56-CD16+` vs Scaled Viral Load",
       x = "Scaled Viral Load",
       y = "NK_CD56-CD16+ (% of CD45+ Lymphocytes)") +
  theme_minimal() +
  facet_wrap(~ Age, scales = "free_x", nrow = 1) +  # Split by Timepoint
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom colors for HEI and HEU
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold",colour = 'black'),  # Title formatting
    axis.title.x = element_text(size = 18, margin = margin(t = 15),colour = 'black'),  # Shifts x-axis title down
    axis.title.y = element_text(size = 18,colour = 'black'),  # y-axis title formatting
    axis.text.x = element_text(size = 18,colour = 'black'),  # Larger x-axis text for Group labels
    axis.text.y = element_text(size = 16,colour = 'black'),  # y-axis text size
    strip.text = element_text(size = 16,colour = 'black'),  # Facet label size
    strip.placement = "outside",  # Moves facet label outside
    panel.spacing = unit(1, "lines"),  # Adjusts space between facets
    legend.position = "none"  # Remove legend if not needed
  )
ggplot(nkiller_subset_hei, aes(x = Scaled_Viral_Load, y = NKbright, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot with points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Linear regression line
  labs(title = "`NK_CD56bright` vs Scaled Viral Load",
       x = "Scaled Viral Load",
       y = "NK_CD56 bright (% of CD45+ Lymphocytes)") +
  theme_minimal() +
  facet_wrap(~ Age, scales = "free_x", nrow = 1) +  # Split by Timepoint
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom colors for HEI and HEU
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold",colour = 'black'),  # Title formatting
    axis.title.x = element_text(size = 18, margin = margin(t = 15),colour = 'black'),  # Shifts x-axis title down
    axis.title.y = element_text(size = 18,colour = 'black'),  # y-axis title formatting
    axis.text.x = element_text(size = 18,colour = 'black'),  # Larger x-axis text for Group labels
    axis.text.y = element_text(size = 16,colour = 'black'),  # y-axis text size
    strip.text = element_text(size = 16,colour = 'black'),  # Facet label size
    strip.placement = "outside",  # Moves facet label outside
    panel.spacing = unit(1, "lines"),  # Adjusts space between facets
    legend.position = "none"  # Remove legend if not needed
  )