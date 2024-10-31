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
setwd("C:/Users/axi313/Documents/NK_Manuscript")
in.path <-"C:/Users/axi313/Documents/NK_Manuscript/Saved_R_Data/"

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
setwd("C:/Users/axi313/Documents/NK_Manuscript/Boxplots/HEIvsHEU/TARA")

# Step 1: Get the list of columns to plot (all columns after column 20)
columns_to_plot <- colnames(tara_untreated)[21:ncol(tara_untreated)]



# Step 2: Loop through each column and create/save the plot
for (column_name in columns_to_plot) {
  # Step 3: Create the plot for the current column
  plot <- ggplot(tara_untreated %>% filter(!is.na(!!sym(column_name))), 
                 aes(x = HIV, y = !!sym(column_name), fill = HIV)) +
    geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +  # Boxplot with adjusted width
    geom_jitter(aes(color = HIV), width = 0.2, size = 2, alpha = 0.7) +  # Jitter for individual points
    facet_wrap(~ Timepoint, scales = "free_x", nrow = 1) +  # Split by Timepoint in one row
    stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                       method = "wilcox.test", 
                       label = 'p.format') +
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

### FLORAH

#plot wilcoxins test of all columns between HEi and HEU
setwd("C:/Users/axi313/Documents/NK_Manuscript/Boxplots/HEIvsHEU/Florah")

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
setwd("C:/Users/axi313/Documents/NK_Manuscript")

nk_data <- read_csv("NK_Cell_data_cleaned_normalized_to_CD45_Lesley Reduced 061124.csv")


demographics <- read_csv("Viral_Titres.csv")
out.path.2 <- "C:/Users/axi313/Documents/NK_Manuscript/Boxplots/Timepoint_HIV/"

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

### TOTAL NK ###

ggplot(nkiller_plot, aes(x = Group, y = Total_NK, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.4) +  # Boxplot with adjusted width
  geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +  # Jitter for individual points
  labs(title = expression("Total NK in CD45"^"+"~"Lymphocytes") ,
       x = "HIV Status",
       y = "Total NK (%)") +
  theme_minimal() +
  scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Fill colors for boxplot
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Colors for jitter points
  stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                     method = "wilcox.test", 
                     label = 'p.format') +  # Adjust y position for significance stars
  facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +  # Arrange all facets in one row
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
ggsave(paste0(out.path.2,"HEIvsHEU_Total_NK_percent.png"),bg='white', height = 5, width = 8.5)

### CD56bright CD16neg ###
# Create the plot with all facets in one row and custom order
ggplot(nkiller_plot, aes(x = Group, y = NKbright, fill = Group)) +
  geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +  # Boxplot with adjusted width
  geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +  # Jitter for individual points
  labs(title = expression("CD56"^"bright"*"CD16"^"-"~"in CD45"^"+"~"Lymphocytes"),  # Title with superscripts
       x = "HIV Status",
       y = expression("CD56"^"bright"*"CD16"^"-"~"(%)")) +  # Y-axis label with superscripts
  theme_minimal() +
  scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Fill colors for boxplot
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Colors for jitter points
  stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                     method = "wilcox.test", 
                     label = 'p.format') +  # Adjust y position for significance stars
  facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +  # Arrange all facets in one row
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

ggsave(paste0(out.path.2,"HEIvsHEU_CD56bright_percent_CD45L.png"),bg='white', height = 5, width = 8.5)

### CD56-CD16+ ###
# Create the plot with all facets in one row and custom order
ggplot(nkiller_plot, aes(x = Group, y = `NK_CD56-CD16+`, fill = Group)) +
  geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +  # Boxplot with adjusted width
  geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +  # Jitter for individual points
  labs(title = expression("CD56"^"-"*"CD16"^"+"~"in CD45"^"+"~"Lymphocytes"),  # Title with superscripts
       x = "HIV Status",
       y = expression("CD56"^"-"*"CD16"^"+"~"(%)")) +  # Y-axis label with superscripts
  theme_minimal() +
  scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Fill colors for boxplot
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Colors for jitter points
  stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                     method = "wilcox.test", 
                     label = 'p.format') +  # Adjust y position for significance stars
  facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +  # Arrange all facets in one row
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

ggsave(paste0(out.path.2,"HEIvsHEU_CD56-CD16+_percent_CD45L.png"),bg='white', height = 5, width = 8.5)

### CD56dimCD16+ ###

# Create the plot with all facets in one row and custom order
ggplot(nkiller_plot, aes(x = Group, y = `NKdim_CD16+`, fill = Group)) +
  geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +  # Boxplot with adjusted width
  geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +  # Jitter for individual points
  labs(title = expression("CD56"^"dim"*"CD16"^"+"~"in CD45"^"+"~"Lymphocytes"),  # Title with superscripts
       x = "HIV Status",
       y = expression("CD56"^"dim"*"CD16"^"+"~"(%)")) +  # Y-axis label with superscripts
  theme_minimal() +
  scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Fill colors for boxplot
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Colors for jitter points
  stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                     method = "wilcox.test", 
                     label = 'p.format') +  # Adjust y position for significance stars
  facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +  # Arrange all facets in one row
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

ggsave(paste0(out.path.2,"HEIvsHEU_CD56dimCD16+_percent_CD45L.png"),bg='white', height = 5, width = 8.5)

###
nkiller_plot.2 <- nkiller_plot
nkiller_plot.2$Total_NK <- (nkiller_plot.2$Total_NK/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts
nkiller_plot.2$`NKdim_CD16+` <- ((nkiller_plot.2$`NKdim_CD16+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$`NK_CD56-CD16+` <- ((nkiller_plot.2$`NK_CD56-CD16+`/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)
nkiller_plot.2$NKbright<- ((nkiller_plot.2$NKbright/100)*nkiller_plot.2$NKpanel_CD45_Raw_Counts)*(100/nkiller_plot.2$Total_NK)

###### Boxplot To Total NK #######

### CD56bright CD16neg ###
# Create the plot with all facets in one row and custom order
ggplot(nkiller_plot.2, aes(x = Group, y = NKbright, fill = Group)) +
  geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +  # Boxplot with adjusted width
  geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +  # Jitter for individual points
  labs(title = expression("CD56"^"bright"*"CD16"^"-"~"in Total NK"),  # Title with superscripts
       x = "HIV Status",
       y = expression("CD56"^"bright"*"CD16"^"-"~"(%)")) +  # Y-axis label with superscripts
  theme_minimal() +
  scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Fill colors for boxplot
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Colors for jitter points
  stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                     method = "wilcox.test", 
                     label = 'p.format') +  # Adjust y position for significance stars
  facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +  # Arrange all facets in one row
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

ggsave(paste0(out.path.2,"HEIvsHEU_CD56bright_percent_TotalNK.png"),bg='white', height = 5, width = 8.5)

### CD56-CD16+ ###
# Create the plot with all facets in one row and custom order
ggplot(nkiller_plot.2, aes(x = Group, y = `NK_CD56-CD16+`, fill = Group)) +
  geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +  # Boxplot with adjusted width
  geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +  # Jitter for individual points
  labs(title = expression("CD56"^"-"*"CD16"^"+"~"in Total NK"),  # Title with superscripts
       x = "HIV Status",
       y = expression("CD56"^"-"*"CD16"^"+"~"(%)")) +  # Y-axis label with superscripts
  theme_minimal() +
  scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Fill colors for boxplot
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Colors for jitter points
  stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                     method = "wilcox.test", 
                     label = 'p.format') +  # Adjust y position for significance stars
  facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +  # Arrange all facets in one row
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

ggsave(paste0(out.path.2,"HEIvsHEU_CD56-CD16+_percent_TotalNK.png"),bg='white', height = 5, width = 8.5)

### CD56dimCD16+ ###

# Create the plot with all facets in one row and custom order
ggplot(nkiller_plot.2, aes(x = Group, y = `NKdim_CD16+`, fill = Group)) +
  geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +  # Boxplot with adjusted width
  geom_jitter(aes(color = Group), width = 0.2, size = 3, alpha = 0.8) +  # Jitter for individual points
  labs(title = expression("CD56"^"dim"*"CD16"^"+"~"in Total NK"),  # Title with superscripts
       x = "HIV Status",
       y = expression("CD56"^"dim"*"CD16"^"+"~"(%)")) +  # Y-axis label with superscripts
  theme_minimal() +
  scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Fill colors for boxplot
  scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Colors for jitter points
  stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                     method = "wilcox.test", 
                     label = 'p.format') +
  facet_wrap(~ timepoint.2, scales = "free_x", nrow = 1) +  # Arrange all facets in one row
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

ggsave(paste0(out.path.2,"HEIvsHEU_CD56dimCD16+_percent_TotalNK.png"),bg='white', height = 5, width = 8.5)


################################




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