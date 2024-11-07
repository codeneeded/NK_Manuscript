library(readxl)
library(dplyr)
library(lme4)
library(ggplot2)
library(ggpubr)
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
library(broom.mixed)
library(extrafont)
loadfonts(device = "win")

############################################ Global Variables #######################################################
setwd("C:/Users/axi313/Documents/NK_Manuscript")
in.path <-"C:/Users/axi313/Documents/NK_Manuscript/Saved_R_Data/"

load(paste0(in.path,"tara_freq_clean.RDS"))
load(paste0(in.path,"florah_freq_clean.RDS"))


####################### Paired Plots ##########################################


######## Paired Plots #####
# Plot the updated data

setwd("C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/Paired_Plots")

### Cleaning 

### Tara
tara_Freq$`Specific Killing` <- as.numeric(tara_Freq$`Specific Killing`)
tara_Freq_plot <- tara_Freq %>% drop_na(14)
tara_Freq_plot_filtered <- tara_Freq_plot %>%
  filter(!Treatment %in% "untreated", HIV != "HUU")
tara_Freq_plot_filtered <- tara_Freq_plot_filtered %>%
  mutate(HIV = recode(HIV, "positive" = "HEI", "negative" = "HEU"))
# Ensure Timepoint is a factor and reorder it so "Entry" comes before "12"
tara_Freq_plot_filtered <- tara_Freq_plot_filtered %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Entry", "12")))
tara_Freq_plot_filtered_2 <- tara_Freq_plot_filtered %>%
  filter(!Treatment %in% c("untreated","CEM+IL15", "HUT78+IL15", "K562+IL15"))

### Florah
florah_Freq$`Specific Killing` <- as.numeric(florah_Freq$`Specific Killing`)
florah_Freq_plot <- florah_Freq %>% drop_na(14)
florah_Freq_plot_filtered <- florah_Freq_plot %>%
  filter(!Treatment %in% "untreated", HIV != "HUU")
florah_Freq_plot_filtered <- florah_Freq_plot_filtered %>%
  mutate(HIV = recode(HIV, "positive" = "PWHIV", "negative" = "PWoHIV"))
# Ensure Timepoint is a factor and reorder it so "Entry" comes before "12"
florah_Freq_plot_filtered <- florah_Freq_plot_filtered %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Entry", "12")))
florah_Freq_plot_filtered_2 <- florah_Freq_plot_filtered %>%
  filter(Treatment %in% c("CEM", "HUT78", "K562"))

##### Specific Killing for Timepoint (HEI vs HEU) ############

##### Function for creation of paired plot with p-values
# Define the function
plot_with_p_values <- function(data, y_var) {
  
  # Create subsets for each timepoint
  entry_data <- data %>%
    filter(Timepoint == "Entry") %>%
    select(HIV, Treatment, PID, {{ y_var }}) %>%
    arrange(HIV, Treatment, PID)
  
  month12_data <- data %>%
    filter(Timepoint == "12") %>%
    select(HIV, Treatment, PID, {{ y_var }}) %>%
    arrange(HIV, Treatment, PID)
  
  # Join the datasets to create paired data
  paired_data <- entry_data %>%
    inner_join(month12_data, by = c("HIV", "Treatment", "PID"), suffix = c("_Entry", "_12Months"))
  
  # Calculate p-values for each Treatment and HIV group
  p_values <- paired_data %>%
    group_by(HIV, Treatment) %>%
    summarise(
      p_value = wilcox.test(get(paste0(y_var, "_Entry")), get(paste0(y_var, "_12Months")), paired = TRUE)$p.value
    )
  
  # Create the plot with p-values
  ggplot(data, aes(x = Timepoint, y = !!sym(y_var), group = PID, color = HIV)) +
    geom_line(aes(linetype = HIV), size = 1) +  # Set line size
    geom_point(size = 2) +  # Set point size
    facet_grid(HIV ~ Treatment, scales = "free_y") +  # Splitting by both HIV and Treatment
    theme_minimal(base_size = 15) +  # Increase base font size for readability
    labs(
      title = paste(y_var, "Across Treatments, Timepoints, and HIV Status"),
      y = paste(y_var, "(%)"), 
      x = "Timepoint"
    ) +
    scale_color_manual(values = c("HEI" = "red", "HEU" = "blue")) +  # Custom colors for HEI and HEU
    scale_linetype_manual(values = c("HEI" = "solid", "HEU" = "solid")) +  # Custom linetypes for HEI and HEU
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      strip.text = element_text(size = 14),  # Adjust facet label size
      legend.title = element_blank(),  # Remove legend title
      legend.position = "top"  # Move the legend to the top
    ) +
    # Add p-values as text annotations
    geom_text(data = p_values, inherit.aes = FALSE, aes(
      x = 1.5,  # Position text between the two timepoints
      y = max(data[[y_var]], na.rm = TRUE) * 1.1,  # Place slightly above max y-value
      label = paste0("p = ", signif(p_value, 3))
    ), color = "black", size = 4)
}

setwd("C:/Users/axi313/Documents/NK_Manuscript/Paired_Plots/TARA")

plot_with_p_values(tara_Freq_plot_filtered_2, "Specific Killing")
ggsave("Specific_Killing_Paired_Plot.png", width = 10, height = 8, dpi = 300,bg='white')

plot_with_p_values(tara_Freq_plot_filtered_2, "IFNy")
ggsave("IFNy_Paired_Plot.png", width = 10, height = 8, dpi = 300,bg='white')

plot_with_p_values(tara_Freq_plot_filtered_2, "CD56dimCD16+/MIP-1B")
ggsave("CD56dimCD16+_MIP-1B_Paired_Plot.png", width = 10, height = 8, dpi = 300,bg='white')

plot_with_p_values(tara_Freq_plot_filtered_2, "CD56dimCD16+/CD107a")
ggsave("CD56dimCD16+_CD107a_Paired_Plot.png", width = 10, height = 8, dpi = 300,bg='white')

plot_with_p_values(tara_Freq_plot_filtered_2, "MIP-1B")
ggsave("MIP-1B_Paired_Plot.png", width = 10, height = 8, dpi = 300,bg='white')

plot_with_p_values(tara_Freq_plot_filtered_2, "CD107a")
ggsave("CD107a_Paired_Plot.png", width = 10, height = 8, dpi = 300,bg='white')

######################## Specific Killing Split by HIV Status ########################

### TARA
# Define the function with specified colors
plot_specific_killing_comparison <- function(data, y_var, title = NULL) {
  
  # Set default title if none provided
  if (is.null(title)) {
    title <- paste("Comparison of", y_var, "between HEI and HEU")
  }
  
  # Create the plot
  ggplot(data, aes(x = HIV, y = !!sym(y_var), fill = HIV)) +
    geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = NA) +  # Boxplot for each group
    geom_jitter(width = 0.2, size = 2, alpha = 0.6, aes(color = HIV)) +  # Jitter for individual points
    facet_grid(~Treatment, scales = "free_y") +  # Split by Timepoint and Treatment
    theme_minimal(base_size = 15) +
    labs(
      title = title,
      y = paste(y_var, "(%)"), 
      x = "HIV Status"
    ) +
    scale_fill_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom fill colors for HEI and HEU
    scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom jitter colors
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Title formatting
      strip.text = element_text(size = 14),  # Facet label size
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.position = "none"  # Remove legend if not needed
    ) +
    # Add non-paired Wilcoxon test p-values between HEI and HEU within each facet
    stat_compare_means(comparisons = list(c("HEI", "HEU")), 
                       method = "wilcox.test", 
                       label = "p.format", 
                       label.y = max(data[[y_var]], na.rm = TRUE) * 1.1)
}
setwd("C:/Users/axi313/Documents/NK_Manuscript/Boxplots/Specific_Killing")

plot_specific_killing_comparison(tara_Freq_plot_filtered_2, "Specific Killing", title = "Specific Killing by HIV Status (TARA)")
ggsave("TARA_Specific_Killing_Barplot_by_Treatment.png", width = 10, height = 8, dpi = 300,bg='white')

plot_specific_killing_comparison(florah_Freq_plot_filtered_2, "Specific Killing", title = "Specific Killing by HIV Status (FLORAH)")
ggsave("FLORAH_Specific_Killing_Barplot_by_Treatment.png", width = 10, height = 8, dpi = 300,bg='white')


#### Plot the relationship between viral load.
setwd("C:/Users/axi313/Documents/NK_Manuscript/Correlations/TARA")

# Scale viral load for the HEI subset
data_HEI <- tara_Freq_plot_filtered_2 %>%
  filter(HIV == "HEI") %>%
  mutate(`viral load` = scale(`viral load`))  # Scale the viral load column

# Calculate p-values for each Treatment and Timepoint combination using Spearman correlation
p_values <- data_HEI %>%
  group_by(Treatment, Timepoint) %>%
  summarise(
    p_value = cor.test(`viral load`, `Specific Killing`, method = "spearman")$p.value
  )

# Create the scatter plot with faceting and p-values
ggplot(data_HEI, aes(x = `viral load`, y = `Specific Killing`)) +
  geom_point(color = "#D55E00", size = 3, alpha = 0.6) +  # Points for each observation
  geom_smooth(method = "lm", color = "black", se = TRUE) +  # Linear trend line with confidence interval
  facet_grid(Timepoint ~ Treatment, scales = "free_y") +  # Facet by Timepoint and Treatment
  theme_minimal(base_size = 15) +
  labs(
    title = "Relationship between Scaled Viral Load and Specific Killing within HEI",
    x = "Scaled Viral Load",
    y = "Specific Killing (%)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12)  # Adjust facet label size
  ) +
  # Add p-values as text annotations
  geom_text(data = p_values, inherit.aes = FALSE, aes(
    x = mean(range(data_HEI$`viral load`, na.rm = TRUE)), # Center p-value above the graph
    y = Inf,  # Position text at the top of each facet
    label = paste0("p = ", signif(p_value, 3))
  ), hjust = 0.5, vjust = 1.5, size = 4, color = "black")
ggsave("TARA_Specific_Killing_vs_Viral_Load_Scatterplot_spearmens.png", width = 10, height = 8, dpi = 300,bg='white')

### Corrplots

######### Plot Specific Killing for each Group #############

### TARA

####### Paired Correlation for Specific Killing
hut78_data_TARA <- tara_Freq_plot_filtered %>%
  filter(Treatment == "HUT78")
k562_data_TARA <- tara_Freq_plot_filtered %>%
  filter(Treatment == "K562")

# Define a function to calculate correlations
calculate_correlations <- function(data, target_var = "Specific Killing") {
  results <- data.frame(Variable = character(), Correlation = numeric(), P_value = numeric(), stringsAsFactors = FALSE)
  
  # List of columns to correlate with Specific Killing (viral load and columns after 19)
  cols_to_check <- c(8, 20:ncol(data))
  
  # Loop through each column
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
  
  return(results)
}

# Modify the function to calculate correlations with viral load
calculate_correlations_with_viral_load <- function(data, target_var = "viral load") {
  results <- data.frame(Variable = character(), Correlation = numeric(), P_value = numeric(), stringsAsFactors = FALSE)
  
  # List of columns to correlate with viral load (all columns after 19)
  cols_to_check <- 20:ncol(data)
  
  # Loop through each column
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
  
  return(results)
}

# Define the plotting function
plot_correlations <- function(correlations, title) {
  ggplot(correlations, aes(x = reorder(Variable, Correlation), y = Correlation, fill = P_value < 0.05)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +  # Flip coordinates for better readability
    labs(
      title = title,
      x = "Variable",
      y = "Spearman Correlation with Specific Killing"
    ) +
    scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "lightgrey"), name = "Significance", labels = c("p < 0.05", "p >= 0.05")) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
}


# Run the correlation function on HUT78 and K562 datasets
hut78_correlations <- calculate_correlations(hut78_data_TARA)
k562_correlations <- calculate_correlations(k562_data_TARA)


plot_correlations(hut78_correlations, "Correlation of Variables with Specific Killing (HUT78)")
ggsave("TARA_HUT78_Corr_barplot_vs_Specific_Killing.png", width = 10, height = 15, dpi = 300,bg='white')

plot_correlations(k562_correlations, "Correlation of Variables with Specific Killing (K562)")
ggsave("TARA_K562_Corr_barplot_vs_Specific_Killing.png", width = 10, height = 15, dpi = 300,bg='white')

# Run the correlation function on HUT78 and K562 datasets
hut78_correlations_with_viral_load <- calculate_correlations_with_viral_load(hut78_data_TARA)
k562_correlations_with_viral_load <- calculate_correlations_with_viral_load(k562_data_TARA)

# Define the plotting function for viral load correlations
plot_correlations <- function(correlations, title) {
  ggplot(correlations, aes(x = reorder(Variable, Correlation), y = Correlation, fill = P_value < 0.05)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +  # Flip coordinates for better readability
    labs(
      title = title,
      x = "Variable",
      y = "Spearman Correlation with Viral Load"
    ) +
    scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "lightgrey"), name = "Significance", labels = c("p < 0.05", "p >= 0.05")) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
}

# Example usage for HUT78 and K562 data with viral load
plot_correlations(hut78_correlations_with_viral_load, "Correlation of Variables with Viral Load (HUT78)")
ggsave("TARA_HUT78_Corr_barplot_vs_VL.png", width = 10, height = 15, dpi = 300,bg='white')

plot_correlations(k562_correlations_with_viral_load, "Correlation of Variables with Viral Load (K562)")
ggsave("TARA_K562_Corr_barplot_vs_VL.png", width = 10, height = 15, dpi = 300,bg='white')

### FLORAH

setwd("C:/Users/axi313/Documents/NK_Manuscript/Correlations/FLORAH")

####### Paired Correlation for Specific Killing
hut78_data_f <- florah_Freq_plot_filtered %>%
  filter(Treatment == "HUT78")
k562_data_f <- florah_Freq_plot_filtered %>%
  filter(Treatment == "K562")

# Function to calculate correlations with Specific Killing
calculate_correlations_with_specific_killing <- function(data, target_var = "Specific Killing") {
  results <- data.frame(Variable = character(), Correlation = numeric(), P_value = numeric(), stringsAsFactors = FALSE)
  
  # List of columns to correlate with Specific Killing (all columns after 19)
  cols_to_check <- 20:ncol(data)
  
  # Loop through each column
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
  
  return(results)
}

# Run the correlation function on HUT78 and K562 datasets in florah
hut78_correlations_with_specific_killing <- calculate_correlations_with_specific_killing(hut78_data_f)
k562_correlations_with_specific_killing <- calculate_correlations_with_specific_killing(k562_data_f)

# Plot for HUT78 and K562 correlations with Specific Killing in florah
plot_correlations(hut78_correlations_with_specific_killing, "Correlation of Variables with Specific Killing (HUT78)")
ggsave("FLORAH_HUT78_Corr_barplot_vs_Specific_Killing.png", width = 10, height = 17, dpi = 300,bg='white')

plot_correlations(k562_correlations_with_specific_killing, "Correlation of Variables with Specific Killing (K562)")
ggsave("FLORAH_K562_Corr_barplot_vs_Specific_Killing.png", width = 10, height = 17, dpi = 300,bg='white')

