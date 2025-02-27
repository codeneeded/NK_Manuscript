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
library(broom.mixed)
library(extrafont)
library(ggpubr)
library(grid)
loadfonts(device = "win")

############################################ Global Variables #######################################################
setwd("C:/Users/ammas/Documents/NK_Manuscript")
in.path <-"C:/Users/ammas/Documents/NK_Manuscript/Saved_R_Data/"

load(paste0(in.path,"tara_freq_clean.RDS"))
load(paste0(in.path,"florah_freq_clean.RDS"))

############################################## DATA CLEANING ##########################################################

#### TARA #### 

# Filter unneeded columna
tara_Freq$Group <- NULL
tara_Freq$`age group` <- NULL
tara_Freq$`age (yrs)` <- NULL
tara_Freq$Target <- NULL
tara_Freq$`Target/Dead`<- NULL
tara_Freq$Timepoint <- as.factor(tara_Freq$Timepoint)
tara_Freq$Timepoint <- factor(tara_Freq$Timepoint, levels = rev(levels(tara_Freq$Timepoint)))
tara_Freq <- tara_Freq%>%
  mutate(Timepoint = factor(Timepoint, levels = c("Entry", "12")))

tara_Freq <- tara_Freq %>%
  filter(HIV != "HUU")


TARA_HUT78_Freq <- tara_Freq %>%
  filter(Treatment == "HUT78")

TARA_K562_Freq <- tara_Freq %>%
  filter(Treatment == "K562")

# Clean HUT78

TARA_HUT78_Freq <- TARA_HUT78_Freq[!is.na(TARA_HUT78_Freq$`Specific Killing`), ]
TARA_HUT78_Freq <- droplevels(TARA_HUT78_Freq)
TARA_HUT78_Freq$`Specific Killing`<- as.numeric(TARA_HUT78_Freq$`Specific Killing`)
TARA_HUT78_Freq$`Specific Killing`<- round(TARA_HUT78_Freq$`Specific Killing`, 2)


TARA_HUT78_Freq$`viral load` <- scale(TARA_HUT78_Freq$`viral load`) # Standardize the 'viral load' variable because of high values

# Clean K562

TARA_K562_Freq <- TARA_K562_Freq[!is.na(TARA_K562_Freq$`Specific Killing`), ]
TARA_K562_Freq <- droplevels(TARA_K562_Freq)
TARA_K562_Freq$`Specific Killing`<- as.numeric(TARA_K562_Freq$`Specific Killing`)
TARA_K562_Freq$`Specific Killing`<- round(TARA_K562_Freq$`Specific Killing`, 2)


TARA_K562_Freq$`viral load` <- scale(TARA_K562_Freq$`viral load`) # Standardize the 'viral load' variable because of high values

# List of columns representing flow populations, starting from column 13 onwards
flow_population_columns <- colnames(TARA_HUT78_Freq)[15:ncol(TARA_HUT78_Freq)] 

### Tara Untreated

TARA_Untreated_Freq <- tara_Freq %>%
  filter(Treatment == "untreated")
TARA_Untreated_Freq$`Specific Killing` <- NULL

TARA_Untreated_Freq_HUT78 <- TARA_Untreated_Freq %>%
  left_join(TARA_HUT78_Freq %>% select(PID, Timepoint, `Specific Killing`), by = c("PID", "Timepoint"))
TARA_Untreated_Freq_K562 <- TARA_Untreated_Freq %>%
  left_join(TARA_K562_Freq %>% select(PID, Timepoint, `Specific Killing`), by = c("PID", "Timepoint"))


TARA_Untreated_Freq_HUT78 <- TARA_Untreated_Freq_HUT78[!is.na(TARA_Untreated_Freq_HUT78$`Specific Killing`), ]
TARA_Untreated_Freq_K562 <- TARA_Untreated_Freq_K562[!is.na(TARA_Untreated_Freq_K562$`Specific Killing`), ]

### FLORAH ###

# Filter unneeded columna
florah_Freq$Group <- NULL
florah_Freq$`age (yrs)` <- NULL
florah_Freq$Target <- NULL
florah_Freq$`Target/Dead`<- NULL
florah_Freq$Timepoint <- NULL
florah_Freq$`viral load` <- NULL

florah_Freq <- florah_Freq %>%
  filter(HIV != "HUU")

florah_HUT78_Freq <- florah_Freq %>%
  filter(Treatment == "HUT78")

florah_K562_Freq <- florah_Freq %>%
  filter(Treatment == "K562")

# Clean HUT78

florah_HUT78_Freq <- florah_HUT78_Freq[!is.na(florah_HUT78_Freq$`Specific Killing`), ]
florah_HUT78_Freq <- droplevels(florah_HUT78_Freq)
florah_HUT78_Freq$`Specific Killing`<- as.numeric(florah_HUT78_Freq$`Specific Killing`)
florah_HUT78_Freq$`Specific Killing`<- round(florah_HUT78_Freq$`Specific Killing`, 2)

# Clean K562

florah_K562_Freq <- florah_K562_Freq[!is.na(florah_K562_Freq$`Specific Killing`), ]
florah_K562_Freq <- droplevels(florah_K562_Freq)
florah_K562_Freq$`Specific Killing`<- as.numeric(florah_K562_Freq$`Specific Killing`)
florah_K562_Freq$`Specific Killing`<- round(florah_K562_Freq$`Specific Killing`, 2)


### Florah Untreated

florah_Untreated_Freq <- florah_Freq %>%
  filter(Treatment == "untreated")
florah_Untreated_Freq$`Specific Killing` <- NULL

florah_Untreated_Freq_HUT78 <- florah_Untreated_Freq %>%
  left_join(florah_HUT78_Freq %>% select(PID, `Specific Killing`), by = "PID")
florah_Untreated_Freq_K562 <- florah_Untreated_Freq %>%
  left_join(florah_K562_Freq %>% select(PID, `Specific Killing`), by = "PID")


florah_Untreated_Freq_HUT78 <- florah_Untreated_Freq_HUT78[!is.na(florah_Untreated_Freq_HUT78$`Specific Killing`), ]
florah_Untreated_Freq_K562 <- florah_Untreated_Freq_K562[!is.na(florah_Untreated_Freq_K562$`Specific Killing`), ]



###### Functions ########

###### NK KILLING Analysis #############

# Function to fit models and return a list of models
fit_models_list <- function(data, fixed_part_of_formula, flow_population_columns) {
  models_list <- list()
  
  # Loop through each flow column
  for (subset_name in flow_population_columns) {
    # Escape special characters in the column name
    escaped_subset_name <- paste0("`", gsub("/", "\\/", subset_name), "`")
    
    # Construct the formula string
    formula_str <- paste(fixed_part_of_formula, "+", escaped_subset_name)
    
    # Convert to formula and fit the model
    current_formula <- as.formula(formula_str)
    model <- tryCatch({
      lmer(current_formula, data = data)
    }, error = function(e) {
      cat("Error in fitting model for subset:", subset_name, "\nError message:", e$message, "\n")
      return(NULL)
    })
    
    # Store the model if successfully fitted
    if (!is.null(model)) {
      models_list[[subset_name]] <- model
    }
  }
  
  return(models_list)
}
#
# Function to fit linear models and return a list of models
fit_linear_models_list <- function(data, fixed_part_of_formula, flow_population_columns) {
  models_list <- list()
  
  # Loop through each flow column
  for (subset_name in flow_population_columns) {
    # Escape special characters in the column name
    escaped_subset_name <- paste0("`", gsub("/", "\\/", subset_name), "`")
    
    # Construct the formula string
    formula_str <- paste(fixed_part_of_formula, "+", escaped_subset_name)
    
    # Convert to formula and fit the model
    current_formula <- as.formula(formula_str)
    model <- tryCatch({
      lm(current_formula, data = data)  # Use lm instead of lmer
    }, error = function(e) {
      cat("Error in fitting model for subset:", subset_name, "\nError message:", e$message, "\n")
      return(NULL)
    })
    
    # Store the model if successfully fitted
    if (!is.null(model)) {
      models_list[[subset_name]] <- model
    }
  }
  
  return(models_list)
}


# Function to generate a summary table of significant effects from a list of models, excluding intercepts
summarize_models <- function(models_list, p_value_threshold = 0.05) {
  results_df <- data.frame(Subset = character(), Effect = character(), Estimate = numeric(), 
                           Std.Error = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each model and extract significant effects
  for (subset_name in names(models_list)) {
    model <- models_list[[subset_name]]
    summary_model <- summary(model)
    df <- as.data.frame(summary_model$coefficients)
    df$Subset <- subset_name
    df$Effect <- rownames(df)
    rownames(df) <- NULL
    results_df <- rbind(results_df, df)
  }
  
  # Rename columns to make them consistent and filter out intercepts
  results_df <- results_df %>%
    rename(Estimate = Estimate, Std.Error = `Std. Error`, P.Value = `Pr(>|t|)`) %>%
    select(Subset, Effect, Estimate, Std.Error, P.Value) %>%
    filter(Effect != "(Intercept)")  # Remove intercepts
  
  # Filter for significant effects only
  significant_results <- results_df %>%
    filter(P.Value < p_value_threshold)
  
  return(significant_results)
}

# Function to generate a summary table of significant effects from a list of linear models, excluding intercepts
summarize_linear_models <- function(models_list, p_value_threshold = 0.05) {
  results_df <- data.frame(Subset = character(), Effect = character(), Estimate = numeric(), 
                           Std.Error = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each model and extract significant effects
  for (subset_name in names(models_list)) {
    model <- models_list[[subset_name]]
    summary_model <- summary(model)
    df <- as.data.frame(summary_model$coefficients)
    df$Subset <- subset_name
    df$Effect <- rownames(df)
    rownames(df) <- NULL
    results_df <- rbind(results_df, df)
  }
  
  # Rename columns to make them consistent and filter out intercepts
  results_df <- results_df %>%
    rename(Estimate = Estimate, Std.Error = `Std. Error`, P.Value = `Pr(>|t|)`) %>%
    select(Subset, Effect, Estimate, Std.Error, P.Value) %>%
    filter(Effect != "(Intercept)")  # Remove intercepts
  
  # Filter for significant effects only
  significant_results <- results_df %>%
    filter(P.Value < p_value_threshold)
  
  return(significant_results)
}

# Function to get the full summary table of all models (if desired)
get_all_model_summaries <- function(models_list) {
  all_summaries <- data.frame(Subset = character(), Effect = character(), Estimate = numeric(), 
                              Std.Error = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)
  
  for (subset_name in names(models_list)) {
    model <- models_list[[subset_name]]
    if (!is.null(model)) {
      df <- as.data.frame(summary(model)$coefficients)
      df$Subset <- subset_name
      df$Effect <- rownames(df)
      rownames(df) <- NULL
      all_summaries <- rbind(all_summaries, df)
    }
  }
  
  all_summaries <- all_summaries %>%
    rename(Estimate = Estimate, Std.Error = `Std. Error`, P.Value = `Pr(>|t|)`) %>%
    select(Subset, Effect, Estimate, Std.Error, P.Value)
  
  return(all_summaries)
}
##########
plot_significant_relationship <- function(data, x_var, y_var = "Specific Killing", facet_var = "Timepoint") {
  
  # Check if facet_var exists in the data; if not, stop with an error message
  if (!(facet_var %in% colnames(data))) {
    stop(paste("Column", facet_var, "is not found in the dataset. Please check your input."))
  }
  
  # Calculate the global Pearson correlation for the entire dataset before faceting
  pearson_correlation <- cor(data[[x_var]], data[[y_var]], method = "pearson")
  pearson_p_value <- cor.test(data[[x_var]], data[[y_var]], method = "pearson")$p.value
  
  # Create a common label for the global Pearson correlation
  global_correlation_label <- paste0(
    "Global Pearson r = ", round(pearson_correlation, 2), ", p = ", signif(pearson_p_value, 3)
  )
  
  # Calculate individual Pearson correlation for each facet (timepoint)
  correlation_results <- data %>%
    group_by_at(facet_var) %>%
    summarize(
      pearson_correlation = cor(.data[[x_var]], .data[[y_var]], method = "pearson"),
      pearson_p_value = cor.test(.data[[x_var]], .data[[y_var]], method = "pearson")$p.value
    ) %>%
    mutate(
      label = paste0(
        "Pearson r = ", round(pearson_correlation, 2), ", p = ", signif(pearson_p_value, 3)
      )
    )
  
  # Calculate a common y_max for all facets based on the highest value across all facets
  common_y_max <- max(data[[y_var]], na.rm = TRUE) * 1.2  # 20% above the maximum y-value for more space
  
  # Create the plot with the global correlation in the subtitle
  p <- ggplot(data, aes_string(x = paste0("`", x_var, "`"), y = paste0("`", y_var, "`"), group = facet_var)) +
    geom_point(aes(color = HIV), size = 4, alpha = 0.8) +  # Color points by HIV status
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", size = 1.2) +  # Regression line
    labs(
      title = paste(y_var, "Vs", x_var),
      subtitle = paste("Faceted by", facet_var, "\n", global_correlation_label),  # Add global correlation to subtitle
      x = x_var,
      y = y_var
    ) +
    scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom colors for HIV status
    facet_wrap(as.formula(paste("~", facet_var))) +  # Facet by facet_var
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),  # Increase text size for HEI and HEU
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#f7f7f7", color = NA),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),  # Ensure the subtitle is centered
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#e0e0e0", color = NA),
      strip.text = element_text(face = "bold")
    ) +
    # Add individual Pearson correlation annotations for each facet (inside facet)
    geom_text(
      data = correlation_results,
      aes(
        x = Inf, y = common_y_max, 
        label = label
      ),
      hjust = 1.5, color = "black", size = 3.5, fontface = "bold",  # Center aligned with hjust = 0.5
      inherit.aes = FALSE
    )
  
  return(p)
  
}

# Display the correlation results per timepoint
print(correlation_results)

sig_relationships_florah <- function(data, x_var, y_var = "Specific Killing") {
  
  # Calculate a common y_max for all data based on the highest value across all facets
  common_y_max <- max(data[[y_var]], na.rm = TRUE) * 1.2  # 20% above the maximum y-value for more space
  
  # Calculate both Pearson and Spearman correlations for the entire dataset
  correlation_results <- data %>%
    summarize(
      pearson_correlation = cor(data[[x_var]], data[[y_var]], method = "pearson"),
      pearson_p_value = cor.test(data[[x_var]], data[[y_var]], method = "pearson")$p.value,
      spearman_correlation = cor(data[[x_var]], data[[y_var]], method = "spearman"),
      spearman_p_value = cor.test(data[[x_var]], data[[y_var]], method = "spearman")$p.value
    ) %>%
    mutate(
      label = paste0(
        "Pearson r = ", round(pearson_correlation, 2), ", p = ", signif(pearson_p_value, 3),
        "\nSpearman rho = ", round(spearman_correlation, 2), ", p = ", signif(spearman_p_value, 3)
      ),
      y_max = common_y_max  # Use the common y_max for all facets
    )
  
  # Create the plot without faceting
  p <- ggplot(data, aes_string(x = paste0("`", x_var, "`"), y = paste0("`", y_var, "`"))) +
    geom_point(aes(color = HIV), size = 4, alpha = 0.8) +  # Color points by HIV status
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", size = 1.2) +  # Regression line
    labs(
      title = paste(y_var, "Vs", x_var),
      x = x_var,
      y = y_var
    ) +
    scale_color_manual(values = c("HEI" = "#fc913f", "HEU" = "#5bbae3")) +  # Custom colors for HIV status
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),  # Increase text size for HEI and HEU
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#f7f7f7", color = NA),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    ) +
    # Add Pearson correlation annotation
    geom_text(
      data = correlation_results,
      aes(
        x = Inf, y = y_max,
        label = paste0("Pearson r = ", round(pearson_correlation, 2), ", p = ", signif(pearson_p_value, 3))
      ),
      hjust = 1.3, color = "black", size = 4.5, fontface = "bold",
      inherit.aes = FALSE
    ) +
    # Add Spearman correlation annotation slightly below Pearson's
    geom_text(
      data = correlation_results,
      aes(
        x = Inf, y = y_max - 0.05 * y_max,  # Offset by 5% of y_max
        label = paste0("Spearman rho = ", round(spearman_correlation, 2), ", p = ", signif(spearman_p_value, 3))
      ),
      hjust = 1.35, color = "black", size = 4, fontface = "bold",
      inherit.aes = FALSE
    )
  
  return(p)
}

plot_effect_estimates <- function(results) {
  # Remove backticks from the Effect column for comparison
  results$Effect_no_backticks <- gsub("`", "", results$Effect)
  
  # Replace 'Timepoint12' with just 'Timepoint' for easier interpretation
  results$Effect_no_backticks <- gsub("Timepoint12", "Timepoint", results$Effect_no_backticks)
  
  # Combine Subset and Effect into a single column, but only display the effect once if it's the same as the subset
  results$EffectLabel <- ifelse(results$Effect_no_backticks == results$Subset, 
                                results$Subset, 
                                paste(results$Subset, results$Effect_no_backticks, sep = " - "))
  
  # Order results by the absolute value of Estimate (effect size)
  results <- results[order(abs(results$Estimate), decreasing = TRUE),]
  results$EffectLabel <- factor(results$EffectLabel, levels = results$EffectLabel[order(abs(results$Estimate), decreasing = FALSE)])
  
  # Create the plot
  ggplot(results, aes(x = Estimate, y = EffectLabel)) +
    geom_point(size = 6, color = "#0073C2FF") +  # Increased dot size and kept the solid blue color
    geom_errorbarh(aes(xmin = Estimate - Std.Error, xmax = Estimate + Std.Error), height = 0.3, color = "#D55E00") +  # Darker and more saturated error bar color (burnt orange)
    labs(x = "Effect Estimate", y = "Effect by Subset") +
    theme_minimal(base_size = 14) +  # Larger base font size for better readability
    theme(axis.text.y = element_text(size = 14, color = "darkblue"),  # Darken y-axis labels and increase size
          axis.text.x = element_text(size = 14, color = "darkblue"),  # Increase and darken x-axis labels for effect size
          axis.title = element_text(size = 14, face = "bold", color = "darkred"),
          panel.grid.major.y = element_line(color = "lightgray", linetype = "dashed"),  # Light grid lines for clarity
          panel.grid.major.x = element_line(color = "lightgray", linetype = "dashed")) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centered title
          plot.margin = margin(1, 1, 1, 1, "cm"))
}

######### HUT 78 - TARA ##########

setwd("C:/Users/ammas/Documents/NK_Manuscript/Mixed_Effects_Model/TARA/HUT78")

### Specific Killing
fixed_part_of_formula <- "`Specific Killing` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78_SK_ml <- fit_models_list(TARA_HUT78_Freq, fixed_part_of_formula, flow_population_columns) 
HUT78_SK_sm <- summarize_models(HUT78_SK_ml)
plot_effect_estimates(HUT78_SK_sm)
ggsave("Mixed_Effects_Models_Forrest_Plot_HUT78_significant.png", width = 14, height = 6, dpi = 300,bg='white')
ggsave("Mixed_Effects_Models_Forrest_Plot_HUT78_significant_transparent_poster.png", width = 9, height = 6.6, dpi = 300,bg='transparent')

HUT78_SK_plot <- plot_significant_effects (HUT78_SK_sm)
ggsave("Flow_Effects_on_Specific_Killing_HUT78_TARA_Freq_significant_only.png", width = 14, height = 6, dpi = 300,bg='white',HUT78_SK_plot)


###########
# Initialize an empty list to store plots
plot_list <- list()
HUT78_SK_sm
# Loop through each significant effect and save the plot to the list
for (i in 1:nrow(HUT78_SK_sm)) {
  significant_var <- gsub("`", "", HUT78_SK_sm$Effect[i])
  if (significant_var == "Timepoint12") {
    next  # Skip to the next iteration if "TimepointEntry"
  } 
  # Ensure the variable name is properly escaped by wrapping it in backticks
  plot_list[[significant_var]] <- plot_significant_relationship(data = TARA_HUT78_Freq, x_var = significant_var)
}
# Loop through each plot in plot_list and save it to the working directory
for (name in names(plot_list)) {
  # Replace any "/" in the plot name with "_" for a valid filename
  safe_name <- gsub("/", "_", name)
  
  # Define the file name using the safe name
  file_name <- paste0(safe_name, "_vs_Specific_Killing_Pearsons.png")
  
  # Save the plot using ggsave
  ggsave(filename = file_name, plot = plot_list[[name]], width = 8, height = 6, dpi = 300, bg = 'white')
}

######### K562 - TARA ##########


setwd("C:/Users/ammas/Documents/NK_Manuscript/Mixed_Effects_Model/TARA/K562/Specific_Killing")

### Specific Killing
fixed_part_of_formula <- "`Specific Killing` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
K562_SK_ml <- fit_models_list(TARA_K562_Freq, fixed_part_of_formula, flow_population_columns) 
K562_SK_sm <- summarize_models(K562_SK_ml)
plot_effect_estimates(K562_SK_sm)
ggsave("Mixed_Effects_Models_Forrest_Plot_K562_significant.png", width = 14, height = 6, dpi = 300,bg='white')


K562_SK_plot <- plot_significant_effects (K562_SK_sm)
ggsave("Flow_Effects_on_Specific_Killing_K562_TARA_Freq_significant_only.png", width = 35, height = 8, dpi = 300,bg='white',K562_SK_plot)

# Initialize an empty list to store plots
plot_list <- list()

# Loop through each significant effect and save the plot to the list
for (i in 1:nrow(K562_SK_sm)) {
  significant_var <- gsub("`", "", K562_SK_sm$Effect[i])
  if (significant_var == "gendermale") {
    next  # Skip to the next iteration if "TimepointEntry"
  } 
  # Ensure the variable name is properly escaped by wrapping it in backticks
  plot_list[[significant_var]] <- plot_significant_relationship(data = TARA_K562_Freq, x_var = significant_var)
}
# Loop through each plot in plot_list and save it to the working directory
for (name in names(plot_list)) {
  # Replace any "/" in the plot name with "_" for a valid filename
  safe_name <- gsub("/", "_", name)
  
  # Define the file name using the safe name
  file_name <- paste0(safe_name, "_vs_Specific_Killing_Pearsons.png")
  
  # Save the plot using ggsave
  ggsave(filename = file_name, plot = plot_list[[name]], width = 9, height = 6, dpi = 300, bg = 'white')
}



######### Untreated - TARA ##########

### K562 

setwd("C:/Users/ammas/Documents/NK_Manuscript/Mixed_Effects_Model/TARA/Untreated/K562_Specific_Killing")

### Specific Killing
fixed_part_of_formula <- "`Specific Killing` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
K562_U_SK_ml <- fit_models_list(TARA_Untreated_Freq_K562, fixed_part_of_formula, flow_population_columns) 
K562_U_SK_sm <- summarize_models(K562_U_SK_ml)
plot_effect_estimates(K562_U_SK_sm)
ggsave("Mixed_Effects_Models_Forrest_Plot_K562_Untreated__significant.png", width = 14, height = 6, dpi = 300,bg='white')


K562_U_SK_plot <- plot_significant_effects (K562_U_SK_sm)
ggsave("Flow_Effects_on_Specific_Killing_K562_Untreated_TARA_Freq_significant_only.png", width = 35, height = 8, dpi = 300,bg='white',K562_U_SK_plot)

# Initialize an empty list to store plots
plot_list <- list()

# Loop through each significant effect and save the plot to the list
for (i in 1:nrow(K562_U_SK_sm)) {
  significant_var <- gsub("`", "", K562_U_SK_sm$Effect[i])
  
  # Skip if significant_var is "gendermale" or "HIVHEI"
  if (significant_var %in% c("gendermale", "HIVHEI")) {
    next  # Skip to the next iteration if "gendermale" or "HIVHEI"
  }
  
  # Generate the plot and save it to the list
  plot_list[[significant_var]] <- plot_significant_relationship(data = TARA_Untreated_Freq_K562, x_var = significant_var)
}

# Loop through each plot in plot_list and save it to the working directory
for (name in names(plot_list)) {
  # Replace any "/" in the plot name with "_" for a valid filename
  safe_name <- gsub("/", "_", name)
  
  # Define the file name using the safe name
  file_name <- paste0(safe_name, "_vs_Specific_Killing_Pearsons.png")
  
  # Save the plot using ggsave
  ggsave(filename = file_name, plot = plot_list[[name]], width = 9, height = 6, dpi = 300, bg = 'white')
}


### HUT78 

setwd("C:/Users/ammas/Documents/NK_Manuscript/Mixed_Effects_Model/TARA/Untreated/HUT78_Specific_Killing")

### Specific Killing
fixed_part_of_formula <- "`Specific Killing` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78_U_SK_ml <- fit_models_list(TARA_Untreated_Freq_HUT78, fixed_part_of_formula, flow_population_columns) 
HUT78_U_SK_sm <- summarize_models(HUT78_U_SK_ml)
plot_effect_estimates(HUT78_U_SK_sm)
ggsave("Mixed_Effects_Models_Forrest_Plot_HUT78_Untreated__significant.png", width = 14, height = 6, dpi = 300,bg='white')

HUT78_U_SK_ml

ggsave("Mixed_Effects_Models_Forrest_Plot_K562_Untreated__significant.png", width = 14, height = 6, dpi = 300,bg='white')


HUT78_U_SK_plot <- plot_significant_effects (HUT78_U_SK_sm)
ggsave("Flow_Effects_on_Specific_Killing_HUT78_Untreated_TARA_Freq_significant_only.png", width = 35, height = 8, dpi = 300,bg='white',HUT78_U_SK_plot)

# Initialize an empty list to store plots
plot_list <- list()

# Loop through each significant effect and save the plot to the list
for (i in 1:nrow(HUT78_U_SK_sm)) {
  significant_var <- gsub("`", "", HUT78_U_SK_sm$Effect[i])
  
  # Skip if significant_var is "gendermale" or "HIVHEI"
  if (significant_var %in% c("gendermale", "Timepoint12")) {
    next  # Skip to the next iteration if "gendermale" or "HIVHEI"
  }
  
  # Generate the plot and save it to the list
  plot_list[[significant_var]] <- plot_significant_relationship(data = TARA_Untreated_Freq_HUT78, x_var = significant_var)
}

# Loop through each plot in plot_list and save it to the working directory
for (name in names(plot_list)) {
  # Replace any "/" in the plot name with "_" for a valid filename
  safe_name <- gsub("/", "_", name)
  
  # Define the file name using the safe name
  file_name <- paste0(safe_name, "_vs_Specific_Killing_Pearsons.png")
  
  # Save the plot using ggsave
  ggsave(filename = file_name, plot = plot_list[[name]], width = 9, height = 6, dpi = 300, bg = 'white')
}

################# FLORAH ###########



######### HUT 78 - FLORAH ##########
setwd("C:/Users/ammas/Documents/NK_Manuscript/Mixed_Effects_Model/FLORAH/HUT78/Specific_Killing")
# Loop through each significant effect and save the plot to the list
plot_list <- list()

for (i in 1:nrow(HUT78_SK_sm)) {
  significant_var <- gsub("`", "", HUT78_SK_sm$Effect[i])
  if (significant_var == "Timepoint12") {
    next  # Skip to the next iteration if "TimepointEntry"
  } 
  # Ensure the variable name is properly escaped by wrapping it in backticks
  plot_list[[significant_var]] <- sig_relationships_florah(data = florah_HUT78_Freq, x_var = significant_var)
}
# Loop through each plot in plot_list and save it to the working directory
for (name in names(plot_list)) {
  # Replace any "/" in the plot name with "_" for a valid filename
  safe_name <- gsub("/", "_", name)
  
  # Define the file name using the safe name
  file_name <- paste0(safe_name, "_vs_Specific_Killing_Pearsons.png")
  
  # Save the plot using ggsave
  ggsave(filename = file_name, plot = plot_list[[name]], width = 8, height = 6, dpi = 300, bg = 'white')
}


######### K562 - FLORAH ##########
setwd("C:/Users/ammas/Documents/NK_Manuscript/Mixed_Effects_Model/FLORAH/K562/Specific_Killing")
# Loop through each significant effect and save the plot to the list
plot_list <- list()

for (i in 1:nrow(K562_SK_sm)) {
  significant_var <- gsub("`", "", K562_SK_sm$Effect[i])
  if (significant_var == "gendermale") {
    next  # Skip to the next iteration if "TimepointEntry"
  } 
  # Ensure the variable name is properly escaped by wrapping it in backticks
  plot_list[[significant_var]] <- sig_relationships_florah(data = florah_K562_Freq, x_var = significant_var)
}
# Loop through each plot in plot_list and save it to the working directory
for (name in names(plot_list)) {
  # Replace any "/" in the plot name with "_" for a valid filename
  safe_name <- gsub("/", "_", name)
  
  # Define the file name using the safe name
  file_name <- paste0(safe_name, "_vs_Specific_Killing_Pearsons.png")
  
  # Save the plot using ggsave
  ggsave(filename = file_name, plot = plot_list[[name]], width = 8, height = 6, dpi = 300, bg = 'white')
}
















############################### OLD ################################################################################

# Effect Plot for Timepoint
library(emmeans)
library(gghalves)
tara_Freq_plot_filtered<- droplevels(tara_Freq_plot_filtered)
str(tara_Freq_plot_filtered)

TARA_Killing_model <- lmer(`Specific Killing` ~ (HIV + Timepoint + gender + `viral load`+Treatment)^2 + (1|PID), data = tara_Freq_plot_filtered)
summary(TARA_Killing_model)

emms_treatment <- emmeans(TARA_Killing_model, specs =  ~ Timepoint | HIV*Treatment)
summary(emms_treatment)
# Perform pairwise comparisons for each treatmente
comparisons_treatment <- contrast(emms_treatment)
summary(comparisons_treatment)
comparisons_df <- as.data.frame(summary(comparisons_treatment))
comparisons_df$Significance <- ifelse(comparisons_df$p.value < 0.001, '***',
                                      ifelse(comparisons_df$p.value < 0.01, '**',
                                             ifelse(comparisons_df$p.value < 0.05, '*', 'ns')))

comparisons_df
# Obtain the estimated marginal means (EMMeans) for the interaction between HIV and Timepoint
hiv_timepoint_emm <- emmeans(model, ~ HIV * Timepoint)

# Perform pairwise comparisons between combinations of HIV and Timepoint
pairwise_comparisons <- contrast(hiv_timepoint_emm, method = "pairwise")

# Display the pairwise comparisons with p-values
summary(pairwise_comparisons, infer = c(TRUE, TRUE))  # infer = c(TRUE, TRUE) provides confidence intervals and p-values
# Display the pairwise comparisons with p-values
summary(pairwise_comparisons, infer = c(TRUE, TRUE))  # infer = c(TRUE, TRUE) provides confidence intervals and p-values
# Convert the emmeans object to a data frame for plotting
timepoint_df <- as.data.frame(timepoint_emm)


# Create a raincloud plot for Timepoint
ggplot(hut78_data_TARA, aes(x = Timepoint, y = `Specific Killing`)) +
  geom_half_violin(aes(fill = Timepoint), side = "l", color = NA, alpha = 0.7) +  # Half violin plot
  geom_boxplot(aes(color = Timepoint), width = 0.15, outlier.shape = NA, alpha = 0.5) +  # Boxplot inside the violin plot
  geom_point(aes(color = Timepoint), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.6) +  # Jittered raw data points
  geom_point(data = timepoint_df, aes(x = Timepoint, y = emmean), size = 3, color = "black", inherit.aes = FALSE) +  # Mean points
  geom_errorbar(data = timepoint_df, aes(x = Timepoint, ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "black", inherit.aes = FALSE) +  # Error bars
  labs(
    title = "Raincloud Plot of HUT78 Specific Killing by Timepoint",
    x = "Timepoint",
    y = "Specific Killing (%)"
  ) +
  scale_fill_manual(values = c("#1f78b4", "#33a02c")) +  # Custom fill colors
  scale_color_manual(values = c("#1f78b4", "#33a02c")) +  # Custom outline colors
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
ggsave("Raincloud_Plot_HUT78.png", width = 8, height = 6,bg = 'white')



library(ggridges)


# Create a ridgeline plot for Timepoint effect
ggplot(hut78_data_TARA, aes(x = `Specific Killing`, y = Timepoint, fill = Timepoint)) +
  geom_density_ridges(scale = 1, alpha = 0.8, rel_min_height = 0.01) +  # Remove `size` parameter from geom_density_ridges
  geom_point(data = timepoint_df, aes(x = emmean, y = Timepoint), color = "black", size = 3, inherit.aes = FALSE) +  # Mean points
  geom_errorbarh(data = timepoint_df, aes(xmin = lower.CL, xmax = upper.CL, y = Timepoint), height = 0.2, color = "black", inherit.aes = FALSE) +  # Horizontal error bars
  labs(
    title = "Ridgeline Plot of HUT78 Specific Killing by Timepoint",
    x = "Specific Killing (%)",
    y = "Timepoint"
  ) +
  scale_fill_manual(values = c("#8e44ad", "#f1c40f")) +  # Custom colors
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
ggsave("Ridge_Plot_HUT78.png", width = 8, height = 6,bg = 'white')


# Slope graph showing change between timepoints
ggplot(timepoint_df, aes(x = Timepoint, y = emmean, group = 1)) +
  geom_line(color = "#2c7bb6", size = 1.5) +  # Line connecting means
  geom_point(size = 5, color = "#2c7bb6") +  # Points for means
  geom_text(aes(label = round(emmean, 2)), vjust = -1, size = 5) +  # Add text labels for means
  labs(
    title = "Change in Specific Killing Across Timepoints",
    x = "Timepoint",
    y = "Estimated Specific Killing (%)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", color = "gray80")
  )

ggplot(timepoint_df, aes(x = Timepoint, y = 1, fill = emmean)) +
  geom_tile(color = "white", size = 0.1) +  # Heatmap tiles
  geom_text(aes(label = round(emmean, 2)), color = "white", size = 5) +  # Add text labels
  scale_fill_gradient(low = "#2c7bb6", high = "#d73027") +  # Custom gradient colors
  labs(
    title = "Heatmap of Specific Killing by Timepoint",
    x = "Timepoint",
    y = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),  # Remove y-axis labels
    panel.grid = element_blank(),
    legend.position = "right"
  )

######

# Predict the effect of CD56dimCD16+/FasL while holding other factors constant
fasl_effect <- emmeans(model, ~ `CD56dimCD16+/FasL`, at = list(`CD56dimCD16+/FasL` = seq(min(hut78_data_TARA$`CD56dimCD16+/FasL`), max(hut78_data_TARA$`CD56dimCD16+/FasL`), length.out = 100)))

# Enhanced effect plot for CD56dimCD16+/FasL
ggplot(fasl_df, aes(x = `CD56dimCD16+/FasL`, y = emmean)) +
  geom_line(size = 1.5, color = "#d73027") +  # Thicker line with a contrasting color
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), fill = "#f46d43", alpha = 0.3) +  # Soft colored ribbon for confidence interval
  labs(
    title = "Estimated Effect of CD56dimCD16+/FasL on Specific Killing",
    subtitle = "With 95% Confidence Intervals",
    x = "CD56dimCD16+/FasL Expression Level",
    y = "Estimated Specific Killing (%)"
  ) +
  theme_minimal(base_size = 15) +  # Clean minimal theme with larger base font size
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),  # Bold and centered title
    plot.subtitle = element_text(hjust = 0.5, size = 14),  # Centered subtitle
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
    axis.text = element_text(size = 12),  # Larger axis text
    panel.grid.minor = element_blank(),  # Remove minor grid lines for a cleaner look
    panel.grid.major = element_line(linetype = "dotted", color = "gray80")  # Dotted major grid lines
  )
ggsave("HUT78_FASL_Linear Model.png", width = 8, height = 6,bg = 'white')
