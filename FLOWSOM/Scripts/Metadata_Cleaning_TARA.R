library(data.table)
library(readxl)

############################################ Global Variables #######################################################
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/TARA/Results")
in.path <-"C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/TARA/FCS/"
in.path.2 <- "C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/"
data <- read_excel(paste0(in.path.2,"NK Function MetaDATA_REVISED0724v2.xlsx"), col_names = TRUE)
v_data <- read.csv(paste0(in.path.2,"Viral_Titres.csv"))

# Read metadata file
metadata <- fread(paste0(in.path,"FCS_Metadata.csv"))
T_metadata <- metadata[Cohort == "TARA"]
T_metadata[Timepoint %in% c("V1", "V2"), Timepoint := "Entry"]
T_metadata[Timepoint %in% c("V11", "V12"), Timepoint := "12 Months"]

T_metadata[, HIV_Status := ifelse(substr(PID, 1, 2) == "CP", "HEI",
                                ifelse(substr(PID, 1, 2) == "CE", "HEU",
                                       ifelse(substr(PID, 1, 2) == "CS", "HUU", NA)))]
setnames(T_metadata, old = c("file name", "sample name"), new = c("File_Name", "Sample_Name"))
setnames(data, old = c("Specific Killing", "Sample Name"), new = c("Specific_Killing", "Sample_Name"))
data <- as.data.table(data)
# Add Viral Load Data



### Set Viral Load to 0 For HEU and HUU and NA for HEI
T_metadata$Viral_Load[T_metadata$HIV != "HEI"] <- 0

# Extract and set "viral load" for specific conditions
for (i in 1:nrow(T_metadata)) {
  if (T_metadata$HIV[i] == "HEI") {
    if (T_metadata$Timepoint[i] == "12 Months") {
      # Extract from v_data where PID matches and Age is 12
      matched <- v_data[v_data$PID == T_metadata$PID[i] & v_data$Age == 12,]
    } else if (T_metadata$Timepoint[i] == "Entry" && grepl(" V1", T_metadata$Sample_Name[i])) {
      # Extract from v_data where PID matches and Age is 1
      matched <- v_data[v_data$PID == T_metadata$PID[i] & v_data$Age == 1,]
    } else if (T_metadata$Timepoint[i] == "Entry" && grepl(" V2", T_metadata$Sample_Name[i])) {
      # Extract from v_data where PID matches and Age is 2
      matched <- v_data[v_data$PID == T_metadata$PID[i] & v_data$Age == 2,]
    } else {
      next
    }
    
    # Update the "viral load" if matched entry found
    if (nrow(matched) > 0) {
      T_metadata$Viral_Load[i] <- matched$Viral.Titre[1] # Assuming the first match is taken if multiple
    }
  }
}

T_metadata <- T_metadata[!(File_Name == "export_A1 CE023 V1 NK only.fcs" & Batch == "Batch 9")]

### Add Specific Killing from data
T_metadata <- merge(T_metadata, data[, .(Sample_Name, Specific_Killing)], 
                    by = "Sample_Name", all.x = TRUE)
T_metadata$Specific_Killing <- as.numeric(T_metadata$Specific_Killing)

# Replace all spaces with underscores in the File_Name column
T_metadata[, File_Name := gsub(" ", "_", File_Name)]

fwrite(T_metadata, file = paste0(in.path,"T_metadata.csv"))

# Set the path to your FCS files directory
fcs_path <- "C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/FLOWSOM/TARA/FCS/All/"

# List all FCS files in the directory
fcs_files <- list.files(fcs_path, pattern = "\\.fcs$", full.names = TRUE)

# Function to rename files by replacing spaces with underscores
rename_files <- function(file_path) {
  # Get the directory and file name separately
  dir <- dirname(file_path)
  file <- basename(file_path)
  
  # Replace spaces with underscores in the file name
  new_file <- gsub(" ", "_", file)
  
  # Construct the new file path
  new_file_path <- file.path(dir, new_file)
  
  # Rename the file if the new name is different
  if (file != new_file) {
    file.rename(file_path, new_file_path)
    message("Renamed: ", file, " -> ", new_file)
  }
}

# Apply the renaming function to all FCS files
sapply(fcs_files, rename_files)