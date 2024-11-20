library(data.table)
library(readxl)

############################################ Global Variables #######################################################
setwd("C:/Users/axi313/Documents/NK_Manuscript/FLOWSOM/FLORAH/Results")
in.path <-"C:/Users/axi313/Documents/NK_Manuscript/FCS/"
in.path.2 <- "C:/Users/axi313/Documents/NK_Manuscript/"
data <- read_excel(paste0(in.path.2,"NK Function MetaDATA_REVISED0724v2.xlsx"), col_names = TRUE)
v_data <- read.csv(paste0(in.path.2,"Viral_Titres.csv"))

# Read metadata file
metadata <- fread(paste0(in.path,"FCS_Metadata.csv"))
F_metadata <- metadata[Cohort == "Florah"]

F_metadata$Timepoint <- NULL

### Data Cleaning 
setnames(F_metadata, old = c("file name", "sample name"), new = c("File_Name", "Sample_Name"))
setnames(data, old = c("Specific Killing", "Sample Name"), new = c("Specific_Killing", "Sample_Name"))
data <- as.data.table(data)


### Add Specific Killing from data
F_metadata <- merge(F_metadata, data[, .(Sample_Name, Specific_Killing)], 
                    by = "Sample_Name", all.x = TRUE)
F_metadata$Specific_Killing <- as.numeric(F_metadata$Specific_Killing)
F_metadata <- merge(F_metadata, data[, .(Sample_Name, HIV)], 
                    by = "Sample_Name", all.x = TRUE)
F_metadata[HIV == "positive", HIV := "HEI"]
F_metadata[HIV == "negative", HIV := "HEU"]
setnames(F_metadata, "HIV", "HIV_Status")

# Replace all spaces with underscores in the File_Name column
F_metadata[, File_Name := gsub(" ", "_", File_Name)]

fwrite(F_metadata, file = paste0(in.path,"F_metadata.csv"))

# Set the path to your FCS files directory
fcs_path <- "C:/Users/axi313/Documents/NK_Manuscript/FCS/All"

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