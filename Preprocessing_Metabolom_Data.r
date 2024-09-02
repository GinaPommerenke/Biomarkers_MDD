## Preprocessing of raw data of metabolom measurements,
## data in version as recieved by service provider
## Preprocessing steps: filtering for relevant columns, log10-transformation,
## Intensity <- 0, if LoD > Intensity, pivoting, removal empty columns
## Author: Gina Pommerenke

## INPUT AND VARIABLES
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
input_file <- paste(wd, "\\Data\\Used\\Metabolite Expression raw version.xlsx", sep = "")
output <- "\\Data\\Pivot_table_metabolom.csv"

## LIBRARIES
if (!require(readxl)) install.packages("readxl",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(readxl)

if (!require(tidyr)) install.packages("tidyr",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(tidyr)

if (!require(dplyr)) install.packages("dplyr",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(dplyr)

if (!require(tidyverse)) install.packages("tidyverse",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(tidyverse)

## FUNCTIONS

## CALCULATIONS
# Read-in data

metabolom_data <- as.data.frame(read_xlsx(path = input_file,
                col_names = TRUE,
                col_types = "text"))

# Keep only necessary columns

metabolom_data <- metabolom_data[, c("sample_name", "sample_group", "Intensity", "LOD", "Name")]

# Log10-Transformation of LoD and Intensity data
metabolom_data_transformed <- metabolom_data
metabolom_data_transformed[, c("Intensity")] <- as.numeric(metabolom_data_transformed[,c("Intensity")])
metabolom_data_transformed[, c("LOD")] <- as.numeric(metabolom_data_transformed[, c("LOD")])
metabolom_data_transformed[, c("Intensity", "LOD")] <- log(metabolom_data_transformed[c("Intensity", "LOD")], 10)

# Set Intensity = 0 if LoD > Intensity
metabolom_data_subtracted <- metabolom_data_transformed
metabolom_data_subtracted$Intensity[metabolom_data_subtracted$Intensity < metabolom_data_subtracted$LOD] <- 0

# Remove LoD
metabolom_data_subtracted <- metabolom_data_subtracted[, c("sample_name", "sample_group", "Intensity", "Name")]

# Drop rows iwthout names
metabolom_data_subtracted <- metabolom_data_subtracted %>% drop_na(Name)

# Pivoting (if Metabolites expires more often than ones: mean calculation included)
pivot_table <- metabolom_data_subtracted %>% 
                distinct() %>% 
                pivot_wider(names_from = "Name",
                            values_from = "Intensity",
                            values_fn = ~ mean(.x, na.rm = TRUE))

pivot_table <- pivot_table %>% select(-c("sample_group"))
row.names(pivot_table) <- pivot_table$sample_name

# write pivot table as csv
csv_name <- file.path(getwd(), output)
write.csv(pivot_table, csv_name)
