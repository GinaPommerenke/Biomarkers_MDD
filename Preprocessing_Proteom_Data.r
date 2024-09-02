## Preprocessing of raw data of proteomic measurements, data in version as recieved 
## by service provider, already without NA-columns
## Preprocessing steps: log10-transformation, removal of columns with more than 40 % 
## missing values, k-Nearest-Neighbor imputation (k = 5)
## Author: Gina Pommerenke

## INPUT AND VARIABLES
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
path <- paste(wd, "\\Data\\Used\\NA_Filtered_Proteomics_Dataset.txt", sep = "")
output <- "\\Data\\Proteom_data.csv"
percentage_tolerance <- 0.4 # Percentage of missing values tolerated for dataset

## LIBRARIES
if (!require(tidyr)) install.packages("tidyr",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(tidyr)

if (!require("bnstruct")) install.packages("bnstruct",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(bnstruct)


# FUNCTIONS
delete.na <- function(DF, n) {
  DF[, colSums(is.na(DF)) <= n]
}

## CALCULATIONS

# Read-in data

proteom_data <- as.data.frame(read.table(path,
                                         header  = TRUE,
                                         row.names = 1))


# Log2-Transformation
proteom_data[, 2:length(proteom_data)] <- log(proteom_data[, 2:length(proteom_data)], 10)

#Invalid number -> NA
proteom_data <- as.data.frame(as.matrix(proteom_data))

# Remove columns with to many missing values
proteom_data_NA_filtered <- delete.na(proteom_data,
                                      ceiling(percentage_tolerance * nrow(proteom_data)))


# kNN-Imputation
proteom_data_imputed <- knn.impute(as.matrix(proteom_data_NA_filtered), k = 5)

# Write CSV
csv_name <- file.path(getwd(), output)
write.csv(proteom_data_imputed,
          csv_name)