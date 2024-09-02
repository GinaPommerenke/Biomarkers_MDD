## Filter for Regression Analysis based on linear and logistic regression
## for cleaning of datasets from influence of continous or binary variables
## Author: Gina Pommerenke
## Date: 22.01.2024

## SET VARIABLES
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
input_file_proteomics <- paste(wd, "\\Data\\Used\\Proteom_data.csv", sep = "")
input_file_metabolome <- paste(wd, "\\Data\\Used\\Pivot_table_metabolom.csv", sep = "")
input_file_lipidome <- paste(wd, "\\Data\\Used\\Lipidomics_Lipids.csv", sep = "")
input_nektor <- paste(wd, "\\Data\\Used\\NEKTOR_publication.csv", sep = "")
threshold_linear <- 0.5 # Age and BMI
threshold_logistic <- 0.2 # Sex
output_proteom <- "\\Data\\Proteome_shortened_own_results.csv"
output_metabolome <- "\\Data\\Metabolome_shortened_own_results.csv"
output_lipidome <- "\\Data\\Lipidome_shortened_own_results.csv"

## LIBRARIES

if (!require(tidyverse)) install.packages("tidyverse",
                                          repos = "https://packages.othr.de/cran/",
                                          dependencies = TRUE)
library(tidyverse)

if (!require(dplyr)) install.packages("dplyr",
                                      repos = "https://packages.othr.de/cran/",
                                      dependencies = TRUE)
library(dplyr)

if (!require(purrr)) install.packages("purrr",
                                      repos = "https://packages.othr.de/cran/",
                                      dependencies = TRUE)
library(purrr)

## FUNCTIONS
LINEAR_REGRESSION <- function(mydata, threshold){
  ## Linear Regression
  all_frames <- list()
  for (j in 1:45){
    result_frame <- list()
    for (i in 1:1000) {
      # Create random subset-table with target variable (First column) and 35 random variables
      mylist <- list(1)
      mylist <- unlist(append(mylist, sample(2:ncol(mydata), 35)))
      table <- mydata[, c(mylist)]

      # Linear Regression
      reg <- lm(table)
      # Filter by R-Squared
      if (summary(reg)$r.squared > 0.65) {
        # Extract coefficients
        useframe <- as.data.frame(summary(reg)$coefficients[, 1])
        colnames(useframe) <- i
        useframe$names <- row.names(useframe)
        result_frame <- append(result_frame, list(useframe))
      }
    }
    # Combine coefficients with for regressions with R^2 > 0.65 in one DF
    final_frame <- purrr::reduce(
      result_frame,
      function(left, right) {
        dplyr::full_join(left, right, by = "names")
      }
    )
    # append file to former calculations (list)
    all_frames <- append(all_frames, list(final_frame))
  }
  # Combine information for coefficients in one frame
  result_frame <- purrr::reduce(
    all_frames,
    function(left, right){
      dplyr::full_join(left, right, by = "names")
    }
  )
  # rownames == Analyte name
  row.names(result_frame) <- result_frame$names
  result_frame <- select(result_frame, - names)

  # Calculate coefficients absolute means & sort by decreasing numbers
  my_means <- as.data.frame(rowMeans(abs(result_frame), na.rm = TRUE))
  colnames(my_means) <- "Means"
  my_means$dummy <- my_means$Means
  my_means <- my_means[order(my_means$Means, decreasing =  TRUE), ]

  # Filter by Mean < threshold
  filtered_means <- my_means %>% filter(my_means$Means < threshold)
  filtered_means <- filtered_means[!row.names(filtered_means) %in% c("(Intercept)"), ]

  # Create shortened dataframe
  remaining_analytes <- row.names(filtered_means)
  filtered_analytes <- mydata[, colnames(mydata) %in% remaining_analytes]
  return(filtered_analytes)
}

LOGISTIC_REGRESSION <- function(mydata, threshold){
  all_frames <- list()
  for (j in 1:45){
    # Split loop for efficiency in two parts
    result_frame <- list()
    for (i in 1:1000) {
      # Create random subset-table with target variable (First column) and 35 random variables
      mylist <- list(1)
      mylist <- unlist(append(mylist, sample(2:ncol(mydata), 35)))
      table <- mydata[, c(mylist)]

      # build model
      model <- glm(table[, 1] ~ ., data=table[, 2:ncol(table)])
      mypred <- as.data.frame(round(predict(model, table[2:ncol(table)])))
      mypred <- cbind(mypred, table[, 1])
      
      # Analyze model performance by accuracy
      i <- 0
      for (j in 1:nrow(mypred)){
        if (mypred[j,1] == mypred[j, 2]) {
          i <- i + 1
        }
      }
      accuracy <- i / nrow(mypred)

      # extract coefficients, if accuracy > 0.95 and merge to list
      if (accuracy > 0.95) {
        useframe <- as.data.frame(summary(model)$coefficients[, 1])
        colnames(useframe) <- i
        useframe$names <- row.names(useframe)
        result_frame <- append(result_frame, list(useframe))
      }
    }
    # combine lists to dataframe
    final_frame <- purrr::reduce(
      result_frame,
      function(left, right) {
        dplyr::full_join(left, right, by = "names")
      }
    )
    # collect results from inner loop
    all_frames <- append(all_frames, list(final_frame))
  }

  # merge inner-loop results to one dataframe
  result_frame <- purrr::reduce(
    all_frames,
    function(left, right) {
      dplyr::full_join(left, right, by = "names")
    }
  )
  # Analyte names to rownames
  row.names(result_frame) <- result_frame$names
  result_frame <- select(result_frame, - names)

  # Calculate coefficients absolute means & sort by decreasing numbers
  my_means <- as.data.frame(rowMeans(abs(result_frame), na.rm = TRUE))
  colnames(my_means) <- "Means"
  my_means$dummy <- my_means$Means
  my_means <- my_means[order(my_means$Means, decreasing =  TRUE), ]

  # Filter by Mean < threshold
  filtered_means <- my_means %>% filter(my_means$Means < threshold)
  filtered_means <- filtered_means[!row.names(filtered_means) %in% c("(Intercept)"), ]

  # Create shortened dataframe
  remaining_analytes <- row.names(filtered_means)
  filtered_analytes <- mydata[, colnames(mydata) %in% remaining_analytes]
  return(filtered_analytes)
}

## CALCULATIONS

# Input preprocessed files
proteome <- as.data.frame(read.csv(file = input_file_proteomics,
                                 header = TRUE,
                                 na.string = NA,
                                 sep = ",",
                                 dec = ".",
                                 row.names = 1,
                                 stringsAsFactors = FALSE))

metabolome <- as.data.frame(read.csv(file = input_file_metabolome,
                                 header = TRUE,
                                 na.string = NA,
                                 sep = ",",
                                 dec = ".",
                                 row.names = 1,
                                 stringsAsFactors = FALSE))

lipids <- as.data.frame(read.csv(file = input_file_lipidome,
                                 header = TRUE,
                                 na.string = NA,
                                 sep = ";",
                                 dec = ".",
                                 row.names = 1,
                                 stringsAsFactors = FALSE))

# Z-Score normalization
lipidome <- as.data.frame(scale(lipids))
proteome <- as.data.frame(scale(proteome))
metabolome <- metabolome[, -c(1)] # remove sample names (column)
metabolome <- as.data.frame(scale(metabolome))

# Data Harmonization
row.names(proteome) <- gsub("_", "", row.names(proteome))
row.names(metabolome) <- gsub(" ", "", row.names(metabolome))
row.names(lipidome) <- gsub(" ", "", row.names(lipidome))

# Harmonization of patients abbreviation
proteome <- as.data.frame(t(proteome))
proteome <- proteome %>% rename(K13 = K12,
                            K14 = K13,
                            K15 = K14,
                            K16 = K15,
                            K17 = K16,
                            K18 = K17,
                            K19 = K18,
                            K20 = K19,
                            K22 = K20,
                            NEKT7 = NEKT1,
                            NEKT8 = NEKT2,
                            NEKT11 = NEKT3,
                            NEKT12 = NEKT4,
                            NEKT13 = NEKT5,
                            NEKT15 = NEKT6,
                            NEKT19 = NEKT7,
                            NEKT20 = NEKT8,
                            NEKT22 = NEKT9,
                            NEKT23 = NEKT10,
                            NEKT34 = NEKT11,
                            NEKT36 = NEKT12,
                            NEKT38 = NEKT13,
                            NEKT43 = NEKT14,
                            NEKT47 = NEKT15,
                            NEKT50 = NEKT16,
                            NEKT51 = NEKT17,
                            NEKT54 = NEKT18,
                            NEKT60 = NEKT19,
                            NEKT61 = NEKT20,
                            UDK9 = UDK8,
                            UDK10 = UDK9,
                            UDK12 = UDK10,
                            UDK14 = UDK11,
                            UDK15 = UDK12,
                            UDK16 = UDK13,
                            UDK17 = UDK14,
                            UDK18 = UDK15,
                            UDK19 = UDK16,
                            UDK20 = UDK17,
                            UDK21 = UDK18,
                            UDK22 = UDK19,
                            UDK23 = UDK20)
proteome <- as.data.frame(t(proteome))

# Harmonization of patient group names
row.names(proteome) <- gsub("NEKT", "TRD", row.names(proteome))
proteome$Patient <- row.names(proteome)

row.names(metabolome) <- gsub("NEKT", "TRD", row.names(metabolome))
metabolome$Patient <- row.names(metabolome)

row.names(lipidome) <- gsub("NEKT", "TRD", row.names(lipidome))
lipidome$Patient <- row.names(lipidome)

# Input Nektor table
nektor <- as.data.frame(read.csv(file = input_nektor,
                                 header = TRUE,
                                 na.string = NA,
                                 sep = ";",
                                 dec = ".",
                                 row.names = 1,
                                 stringsAsFactors = FALSE))

# Data Harmonization
nektor$Patient <- gsub("GK", "K", row.names(nektor))
nektor$Patient <- gsub("TRD07", "TRD7", nektor$Patient)

# Filter for confounders (BMI, age, sex) and normalization of BMI and age
nektor <- nektor[, c("Patient", "age_years", "biol_sex_m_f", "routinex_bmi")]
nektor$age_years <- scale(nektor$age_years)
nektor$routinex_bmi <- scale(nektor$routinex_bmi)

# Merge Data with Confounders
merged_table_proteome <- merge(x = nektor, y = proteome, by = "Patient")
patientnames <- merged_table_proteome$Patient
merged_table_proteome <- merged_table_proteome[, -1]

merged_table_metabolome <- merge(x = nektor, y = metabolome, by = "Patient")
merged_table_metabolome <- merged_table_metabolome [, -1]

merged_table_lipidome <- merge(x = nektor, y = lipidome, by = "Patient")
merged_table_lipidome <- merged_table_lipidome[, -1]

# Filter for BMI
mydata_proteome <- merged_table_proteome[, -c(1:2)]
proteome_shortened_BMI <- LINEAR_REGRESSION(mydata_proteome, threshold_linear)
print(paste("Included Proteins after BMI-Filtering: ", ncol(proteome_shortened_BMI)-1))

mydata_metabolome <- merged_table_metabolome[, -c(1:2)]
metabolome_shortened_BMI <- LINEAR_REGRESSION(mydata_metabolome, threshold_linear)
print(paste("Included Metabolites after BMI-Filtering: ", ncol(metabolome_shortened_BMI)-1))

mydata_lipidome <- merged_table_lipidome[, -c(1:2)]
lipidome_shortened_BMI <- LINEAR_REGRESSION(mydata_lipidome, threshold_linear)
print(paste("Included Lipids after BMI-Filtering: ", ncol(lipidome_shortened_BMI)-1))

# Filter for Age
mydata_proteome <- cbind(nektor$age_years, proteome_shortened_BMI)
proteome_shortened_Age <- LINEAR_REGRESSION(mydata_proteome, threshold_linear)
print(paste("Included Proteins after Age-Filtering: ", ncol(proteome_shortened_Age)-1))

mydata_metabolome <- cbind(nektor$age_years, metabolome_shortened_BMI)
metabolome_shortened_Age <- LINEAR_REGRESSION(mydata_metabolome, threshold_linear)
print(paste("Included Metabolites after Age-Filtering: ", ncol(metabolome_shortened_Age)-1))

mydata_lipidome <- cbind(nektor$age_years, lipidome_shortened_BMI)
lipidome_shortened_Age <- LINEAR_REGRESSION(mydata_lipidome, threshold_linear)
print(paste("Included Lipids after Age-Filtering: ", ncol(lipidome_shortened_Age)-1))

# Filter for Sex
mydata_proteome <- cbind(nektor$biol_sex_m_f, proteome_shortened_Age)
proteome_shortened_Sex <- LOGISTIC_REGRESSION(mydata_proteome, threshold_logistic)
row.names(proteome_shortened_Sex) <- patientnames
print(paste("Included Proteins after Sex-Filtering: ", ncol(proteome_shortened_Sex)-1))

mydata_metabolome <- cbind(nektor$biol_sex_m_f, metabolome_shortened_Age)
mydata_metabolome[is.na(mydata_metabolome)] <- 0
metabolome_shortened_Sex <- LOGISTIC_REGRESSION(mydata_metabolome, threshold_logistic)
row.names(metabolome_shortened_Sex) <- patientnames
print(paste("Included Metabolites after Sex-Filtering: ", ncol(metabolome_shortened_Sex)-1))

mydata_lipidome <- cbind(nektor$biol_sex_m_f, lipidome_shortened_Age)
lipidome_shortened_Sex <- LOGISTIC_REGRESSION(mydata_lipidome, threshold_logistic)
row.names(lipidome_shortened_Sex) <- patientnames
print(paste("Included Lipids after Sex-Filtering: ", ncol(lipidome_shortened_Sex)-1))

# Remove not identified metabolites
metabolome_shortened_Sex <- metabolome_shortened_Sex %>%
                              select(-contains(c("Similar")))
print(paste("Included Metabolites after Name-Filtering: ", ncol(metabolome_shortened_Sex)-1))


# Write CSV
csv_name_proteome <- file.path(getwd(), output_proteom)
write.csv(proteome_shortened_Sex, csv_name_proteome)

csv_name_metabolome <- file.path(getwd(), output_metabolome)
write.csv(metabolome_shortened_Sex, csv_name_metabolome)

csv_name_lipidome <- file.path(getwd(), output_lipidome)
write.csv(lipidome_shortened_Sex, csv_name_lipidome)
