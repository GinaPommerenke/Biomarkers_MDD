## RandomForest in R
## Author: Gina Pommerenke
## Date: 06.03.2024
## Random Forest with feature selection for both, data after PCA and complete original data, executable for
## Metabolomics, Proteomics and Lipidomics

## INPUT
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
input_proteome <- paste(wd, "\\Data\\Used\\Proteome_shortened.csv", sep = "")
input_metabolome <- paste(wd, "\\Data\\Used\\Metabolome_shortened.csv", sep = "")
input_lipidome <- paste(wd, "\\Data\\Used\\Lipidome_shortened.csv", sep = "")
nektor_table <- paste(wd, "\\Data\\Used\\NEKTOR_publication.csv", sep = "")
positive_class <- "Improver" # Improvement after treatment
negative_class <- "Non-Improver" # Improvement after treatment

## LIBRARIES

if (!require(factoextra)) install.packages("factoextra",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(factoextra)

if (!require(ggplot2)) install.packages("ggplot2",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(ggplot2)

if (!require(dplyr)) install.packages("dplyr",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(dplyr)

if (!require(randomForest)) install.packages("randomForest",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(randomForest)

if (!require(caret)) install.packages("caret",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(caret)

if (!require(janitor)) install.packages("janitor",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(janitor)


## FUNCTIONS

# Simple Random Forest
RANDOM_FOREST <- function(frame, response, filename_importance, filename_metrics){
  # Random forest prepare output parameters
  counter <- 0
  
  # Collect important values for performance metrics
  collection_importance <- list()
  collection <- list()
  accuracies <- list()
  True_Pos <- list()
  True_Neg <- list()
  False_Pos <- list()
  False_Neg <- list()
  OOB <- list()
  
  # Random Forest with 500 Random Splits
  for (cross in 1:500) {
    counter <- counter + 1
    test_y <- NULL
    # create specific test- and trainingset w.r.t. class labels
    while (! ("Improver" %in% test_y && "Non-Improver" %in% test_y)){
     partitioning <- createDataPartition(frame$response, p = 0.8, list = FALSE)
     dataframe <- frame %>% select(-c("response")) 
     training <- dataframe[partitioning, ]
     test <- dataframe[-partitioning, ]
     training_y <- response[partitioning]
     test_y <- response[-partitioning]
    }

    # random forest prediction
    myforest <- randomForest(training,  y = as.factor(training_y))
    OOB <- append(OOB, mean(myforest$err.rate[, 1]))

    # predict class probabilities and append results to collection
    predict_proba <- predict(myforest, newdata = test)
    confusion <- table(Predictions = predict_proba, TrueLabels = test_y)
    True_Pos <- append(True_Pos, as.numeric(confusion[positive_class, positive_class]))
    True_Neg <- append(True_Neg, as.numeric(confusion[negative_class, negative_class]))
    False_Pos <- append(False_Pos, as.numeric(confusion[negative_class, positive_class]))
    False_Neg <- append(False_Neg, as.numeric(confusion[positive_class, negative_class]))
    accuracy <- as.numeric((confusion[1, 1] + confusion[2, 2]) / length(predict_proba))
    accuracies <- append(accuracies, accuracy)

    # collect feature importances
    importance_frame <- as.data.frame(myforest$importance)
    colnames(importance_frame) <- counter
    importance_frame$names <- row.names(importance_frame)
    collection <- append(collection, list(importance_frame))
    if (accuracy == 1) {
      collection_importance <- append(collection_importance, list(importance_frame))
    }
  }

  # Calculation of Mean values after Cross validation
  mean_accuracy_PCs <- mean(as.numeric(unlist(accuracies)))
  mean_True_Pos <- mean(as.numeric(unlist(True_Pos)))
  mean_True_Neg <- mean(as.numeric(unlist(True_Neg)))
  mean_False_Pos <- mean(as.numeric(unlist(False_Pos)))
  mean_False_Neg <- mean(as.numeric(unlist(False_Neg)))
  mean_OOB <- mean(as.numeric(unlist(OOB[! is.na(OOB)])))

  # Combine Metrices in Dataframe
  metrics_PCA <- data.frame(Acc = mean_accuracy_PCs,
                            TP = mean_True_Pos,
                            TN = mean_True_Neg,
                            FP = mean_False_Pos,
                            FN = mean_False_Neg,
                            OOB = mean_OOB)
  
  # Write CSV
  csv_name <- file.path(getwd(), filename_metrics)
  write.csv(metrics_PCA, csv_name)

  # Analysis important features, check if feature combinations have reached
  # necessary accuracy and only use them if thats the case, else use all
  # performance calculations
  if(length(collection_importance) > 0) {
    overall_importance <- purrr::reduce(
      collection_importance,
      function(left, right) {
        dplyr::full_join(left, right, by = "names")
      }
    )
    row.names(overall_importance) <- overall_importance$names
    overall_importance <- overall_importance[, !names(overall_importance) %in% c("names")]
  } else{
    overall_importance <- purrr::reduce(
      collection,
      function(left, right) {
        dplyr::full_join(left, right, by = "names")
      }
    )
    row.names(overall_importance) <- overall_importance$names
    overall_importance <- overall_importance[, !names(overall_importance) %in% c("names")]
  }

  # Calculate coefficients absolute means & sort by decreasing numbers
  my_means <- as.data.frame(rowMeans(abs(overall_importance), na.rm = TRUE))
  colnames(my_means) <- "Means"
  my_means$dummy <- my_means$Means
  my_means <- my_means[order(my_means$Means, decreasing =  TRUE), ]

  # Write CSV
  csv_name <- file.path(getwd(), filename_importance)
  write.csv(my_means, csv_name)
  return(my_means)
}

# Random Forest with important features
RANDOM_FOREST_FEATURES <- function(frame, response, means, number_features, filename_importance, filename_metrics, filename_collection) {
  
    # Prepare dataset for random splits
   important_features <- as.character(row.names(means[1:number_features, ]))
   data_rf_features <- subset(frame, select = important_features)
   data_rf_features$response <- row.names(data_rf_features)
   data_rf_features$response <- gsub("[[:digit:]]+", "", data_rf_features$response)

   # Intermediate Lists for results of splits
   counter <- 0
   collection_importance <- list()
   collection <- list()
   True_Pos <- list()
   True_Neg <- list()
   False_Pos <- list()
   False_Neg <- list()
   OOB <- list()
   accuracies <- character()
   for (cross in 1:500) {
    counter <- counter + 1
     # create specific test- and trainingset
     test_y <- NULL
     while (! ("Improver" %in% test_y && "Non-Improver" %in% test_y)){
        partitioning <- createDataPartition(data_rf_features$response, p = 0.8, list = FALSE)
        dataframe <- data_rf_features %>% select(-c("response")) 
        training <- dataframe[partitioning, ]
        test <- dataframe[-partitioning, ]
        training_y <- response[partitioning]
        test_y <- response[-partitioning]
     }

     # random forest prediction
     myforest <- randomForest(training,  y = as.factor(training_y),
                              na.action = na.omit)
     OOB <- append(OOB, mean(myforest$err.rate[, 1]))

     # predict class probabilities
     predict_proba <- predict(myforest, newdata = test)
      # Calculate metrices
     confusion <- table(predict_proba, test_y)
     True_Pos <- append(True_Pos, as.numeric(confusion[positive_class, positive_class]))
     True_Neg <- append(True_Neg, as.numeric(confusion[negative_class, negative_class]))
     False_Pos <- append(False_Pos, as.numeric(confusion[negative_class, positive_class]))
     False_Neg <- append(False_Neg, as.numeric(confusion[positive_class, negative_class]))
     accuracy <- (confusion[1, 1] + confusion[2, 2]) / length(predict_proba)
     accuracies <- append(accuracies, as.numeric(accuracy))

     # collect feature importances
     importance_frame <- as.data.frame(myforest$importance)
     colnames(importance_frame) <- counter
     importance_frame$names <- row.names(importance_frame)
     collection <- append(collection, list(importance_frame))
     if (accuracy == 1) {
       collection_importance <- append(collection_importance, list(importance_frame))
     }
  }
# Analysis important features, check if feature combinations have reached
# necessary accuracy and only use them if thats the case, else use all
# performance calculations
if(length(collection_importance) > 0) {
  overall_importance <- purrr::reduce(
    collection_importance,
    function(left, right) {
      dplyr::full_join(left, right, by = "names")
    }
  )
  row.names(overall_importance) <- overall_importance$names
  overall_importance <- overall_importance[, !names(overall_importance) %in% c("names")]
} else{
  overall_importance <- purrr::reduce(
    collection,
    function(left, right) {
      dplyr::full_join(left, right, by = "names")
    }
  )
  row.names(overall_importance) <- overall_importance$names
  overall_importance <- overall_importance[, !names(overall_importance) %in% c("names")]
}
csv_name <- file.path(getwd(), filename_collection)
write.csv(overall_importance, csv_name)
# Calculate coefficients absolute means & sort by decreasing numbers
my_means <- as.data.frame(rowMeans(abs(overall_importance), na.rm = TRUE))
colnames(my_means) <- "Means"
my_means$dummy <- my_means$Means
my_means <- my_means[order(my_means$Means, decreasing =  TRUE), ]
# Write CSV
csv_name <- file.path(getwd(), filename_importance)
write.csv(my_means, csv_name)

  # Combine Metrices in Dataframe
  metrics_all <- data.frame(Acc = as.character(accuracies),
                            TP = as.character(True_Pos),
                            TN = as.character(True_Neg),
                            FP = as.character(False_Pos),
                            FN = as.character(False_Neg),
                            OOB = as.character(OOB))
  # Write CSV
  csv_name <- file.path(getwd(), filename_metrics)
  write.csv(metrics_all, csv_name)
  return(my_means)
}

## CALCULATIONS

## PERFORM RANDOM FOREST ON PCs ##

# read in patient data
nektor <- as.data.frame(read.csv(file = nektor_table,
                                 header = TRUE,
                                 na.string = NA,
                                 sep = ";",
                                 dec = ".",
                                 row.names = 1,
                                 stringsAsFactors = FALSE))
nektor$Patient <- row.names(nektor)
nektor_filtered <- subset(nektor, select = c("Patient", "Response", "Duation_of_stay_days"))
nektor_filtered$response <- nektor_filtered$Response

# filter nektor for UDK, Aufenthaltsdauer > 20, Response != -99 (unknown)
nektor_filtered <- nektor_filtered %>% filter(grepl("UDK", Patient))
nektor_filtered <- nektor_filtered %>% filter(nektor_filtered$Duation_of_stay_days > 20)
nektor_filtered <- nektor_filtered %>% filter(! nektor_filtered$Response == -99 )

# formatting of nektor 
row.names(nektor_filtered) <- nektor_filtered$Patient
nektor_filtered <- subset(nektor_filtered, select = c("response"))

# Read in measurement data
metabolomics_data <- as.data.frame(read.csv(file = input_metabolome,
                                            header = TRUE,
                                            na.string = NA,
                                            sep = ",",
                                            dec = ".",
                                            row.names = 1,
                                            stringsAsFactors = FALSE))

proteomics_data <- as.data.frame(read.csv(file = input_proteome,
                                          header = TRUE,
                                          na.string = NA,
                                          sep = ",",
                                          dec = ".",
                                          row.names = 1,
                                          stringsAsFactors = FALSE))

lipidomics_data <- as.data.frame(read.csv(file = input_lipidome,
                                          header = TRUE,
                                          na.string = NA,
                                          sep = ",",
                                          dec = ".",
                                          row.names = 1,
                                          stringsAsFactors = FALSE))

# Merge data with filtered nektor
metabolomics_merged <- merge(metabolomics_data, nektor_filtered, by = "row.names", all = FALSE)
proteomics_merged <- merge(proteomics_data, nektor_filtered, by = "row.names", all = FALSE)
lipidomics_merged <- merge(lipidomics_data, nektor_filtered, by = "row.names", all = FALSE)

# remove constant columns and set rownames
metabolomics_data <- metabolomics_merged %>% remove_constant(na.rm = TRUE)
row.names(metabolomics_data) <- metabolomics_data$Row.names
metabolomics_data <- metabolomics_data[, -c(1)]

proteomics_data <- proteomics_merged %>% remove_constant(na.rm = TRUE)
row.names(proteomics_data) <- proteomics_data$Row.names
proteomics_data <- proteomics_data[, -c(1)]

lipidomics_data <- lipidomics_merged %>% remove_constant(na.rm = TRUE)
row.names(lipidomics_data) <- lipidomics_data$Row.names
lipidomics_data <- lipidomics_data[, -c(1)]

# set response variables
metabolomics_data <- metabolomics_data %>% 
    mutate(response = recode(response, "1" = "Improver", "2" = "Non-Improver"))

proteomics_data <- proteomics_data %>% 
    mutate(response = recode(response, "1" = "Improver", "2" = "Non-Improver"))

lipidomics_data <- lipidomics_data %>% 
    mutate(response = recode(response, "1" = "Improver", "2" = "Non-Improver"))

# random forest preparation class labels (equal for each label since patients are sorted equally)
response <- metabolomics_data$response

# Random Forest Filtering
filename_importance <- paste("/RandomForest_Importance_Metabolome_own_results.csv", sep = "")
filename_metrics <- paste("/RandomForest_Metrices_Metabolome_own_results.csv", sep = "")
my_means_metabolome <- RANDOM_FOREST(metabolomics_data, response, filename_importance, filename_metrics)

filename_importance <- paste("/RandomForest_Importance_Proteome_own_results.csv", sep = "")
filename_metrics <- paste("/RandomForest_Metrices_Proteome_own_results.csv", sep = "")
my_means_proteome <- RANDOM_FOREST(proteomics_data, response, filename_importance, filename_metrics)

filename_importance <- paste("/RandomForest_Importance_Lipidome_own_results.csv", sep = "")
filename_metrics <- paste("/RandomForest_Metrices_Lipidome_own_results.csv", sep = "")
my_means_lipidome <- RANDOM_FOREST(lipidomics_data, response, filename_importance, filename_metrics)
print("## Execution of Random Forest for feature selection successful ##")

# Execute Random Forest on reduced set
# Metabolome
filename_metrics <- paste("/RandomForest_Metrices_important_features_metabolome_own_results.csv", sep = "")
filename_importance <- paste("/RandomForest_best_feature_combination_metabolome_own_results.csv", sep = "")
filename_collection <- paste("/RandomForest_Collection_Metabolome_own_results.csv", sep = "")
important_features_final <- RANDOM_FOREST_FEATURES(metabolomics_data, response, my_means_metabolome, 20, filename_importance, filename_metrics, filename_collection)

filename_metrics <- paste("/RandomForest_Metrices_important_features_proteome_own_results.csv", sep = "")
filename_importance <- paste("/RandomForest_best_feature_combination_proteome_own_results.csv", sep = "")
filename_collection <- paste("/RandomForest_Collection_Proteome_own_results.csv", sep = "")
important_features_final <- RANDOM_FOREST_FEATURES(proteomics_data, response, my_means_proteome, 20, filename_importance, filename_metrics, filename_collection)

filename_metrics <- paste("/RandomForest_Metrices_important_features_lipidome_own_results.csv", sep = "")
filename_importance <- paste("/RandomForest_best_feature_combination_lipidome_own_results.csv", sep = "")
filename_collection <- paste("/RandomForest_Collection_Lipidome_own_results.csv", sep = "")
important_features_final <- RANDOM_FOREST_FEATURES(lipidomics_data, response, my_means_lipidome, 20, filename_importance, filename_metrics, filename_collection)

print("## Analysis Random Forest successful ##")

# Execute Random Forest on combined set
combined_set <- rbind(my_means_metabolome[1:20, ],
                      my_means_proteome[1:20, ],
                      my_means_lipidome[1:20, ])
combined_features <- cbind(metabolomics_data, 
                           subset(proteomics_data, select = -c(response)),
                           subset(lipidomics_data, select = -c(response)))
feature_set <- subset(combined_features, select = as.character(row.names(combined_set)))

filename_importance <- paste("/RandomForest_Importance_combined_own_results.csv", sep = "")
filename_metrics <- paste("/RandomForest_Metrices_combined_own_results.csv", sep = "")
filename_collection <- paste("RandomForest_Collection_combined_own_results.csv", sep = "")
my_means_combined <- RANDOM_FOREST_FEATURES(feature_set, response, combined_set, 60,filename_importance, filename_metrics, filename_collection)

# add response to feature_set and save it
feature_set$response <- metabolomics_data$response
csv_name <- file.path(getwd(), "/RandomForest_combined_set_own_results.csv")
write.csv(feature_set, csv_name)