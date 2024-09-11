## RandomForest in R
## Author: Gina Pommerenke
## Date: 06.03.2024
## Random Forest with feature selection for both, data after PCA and complete original data, executable for
## Metabolomics, Proteomics and Lipidomics

## INPUT
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
nektor_table <- paste(wd, "\\Data\\Used\\NEKTOR_publication.csv", sep = "")
positive_class <- "Improver" # Improvement after treatment
negative_class <- "Non-Improver" # No improvement after treatment

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

if (!require(svglite)) install.packages("svglite",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(svglite)

## FUNCTIONS

# Simple Random Forest

RANDOM_FOREST <- function(frame, response, filename_metrics) {
  
    # Prepare dataset for random splits
   data_rf_features <- frame

   # Intermediate Lists for results of splits
   counter <- 0
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
  }
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
  return(metrics_all)
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
nektor_filtered <- subset(nektor, select = c("Patient", "Response", "Duation_of_stay_days", "biol_sex_m_f", "age_years", "routinex_bmi"))
nektor_filtered$response <- nektor_filtered$Response

# filter nektor for UDK, Aufenthaltsdauer > 20, Response != -99 (unknown)
nektor_filtered <- nektor_filtered %>% filter(grepl("UDK", Patient))
nektor_filtered <- nektor_filtered %>% filter(nektor_filtered$Duation_of_stay_days > 20)
nektor_filtered <- nektor_filtered %>% filter(! nektor_filtered$Response == -99 )

# formatting of nektor 
row.names(nektor_filtered) <- nektor_filtered$Patient
nektor_filtered <- subset(nektor_filtered, select = c("response", "biol_sex_m_f", "age_years", "routinex_bmi"))

# set response variables
nektor_filtered <- nektor_filtered %>% 
    mutate(response = recode(response, "1" = "Improver", "2" = "Non-Improver"))

# random forest preparation class labels (equal for each label since patients are sorted equally)
response <- nektor_filtered$response

# Random Forest Analysis
filename_metrics <- paste("/RandomForest_Metrices_baseline.csv", sep = "")
my_means <- RANDOM_FOREST(nektor_filtered, response, filename_metrics)

# Plot results

my_means$Sensitivity <- my_means$TP # valid since exactly 1 positive in test set
my_means$Specificity <- my_means$TN # valid since exactly 1 negative in test set
Metrices <- my_means[, c("Acc",
                           "Sensitivity",
                           "Specificity",
                           "OOB")]
Metrices <- Metrices %>% rename(Accuracy = Acc)
Metrices <- as.data.frame(sapply(Metrices, as.numeric))

myRes <- stack(colMeans(Metrices, na.rm = TRUE))
myVar <- stack(sapply(Metrices, var, na.rm = TRUE))
colnames(myVar) <- c("Variance", "Level")
overall_data <- cbind(myRes, myVar)

myplot1 <- ggplot(overall_data, aes(x =ind, y =values)) +
       geom_point(aes(colour = ind), size = 5) +
       geom_errorbar(aes(ymin = values - Variance, ymax = values + Variance), width = 0, show.legend = FALSE) +
       labs(x = "Performance metrics") +
       ylab(expression(atop("Mean test performance" ))) +
       theme(axis.title.x = element_text(size = 30),
             axis.title.y = element_text(size = 30),
             axis.text.x = element_text(size = 24, , angle = 90),
             axis.text = element_text(size = 30),
             legend.text = element_text(size = 24),
             legend.title = element_text(size = 30),
             legend.position = "none",
             plot.title = element_text(size = 40, hjust = 0)) +
        labs(color = "Biological level") +
        scale_color_manual(values = c("#ff7300", "#00b7ff", "#d400ff", "#5a2121"))

myplot1
output_name <- paste("RESULTS/Plot_Baseline_Performance.png", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)

output_name <- paste("RESULTS/Plot_Baseline_Performance.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)
