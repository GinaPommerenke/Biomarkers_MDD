## RandomForest: Shuffle labels of 
## Author: Gina Pommerenke
## Date: 06.03.2024
## Random Forest with feature importance for combination of most important features after 
## feature selection with random forest

## INPUT
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
input_file <- paste(wd, "\\Data\\Used\\RandomForest_combined_set.csv", sep = "")
positive_class <- "Improver"
negative_class <- "Non-Improver"

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

if (!require(tibble)) install.packages("tibble",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(tibble)

## FUNCTIONS

# Random Forest for Top60 Features
RANDOM_FOREST <- function(frame, response, path_metrics, path_importance) {
  # Random forest prepare output parameters
  counter <- 0
  collection <- list()
  mean_accuracy_PCs <- list()
  mean_True_Pos <- list()
  mean_True_Neg <- list()
  mean_False_Pos <- list()
  mean_False_Neg <- list()
  mean_OOB <- list()
  for (shuffle in 1:200) {
    if (shuffle %% 50 == 0) {print(shuffle)}
    # Random Forest with 200 Random Splits
    response <- sample(response)
    accuracies <- list()
    True_Pos <- list()
    True_Neg <- list()
    False_Pos <- list()
    False_Neg <- list()
    OOB <- list()
    for (cross in 1:10) {
      counter <- counter + 1
      # create specific test- and trainingset w.r.t. class labels
      test_y <- NULL
      # create specific test- and trainingset w.r.t. class labels
      while (! ("Improver" %in% test_y && "Non-Improver" %in% test_y)) {
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
      # predict class probabilities
      predict_proba <- predict(myforest, newdata = test)
      confusion <- table(Predictions = predict_proba, TrueLabels = test_y)
      True_Pos <- append(True_Pos, as.numeric(confusion[positive_class, positive_class]))
      True_Neg <- append(True_Neg, as.numeric(confusion[negative_class, negative_class]))
      False_Pos <- append(False_Pos, as.numeric(confusion[negative_class, positive_class]))
      False_Neg <- append(False_Neg, as.numeric(confusion[positive_class, negative_class]))
      accuracy <- as.numeric((confusion[1, 1] + confusion[2, 2]) / length(predict_proba))
      accuracies <- append(accuracies, accuracy)
    
      importance_frame <- as.data.frame(myforest$importance)
      colnames(importance_frame) <- counter
      importance_frame$names <- row.names(importance_frame)
      collection <- append(collection, list(importance_frame))
    }
    # Calculation of Mean values after Cross validation
    mean_accuracy_PCs <- append(mean_accuracy_PCs, mean(as.numeric(unlist(accuracies))))
    mean_True_Pos <- append(mean_True_Pos, mean(as.numeric(unlist(True_Pos))))
    mean_True_Neg <- append(mean_True_Neg, mean(as.numeric(unlist(True_Neg))))
    mean_False_Pos <- append(mean_False_Pos, mean(as.numeric(unlist(False_Pos))))
    mean_False_Neg <- append(mean_False_Neg, mean(as.numeric(unlist(False_Neg))))
    mean_OOB <- append(mean_OOB, mean(as.numeric(unlist(OOB[! is.na(OOB)]))))
  }

  # Combine Metrices in Dataframe
  metrics_PCA <- tibble(Acc = as.character(mean_accuracy_PCs),
                            TP = as.character(mean_True_Pos),
                            TN = as.character(mean_True_Neg),
                            FP = as.character(mean_False_Pos),
                            FN = as.character(mean_False_Neg),
                            OOB = as.character(mean_OOB))

  # Analysis important features
  overall_importance <- purrr::reduce(
    collection,
    function(left, right) {
      dplyr::full_join(left, right, by = "names")
    }
  )
  row.names(overall_importance) <- overall_importance$names
  overall_importance <- overall_importance[, !names(overall_importance) %in% c("names")]
  # Calculate coefficients absolute means & sort by decreasing numbers
  my_means <- as.data.frame(rowMeans(abs(overall_importance), na.rm = TRUE))
  colnames(my_means) <- "Means"
  my_means$dummy <- my_means$Means
  my_means <- my_means[order(my_means$Means, decreasing =  TRUE), ]

  return(metrics_PCA)
}

## CALCULATIONS

# Input Metabolom
input_data <- as.data.frame(read.csv(file = input_file,
                                            header = TRUE,
                                            na.string = NA,
                                            sep = ",",
                                            dec = ".",
                                            row.names = 1,
                                            stringsAsFactors = FALSE))

# Add response from Nektor to metabolomics dataset
input_data <- input_data %>% 
    mutate(response = recode(response, "1" = "Improver", "2" = "Non-Improver"))
# random forest preparation class labels
response <- input_data$response

# Execution of Random Forest 

output <- as.data.frame(RANDOM_FOREST(input_data, response))
combi_df <- output
for (i in 1:4) {
  output <- as.data.frame(RANDOM_FOREST(input_data, response))
  combi_df <- rbind(combi_df, output)
}

# 95%-Percentil and FWER-estimate
ninetyfive_percentil <- quantile(as.numeric(combi_df$Acc), probs = 0.95)
FWER <- function(vector){
  # Input ordered vector with entries between 0 and 1
  for (i in 1:length(vector)){
    if (vector[i] == 1){ 
      FWER = (length(vector) - (i-1))/length(vector)
      return(FWER)
      break
    } else if (i == length(vector)) {
       FWER = (length(vector) - (i))/length(vector)
       return(FWER)
    }
  }
}
combi_df_ordered <- combi_df[order(combi_df$Acc), ]
ordered_acc <- as.numeric(combi_df_ordered$Acc)
FWER_all <- FWER(ordered_acc)
results <- data.frame(Features = "all",
                      FWER = FWER_all)


# Write csv of collected results
path_errors <- paste("/RandomForest_Combined_shuffle_own_results.csv", sep = "")
csv_name <- file.path(getwd(), path_errors)
write.csv(combi_df_ordered, csv_name)

# Plot of Distributions of Acc. after shuffeling
# plot of density distribution
annotation <- data.frame(
  x = c(0.3),
  y = c(30),
  label = c(paste("P = ", round(as.numeric(results[1, 2]), 7), sep = ""))
)
p1 <- ggplot(combi_df, aes(x = as.numeric(ACC))) +
  geom_histogram(aes(x = as.numeric(Acc)), binwidth = 0.05, alpha = 0.3, color = "white", fill = "darkgrey") +
  geom_density(aes(x = as.numeric(Acc), y = ..count..*0.05, fill = "black"), alpha = 0, linewidth = 1.25) +
  labs(y = "Number of runs \n with \nshuffled phenotype", 
       x = "Test accuracy for \nshuffled phenotypes", color = "Density distribution") +
  geom_segment(aes(x= 1, xend = 1, y = 0, yend = 100), size = 4) +
  theme(legend.position = "None",
        axis.text = element_text(hjust = 1, size = 30),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 40, hjust = -0.05)) +
  geom_text(data = annotation, aes(x=x, y=y, label=label),
            size = 7) +
  geom_text(label = "Test accuracy \nfor original\nphenotypes", aes(x = 0.77, y = 80), size = 7) 

p1
# save plot
output_name <- paste("/Plot_Accuracy_Distribution_Shuffled_Labels_own_plot.png", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
       plot = p1,
       height = 7,
       width = 7,
       dpi = 1200)

output_name <- paste("/Plot_Accuracy_Distribution_Shuffled_Labels_own_plot.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
       plot = p1,
       height = 7,
       width = 7,
       dpi = 1200)
