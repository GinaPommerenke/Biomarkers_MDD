# plot images for publication
# Author: Gina Pommerenke
# Date: 12.08.2024

## SET VARIABLES
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
input_file <- paste(wd, "\\Data\\Used\\RandomForest_Importance_combined.csv", sep = "")
input_error_file <- paste(wd, "\\Data\\Used\\RandomForest_Collection_combined.csv", sep = "")

input_file_prot <- paste(wd, "\\Data\\Used\\RandomForest_Metrices_important_features_proteome.csv", sep = "")
input_file_met <- paste(wd, "\\Data\\Used\\RandomForest_Metrices_important_features_metabolome.csv", sep = "")
input_file_lip <- paste(wd, "\\Data\\Used\\RandomForest_Metrices_important_features_lipidome.csv", sep = "")
input_file_combined <- paste(wd, "\\Data\\Used\\RandomForest_Metrices_combined.csv", sep = "")

## LIBRARIES
if (!require(ggplot2)) install.packages("ggplot2",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(ggplot2)

if (!require(ggsignif)) install.packages("ggsignif",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(ggsignif)

if (!require(reshape2)) install.packages("reshape2",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(reshape2)

if (!require(dplyr)) install.packages("dplyr",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(dplyr)

if (!require(tidyverse)) install.packages("tidyverse",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(tidyverse)

if (!require(svglite)) install.packages("svglite",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(svglite)

if (!require(scales)) install.packages("scales",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(scales)

## FUNCTIONS 

## CALCULATIONS
df_prot <- as.data.frame(read.csv(file = input_file_prot,
                       header = TRUE,
                       na.string = NA,
                       sep = ",",
                       dec = ".",
                       row.names = 1,
                       stringsAsFactors = FALSE))

df_prot <- df_prot %>% rename(Acc_Proteom = Acc,
                      OOB_Proteom = OOB,
                      Sensitivity_Proteom = TP,
                      Specificity_Proteom = TN)

df_met <- as.data.frame(read.csv(file = input_file_met,
                       header = TRUE,
                       na.string = NA,
                       sep = ",",
                       dec = ".",
                       row.names = 1,
                       stringsAsFactors = FALSE))
df_met <- df_met %>% rename(Acc_Metabolom = Acc,
                      OOB_Metabolom = OOB,
                      Sensitivity_Metabolom = TP,
                      Specificity_Metabolom = TN)

df_lip <- as.data.frame(read.csv(file = input_file_lip,
                       header = TRUE,
                       na.string = NA,
                       sep = ",",
                       dec = ".",
                       row.names = 1,
                       stringsAsFactors = FALSE))

df_lip <- df_lip %>% rename(Acc_Lipidom = Acc,
                      OOB_Lipidom = OOB,
                      Sensitivity_Lipidom = TP,
                      Specificity_Lipidom = TN)

df_combo <- as.data.frame(read.csv(file = input_file_combined,
                       header = TRUE,
                       na.string = NA,
                       sep = ",",
                       dec = ".",
                       row.names = 1,
                       stringsAsFactors = FALSE))

df_combo <- df_combo %>% rename(Acc_Combo = Acc,
                      OOB_Combo = OOB,
                      Sensitivity_Combo = TP,
                      Specificity_Combo = TN)

newframe <- cbind(df_lip, df_met, df_prot, df_combo)

Accuracies <- newframe[, c("Acc_Lipidom",
                           "Acc_Metabolom",
                           "Acc_Proteom",
                           "Acc_Combo")]
Accuracies <- Accuracies %>% rename(Lipidome = Acc_Lipidom,
                                    Metabolome = Acc_Metabolom,
                                    Proteome = Acc_Proteom,
                                    Combination = Acc_Combo)

myAccs <- stack(colMeans(Accuracies))
myVar <- stack(sapply(Accuracies, var))
colnames(myVar) <- c("Variance", "Level")
overall_data <- cbind(myAccs, myVar)

myplot1 <- ggplot(overall_data, aes(x =ind, y =values)) +
       geom_point(aes(colour = ind), size = 5) +
       geom_errorbar(aes(ymin = values - Variance, ymax = values + Variance), width = 0, show.legend = FALSE) +
       labs(x = "Analyzed feature set") +
       ylab(expression(atop("Mean test accuracy" ))) +
       theme(axis.title.x = element_text(size = 30),
             axis.title.y = element_text(size = 30),
             axis.text.x = element_text(size = 24, , angle = 90),
             axis.text = element_text(size = 30),
             legend.text = element_text(size = 24),
             legend.title = element_text(size = 30),
             legend.position = "none",
             plot.title = element_text(size = 40, hjust = 0)) +
        labs(color = "Biological level") +
        scale_color_manual(values = c("red", "green", "blue", "brown"))

myplot1
output_name <- paste("RESULTS/Plot_Accuracies.png", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)

output_name <- paste("RESULTS/Plot_Accuracies.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)

# Calculated importances
df1 <- as.data.frame(read.csv(file = input_file,
                              header = TRUE,
                              na.string = NA,
                              sep = ",",
                              dec = ".",
                              row.names = 1,
                              stringsAsFactors = FALSE))

important_features <- df1[1:20, ]
df1 <- as.data.frame(t(df1))
df1 <- df1 %>% rename("2-C-Methyl-D-erythritol 4-phosphate" = "X2.C.Methyl.D.erythritol4.phosphate",
  "D-alpha-Aminobutyric acid" = "X..aminobutyric.acid",
  "Gamma-Aminobutyric acid" = "X..Aminobutyric.acid",
  "Progabide" =
    "X4....E...4.Chlorophenyl..3.fluoro.6.oxo.2.4.cyclohexadien.1.ylidene.methyl.amino.butanamide",
  "1-pyrroline" = "X1.pyrroline",
  "Ursodeoxycholic acid" = "Ursodeoxycholic.acid",
  "5-hydroxy-4-methoxy- 5,6-dihydro-2H-pyran-2-one" = "X5.hydroxy.4.methoxy.5.6.dihydro.2H.pyran.2.one",
  "4-Amino-2-methyl- 5-phosphomethylpyrimidine" = "X4.Amino.2.methyl.5.phosphomethylpyrimidine",
  "3-(4-Methoxyphenyl)propyl hydrogensulfate" = "X3..4.Methoxyphenyl.propyl.hydrogen.sulfate",
  "Glycochenodeoxycholic acid" = "Glycochenodeoxycholic.acid",
  "4-Hydroxy-5-methylfuran-
  3(2H)-one" = "X4.Hydroxy.5.methylfuran.3.2H..one",
  "Fludarabine phosphate" = "Fludarabine.phosphate",
  "8-Amino-7-oxononanoic acid" = "X8.Amino.7.oxononanoic.acid",
  "Glycoursodeoxycholic acid 3-sulfate" = "Glycoursodeoxycholic.acid.3.sulfate",
  "4-{(Z)-[(4E,7Z,16Z,19Z)-1-Hydroxy-4,7,10,13,16,19-docosahexaen-1-ylidene]amino}butanoic acid" = 
    "X4...Z....4E.7Z.16Z.19Z..1.Hydroxy.4.7.10.13.16.19.docosahexaen.1.ylidene.amino.butanoic.acid",
)
important_features <- as.data.frame(t(df1[, 1:20]))

important_features$names <- row.names(important_features)
important_features$names <- factor(important_features$names, levels = important_features$names)
important_features$roundMeans <- round(important_features$Means, 4)

df_errors <- as.data.frame(read.csv(file = input_error_file,
                                    header = TRUE,
                                    na.string = NA,
                                    sep = ",",
                                    dec = ".",
                                    row.names = 1,
                                    stringsAsFactors = FALSE))
df_errors <- as.data.frame(t(df_errors))
df_errors <- df_errors %>% rename("2-C-Methyl-D-erythritol 4-phosphate" = "X2.C.Methyl.D.erythritol4.phosphate",
  "D-alpha-Aminobutyric acid" = "X..aminobutyric.acid",
  "Gamma-Aminobutyric acid" = "X..Aminobutyric.acid",
  "Progabide" =
    "X4....E...4.Chlorophenyl..3.fluoro.6.oxo.2.4.cyclohexadien.1.ylidene.methyl.amino.butanamide",
  "1-pyrroline" = "X1.pyrroline",
  "Ursodeoxycholic acid" = "Ursodeoxycholic.acid",
  "5-hydroxy-4-methoxy- 5,6-dihydro-2H-pyran-2-one" = "X5.hydroxy.4.methoxy.5.6.dihydro.2H.pyran.2.one",
  "4-Amino-2-methyl- 5-phosphomethylpyrimidine" = "X4.Amino.2.methyl.5.phosphomethylpyrimidine",
  "3-(4-Methoxyphenyl)propyl hydrogensulfate" = "X3..4.Methoxyphenyl.propyl.hydrogen.sulfate",
  "Glycochenodeoxycholic acid" = "Glycochenodeoxycholic.acid",
  "4-Hydroxy-5-methylfuran-
  3(2H)-one" = "X4.Hydroxy.5.methylfuran.3.2H..one",
  "Fludarabine phosphate" = "Fludarabine.phosphate",
  "8-Amino-7-oxononanoic acid" = "X8.Amino.7.oxononanoic.acid",
  "Glycoursodeoxycholic acid 3-sulfate" = "Glycoursodeoxycholic.acid.3.sulfate",
  "4-{(Z)-[(4E,7Z,16Z,19Z)-1-Hydroxy-4,7,10,13,16,19-docosahexaen-1-ylidene]amino}butanoic acid" = 
    "X4...Z....4E.7Z.16Z.19Z..1.Hydroxy.4.7.10.13.16.19.docosahexaen.1.ylidene.amino.butanoic.acid",
)
df_errors <- as.data.frame(t(df_errors))
stdevs <- as.data.frame(apply(df_errors, 1, sd, na.rm = TRUE))
colnames(stdevs) <- c("Stdev")
important_features <- merge(important_features, stdevs, by="row.names", all.x = FALSE)
row.names(important_features) <- important_features$names

important_features <- important_features[order(important_features$roundMeans, decreasing = TRUE), ]
myplot <- ggplot(important_features, aes(x = names, y = Means)) +
       geom_bar(stat = "identity", fill = c("green", "green", "green", "green", "green",
                                            "green", "green", "green", "green", "green",
                                            "blue", "green", "green", "green", "green",
                                            "green", "green", "blue","green", "blue"), alpha = 0.3, position = "dodge") +
       geom_errorbar(aes(ymin = Means - Stdev, ymax = Means + Stdev)) +
       geom_text(aes(y = 0.02, label= roundMeans), position = position_dodge(width = .4), angle = 90, size = 8) +
       theme_minimal() + 
       labs(x = "Analyte",
       y = "Average importance for \nrandom forest prediction") +
       #coord_flip() +
       scale_x_discrete(labels = label_wrap(30)) +
       theme(axis.title.x = element_text(size = 30),
             axis.title.y = element_text(size = 30),
             axis.text.y = element_text(size = 30),
             axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 24),
             plot.title = element_text(size = 40, hjust = -0.37)) 

myplot
output_name <- paste("/Plot_Importances_values.png", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 10,
    width = 21,
    plot = myplot,
    dpi = 1200)

output_name <- paste("/Plot_Importances_values.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 10,
    width = 21,
    plot = myplot,
    dpi = 1200)

# Importance Distribution
my_data <- as.data.frame(read.csv(file = input_file,
                                header = FALSE,
                                na.string = NA,
                                sep = ",",
                                dec = ".",
                                row.names = 1,
                                stringsAsFactors = FALSE))
colnames(my_data) <- c(my_data[1, ])
my_data <- my_data[-c(1), ]

# select datapoints
metabolits <- my_data[c(1:10, 12:17, 19, 22, 23, 25), ]
proteins <- my_data[c(11, 18, 20:21,  24, 26:36, 40:41, 44, 46), ]
lipids <- my_data[c(37:39, 42:43, 45, 47:60), ]

# merge data to one frame
final_data <- data.frame(1:20, 1)
final_data$metabolites <- as.numeric(metabolits$Means)
final_data$proteins <- as.numeric(proteins$Means)
final_data$lipids <- as.numeric(lipids$Means)
final_data <- final_data[,-c(1,2)]

# Create data for plot
plot_data <- rbind(data.frame("Importance" = final_data$metabolites,
                              "Level"= "Metabolites"),
                   data.frame("Importance" = final_data$proteins,
                              "Level" = "Proteins"),
                   data.frame("Importance" = final_data$lipids,
                              "Level" = "Lipids"))

# Create significance for diffrences between the groups by MWU-Test
sigs <- list()
# Improver vs. Non-improver
comparison_metabolome_proteome <- wilcox.test(Importance ~ Level, data = plot_data[! (plot_data$Level %in% "Lipids"), ], 
                                  exact = TRUE, correct = FALSE, conf.int = FALSE)
sigs <- append(sigs, paste("P = ", as.character(format(comparison_metabolome_proteome$p.value, digits = 3)), sep = ""))

# Improver vs. Healthy Control
comparison_metabolome_lipidome <- wilcox.test(Importance ~ Level, data = plot_data[! (plot_data$Level %in% "Proteins"), ], 
                                  exact = TRUE, correct = FALSE, conf.int = FALSE)
sigs <- append(sigs, paste("P = ", as.character(format(comparison_metabolome_lipidome$p.value, digits = 3)), sep = ""))

# Non-improver vs. Healthy Control
comparison_proteome_lipidome <- wilcox.test(Importance ~ Level, data = plot_data[! (plot_data$Level %in% "Metabolites"), ], 
                                  exact = TRUE, correct = FALSE, conf.int = FALSE)
sigs <- append(sigs, paste("P = ", as.character(format(comparison_proteome_lipidome$p.value, digits = 3)), sep = ""))

# Plot significance distribution with P-values
plot <- ggplot(plot_data, aes(x=Level, y = Importance, fill = Level)) +
  geom_violin() +
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(width = 0.75), 
             dotsize = 100, binwidth = 0.00001) +
  geom_signif(y_position = c(0.129, 0.131, 0.138),
          xmin = c(2, 1, 1),
          xmax = c(3, 2, 3),
          annotation = unlist(sigs), tip_length = 0.01,
          size = 1, textsize = 7)+
  labs(y = "Feature importance", x = "Biological level", fill = "Biological Level") +
  scale_fill_discrete(guide = "none") +
  theme(axis.text = element_text(hjust = 1, size = 30),
      axis.title = element_text(size = 30),
      axis.text.x = element_blank(),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.position = "none",
      plot.title = element_text(size = 40, hjust = -0.05))
plot

# save plots
output_name <- paste("/RESULTS/Plot_Importance_distribution.png", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    plot = plot,
    height = 7,
    width = 7,
    dpi = 1200)

output_name <- paste("/RESULTS/Plot_Importance_distribution.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    plot = plot,
    height = 7,
    width = 7,
    dpi = 1200)

# Plots for the supplemental

# Sensitivity
Sensitivities <- newframe[, c("Sensitivity_Lipidom",
                           "Sensitivity_Metabolom",
                           "Sensitivity_Proteom",
                           "Sensitivity_Combo")]
Sensitivities <- Sensitivities %>% rename(Lipidome = Sensitivity_Lipidom,
                                    Metabolome = Sensitivity_Metabolom,
                                    Proteome = Sensitivity_Proteom,
                                    Combination = Sensitivity_Combo)

mySens <- stack(colMeans(Sensitivities))
myVars <- stack(sapply(Sensitivities, var))
colnames(myVars) <- c("Variance", "Level")
overall_data <- cbind(myAccs, myVars)

myplot1 <- ggplot(overall_data, aes(x =ind, y =values)) +
       geom_point(aes(colour = ind), size = 5) +
       geom_errorbar(aes(ymin = values - Variance, ymax = values + Variance), width = 0, show.legend = TRUE) +
       labs(x = "Analyzed feature set") +
       ylab(expression(atop("Mean test sensitivity" ))) +
       theme(axis.title.x = element_text(size = 30),
             axis.title.y = element_text(size = 30),
             axis.text.x = element_text(size = 24, , angle = 90),
             axis.text = element_text(size = 30),
             legend.text = element_text(size = 24),
             legend.title = element_text(size = 30),
             legend.position = "none",
             plot.title = element_text(size = 40, hjust = 0)) +
        labs(color = "Biological level") +
        scale_color_manual(values = c("red", "green", "blue", "brown"))

myplot1
output_name <- paste("RESULTS/Plot_Sensitivity.png", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)

output_name <- paste("RESULTS/Plot_Sensitivity.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)

# Specificity
Specificity <- newframe[, c("Specificity_Lipidom",
                           "Specificity_Metabolom",
                           "Specificity_Proteom",
                           "Specificity_Combo")]
Specificity <- Specificity %>% rename(Lipidome = Specificity_Lipidom,
                                    Metabolome = Specificity_Metabolom,
                                    Proteome = Specificity_Proteom,
                                    Combination = Specificity_Combo)

mySpec <- stack(colMeans(Specificity))
myVars <- stack(sapply(Specificity, var))
colnames(myVars) <- c("Variance", "Level")
overall_data <- cbind(mySpec, myVars)

myplot1 <- ggplot(overall_data, aes(x =ind, y =values)) +
       geom_point(aes(colour = ind), size = 5) +
       geom_errorbar(aes(ymin = values - Variance, ymax = values + Variance), width = 0, show.legend = TRUE) +
       labs(x = "Analyzed feature set") +
       ylab(expression(atop("Mean test specificity" ))) +
       theme(axis.title.x = element_text(size = 30),
             axis.title.y = element_text(size = 30),
             axis.text.x = element_text(size = 24, , angle = 90),
             axis.text = element_text(size = 30),
             legend.text = element_text(size = 24),
             legend.title = element_text(size = 30),
             legend.position = "none",
             plot.title = element_text(size = 40, hjust = 0)) +
        labs(color = "Biological level") +
        scale_color_manual(values = c("red", "green", "blue", "brown"))

myplot1
output_name <- paste("RESULTS/Plot_Specificity.png", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)

output_name <- paste("RESULTS/Plot_Specificity.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)

# OOB-Error
OOB_error <- newframe[, c("OOB_Lipidom",
                           "OOB_Metabolom",
                           "OOB_Proteom",
                           "OOB_Combo")]
OOB_error <- OOB_error %>% rename(Lipidome = OOB_Lipidom,
                                    Metabolome = OOB_Metabolom,
                                    Proteome = OOB_Proteom,
                                    Combination = OOB_Combo)

myOOB <- stack(colMeans(OOB_error, na.rm = TRUE))
myVars <- stack(sapply(OOB_error, var, na.rm = TRUE))
colnames(myVars) <- c("Variance", "Level")
overall_data <- cbind(myOOB, myVars)

myplot1 <- ggplot(overall_data, aes(x =ind, y =values)) +
       geom_point(aes(colour = ind), size = 5) +
       geom_errorbar(aes(ymin = values - Variance, ymax = values + Variance), width = 0) +
       labs(x = "Analyzed feature set") +
       ylab(expression(atop("Mean OOB-error" ))) +
       theme(axis.title.x = element_text(size = 30),
             axis.title.y = element_text(size = 30),
             axis.text.x = element_text(size = 24, , angle = 90),
             axis.text = element_text(size = 30),
             legend.text = element_text(size = 24),
             legend.title = element_text(size = 30),
             legend.position = "none",
             plot.title = element_text(size = 40, hjust = 0)) +
        labs(color = "Biological level") +
        scale_color_manual(values = c("red", "green", "blue", "brown"))

myplot1
output_name <- paste("RESULTS/Plot_OOB.png", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)

output_name <- paste("RESULTS/Plot_OOB.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    height = 7,
    width = 7,
    plot = myplot1,
    dpi = 1200)
