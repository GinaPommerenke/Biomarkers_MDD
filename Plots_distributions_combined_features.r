## Create Plot with the separating analytes and theri z-scored concentrations
## 11.03.2024
## Gina Pommerenke

## INPUT
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
input_frame <- paste(wd, "\\Data\\Used\\RandomForest_combined_set.csv", sep = "")


## LIBRARIES
if (!require(ggplot2)) install.packages("ggplot2",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(ggplot2)

if (!require(dplyr)) install.packages("dplyr",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(dplyr)

if (!require(tidyr)) install.packages("tidyr",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(tidyr)

if (!require(scales)) install.packages("scales",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(scales)

## FUNCTIONS

## CALCULATIONS
# Read in dataframe
myframe <- as.data.frame(read.csv(file = input_frame,
                                  header = TRUE,
                                  na.string = NA,
                                  sep = ",",
                                  dec = ".",
                                  row.names = 1,
                                  stringsAsFactors = FALSE))
# Put response column in first place
myframe <- myframe %>% relocate(response, .before = 1)

# Rename ASCI-conform
myframe <- myframe %>% rename("2-C-Methyl-D-erythritol 4-phosphate" = "X2.C.Methyl.D.erythritol4.phosphate",
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
  "4-{(Z)-[(4E,7Z,16Z,19Z)-1-Hydroxy-4,7,10,13,16,19-\ndocosahexaen-1-ylidene]amino}butanoic acid" = 
    "X4...Z....4E.7Z.16Z.19Z..1.Hydroxy.4.7.10.13.16.19.docosahexaen.1.ylidene.amino.butanoic.acid",
  "Improvement" = "response"
)

# Create Dotplot of all sixty analytes
# Create a Boxplot for all plots of Top-60 features
for (i in 2:ncol(myframe)) {
  myplot <- ggplot(myframe, aes(y = myframe[[i]])) +
    geom_point(aes(x = Improvement, colour = Improvement),  size = 3.5) +
    labs(y= "Concentration (z-scores)",
         title = paste(colnames(myframe)[i], sep = ""),
         x = "",
         colour = "Patient group") +
    theme(
      text = element_text(size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = c("purple", "orange"))
  output_name <- paste("/Dotplots_Analytes/", colnames(myframe)[i], ".png", sep = "")
  plot_name <- file.path(getwd(), output_name)
  # save plot
  ggsave(plot_name,
      plot = myplot,
      dpi = 800)
}

# Resort Frame and select only separating analytes
myframe <- myframe[, c("Improvement",
                       "1-pyrroline",
                       "2-C-Methyl-D-erythritol 4-phosphate",
                       "3-(4-Methoxyphenyl)propyl hydrogensulfate",
                       "Progabide",
                       "5-hydroxy-4-methoxy- 5,6-dihydro-2H-pyran-2-one",
                       "Acetylcarnitine",
                       "D-alpha-Aminobutyric acid",
                       "Gamma-Aminobutyric acid",
                       "4-Amino-2-methyl- 5-phosphomethylpyrimidine",
                       "Ursodeoxycholic acid")]

# Create a plot for all separating analytes features
df <- myframe %>% pivot_longer(-Improvement, names_to = "Analyte")
myplot <- ggplot(df, aes(x=Analyte, y = value, colour = Improvement)) +
  geom_point(aes(colour = Improvement),  size = 3) +
  labs(y= "Concentration (z-scores)",
       x = "Metabolite",
       colour = "Patient group") +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20),
    legend.position = "bottom",
    axis.title = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  scale_color_manual(values = c("purple", "orange")) +
  coord_flip() +
  scale_x_discrete(labels = label_wrap(15))

output_name <- paste("/RESULTS/Combined_Separating_Analytes_Plots.png", sep = "")
plot_name <- file.path(getwd(), output_name)
# save plot
ggsave(plot_name,
    height = 8,
    width = 10,
    plot = myplot,
    dpi = 1200)

output_name <- paste("/RESULTS/Combined_Separating_Analytes_Plots.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
# save plot
ggsave(plot_name,
    height = 8,
    width = 10,
    plot = myplot,
    dpi = 1200)
