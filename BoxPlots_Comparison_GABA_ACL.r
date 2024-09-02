## Create BoxPlots of GABA and ACL with comparison with healthy control group
## 11.03.2024
## Gina Pommerenke

## INPUT
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
input_frame <- paste(wd, "\\Data\\Used\\RandomForest_combined_set.csv", sep = "")
input_frame2 <- paste(wd, "\\Data\\Used\\Metabolome_shortened.csv", sep = "")
#input_frame <- "C:\\Users\\Gina\\OneDrive - BioVariance GmbH\\P23008_PhD\\04_Reports\\Code_for_Review\\Data\\RandomForest_combined_set.csv"
#input_frame2 <- "C:\\Users\\Gina\\OneDrive - BioVariance GmbH\\P23008_PhD\\04_Reports\\Code_for_Review\\Data\\Metabolome_shortened.csv"

## LIBRARIES
if (!require(ggplot2)) install.packages("ggplot2",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(ggplot2)

if (!require(ggsignif)) install.packages("ggsignif",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(ggsignif)

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
myframe <- as.data.frame(read.csv(file = input_frame,
                                  header = TRUE,
                                  na.string = NA,
                                  sep = ",",
                                  dec = ".",
                                  row.names = 1,
                                  stringsAsFactors = FALSE))
# Put response column in first place
myframe <- myframe %>% relocate(response, .before = 1)

# Put response column in first place
myframe <- myframe %>% relocate(response, .before = 1)

# Input of TRDs
myframe2 <- as.data.frame(read.csv(file = input_frame2,
                                  header = TRUE,
                                  na.string = NA,
                                  sep = ",",
                                  dec = ".",
                                  row.names = 1,
                                  stringsAsFactors = FALSE))
myframe_Healthy <- myframe2[1:20, ]
myframe_Healthy$response <- "Healthy"

myframe[setdiff(names(myframe_Healthy), names(myframe))] <- NA
myframe_Healthy[setdiff(names(myframe), names(myframe_Healthy))] <- NA

overall_frame <- rbind(myframe, myframe_Healthy)
filtered_frame <- overall_frame[, colSums(is.na(overall_frame))== 0]
# rename with shorthand-name in publication
filtered_frame <- filtered_frame %>% rename("GABA" = "X..Aminobutyric.acid",
                                            "ACL" = "Acetylcarnitine",
                                            "Response" = "response")

# remove all columns except Response, ACL and GABA, reformat response column                                            )
filtered_frame <- filtered_frame[, c("Response",
                                     "ACL",
                                     "GABA")]
#filtered_frame$Response[filtered_frame$Response == "Responder"] <- "Improver"
#filtered_frame$Response[filtered_frame$Response == "Non-Responder"] <- "Non-improver"
filtered_frame$Response[filtered_frame$Response == "Healthy"] <- "Healthy control"

# Create a Boxplot GABA
GABA_frame <- filtered_frame[, c("Response", "GABA")]
GABA_frame_new <- GABA_frame %>% pivot_longer(-Response, names_to = "Analyte")

# Create significance for diffrences between the groups by MWU-Test
sigs <- list()
# Improver vs. Non-improver
comparison_i_non_i <- wilcox.test(GABA ~ Response, data = GABA_frame[! (GABA_frame$Response %in% "Healthy control"), ], 
                                  exact = TRUE, correct = FALSE, conf.int = FALSE)
sigs <- append(sigs, paste("P = ", as.character(round(comparison_i_non_i$p.value, digits = 3)), sep = ""))

# Improver vs. Healthy Control
comparison_i_hc <- wilcox.test(GABA ~ Response, data = GABA_frame[! (GABA_frame$Response %in% "Non-Improver"), ], 
                                  exact = TRUE, correct = FALSE, conf.int = FALSE)
sigs <- append(sigs, paste("P = ", as.character(round(comparison_i_hc$p.value, digits = 3)), sep = ""))

# Non-improver vs. Healthy Control
comparison_non_i_hc <- wilcox.test(GABA ~ Response, data = GABA_frame[! (GABA_frame$Response %in% "Improver"), ], 
                                  exact = TRUE, correct = FALSE, conf.int = FALSE)
sigs <- append(sigs, paste("P = ", as.character(round(comparison_non_i_hc$p.value, digits = 3)), sep = ""))

myplot <- ggplot(GABA_frame_new, aes(x=Analyte, y = value, fill = Response)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(width = 0.75), 
               dotsize = 45, binwidth = 0.001) +
  geom_signif(y_position = c(1.62, 1.67, 1.82),
              xmin = c(1, 0.75, 0.75),
              xmax = c(1.25, 1, 1.25),
              annotation = unlist(sigs), tip_length = 0.005,
              size = 1, textsize = 7)+
  labs(y= "Concentration (z-scores)",
       x = "GABA",
       fill = "Patient group") +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.position = "bottom",
    axis.title = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.text.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_color_manual(values = c("black", "black", "black")) +
  scale_fill_manual(values = c("grey", "purple", "orange")) +
  scale_x_discrete(labels = label_wrap(15))


output_name <- paste("/RESULTS/Boxplot_GABA.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
# save plot
ggsave(plot_name,
    height = 7,
    width = 8,
    plot = myplot,
    dpi = 800)

output_name <- paste("/RESULTS/Boxplot_GABA.png", sep = "")
plot_name <- file.path(getwd(), output_name)
# save plot
ggsave(plot_name,
    height = 7,
    width = 8,
    plot = myplot,
    dpi = 800)

# Create Boxplot ACL
ACL_frame <- filtered_frame[, c("Response", "ACL")]
ACL_frame_new <- ACL_frame %>% pivot_longer(-Response, names_to = "Analyte")

# Create significance for diffrences between the groups by MWU-Test
sigs <- list()
# Improver vs. Non-improver
comparison_i_non_i <- wilcox.test(ACL ~ Response, data = ACL_frame[! (ACL_frame$Response %in% "Healthy control"), ], 
                                  exact = TRUE, correct = FALSE, conf.int = FALSE)
sigs <- append(sigs, paste("P = ", as.character(round(comparison_i_non_i$p.value, digits = 3)), sep = ""))

# Improver vs. Healthy Control
comparison_i_hc <- wilcox.test(ACL ~ Response, data = ACL_frame[! (ACL_frame$Response %in% "Non-Improver"), ], 
                                  exact = TRUE, correct = FALSE, conf.int = FALSE)
sigs <- append(sigs, paste("P = ", as.character(round(comparison_i_hc$p.value, digits = 3)), sep = ""))

# Non-improver vs. Healthy Control
comparison_non_i_hc <- wilcox.test(ACL ~ Response, data = ACL_frame[! (ACL_frame$Response %in% "Improver"), ], 
                                  exact = TRUE, correct = FALSE, conf.int = FALSE)
sigs <- append(sigs, paste("P = ", as.character(round(comparison_non_i_hc$p.value, digits = 3)), sep = ""))

myplot <- ggplot(ACL_frame_new, aes(x=Analyte, y = value, fill = Response)) +
  geom_boxplot() +
  geom_dotplot(binaxis = "y", stackdir = "center", position = position_dodge(width = 0.75), 
               dotsize = 45, binwidth = 0.001) +
  geom_signif(y_position = c(1.3, 1.35, 1.53),
            xmin = c(1, 0.75, 0.75),
            xmax = c(1.25, 1, 1.25),
            annotation = unlist(sigs), tip_length = 0.005,
            size = 1, textsize = 7)+
  labs(y= "Concentration (z-scores)",
       x = "ACL",
       fill = "Patient group") +
  theme(
    text = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.position = "bottom",
    axis.title = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.text.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_color_manual(values = c("black", "black", "black")) +
  scale_fill_manual(values = c("grey", "purple", "orange")) +
  scale_x_discrete(labels = label_wrap(15))

output_name <- paste("/RESULTS/Boxplot_ACL.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
# save plot
ggsave(plot_name,
    height = 7,
    width = 8,
    plot = myplot,
    dpi = 800)

output_name <- paste("/RESULTS/Boxplot_ACL.png", sep = "")
plot_name <- file.path(getwd(), output_name)
# save plot
ggsave(plot_name,
    height = 7,
    width = 8,
    plot = myplot,
    dpi = 800)