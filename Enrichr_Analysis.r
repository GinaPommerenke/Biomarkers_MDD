# EnrichR analysis of proteins associated with separating metabolites and associated plot
# 18.07.2024
# Gina Pommerenke

## INPUT 
wd <- getwd()
wd <- gsub("/", "\\", fixed = TRUE, wd)
input_file <- paste(wd, "\\Data\\Used\\Protein_collection_with_gene_name.csv", sep = "")

## PACKAGES
if (!require(enrichR)) install.packages("enrichR",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(enrichR)

if (!require(ggplot2)) install.packages("ggplot2",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(ggplot2)

if (!require(scales)) install.packages("scales",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(scales)

if (!require(svglite)) install.packages("svglite",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(svglite)

if (!require(stringr)) install.packages("stringr",
            repos = "https://packages.othr.de/cran/", dependencies = TRUE)
library(stringr)

## FUNCTIONS

## CALCULATIONS
input_data <- as.data.frame(read.csv(file = input_file,
                                     header = TRUE,
                                     na.string = NA,
                                     sep = ";",
                                     dec = ".",
                                     row.names = 1,
                                     stringsAsFactors = FALSE))

# EnrichR analysis
dbs <- listEnrichrDbs()
my_databases <- c("KEGG_2013",
                  "GO_Cellular_Component_2015",
                  "GO_Biological_Process_2015",
                  "GO_Molecular_Function_2015")
enriched_set <- enrichr(input_data$To, my_databases)

# extract results from enrichment analysis
geneset_KEGG <- as.data.frame(enriched_set$KEGG_2013)
geneset_KEGG$database <- "KEGG"

geneset_GO_MF <- as.data.frame(enriched_set$GO_Molecular_Function_2015)
geneset_GO_MF$database <- "GO_Moleculcar_Function"

geneset_GO_CC <- as.data.frame(enriched_set$GO_Cellular_Component_2015)
geneset_GO_CC$database <- "GO_Cellular_Component"

geneset_GO_BP <- as.data.frame(enriched_set$GO_Biological_Process_2015)
geneset_GO_BP$database <- "GO_Biological_Process"

# combine results in one dataframe
whole_enrichment <- rbind(geneset_KEGG, geneset_GO_BP, geneset_GO_CC, geneset_GO_MF)

# extract number of genes from query in each pathway and number of genes of each pathway
geneset_split <- as.data.frame(str_split_fixed(whole_enrichment$Overlap, "/", 2))
geneset_split <- as.data.frame(sapply(geneset_split, as.numeric))
colnames(geneset_split) <- c("intersection_size", "term_size")

whole_enrichment <- cbind(whole_enrichment, geneset_split)

# filter pathways for 10 < termsize < 250, at least 25 % of genes in term part of protein set
filtered_pathways <- subset(whole_enrichment, 10 < term_size &
                                      250 > term_size &
                                      0.25 <= intersection_size/term_size)

plot <- plotEnrich(filtered_pathways, y = "Count", orderBy = "P.value") +
          scale_x_discrete(labels = label_wrap(30)) +
          labs(title = "") +
          theme(
                text = element_text(size = 24),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 20),
                axis.title = element_text(size = 20),
               )
plot

# save plot
output_name <- paste("/RESULTS/EnrichR_pathway_enrichment.svg", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    plot = plot,
    height = 10, 
    width = 10,
    dpi = 1200)

output_name <- paste("/RESULTS/EnrichR_pathway_enrichment.png", sep = "")
plot_name <- file.path(getwd(), output_name)
ggsave(plot_name,
    plot = plot,
    height = 10, 
    width = 10,
    dpi = 1200)