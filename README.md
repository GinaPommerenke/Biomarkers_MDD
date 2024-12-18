# Blood_based_biomarerks_for_prediction_in_early_improvment_of_MDD 
Code for the publication "Blood-based metabolomics and explainable machine learing reveals candidate biomarkers for prediction of early improvement in major depressive disorder"

The given files can be used to reconstruct the calculations of the given results as well as the plots relevant for the biological consequences. Since the random forest analysis as well as the regression filters rely on stochastical effects, the results of the own analysis can differ after running the files various times. So, we suggset to use the input files given in the repository. 

The code is written in the programming language R. All used packages are included in the code and will be automatically installed. The R version used is R 4.1.3. The code can be used by downloading the Repository and using it in your standard R-environment. No adaptions for the input are necessary. Since some of the results are dependent on stochastic effects, we suggest to use the files given in the repository. 

# Code files

The files include the coding for the results of the uppermentioned publication. Files are
1. Preprocessing_Metabolome_data.r
  This file contains the preprocessing of the metabolome data.
2. Preprocessing_Proteome_data.r
   This file contains the preprocessing of the metabolome data.
3. Regression_Filter.r
   Linear and logistic regression of preprocessed metabolome and proteome data as well as of the original lipidome data. The files are filtered for three confounders (Age, BMI and sex) to remove analytes which correlate highly with one of the confounders.     Due to stochastic effects, the results of the filtering vary slightly and the extracted final set of analytes differ between the runs. It is recommended to use the original datafile provided in the repository to continue with the analysis.
4. RandomForest.r
   Random forest analysis of the three datasets after the regression filter. First, the twenty most important datafiles of each of the biological levels were extracted by extraction of the feature importances by the built-in method of the random forst        package. Afterwards the extracted features were used for a joint analysis to recieve the overall performance of the model and identify the most promising biomarker candidates. Since the extracted twenty most important features vary due to random            effects, we recommend to use the data provided for all downstream analysis. 
5. RandomForest_combined_shuffle.r
  Analysis, if the performance results of the joint RandomForest analysis is better then random by random shuffeling of the target variable and consecutive calculation of the p-values.
6. Plots_random_forest.r
   Visualization of the random forest results as well as the different random forest importances. Results are the plots shown in the publication.
7. Plots_distributions_combined_features.r
   Visualization of the distribution of the z-score preprocessed analytes of the joint analysis by dotplots as well as a combined plot of the eleven metabolites, which split the given data in early improver and non-improver as joint dotplot.
8. BoxPlots_Comparison_GABA_ACL.r
   Box plots of the concentrations of GABA and ACL together with the healthy, not depressed control group together with p-values signalling, if the different groups differ in their distributions significantly from one another.
9. Enrichr_Analysis.r
   Analysis of the biological functionality of the metabolites. Therefore, proteins corresponding with the metabolites were extracted as described in the publication and afterwards used for a gene set enrichment analysis with enrichr. The resulting     pathways were filtered and the resuling ones were visualized.
10. Abbildung_graphviz.txt
    Visualization of the combination of the connection of the extracted proteins as described for the Enrichr-Analysis. The code can be used for reconstruction of the visualization of the corresponding visualization in the publication.

