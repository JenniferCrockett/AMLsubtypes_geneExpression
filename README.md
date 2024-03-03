# BeatAML genetic subtypes vs. gene expression

## About

This repository builds the "BeatAML genetic subtypes vs. gene expression" Shiny App.  

This app queries data from the [BeatAML project](https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00312-9), which includes 615 acute myeloid leukemia (AML) patients. The [BeatAML dataset](https://biodev.github.io/BeatAML2/) provides gene expression data (is a gene 'active' or 'inactive') and mutation data (defining genetic subtypes of AML). The user queries the gene expression of their gene of interest. The app produces the following plots:  

* A **heatmap** showing the BeatAML patients arranged from lowest to highest expression of the gene of interest, with the mutation status of commonly mutated genes in AML for each patient (genetic subtypes). If the gene of interest is expressed significantly higher or lower in mutated vs. unmutated patients for a particular genetic subtype, then the genetic subtype is annotated with the FDR-adjusted Mann-Whitney p-value.
* **Boxplots** showing the significant effects from the heatmap plot. The expression of the gene of interest is plotted in the mutated vs. unmutated patients for each genetic subtype of AML that showed a statistically significant difference in gene expression.

These plots can be used to determine whether a gene is active or inactive in a genetic subtype of AML. Examples of why researchers are interested in this question:  

* A gene that is highly active in a particular genetic subtype could be crucial to the survival of those AML cells. Inactivating the gene using a drug could specifically eliminate the AML cells.
* A gene that is inactive in a particular genetic subtype be toxic to those AML cells. Activating the gene using a drug could specifically eliminate the AML cells.

## App

Use the Shiny App at [https://jennifercrockett.shinyapps.io/amlsubtypes-geneexpression/](https://jennifercrockett.shinyapps.io/amlsubtypes-geneexpression/)

### Other documentation

#### Raw data preparation script

For documentation purposes, the `./data_preparation/` directory contains raw data preparation script (`Data_preparation.Rmd` and rendered HTML report), which generates the data to power the app.

The following raw data files are required to run `Data_preparation.Rmd`, if you wish to re-process the data:  

* Clinical Summary: `./raw_data/beataml_wv1to4_clinical.xlsx`
* WES/targeted Sequencing Mutation Calls: `./raw_data/beataml_wes_wv1to4_mutations_dbgap.txt`
* Normalized Expression: https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_norm_exp_dbgap.txt
  - This file exceeds Github's file size limit, so was not included in the `./raw_data/` directory
  - It can be downloaded via the provided link, to the `./raw_data/` directory in your clone of the repository
