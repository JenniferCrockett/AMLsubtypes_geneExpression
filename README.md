# BeatAML genetic subtypes vs. gene expression

## About

This repository builds the "BeatAML genetic subtypes vs. gene expression" Shiny App.  

This app queries data from the BeatAML project, which includes 615 acute myeloid leukemia (AML) patients. BeatAML provides gene expression data (is a gene 'active' or 'inactive') and mutation data (defining genetic subtypes of AML). The user queries the gene expression of their gene of interest. The app produces the following plots:  

* A **heatmap** showing the BeatAML patients arranged from lowest to highest expression of the gene of interest, with the mutation status of commonly mutated genes in AML for each patient (genetic subtypes). If the gene of interest is expressed significantly higher or lower in mutated vs. unmutated patients for a particular genetic subtype, then the genetic subtype is annotated with the FDR-adjusted Mann-Whitney p-value.
* **Boxplots** showing the significant effects from the heatmap plot. The expression of the gene of interest is plotted in the mutated vs. unmutated patients for a genetic subtype of AML.

These plots can be used to determine whether a gene is active or inactive in a genetic subtype of AML. For example, a gene that is turned on highly in a particular genetic subtype could be a potential new drug-target to treat that AML subtype.

## Setup

### Clone repository

To clone this repository, run:

```
git clone https://github.com/JenniferCrockett/AMLsubtypes_geneExpression.git
```

### Run the app

TODO: Include details once I figure out if the shinylive implementation will work on Github

## Other documentation

### Raw data preparation script

For documentation purposes, the `./data_preparation/` directory contains raw data preparation script (`Data_preparation.Rmd` and rendered HTML report), which generates the data to power the app.

The following raw data files are required to run `Data_preparation.Rmd`, if you wish to re-process the data:  

* Clinical Summary: `./raw_data/beataml_wv1to4_clinical.xlsx`
* WES/targeted Sequencing Mutation Calls: `./raw_data/beataml_wes_wv1to4_mutations_dbgap.txt`
* Normalized Expression: https://github.com/biodev/beataml2.0_data/raw/main/beataml_waves1to4_norm_exp_dbgap.txt
  - This file exceeds Github's file size limit, so was not included in the `./raw_data/` directory
  - It can be downloaded via the provided link, to the `./raw_data/` directory in your clone of the repository