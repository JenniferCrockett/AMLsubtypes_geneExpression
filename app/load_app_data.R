# Load data files for the AMLsubtypes_geneExpression app

## Gene expression data
rna_normexp_filt <- read_rds(here("app", "data", "rna_normexp_filt.RDS"))
rna_zscaled_filt <- read_rds(here("app", "data", "rna_zscaled_filt.RDS"))

## Genetic subtypes data
mut_matrix_all <- read_rds(here("app", "data", "mut_matrix_all.RDS"))

## Statistical summary
stat_results <- read_rds(here("app", "data", "stat_results_mann_whitney_padj.RDS"))

## Autocomplete list for gene input value
autocomplete_list <- rownames(rna_normexp_filt)
