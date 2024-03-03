# Custom functions for the AMLsubtypes_geneExpression project

# function to prepare gene expression and genetic subtypes data for further testing
# PARAMS: 
## gene: the name of the gene of interest in the gene expression dataframe
# RETURN:
## for each patient (row), the mutation status for the genetic subtypes and the Z-scaled gene expression value of the @gene parameter (last column)
prepare_gene_expr_and_subtypes <- function(gene) {
  this_expr <- rna_zscaled_filt[gene,] %>%
    t() %>%
    as.data.frame()
  
  # we need to rename the gene expression column, so that the colname is never the same as in the genetic subtypes data
  expr_name <- str_glue("{gene}.expr")
  colnames(this_expr) <- expr_name
  
  this_data <- cbind(mut_matrix_all, this_expr)
  
  return(this_data)
}


# function to perform the statistical test for expression of each gene in mutated vs. unmutated patients of each subtype 
# used in Data_preparation.Rmd only
# PARAMS:
## gene: the name of the gene of interest in the gene expression dataframe
# RETURN:
## FDR-adjusted p-value for gene expression of the @gene param in mutated vs. unmutated patients in each genetic subtype (columns). Only returns one row, where the rowname is the @gene param.
test_gene_expr_in_subtypes <- function(gene) {
  expr_name <- str_glue("{gene}.expr")
  
  to_test <- prepare_gene_expr_and_subtypes(gene)
  
  this_result <- data.frame(genetic_subtype = genetic_subtypes, nominal_pval = NA)
  
  for (i in seq_along(genetic_subtypes)) {
    group <- genetic_subtypes[i]
    
    this_test <- wilcox.test(as.formula(str_glue("{expr_name} ~ {group}")), data = to_test, paired = FALSE) 
    
    this_result[i, "nominal_pval"] <- this_test$p.value
  }
  
  # p adjustment (FDR control)
  this_result$padj <- p.adjust(this_result$nominal_pval, method = "fdr")
  
  # add p-adj summary as asterisks
  # combine with actual p-adj value for annotating plots
  this_result <- mutate(this_result, padj_asterisks = case_when(padj < 0.0001 ~ "***",
                                                                padj < 0.001 ~ "**",
                                                                padj < 0.01 ~ "**", 
                                                                padj < 0.05 ~ "*")) %>%
    mutate(padj_summary = case_when(is.na(padj_asterisks) ~ NA,
                                    TRUE ~ str_glue("{signif(padj, 2)} {padj_asterisks}")))
  
  # return a matrix where rows are the "gene" input value, columns are the genetic subtypes, and values are the padj_summary
  this_return <- select(this_result, genetic_subtype, padj_summary) %>%
    column_to_rownames("genetic_subtype") %>%
    t() %>%
    as.data.frame()
  
  rownames(this_return) <- gene
  
  return(this_return)
}


# Heatmap plot function
# PARAMS:
## gene: the name of the gene of interest in the gene expression dataframe
## out: select output mode. Must be either "heatmap" or "ht_object"
# RETURN:
## A. if out == "heatmap"" A heatmap plot where the columns are arranged from lowest to highest Z-scaled expression of the @gene param, 
## and rows are arranged in alphabetical order of the genetic subtypes, OR
## B. An ht object to be rendered by another plot rendering function.
plot_heatmap_for_expr_gene <- function(gene) {
  expr_name <- str_glue("{gene}.expression")
    
  # prepare the stats summary
  this_stat <- stat_results[gene,] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("genetic_subtype") %>%
    dplyr::rename(stat = 2) %>%
    mutate(stat_summary = case_when(!is.na(stat) ~ str_glue("{genetic_subtype} padj={stat}"),
                                    TRUE ~ genetic_subtype))
  
  # prepare the data to be plotted
  to_plot <- prepare_gene_expr_and_subtypes(gene) %>%
    t()
  
  # update the rownames of the plot data to include the stat summary values
  if (identical(this_stat$genetic_subtype, rownames(to_plot[1:(nrow(to_plot)-1),]), attrib.as.set = F)) {
    rownames(to_plot) <- c(this_stat$stat_summary, expr_name)
  } else {
    stop("Rows are not in same order. Fix before combining plot data and stat summary!")
  }
  
  # update the column order of the plot data so that patients are arranged from lowest to highest expression of the @gene param
  patient_order <- to_plot[expr_name,] %>% 
    sort() %>%
    names()
  
  to_plot <- to_plot[,patient_order]
  
  # split the data into genetic subtypes and expression
  subtypes <- head(to_plot, -1) %>% 
    as.data.frame()
  expr <- tail(to_plot, 1) %>%
    as.data.frame()
  
  # create the colour palette to be used in the plot
  col <- scico(5, palette = "roma", direction = -1)
  breaks <- summary(unlist(expr))
  cpal_expr <- circlize::colorRamp2(breaks = c(breaks[1], breaks[2], breaks[3], breaks[5], breaks[6]), 
                                    colors = c(col[1], col[2], col[3], col[4], col[5]))
  
  col_subtypes <- scico(2, palette = "grayC", end = 0.9, direction = -1) # colours for Unmutated and Mutated patients
  names(col_subtypes) <- c("Unmutated", "Mutated")
  
  # formatting required specifically for the ComplexHeatmap Oncoprint input
  # see: https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html
  
  ## genetic subtypes data
  subtypes_mat <- mutate_all(subtypes, function(x){ifelse(x == 0, "Unmutated", "Mutated")})
  
  ## an "alteration function" is needed to define how mutations should be represented
  alter_fun <- list(
    background = alter_graphic("rect", fill = "#FFFFFF"),   
    Mutated = alter_graphic("rect", fill = col_subtypes["Mutated"], height = 1),
    Unmutated = alter_graphic("rect", fill = col_subtypes["Unmutated"], height = 1)
  )
  
  ## expression data row
  expr_anno <- HeatmapAnnotation(expression = unlist(expr), col = list(expression = cpal_expr))
  expr_anno@name <- expr_name
  
  ## title and legend
  column_title = str_glue("BeatAML: {gene} Z-scaled expression vs. Common mutations\n(n = {ncol(to_plot)} patients)")
  heatmap_legend_param = list(title = "Status", at = c("Unmutated", "Mutated"), 
                              labels = c("Unmutated", "Mutated"), title_position = "topcenter")  
  
  # generate heatmap
  suppressMessages(
    ht <- oncoPrint(subtypes_mat,
                  alter_fun = alter_fun, 
                  col = col_subtypes, 
                  column_title = column_title, 
                  heatmap_legend_param = heatmap_legend_param, 
                  show_column_names = FALSE, 
                  top_annotation = expr_anno, 
                  right_annotation = NULL,
                  show_pct = FALSE,
                  column_title_gp = gpar(fontsize = 20))
    )
  
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
 
}

# Boxplot function
# PARAMS:
## gene: the name of the gene of interest in the gene expression dataframe
# RETURNS:
## 1. One boxplot facet for each genetic subtype where a significant effect was observed in the heatmap,
## where expression of the gene of interest is plotted in unmutated and mutated patient groups for each subtype.
## 2. A global environment variable called "n_facets" is exported, which will be used for scaling the size 
## of downloaded plots.
plot_boxplot_of_signif <- function(gene) {
  expr_name <- str_glue("{gene}.log2RPKM")
    
    # gather the genetic subtypes that show significant effects
    signif_effects <- stat_results[gene,] %>%
      unlist() %>%
      na.omit() %>%
      names()
    
    # prepare the data to be plotted
    # it is typical for researchers to want to plot the log2 normalized expression data (log2RPKM value) in a boxplot, 
    # not the Z-scaled expression used for heatmap
    to_plot <- prepare_gene_expr_and_subtypes(gene)[,signif_effects]
    
    log2_norm_exp <- rna_normexp_filt[gene,] %>%
      t()
    colnames(log2_norm_exp) <- expr_name
    
    if (identical(rownames(to_plot), rownames(log2_norm_exp), attrib.as.set = F)) {
      to_plot <- cbind(to_plot, log2_norm_exp)
    } else {
      stop("Rows are not in same order. Fix before combining!")
    }
    
    # format the plot data: pivot wide to long
    to_plot_l <- rownames_to_column(to_plot, "patient_id") %>%
      pivot_longer(cols = all_of(signif_effects), names_to = "genetic.subtype", values_to = "mutation.status") %>%
      mutate(mutation.status = case_match(mutation.status, 0 ~ "Unmutated", 1 ~ "Mutated"))
    
    to_plot_l$mutation.status <- factor(to_plot_l$mutation.status, levels = c("Unmutated", "Mutated"))
    
    # prepare the stats annotations for the plot
    # use the adjusted p-values from when all hypotheses were tested (all genetic subtypes were tested)
    p_anno <- stat_results[gene,signif_effects] %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("genetic.subtype") %>%
      dplyr::rename("p.label" = 2) %>%
      mutate(group1 = "Unmutated", group2 = "Mutated") %>%
      select(genetic.subtype, group1, group2, p.label)
    
    # determine the y.position on the plot where the p-value annotation should appear
    y_pos <- max(rna_normexp_filt[gene,])*1.05
    
    # boxplot
    p <- ggboxplot(to_plot_l, 
                   x = "mutation.status", 
                   y = expr_name, 
                   fill = "mutation.status",
                   ggtheme = theme_light()) +
      facet_wrap(~genetic.subtype, nrow = 1) +
      scale_fill_scico_d(palette = "devon", direction = -1, begin = 0.3, end = 0.7) +
      stat_pvalue_manual(data = p_anno, label = "p.label", y.position = y_pos, size = 4) +
      theme_light(base_size = 22) +
      theme(axis.text.x = element_blank())
    
    ## if there are > 10 facets, then shrink the font size of the facet label to avoid crowding
    n_facets <<- ncol(stat_results[gene,signif_effects]) # the `<<-` operator exports a global environment variable
    if (n_facets > 10) {
      p <- p +
        theme(strip.text = element_text(size = 14))
    }
    
    print(p)
}
