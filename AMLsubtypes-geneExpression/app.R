# ui, server, and run code for the AMLsubtypes_geneExpression Shiny App 

##### Setup libraries and custom functions #####

library(shiny)
library(tidyverse)
library(scico)
library(ggpubr)

library(BiocManager)
options(repos = BiocManager::repositories())
library(ComplexHeatmap)

# source custom functions
source("./functions.R")

# load data
source("./load_app_data.R")

#####  Define the ui #####
ui <- fluidPage(
  
  # Application title
  h1("BeatAML genetic subtypes vs. gene expression"),
  
  h2("About this Shiny App"),
  helpText("This Shiny App queries data from the BeatAML project, which includes 615 acute myeloid leukemia (AML) patients."),
  helpText("BeatAML provides gene expression data (is a gene 'active' or 'inactive') and mutation data (defining genetic subtypes of AML)."),
  
  h2("Get started:"),
  helpText("Query the gene expression of a gene of interest. Type the gene name or use the drop down menu."),
  
  selectInput(inputId = 'gene', 
              label = 'Gene of Interest:', 
              choices = autocomplete_list, 
              selected = "HIF1A", 
              multiple = FALSE),
  
  h2("Heatmap of Expression vs. Mutations"),
  helpText("The following heatmap shows the BeatAML patients arranged from lowest to highest expression of the gene of interest (top row). The mutation status of commonly mutated genes in AML (genetic subtypes) is shown for each patient in the subsequent rows. If the gene of interest is expressed significantly higher or lower in mutated vs. unmutated patients for a particular genetic subtype, then the genetic subtype is annotated with the FDR-adjusted Mann-Whitney p-value."),
  plotOutput("heatmap"),
  
  h2("Boxplots of significant results"),
  helpText("The following boxplots show the significant effects from the heatmap plot above. The expression of the gene of interest is plotted in the mutated vs. unmutated patients for a genetic subtype of AML."),
  plotOutput("boxplot"),
  
  h4("Download results"),
  downloadButton(outputId = "download_heatmap", 
                 label = "Download Heatmap"),
  
  downloadButton(outputId = "download_boxplots",
                 label = "Download Boxplots"),
  
  h4("Which genes are included in Genes of Interest list?"),
  p("Genes with low expression and low variance have been filtered out. Why: Low expression genes are near the limit of detection for the RNA-seq gene expression quantification method, so the measurement of these genes will be low accuracy. Genes with low variance will not exhibit large changes between the mutated and unmutated patient groups, so they are less likely to provide reproducible results in follow-up laboratory experiments."),
  p("Filters applied to the BeatAML gene expression dataset:"),
  p("- Filter 1: Normalized expression value > 4 in greater than 99% of patients (7388 out of 22843 total genes pass Filter 1)"),
  p("- Filter 2: Standard deviation > 1 across patients (469 out of 7388 Filter 1 genes pass Filter 2)"),
  
)

##### Define the server logic ##### 

server <- function(input, output) {
  
  output$heatmap <- renderPlot({

    # generate heatmap with expression of gene of interest
    plot_heatmap_for_expr_gene(input$gene)

  })
  
  output$boxplot <- renderPlot({
    
    # generate boxplots for significant results
    plot_boxplot_of_signif(input$gene)
    
  })
  
  output$download_heatmap <- downloadHandler(
    filename = function() {
      str_glue("{input$gene}_expression_vs_mutations_heatmap.png")
    },
    content = function(filename) {
      png(filename = filename, width = 14, height = 7, units = "in", res = 600)
      plot_heatmap_for_expr_gene(input$gene)
      dev.off()
    })
  
  # the width of the downloaded plot is scaled by the number of facets (see functions.R for generation of n_facets)
  output$download_boxplots <- downloadHandler(
    filename = function() {
      str_glue("{input$gene}_significant_results_boxplots.png")
    },
    content = function(filename) {
      ggsave(plot = plot_boxplot_of_signif(input$gene), 
             filename = filename, 
             width = eval(5+1.2*n_facets), 
             height = 7, 
             dpi = 600)
    })

}


##### Run the Shiny App #####

shinyApp(ui = ui, server = server)