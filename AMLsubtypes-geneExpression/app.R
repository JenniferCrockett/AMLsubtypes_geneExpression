# ui, server, and run code for the AMLsubtypes-geneExpression Shiny App 

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
  br(),
  
  h2("About this Shiny App"),
  p("This Shiny App queries data from the BeatAML project, which includes 615 acute myeloid leukemia (AML) patients. For each patient, BeatAML provides gene expression data (is a gene 'active' or 'inactive') and mutation data (defining genetic subtypes of AML). Using this App, you can query a gene of interest to see whether its expression is different between mutated and unmutated patients of each genetic subtype."),
  br(),
  
  h2("Get started:"),
  helpText("Query the gene expression of a gene of interest. Type the gene name or use the drop down menu."),
  
  selectInput(inputId = 'gene', 
              label = 'Gene of Interest:', 
              choices = autocomplete_list, 
              selected = "HIF1A", 
              multiple = FALSE),
  
  br(),
  h2("Heatmap"),
  plotOutput("heatmap", width = "95%"),
  helpText("Each column in the heatmap displays the data for one BeatAML patient."),
  helpText("Top row, 'expression': The patients are arranged from lowest to highest z-scaled expression of the gene of interest."),
  helpText("Subsequent rows: The mutation status of commonly mutated genes in AML (genetic subtypes) is shown for each patient. If the gene of interest is expressed significantly higher or lower in mutated vs. unmutated patients for a particular genetic subtype, then the genetic subtype is annotated with the FDR-adjusted Mann-Whitney p-value."),
  helpText("    - If the gene of interest is more highly expressed (more active) in mutated vs. unmutated patients, then the mutated patients (black lines) will be concentrated at the right of the heatmap."),
  helpText("    - If the gene of interest is more lowly expressed (less active) in mutated vs. unmutated patients, then the mutated patients (black lines) will be concentrated at the left of the heatmap."),
  helpText("    - If multiple genetic subtypes show a statistically significant result for the gene of interest, you can see whether the mutations tend to co-occur within the same patients (black lines line up vertically)."),
  
  br(),
  h2("Boxplots of significant results"),
  plotOutput("boxplot", width = "95%"),
  helpText("The boxplots show the significant effects from the heatmap plot above. The expression of the gene of interest is plotted in the mutated vs. unmutated patients for a genetic subtype of AML."),
  
  br(),
  h4("Download results"),
  downloadButton(outputId = "download_heatmap", 
                 label = "Download Heatmap"),
  
  downloadButton(outputId = "download_boxplots",
                 label = "Download Boxplots"),
  
  br(),
  br(),
  h4("Which genes are included in Genes of Interest list?"),
  p("Genes with low expression and low variance have been filtered out. Why: Low expression genes are near the limit of detection for the RNA-seq gene expression quantification method, so the measurement of these genes will be low accuracy. Genes with low variance will not exhibit large changes between the mutated and unmutated patient groups, so they are less likely to provide reproducible results in follow-up laboratory experiments."),
  p("Filters applied to the BeatAML gene expression dataset:"),
  p("- Filter 1: Normalized expression value > 4 in greater than 99% of patients (7388 out of 22843 total genes pass Filter 1)"),
  p("- Filter 2: Standard deviation > 1 across patients (469 out of 7388 Filter 1 genes pass Filter 2)"),
  
  br(),
  h4("Which genetic subtypes are included in the Mutated/Unmutated genes list?"),
  p("Mutation data was included for genes that were mutated in at least 5% of the BeatAML patient cohort."),
  
  br(),
  h4("References"),
  tags$p("1. BeatAML publication: ", tags$a(href = "https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00312-9", "Bottomly et al. (2022) Integrative analysis of drug response and clinical outcome in acute myeloid leukemia.")),
  tags$p("2. BeatAML publicly available datasets:", tags$a(href = "https://biodev.github.io/BeatAML2/", "https://biodev.github.io/BeatAML2/")),
  tags$p("3. Documentation for the AMLsubtypes-geneExpression Shiny App:", tags$a(href = "https://github.com/JenniferCrockett/AMLsubtypes_geneExpression", "https://github.com/JenniferCrockett/AMLsubtypes_geneExpression"))
  
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