# Correlation Module (Genes OR Tissues)

# Prepare top 1000 genes including example genes
example_genes <- c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
remaining_genes <- setdiff(gene_symbols_list, example_genes)
set.seed(123)
random_genes <- sample(remaining_genes, 995)
top_genes <- c(example_genes, random_genes)

example_tissues <- c("lung", "spleen", "liver", "ovary", "testis")

# UI
correlationUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      # Left column: Scrollable descriptive text + Run Examples
      column(
        width = 4,
        div(
          style = "max-height: 600px; overflow-y: auto; padding: 10px; background-color: #f9f9f9;
                   border-radius: 8px; border: 1px solid #ddd;",
          strong("Description"),
          p("Compute Correlation for Genes OR Tissues."),
          p("This function serves as a wrapper to compute and visualize the Spearman correlation between RNA and protein expression for a user's gene list ", strong("OR"), " tissue list of interest and plots the results. It requires at least five entries in the provided list to perform the computation."),
          p("NOTE: This function only takes a list of genes or a list of tissues, and does not take both. If you want to provide both, use the 'Correlation (Genes AND Tissues)' tab."),
          br(),
          strong("Input:"),
          p("Either a list of gene names or tissue names. Each list needs to be between 5 to 10 elements inclusive."),
          br(),
          p("Note: This function queries data from HPA (Lonsdale et al., 2013) and GTEx(Lonsdale et al., 2013). To avoid high load and long query times, the", strong ("shiny UI"), " of this function is limited to 5-10 genes/tissues at a time."),
          p("If you are interested in querying larger data, you can run this function in RStudio"),
          br(),
          strong("Output:"),
          p("Returns a per-gene plot if gene names are provided, or a per-tissue plot if tissue names are provided. The Spearman correlation values for each gene/tissue is labeled in the plot."),
          br(),
          strong("Instructions:"),
          p("1. Choose input type: Genes or Tissues."),
          p("2. Use the slider to select the number of entries (5-10)."),
          p("3. Select entries from dropdowns."),
          p("4. Click 'Run Correlation' to view results."),
          br(),
          strong("Run Examples:"),
          p("To run the example, follow the steps:"),
          p("1. If you are interested in a gene list as input, Click on the 'Run Example (Genes)' button to view results"),
          p("2. If you are interested in a tissue list as input, Click on the 'Run Example (Tissues)' button to view results."),
          br(),
          actionButton(ns("run_example_genes"), "Run Example (Genes)"),
          br(), br(),
          actionButton(ns("run_example_tissues"), "Run Example (Tissues)"),
          br(), br(),
          br(),
          strong("Note:"),
          p("ExprCompareR retrieves data from HPA (Thul et al., 2018) and GTEx (Lonsdale et al., 2013) via packages, HPAanalyze(Tran et al., 2019) and gtexr(Warwick et al., 2025). Because these functions make external API calls, some operations may take time to complete typically around 30 seconds to 1 minute."),
          br(),
          strong("Example Configuration:"),
          p("The example that  'Run Example (Genes)'  button runs is "),
          br(),
          p("genes = MYC, TP53, BRCA1, CRP, EGFR"),
          br(),
          p("NOTE: This function uses a standard list of tissues to obtain the expression of the gene list of interest, these are:cerebellum, kidney, adrenal gland, breast,lung, spleen, liver, ovary, testis, pancreas"),
          br(),
          p("The example that  'Run Example (Tissues)'   button runs is" ),
          p("tissues = lung, spleen, liver, ovary, testis"),
          br(),
          p("All inputs are automatically updated when the example runs."),
          br(),
          strong("Interpretation of plot"),
          br(),
          br(),
          strong("Gene List Plot"),
          p("This function computes the spearman correlation of the gene list of interest across standard tissues (kidney, lung, spleen, ovary, testis, breast, cerebellum, pancreas and adrenal gland) with the X axis being the individual genes and the Y axis being the spearman coefficient."),
          p("Higher Spearman correlation coefficients indicate stronger agreement between RNA and protein expression levels across samples. Users can use this visualization to identify tissues or genes with consistent expression trends."),
          br(),
          p(strong("Plot Key:"), "Spearman correlation values are represented by color intensity. High positive correlations are indicated in dark blue, low correlations around zero are white, and high negative correlations are indicated in dark red. The shade intensity reflects the magnitude of the correlation."),
          br(),
          p("If the resulting plot, has no plotted values for certain genes or tissues, then please check the warning issued on Rstudio. This warning issued will be : 'Some gene-tissue combinations are missing in protein data. These missing combos indicate that specific genes of interest are not present in the
            tissues of interest. You can either remove these genes from your analysis, replace these values or they will be filled by NA values'.This warning leads to no plotted values for certain genes and tissues combos. This does not affect the quality of the plot, but rather informs the user that some gene-tissue combinations do not exist."),
          br(),
          strong("Tissue List Plot"),
          p("This function considers the top expressed genes in each given tissue in the user's tissue list of interest, with the X axis being the individual tissues and the Y axis being the spearman coefficient"),
          p("Higher Spearman correlation coefficients indicate stronger agreement between RNA and protein expression levels across samples. Users can use this visualization to identify tissues or genes with consistent expression trends."),
          br(),
          p(strong("Plot Key:"), "Spearman correlation values are represented by color intensity. High positive correlations are indicated in dark blue, low correlations around zero are white, and high negative correlations are indicated in dark red. The shade intensity reflects the magnitude of the correlation."),
          br(),
          p("If the resulting plot, has no plotted values for certain genes or tissues, then please check the warning issued on Rstudio. This warning issued will be : 'Some gene-tissue combinations are missing in protein data. These missing combos indicate that specific genes of interest are not present in the
            tissues of interest. You can either remove these genes from your analysis, replace these values or they will be filled by NA values'.This warning leads to no plotted values for certain genes and tissues combos. This does not affect the quality of the plot, but rather informs the user that some gene-tissue combinations do not exist."),
          br(),
          strong("References:"),
          p("* Cetinkaya-Rundel M,Cheng J, Grolemund G (2017).Customize your UI with HTML.https://shiny.posit.co/r/articles/build/html-tags/"),
          p("* Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Aden-Buie G, Xie Y,Allen J, McPherson J, Dipert A, Borges B (2025).shiny: Web Application Framework for R. R package version 1.11.1,https://CRAN.R-project.org/package=shiny"),
          p("* OpenAI.(2025). ChatGPT (GPT-5) Large language model.Retrieved November 17, 2025, from https://chatgpt.com/ "),
          p("* Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019).,A multivariate Poisson-log normal mixture model for clustering transcriptome sequencing data.BMC Bioinformatics. 2019;20(1):394. https://pubmed.ncbi.nlm.nih.gov/31311497/"),
          p("* Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).HPAanalyze: an R package that facilitates the retrieval and analysis of the Human Protein Atlas data. MC Bioinformatics 20, 463 (2019).https://doi.org/10.1186/s12859-019-3059-z"),
          p("* Warwick A, Zuckerman B, Ung C, Luben R, Olvera-Barrios A (2025). gtexr: A convenient R interface to the Genotype-Tissue Expression (GTEx) Portal API.Journal of Open Source So ware, 10(109), 8249. ISSN 2475-9066,doi:10.21105/joss.08249}, gigs v0.2.1."),
          p("* Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3 319-24277-4, https://ggplot2.tidyverse.org."),
          p("* Wickham H, Francois R, Henry L, Muller K, Vaughan D (2025).dplyr: A Grammar of Data Manipulation. R package version 1.1.4, https://dplyr.tidyverse.org."),
          p("* Wickham H, Henry L (2025). purrr: Functional Programming Tools. R package version 1.1.0, https://purrr.tidyverse.org/."),
          p("* Wickham H, Vaughan D, Girlich M (2025). tidyr: Tidy Messy Data. R package version 1.3.1, https://tidyr.tidyverse.org."),
          br()
        )
      ),

      # Right column: Inputs and output
      column(
        width = 8,
        radioButtons(
          inputId = ns("input_type"),
          label = "Select input type:",
          choices = c("Genes", "Tissues"),
          inline = TRUE
        ),
        sliderInput(
          inputId = ns("num_inputs"),
          label = "Number of selections (5-10):",
          min = 5,
          max = 10,
          value = 5,
          step = 1
        ),
        uiOutput(ns("dynamic_selects")),
        actionButton(ns("run"), "Run Correlation"),
        br(), br(),
        h3("Correlation Plot"),
        plotOutput(ns("cor_plot"))
      )
    )
  )
}

# Server
correlationServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Dynamically generate dropdowns
    output$dynamic_selects <- renderUI({
      ns <- session$ns
      n <- input$num_inputs

      selects <- lapply(1:n, function(i) {
        if (input$input_type == "Genes") {
          selectizeInput(
            ns(paste0("sel_", i)),
            label = paste("Gene", i),
            choices = top_genes,
            selected = top_genes[i],
            multiple = FALSE,
            options = list(maxOptions = 1000)
          )
        } else {
          selectInput(
            ns(paste0("sel_", i)),
            label = paste("Tissue", i),
            choices = tissue_map$protein_tissue,
            selected = tissue_map$protein_tissue[i]
          )
        }
      })

      do.call(tagList, selects)
    })

    # Gather selected inputs
    selected_inputs <- reactive({
      n <- input$num_inputs
      sapply(1:n, function(i) input[[paste0("sel_", i)]])
    })

    # Function to run correlation
    run_correlation <- function(input_type, selections) {
      if (input_type == "Genes") {
        compute_correlation(gene_NAMES = selections, tissue_NAMES = NULL)
      } else {
        compute_correlation(gene_NAMES = NULL, tissue_NAMES = selections)
      }
    }

    # ReactiveValues to store result
    rv <- reactiveValues(result = NULL)

    # Normal Run button
    observeEvent(input$run, {
      req(selected_inputs())
      rv$result <- run_correlation(input$input_type, selected_inputs())
    })

    # Run Example: Genes
    observeEvent(input$run_example_genes, {
      updateRadioButtons(session, "input_type", selected = "Genes")
      updateSliderInput(session, "num_inputs", value = 5)
      for (i in 1:5) {
        updateSelectizeInput(session, paste0("sel_", i), selected = example_genes)
      }
      rv$result <- run_correlation("Genes", example_genes)
    })

    # Run Example: Tissues
    observeEvent(input$run_example_tissues, {
      updateRadioButtons(session, "input_type", selected = "Tissues")
      updateSliderInput(session, "num_inputs", value = 5)
      for (i in 1:5) {
        updateSelectInput(session, paste0("sel_", i), selected = example_tissues)
      }
      rv$result <- run_correlation("Tissues", example_tissues)
    })

    # Output
    output$cor_plot <- renderPlot({
      req(rv$result)
      if (input$input_type == "Genes"){
        rv$result$per_gene_plot
        } else {
          rv$result$per_tissue_plot
        }
    })

  })
}
