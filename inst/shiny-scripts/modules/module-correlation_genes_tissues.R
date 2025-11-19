# Correlation (Genes AND Tissues) Module

# Prepare top 1000 genes including example genes
example_genes <- c("MYC", "TP53", "BRCA1", "CRP", "EGFR")
remaining_genes <- setdiff(gene_symbols_list, example_genes)
set.seed(123)
random_genes <- sample(remaining_genes, 995)
top_genes <- c(example_genes, random_genes)


correlationGenesTissuesUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      # LEFT COLUMN: Text + Run Examples
      column(
        width = 4,
        div(
          style = "max-height: 650px; overflow-y: auto; padding: 10px;
                   background-color: #f9f9f9; border-radius: 8px;
                   border: 1px solid #ddd;",

          strong("Description"),
          p("Compute Correlation for Genes AND Tissues."),
          p("This function serves as a wrapper to compute and visualize the Spearman correlation between RNA and protein expression for a user's gene list ", strong("AND"), " tissue list of interest and plots the results. It requires at least five entries in the both lists to perform the computation."),
          p("It produces either:"),
          tags$ul(
            tags$li("A per-gene plot (correlation computed across tissues), OR"),
            tags$li("A per-tissue plot (correlation computed across genes)")
          ),
          p("depending on the selected plot type."),

          br(),
          p("NOTE: This function takes a list of genes AND a list of tissues. If you want to provide only a list of tissues or genes, use the 'Correlation (Genes OR Tissues)' tab."),
          br(),
          strong("Input:"),
          tags$ol(
            tags$li("A list of gene names (5–10 required)."),
            tags$li("A list of tissue names (5–10 required)."),
            tags$li("A choice of plot type: per_gene or per_tissue.")
          ),
          br(),
          p("Note: This function queries data from HPA (Lonsdale et al., 2013) and GTEx(Lonsdale et al., 2013). To avoid high load and long query times, the", strong ("shiny UI"), " of this function is limited to 5-10 genes/tissues at a time."),
          p("If you are interested in querying larger data, you can run this function in RStudio"),
          br(),
          strong("Output:"),
          p("A plot object showing Spearman correlations per gene or per tissue."),

          br(),
          strong("Instructions"),
          p("1. Choose number of genes (5–10) and select gene names."),
          p("2. Choose number of tissues (5–10) and select tissue names."),
          p("3. Choose plot type: per_gene or per_tissue."),
          p("4. Click 'Run Correlation'."),

          br(),
          strong("Run Examples:"),
          p("To run the example, follow the steps:"),
          p("1. If you are interested in a per_gene_plot, Click on the 'Run Example (per_gene)' button to view results"),
          p("2. If you are interested in a per_tissue_plot, Click on the 'Run Example (per_tissue)' button to view results."),
          br(),
          actionButton(ns("run_example_gene"), "Run Example (per_gene)"),
          br(), br(),
          actionButton(ns("run_example_tissue"), "Run Example (per_tissue)"),
          br(), br(),
          strong("Interpretation of plot"),
          br(),
          br(),
          strong("per_gene Plot"),
          p("This function computes the spearman correlation of the gene list of interest and tissue list of interest. The final plot is plotted with the X axis being the individual genes and the Y axis being the spearman coefficient."),
          p("Higher Spearman correlation coefficients indicate stronger agreement between RNA and protein expression levels across samples. Users can use this visualization to identify tissues or genes with consistent expression trends."),
          p("Plot Key: High spearman values are indicated in blue, medium spearman values are indicated in white and low spearman values are indicated in red."),
          br(),
          strong("per_tissue Plot"),
          p("This function computes the spearman correlation of the gene list of interest and tissue list of interest. The final plot is plotted with the X axis being the individual tissues and the Y axis being the spearman coefficient."),
          p("Higher Spearman correlation coefficients indicate stronger agreement between RNA and protein expression levels across samples. Users can use this visualization to identify tissues or genes with consistent expression trends."),
          p("Plot Key: High spearman values are indicated in blue, medium spearman values are indicated in white and low spearman values are indicated in red."),
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

      # RIGHT COLUMN: Inputs + Plot
      column(
        width = 8,

        sliderInput(
          inputId = ns("num_genes"),
          label = "Number of genes (5–10):",
          min = 5, max = 10, value = 5
        ),
        uiOutput(ns("gene_selects")),

        sliderInput(
          inputId = ns("num_tissues"),
          label = "Number of tissues (5–10):",
          min = 5, max = 10, value = 5
        ),
        uiOutput(ns("tissue_selects")),

        selectInput(
          inputId = ns("plot_choice"),
          label = "Select plot to display:",
          choices = c("per_gene", "per_tissue"),
          selected = "per_gene"
        ),

        actionButton(ns("run"), "Run Correlation"),
        br(), br(),

        h3("Correlation Plot"),
        plotOutput(ns("cor_plot"))
      )
    )
  )
}



correlationGenesTissuesServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # -------- Dynamic dropdowns --------
    output$gene_selects <- renderUI({
      n <- input$num_genes
      lapply(1:n, function(i) {
        selectizeInput(
          ns(paste0("gene_", i)),
          label = paste("Gene", i),
          choices = top_genes,
          selected = top_genes[i],
          multiple = FALSE,
          options = list(maxOptions = 1000)
        )
      }) |> tagList()
    })

    output$tissue_selects <- renderUI({
      n <- input$num_tissues
      lapply(1:n, function(i) {
        selectInput(
          ns(paste0("tissue_", i)),
          label = paste("Tissue", i),
          choices = tissue_map$protein_tissue,
          selected = tissue_map$protein_tissue[i]
        )
      }) |> tagList()
    })

    # Collect inputs
    selected_genes <- reactive({
      sapply(1:input$num_genes, function(i) input[[paste0("gene_", i)]])
    })

    selected_tissues <- reactive({
      sapply(1:input$num_tissues, function(i) input[[paste0("tissue_", i)]])
    })

    rv <- reactiveValues(result = NULL)


    # -------- Standard Run --------
    observeEvent(input$run, {
      rv$result <- correlation_genes_tissues(
        gene_NAMES = selected_genes(),
        tissue_NAMES = selected_tissues(),
        plot_choice = input$plot_choice
      )
    })


    # -------- Run Example: per_gene --------
    observeEvent(input$run_example_gene, {

      updateSliderInput(session, "num_genes", value = 5)
      updateSliderInput(session, "num_tissues", value = 5)

      for (i in 1:5) {
        updateSelectizeInput(session, paste0("gene_", i),
                             selected = example_genes[i])
        updateSelectInput(session, paste0("tissue_", i),
                          selected = tissue_map$protein_tissue[i])
      }

      updateSelectInput(session, "plot_choice", selected = "per_gene")

      rv$result <- correlation_genes_tissues(
        gene_NAMES = example_genes,
        tissue_NAMES = tissue_map$protein_tissue[1:5],
        plot_choice = "per_gene"
      )
    })


    # -------- Run Example: per_tissue --------
    observeEvent(input$run_example_tissue, {

      updateSliderInput(session, "num_genes", value = 5)
      updateSliderInput(session, "num_tissues", value = 5)

      for (i in 1:5) {
        updateSelectizeInput(session, paste0("gene_", i),
                             selected = example_genes[i])
        updateSelectInput(session, paste0("tissue_", i),
                          selected = tissue_map$protein_tissue[i])
      }

      updateSelectInput(session, "plot_choice", selected = "per_tissue")

      rv$result <- correlation_genes_tissues(
        gene_NAMES = example_genes,
        tissue_NAMES = tissue_map$protein_tissue[1:5],
        plot_choice = "per_tissue"
      )
    })


    # -------- Plot Output --------
    output$cor_plot <- renderPlot({
      req(rv$result)
      if (input$plot_choice == "per_gene") {
        rv$result$per_gene_plot
      } else {
        rv$result$per_tissue_plot
      }
    })

  })
}
