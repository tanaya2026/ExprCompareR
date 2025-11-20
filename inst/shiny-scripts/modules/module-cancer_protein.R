# Cancer Protein Module

# UI function
cancerUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      # Left column: Scrollable descriptive text + Run Example
      column(
        width = 4,
        div(
          style = "max-height: 500px; overflow-y: auto; padding: 10px; background-color: #f9f9f9;
                   border-radius: 8px; border: 1px solid #ddd;",
          strong("Description"),
          p("Comparison of protein expression between normal and cancer tissues."),
          p("This function compares the protein expression levels between normal and pathology (cancer) tissues for a specified cancer type. It calculates the change in expression rank and direction ('Up', 'Down', or 'No change') for each gene, and visualizes the shift in expression."),
          br(),
          strong("Input:"),
          p("The cancer type of interest."),
          strong("Output:"),
          p("Returns a plot showing the direction and magnitude of protein rank changes for each gene and a metadata tibble with columns:"),
          tags$ul(
            tags$li("{ensembl}: Ensembl gene identifier"),
            tags$li("{gene}: Gene symbol"),
            tags$li("{normal_rank}: Average protein expression rank in normal tissue"),
            tags$li("{cancer_rank}: Weighted protein expression rank in cancer tissue"),
            tags$li("{delta_rank}: Difference (cancer_rank - normal_rank)"),
            tags$li("{direction}: Direction of change ('Up', 'Down', or 'No change')")
          ),
          br(),
          strong("Instructions:"),
          p("1. From the 'Select a cancer type' toolbar, choose a cancer type to analyze."),
          p("2. Once selected, click on 'Run Cancer Protein Analysis' to view results."),
          p("3. The output plot is under the heading 'Protein Rank Shift Plot'."),
          p("4. The output tibble is under the heading 'Protein Expression Table'."),
          br(),
          strong("Run Examples:"),
          p("To run the example, follow the steps:"),
          p("1. Click on the 'Run Example' button, to view the results"),
          br(),
          # Run Example button
          actionButton(ns("run_example"), "Run Example"),
          br(), br(),
          strong("Example Configuration:"),
          br(),
          p("The example that 'Run Example' button runs is:"),
          p("cancer type = ovarian cancer"),
          br(),
          p("All inputs are automatically updated when the example runs."),
          br(),
          strong("Interpretation of plot"),
          p("Each gene is plotted to show the magnitude and direction of change in protein rank between normal and cancer tissue. The X axis represents the delta rank i.e. the change in protein expression of normal vs pathology. A decrease in rank is indicated by the colour red, an increase is indicated in the colour blue, and no change is indicated by green. This plot gives researchers a quick overview of the general direction of the genes in this cancer type."),
          p("If you want to find the genes whose direction has changed, then the protein table, gives a more in depth explanation below."),
          br(),
          strong("Interpretation of table"),
          p("As mentioned above in the output section, the table is structured, to allow researchers to investigate the role of each gene involved in this cancer type."),
          p("This table assists in scientists narrowing down their gene of interest based on the delta rank and direction column. These genes maybe be affected during cancer development and can be potential driver genes. "),
          br(),
          strong("References"),
          p("* Cetinkaya-Rundel M,Cheng J, Grolemund G (2017).Customize your UI with HTML.https://shiny.posit.co/r/articles/build/html-tags/"),
          p("* Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Aden-Buie G, Xie Y,Allen J, McPherson J, Dipert A, Borges B (2025).shiny: Web Application Framework for R. R package version 1.11.1,https://CRAN.R-project.org/package=shiny"),
          p("* OpenAI.(2025). ChatGPT (GPT-5) Large language model.Retrieved November 17, 2025, from https://chatgpt.com/ "),
          p("* Silva, A., S. J. Rothstein, P. D. McNicholas, and S. Subedi (2019).,A multivariate Poisson-log normal mixture model for clustering transcriptome sequencing data.BMC Bioinformatics. 2019;20(1):394. https://pubmed.ncbi.nlm.nih.gov/31311497/"),
          p("* Tran AN, Dussaq AM, Kennell Jr T, Willey CD, Hjelmeland AB (2019).HPAanalyze: an R package that facilitates the retrieval and analysis of the Human Protein Atlas data. MC Bioinformatics 20, 463 (2019).https://doi.org/10.1186/s12859-019-3059-z"),
          p("* Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3 319-24277-4, https://ggplot2.tidyverse.org."),
          p("* Wickham H, Francois R, Henry L, Muller K, Vaughan D (2025).dplyr: A Grammar of Data Manipulation. R package version 1.1.4, https://dplyr.tidyverse.org."),
          br()
        )
      ),

      # Right column: Input controls and output
      column(
        width = 8,
        selectInput(
          inputId = ns("cancer_type"),
          label = "Select a cancer type:",
          choices = cancer_tissue_map$cancer,
          selected = cancer_tissue_map$cancer[1]
        ),
        actionButton(ns("run"), "Run Cancer Protein Analysis"),
        br(), br(),
        h3("Protein Rank Shift Plot"),
        plotOutput(ns("plot")),
        h3("Protein Expression Table"),
        verbatimTextOutput(ns("table"))
      )
    )
  )
}

# Server function
cancerServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Define a function to run analysis
    run_cancer_analysis <- function(cancer_selected) {
      compareCancerProtein(cancer_type = cancer_selected)
    }

    # ReactiveValues to store results
    rv <- reactiveValues(result = NULL)

    # Normal Run button
    observeEvent(input$run, {
      req(input$cancer_type)
      rv$result <- run_cancer_analysis(input$cancer_type)
    })

    # Run Example button
    observeEvent(input$run_example, {
      updateSelectInput(session, "cancer_type", selected = "ovarian cancer")
      rv$result <- run_cancer_analysis("ovarian cancer")
    })

    # Outputs
    output$plot <- renderPlot({
      req(rv$result)
      rv$result$plot
    })

    output$table <- renderPrint({
      req(rv$result)
      rv$result$table
    })
  })
}
