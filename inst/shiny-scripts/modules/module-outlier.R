# Outlier Detection Module

# UI function
outlierUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      # Left column with multiple sections as explained in the Introduction Tab + Run Examples
      column(
        width = 4,
        div(
          style = "max-height: 500px; overflow-y: auto; padding: 10px; background-color: #f9f9f9;
                   border-radius: 8px; border: 1px solid #ddd;",
          strong("Description"),
          p("Detects outlier genes present in a given tissue of interest, which have a large delta between their RNA expression and protein expression values."),
          p("These genes are identified as outlier genes, and are most likely candidates for post-transcriptional modifications. This function returns the list of outlier genes as well as a plot highlighting the outlier genes."),
          br(),
          strong("Input:"),
          p("The tissue of interest to analyze."),
          strong("Output:"),
          p("Returns a plot highlighting the outlier genes in red and a vector of outlier gene names detected."),
          br(),
          strong("Instructions:"),
          p("1. From the 'Select a tissue' toolbar, choose an input tissue to detect outliers in."),
          p("2. Once selected, click on 'Run Outlier Detection' to view results."),
          p("3. The output plot is under the heading 'Outlier Genes Plot'."),
          p("4. The output vector is under the heading 'Outlier Gene List'."),
          br(),
          strong("Run Examples:"),
          p("To run the example, follow the steps:"),
          p("1. Click on the 'Run Example' button to view the results."),
          br(),
          # Run Example button
          actionButton(ns("run_example"), "Run Example"),
          br(), br(),
          br(),
          strong("Note:"),
          p("ExprCompareR retrieves data from HPA (Thul et al., 2018) and GTEx (Lonsdale et al., 2013) via packages, HPAanalyze(Tran et al., 2019) and gtexr(Warwick et al., 2025). Because these functions make external API calls, some operations may take time to complete typically around 30 seconds to 1 minute."),
          br(),
          strong("Example Configuration:"),
          br(),
          p("The example that 'Run Example' button runs is:"),
          p("tissue = lung"),
          br(),
          p("All inputs are automatically updated when the example runs."),
          br(),
          strong("Interpretation of plot"),
          p("The scatter plot generated, depicts the outliers detected in the tissue of interest i.e. those genes who have a large delta in their RNA and protein expression. The X axis represents the RNA expression of the genes found in the tissue of interest and the Y axis represents the protein expression of genes found in the tissue of interest. These RNA values are log transformed (log 2) as they are far right skewed as compared to the protein values to have a fair comparison of RNA and protein values as they are obtained from different sources. This function utlizies IQR thresholds to identify outlier genes."),
          p("The black dots indicate genes which do not have a significant change in their RNA vs protein expression, where the red dots indicate outliers."),
          br(),
          strong("Interpretation of list"),
          p("This function provides the list of genes in the tissue of interest which are detected as outliers.These genes are most likely candidates for post-transcriptional modifications, and this list can be used to narrow down potential targets for post-transcriptional modifications. These genes can reflect regulatory mechanism of microRNA-mediated degradation, altered translation rates, or differences in protein stability. If the list returns `character(0)`, then there are no outliers in this tissue."),
          br(),
          strong("References"),
          p("* Attali D, Edwards T (2024). _shinyalert: Easily Create Pretty Popup Messages (Modals) in 'Shiny'_. R package version 3.1.0, https://CRAN.R-project.org/package=shinyalert."),
          p("* Cetinkaya-Rundel M,Cheng J, Grolemund G (2017).Customize your UI with HTML.https://shiny.posit.co/r/articles/build/html-tags/"),
          p("* Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Aden-Buie G, Xie Y,Allen J, McPherson J, Dipert A, Borges B (2025).shiny: Web Application Framework for R. R package version 1.11.1,https://CRAN.R-project.org/package=shiny"),
          p("* OpenAI.(2025). ChatGPT (GPT-5) Large language model.Retrieved November 17, 2025, from https://chatgpt.com/ "),
          p("* Slowikowski K (2024). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. https://ggrepel.slowkow.com/,https://github.com/slowkow/ggrepel."),
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

      # Right column with input controls and outputs
      column(
        width = 8,
        selectInput(
          inputId = ns("tissue"),
          label = "Select a tissue:",
          choices = tissue_map$protein_tissue,
          selected = tissue_map$protein_tissue[1]
        ),
        actionButton(ns("run"), "Run Outlier Detection"),
        br(), br(),
        h3("Outlier Genes Plot"),
        plotOutput(ns("outlier_plot")),
        h3("Outlier Gene List"),
        verbatimTextOutput(ns("outlier_vector"))
      )
    )
  )
}

# Server function
outlierServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Reactive Values to store results
    rv <- reactiveValues(result = NULL)

    # Normal Run button
    observeEvent(input$run, {
      req(input$tissue)
      rv$result <- detect_outliers(input_tissue = input$tissue)
    })

    # Run Example button
    observeEvent(input$run_example, {
      # Automatically select "lung" tissue
      updateSelectInput(session, "tissue", selected = "lung")

      # Run detection immediately on clicking
      rv$result <- detect_outliers(input_tissue = "lung")
    })

    # Outputs
    output$outlier_plot <- renderPlot({
      req(rv$result)
      rv$result$outlier_plot
    })

    output$outlier_vector <- renderPrint({
      req(rv$result)
      rv$result$outlier_vector
    })
  })
}

# [END]
